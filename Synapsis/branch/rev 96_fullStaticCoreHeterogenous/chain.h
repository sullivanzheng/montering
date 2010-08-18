#ifndef CHAIN_H
#define CHAIN_H

#define MAXMatrixDet 200

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstring>
#include <sstream>
#include "math/mtrand.h"
#include "chainconst.h"
#include "chaingbvar.h"
#include "math/ezlin.h"
#include "string/strplus.h"
#include "SmallTypes.h"


inline long wrap(long i, long roundnum){
	if (i<-roundnum || i>2*roundnum-1){
		std::cout<<"[In global function 'wrap']"
			"Wrapping error, wrapping more than one round:"<<i<<std::endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	if (i<0)
		return i+roundnum;
	if (i>roundnum-1)
		return i-roundnum;
	else
		return i;
}


class CircularChain; //Finish declaration of Chain for the use of declaration of its nest friend class rigid;

class cls_rigid{
private: 
	cls_rigid();
public:
	CircularChain* target;
	std::vector<long> protect;
	std::vector<std::vector<double> > ref_v;
	std::vector<std::vector<double> > ref_v_xyz;
	cls_rigid(CircularChain * r_target, std::vector<long>r_protect,
		std::vector<std::vector<double> > r_ref_v);
	void update_ref_v_xyz();
};

class allrigid {
private:
	allrigid();
public:
	struct sphere{
		int x0,x1;
		int rg0,rg1;
		double d;
	};
	double E;
	double AxisBeta;
	double RadiusBeta;
	double r;
	double r_siteI;
	std::vector<cls_rigid> R;
	allrigid(char *configfile, CircularChain *target);
	std::vector<long> protect;
	std::vector<sphere> spheres;
	double update_allrigid_and_E();
	int IEV_spheres(long m, long n);
};

class Chain {
public:
	struct stct_stat
    {
		//Statistical variable for the Chain.
		/*Variable names starting with auto_ is 
		  automatically updated after each move
		  by calling 
  			this->dnaChain->auto_updt_stats();

		  Variable names starting with auto_prd 
		  is automatically AND periodically updated
		  after certain number of moves.
		  This number can be decided by defaultSampleCycle and
		  passed to endToEndSampleCycle. */
		counter crk_accepts;
		counter rpt_accepts;
		counter auto_moves;
        statqueue <double> auto_prd_endToEndDistance2;
        statqueue <double> auto_prd_endToEndAngle;
		statqueue <double> gyration_ratio;
		statqueue <double> anglelist[maxa];

        void resetStat(){
            crk_accepts.lap();
			rpt_accepts.lap();
            auto_moves.lap();
            auto_prd_endToEndDistance2.clear();
            auto_prd_endToEndAngle.clear();
			gyration_ratio.clear();
			for (long i=0;i<maxa;i++) anglelist[i].clear();
		}        
	}stats;

	struct segment{
		double x,y,z;
		double dx,dy,dz;
		double l; //segment length
		double bangle;
    }C[maxa];//C stands for "Chain"

protected:
    long totsegnum;
	double contour_length;
	double max_seglength;
    static const long defaultSampleCycle=100;
    long endToEndSampleCycle;
	static const long NORMALIZE_PERIOD=100000;
	long readIniFile(char const *filename);
    void initializeCircle(long num);
	double calAngle(segment &C1, segment &C2);
	void SetRotM_crankshaft(double M[3][3],long m, long n, double a);
	void SetRotM_halfchain(double M[3][3], double rv[3], double a);
	virtual long updateAllBangle() = 0;		
    void updateAllBangle_Ini(bool circular);		
	virtual double updateBangle(long i) = 0;
    void normalize();

public:
	double VolEx_R;
	void auto_updt_stats(){
        stats.auto_moves++;
		if (stats.auto_moves() % endToEndSampleCycle==0) {
            double temp=this->getEndToEndDistance();
            stats.auto_prd_endToEndDistance2.push(temp*temp);
            stats.auto_prd_endToEndAngle.push(this->getEndToEndAngle());
        }
        if (stats.auto_moves() % NORMALIZE_PERIOD == 0){
            this->normalize();
        }
    }
private:
	explicit Chain(void) {};

public:
	explicit Chain(char const *filename, bool circular, long r_totsegnum);
	long dispChainCord();
	virtual double calG_bSum() = 0;
	virtual long crankshaft(long m, long n, double a) = 0;
    inline double getEndToEndAngle(void){
        static segment Ci={0,0,0,0,0,0,0},Cf={0,0,0,0,0,0,0};
        static double endToEndAngleMemo=0.0;
        if (Ci.dx != C[0].dx || Ci.dy != C[0].dy || Ci.dz != C[0].dz 
            || Cf.dx != C[maxnum].dx || Cf.dy != C[maxnum].dy || Cf.dz != C[maxnum].dz )
        {
            Ci=C[0];Cf=C[maxnum];
            endToEndAngleMemo=this->calAngle(C[0],C[maxnum]);
            return endToEndAngleMemo;
        }
        else
        {
            return endToEndAngleMemo;
        }
    }
    inline double getEndToEndDistance(void){
        return modu(C[maxnum].x+C[maxnum].dx-C[0].x,
                    C[maxnum].y+C[maxnum].dy-C[0].y,
                    C[maxnum].z+C[maxnum].dz-C[0].z);
    }
	virtual double dE_TrialCrankshaft(long m, long n, double a) = 0;
	virtual void snapshot(char *filename);
};

class CircularChain: public Chain{

protected:
	void driftProof();
	virtual long updateAllBangle();
	virtual double updateBangle(long i);
public:
	double writhe;
	double topl;

	double E_t; //torsional energy.
	double dLk;

	long productLk(long vertM, long vertN);
	
	CircularChain();
	CircularChain(char const *filename,long totsegnum);
	virtual double G_b(long n);
	virtual double getg(long n);
	virtual double calG_bSum();
	virtual long crankshaft(long m, long n, double a);
	virtual double dE_reptation(long m, long n, long move);
	virtual double dE_TrialCrankshaft(long m, long n, double a);
	int kpoly(long ial[2],long ierr);
	virtual void snapshot(char *filename);
	long IEV_closeboundary(long in, long ik);
	long IEV_Alex_closeboundary(long in, long ik, double info[3]);
	long IEV_with_rigidbody_closeboundary( long in,  long ik, double info[3]);
	double E_t_updateWrithe_E_t(); //Based on _fastWr_topl_update();
	long checkConsistancy();
	long getBranchNumber();
	long scanBranch(char* filename);
private:
	double _bwr(long m, long n);
	long _kndwr_topl_update(double & topl, long & ierr);
	long _kndwr(long &ierr);
	double _fastWr_topl_update();
	double _det(long n, double da[MAXMatrixDet*MAXMatrixDet]);

public:
	double _wrfun(long m, long n);

};
#endif /* CHAIN_H */