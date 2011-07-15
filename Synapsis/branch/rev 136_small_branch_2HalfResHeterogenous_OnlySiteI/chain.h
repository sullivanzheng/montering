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
#include "bias.h"

int printmtrx(long n, double *a);

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
	std::vector<long> ref_vec_basis;
	long anchor;
	cls_rigid(CircularChain * r_target, 
		      std::vector<long>r_protect,
			  std::vector<std::vector<double> > r_ref_v, 
			  std::vector<long>r_ref_vec_basis,
			  long r_anchor);
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
	struct sphere_simple{
		double x;
		double y;
		double z;
		double r;
	};
	double E;
	double unbiasedE;
	double Q;

	//Site II and III comformation parameters
	double AxisBeta;
	double RadiusBeta;
	double r;

	//Site I comformation parameters
	double r_siteI;
	double r_siteI_deviation;
	double siteI_direction;

	std::vector<cls_rigid> R;
	allrigid(char *configfile, CircularChain *target);
	std::vector<long> protect;
	std::vector<sphere> spheres;
	std::vector<sphere_simple> spheres_simple;
	double update_allrigid_and_E();
	int IEV_spheres(long m, long n);
};

struct segment{
		double x,y,z;
		double dx,dy,dz;
		double l; //segment length
		double bangle;
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
		counter crk_counts;
		counter rpt_accepts;
		counter rpt_counts;
		counter rpt_rejection_count;
		counter rpt_rejection_quickrej;
		counter rptsimp_accepts;
		counter rptsimp_counts;
		counter tdm_counts;
		counter tdm_accepts;
		counter auto_moves;
        statqueue <double> auto_prd_endToEndDistance2;
        statqueue <double> auto_prd_endToEndAngle;
		statqueue <double> gyration_ratio;
		statqueue <double> anglelist[maxa];

        void resetStat(){
            crk_accepts.lap();
			crk_counts.lap();
			rpt_accepts.lap();
			rpt_counts.lap();
			rpt_rejection_count.lap();
			rpt_rejection_quickrej.lap();
			rptsimp_accepts.lap();
			rptsimp_counts.lap();
			tdm_counts.lap();
			tdm_accepts.lap();
            auto_moves.lap();
            auto_prd_endToEndDistance2.clear();
            auto_prd_endToEndAngle.clear();
			gyration_ratio.clear();
			for (long i=0;i<maxa;i++) anglelist[i].clear();
		}        
	}stats;

segment C[maxa];//C stands for "Chain"

protected:
    long totsegnum;
	double contour_length;
	double max_seglength;
    static const long defaultSampleCycle=100;
    long endToEndSampleCycle;
	static const long NORMALIZE_PERIOD=1000;
	long readIniFile(char const *filename);
    void initializeCircle(long num);
	double calAngle(segment &C1, segment &C2);
	void SetRotM_crankshaft(double M[3][3],long m, long n, double a);
	void SetRotM_halfchain(double M[3][3], double rv[3], double a);
	virtual long normalizeAllBangle() = 0;		
    void updateAllBangle_Ini(bool circular);		
	virtual double updateBangle(long i) = 0;

public:
	void normalize_X_bangle();

public:
	double VolEx_R;
	void auto_updt_stats(){
        stats.auto_moves++;
		/*if (stats.auto_moves() % endToEndSampleCycle==0) {
            double temp=this->getEndToEndDistance();
            stats.auto_prd_endToEndDistance2.push(temp*temp);
            stats.auto_prd_endToEndAngle.push(this->getEndToEndAngle());
        }*/
		if (stats.auto_moves.getTotCounts() % NORMALIZE_PERIOD == 0){
            this->normalize_X_bangle();
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
	virtual long normalizeAllBangle();
	virtual double updateBangle(long i);
	double _adjustBangle(long m, long dm, double newBangle, 
		segment const C2[maxa], segment Ctemp[maxa]);
	double _bisectBangle(long m, long dm, 
		double a1, double f1, double a2, double f2, double length,
		segment const C2[maxa], segment Ctemp[maxa], double eps);

	int _deformReptSegments_updateInternalBangles(long m, long dm, double length, double Lnow, segment C2[maxa]);
	int _deformReptSegments_updateInternalBangles_Normalize(long m, long dm, double length, double Lnow, segment C2[maxa]);
public:
	double writhe;
	double topl;

	double E_t; //torsional energy.
	double dLk;

	long productLk(long vertM, long vertN);
	long productLk_fast(long vertM, long vertN);
	long productLk2(long vertM, long vertN, long s, long t);
	long overpassing(long vertM, long vertN);
	
	CircularChain();
	CircularChain(char const *filename,long totsegnum);
	virtual double G_b(long n);
	virtual double getg(long n);
	virtual double calG_bSum();
	virtual long crankshaft(long m, long n, double a);
	virtual double dE_reptation_simple(long m, long n, long move);
	virtual double dE_reptation_3_4(long m1, long dm1, long m2, long dm2, int& rejection_sign);
	virtual double dE_TrialCrankshaft(long m, long n, double a);
	virtual double dE_treadmill(double direction);

	int kpoly(long ial[2],long ierr);
	int kpoly2(long ial[2],long ierr);
	int AP(long vertM, long vertN,double s, double t);

	virtual int snapshotseg(char *filename, segment const * Ct, int start, int end);
	virtual void snapshot(char *filename);

	long IEV_closeboundary(long in, long ik);
	long IEV_Alex_closeboundary(long in, long ik, double info[3]);
	long IEV_with_rigidbody_closeboundary( long in,  long ik, double info[3]);
	long IEV_with_rigidbody_closeboundary_fullChain(double info[3]);

	double E_t_updateWrithe_E_t(); //Based on _fastWr_topl_update();
	long checkConsistency(double eps=1e-10);
    long checkBangleConsistency();
	long getBranchNumber();
	long scanBranch(char* filename);
private:
	double _bwr(long m, long n);

public:
	long _kndwr_topl_update(double & topl, long & ierr);
	long _kndwr_topl_update_stable(double & topl, long & ierr);
private:
	long _kndwr(long &ierr);
	double _fastWr_topl_update();
	double _det(long n, double da[MAXMatrixDet*MAXMatrixDet]);

public:
	double _wrfun(long m, long n);

};
#endif /* CHAIN_H */