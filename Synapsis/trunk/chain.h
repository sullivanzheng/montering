#ifndef CHAIN_H
#define CHAIN_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "math\mtrand.h"
#include "chainconst.h"
#include "chaingbvar.h"
#include "math\ezlin.h"
#include "string\strplus.h"
#include "SmallTypes.h"

double G_b(double ang);

class CircularChain; //Finish declaration of Chain for the use of declaration of its nest friend class rigid;

class cls_rigid{
private: 
	cls_rigid();
public:
	CircularChain* target;
	std::vector<int> protect;
	std::vector<std::vector<double> > ref_v;
	std::vector<std::vector<double> > ref_v_xyz;
	cls_rigid(CircularChain * r_target, std::vector<int>r_protect,
		std::vector<std::vector<double> > r_ref_v);
	void update_ref_v_xyz();
};

class allrigid{
private:
	allrigid();
public:
	double E;
	std::vector<cls_rigid> R;
	allrigid(char *configfile, CircularChain *taget);
	std::vector<int> protect;
	double update_allrigid_and_E();
};


class Chain {
public:
	struct stct_stat
    {
		//Statistical variable for the Chain.
		/*Variable names starting with auto_ is 
		  automatically updated after each move.

		  Variable names starting with auto_prd 
		  is automatically AND periodically updated
		  after certain number of moves.
		  This number can be decided by defaultSampleCycle and
		  passed to endToEndSampleCycle. */
		counter accepts;
		counter auto_moves;
        statqueue <double> auto_prd_endToEndDistance2;
        statqueue <double> auto_prd_endToEndAngle;
		statqueue <double> gyration_ratio;
		statqueue <double> anglelist[maxa];

        void resetStat(){
            accepts.lap();
            auto_moves.lap();
            auto_prd_endToEndDistance2.clear();
            auto_prd_endToEndAngle.clear();
			gyration_ratio.clear();
			for (int i=0;i<maxa;i++) anglelist[i].clear();
		}        
	}stats;

	struct segment{
		double x,y,z;
		double dx,dy,dz;
		double bangle;
    }C[maxa];//C stands for "Chain"

protected:
    int length;
    static const long defaultSampleCycle=100;
    long endToEndSampleCycle;
	static const int NORMALIZE_PERIOD=100000;
	int readIniFile(char const *filename);
    void initializeCircle(int num);
	double calAngle(segment &C1, segment &C2);
	void SetRotM_crankshaft(double M[3][3],int m, int n, double a);
	void SetRotM_halfchain(double M[3][3], double rv[3], double a);

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
	Chain::Chain(bool circular,int r_length);
	explicit Chain(char const *filename, bool circular, int r_length);
	int dispChainCord();
	virtual double calG_bSum() = 0;
	virtual int crankshaft(int m, int n, double a) = 0;
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
	virtual double deltaE_TrialCrankshaft_countMove(int m, int n, double a) = 0;
	virtual void snapshot(char *filename);
};

class CircularChain: public Chain{

protected:
	void driftProof();
public:
	CircularChain();
	CircularChain(int length);
	CircularChain(char const *filename,int length);
	virtual double calG_bSum();
	virtual int crankshaft(int m, int n, double a);
	virtual double deltaE_TrialCrankshaft_countMove(int m, int n, double a);
	virtual void snapshot(char *filename);
	int IEV( int in, int ik);
};
#endif /* CHAIN_H */