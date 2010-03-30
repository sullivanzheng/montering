#ifndef CHAIN_H
#define CHAIN_H


#include <iostream>
#include <fstream>
#include <cmath>
#include "math\mtrand.h"
#include "chainconst.h"
#include "chaingbvar.h"
#include "math\ezlin.h"
#include "SmallTypes.h"

double G_b(double ang);
double G_t(int Lk,double L_bp,int kink_num);
double P_t_over_2t(double L_bp,int kinkNum);

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
		counter auto_accepts;
		counter auto_moves;
		signed int kink_num;
        statqueue <double> auto_prd_endToEndDistance2;
        statqueue <double> auto_prd_endToEndAngle;
		statqueue <double> gyration_ratio;
		statqueue <double> anglelist[maxa];

        void resetStat(){
            auto_accepts.lap();
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
	virtual int updateAllBangleKinkNum() = 0;
	void updateAllBangleKinkNum_Ini(bool circular);
	virtual int updateBangleKinkNum(int i) = 0;

    void normalize();

public:
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
	virtual int crankshaft(int m, int n, double a, bool trialMoveFlag=false) = 0;
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
	virtual double deltaE_TrialCrankshaft(int m, int n, double a) = 0;
	virtual void snapshot(char *filename);
};

class CircularChain: public Chain{
public:
	int Lk;
protected:
	virtual int updateAllBangleKinkNum();
	virtual int updateBangleKinkNum(int i);
	void driftProof();
private:
	CircularChain();
public:
	CircularChain(int length);
	CircularChain(char const *filename,int length);
	virtual double calG_bSum();
	virtual int crankshaft(int m, int n, double a,bool trialMoveFlag=false);
	virtual double deltaE_TrialCrankshaft(int m, int n, double a);
	virtual void snapshot(char *filename);
};



#endif /* CHAIN_H */