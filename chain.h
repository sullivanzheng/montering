#ifndef CHAIN_H
#define CHAIN_H


#include <iostream>
#include <fstream>
#include <cmath>
#include "math\mtrand.h"
#include "chainconst.h"
#include "chaingbvar.h"
#include "math\ezlin.h"

double G_b(double ang);
double G_t(int Lk,double L_bp,int kink_num);
double P_t_over_2t(double L_bp,int kinkNum);
class counter{
private:
	unsigned long num;
	unsigned long lastreset;
public:
	counter(void){num=lastreset=0;}
	void clear(void){num=lastreset=0;}
	//inline void count(void){num++;}
	inline void operator++ (int){num++;}
	inline void reset(void){lastreset=num;}
	inline unsigned long operator = (unsigned long reset_sym){
		try{
			if (reset_sym!=0)
				throw ("Only 0 can be assigned to an couter object.");
		reset();
		}
		catch (char * msg) {std::cout << msg;}
        return 0;
	}
    inline unsigned long getNumber(void){return num-lastreset;}
	inline unsigned long operator () (void){return num-lastreset;}
	inline unsigned long getTotCounts(void){return num;}
};

template <typename TypeName> class statqueue{
public:
    long itemNum;
    TypeName squareSum;
    TypeName sum,
        maxItem,minItem;
    statqueue(){
        itemNum=0;
        squareSum=sum=0;
        maxItem=minItem=0;
    }
    void clear(){
        itemNum=0;
        squareSum=sum=0;
        maxItem=minItem=0;
    }
    void push(TypeName r_data){
        itemNum++;
        sum+=r_data;
        squareSum+=(r_data*r_data);
        if (itemNum==1){
            maxItem=minItem=r_data;
        }
        else{
            maxItem=r_data>maxItem?r_data:maxItem;
            minItem=r_data<minItem?r_data:minItem;
        }
    }
    TypeName getMean(){
        return sum/itemNum;
    }
    TypeName getVar(){
        if (itemNum<=1) {
            std::cout<<"statequeue:Could not return variance when itemNum<=1"<<endl;
            getchar();exit(EXIT_FAILURE);
        }
        return (squareSum-(sum*sum)/itemNum)/(itemNum-1);
    }
    TypeName getStdev(){
        return sqrt(getVar());
    }
    TypeName getStdevOfMean(){
        return sqrt(getVar())
            /sqrt(double(itemNum));
    }
    TypeName getmaxItem(){
        if (itemNum<1) {
            std::cout<<"statequeue:Could not return min:max since it is empty"<<endl;
            getchar();exit(EXIT_FAILURE);
        }
        return maxItem;
    }
    TypeName getminItem(){
        if (itemNum<1) {
            std::cout<<"statequeue:Could not return min:max since it is empty"<<endl;
            getchar();exit(EXIT_FAILURE);
        }
        return minItem;
    }
};
class Chain {
public:
	struct
    {
		counter accepts;
		counter moves;
		signed int kink_num;
        statqueue <double> endToEndDistance2;
        statqueue <double> endToEndAngle;
        //For Chain of Conditional Probabilities ONLY
        long ligateAccepts;
        long halfChainSuccess;
        long crankSuccess;

        void resetStat(){
            this->accepts.reset();
            this->moves.reset();
            this->endToEndDistance2.clear();
            this->endToEndAngle.clear();
            this->ligateAccepts=0;
            this->halfChainSuccess=0;
            this->crankSuccess=0;
        }
        /////////
	}stats;

	struct segment{
		double x,y,z;
		double dx,dy,dz;
		double bangle;
    }C[maxa];//C stands for "Chain"

protected:
    int length;
    static const long defaultSampleCycle=100;
    long endToEndDistanceSampleCycle;
	static const int NORMALIZE_PERIOD=100000;
	int readIniFile(char const *filename);
    void readIniFile(int num);
	double calAngle(segment &C1, segment &C2);
	void SetRotM_crankshaft(double M[3][3],int m, int n, double a);
	void SetRotM_halfchain(double M[3][3], double rv[3], double a);
	virtual int updateAllBangleKinkNum() = 0;
	void updateAllBangleKinkNum_Ini(bool circular);
	virtual int updateBangleKinkNum(int i) = 0;

    void normalize();

public:
    void countmove(){
        stats.moves++;
        if (stats.moves() % endToEndDistanceSampleCycle==0) {
            double temp=this->getEndToEndDistance();
            stats.endToEndDistance2.push(temp*temp);
            stats.endToEndAngle.push(this->getEndToEndAngle());
        }
        if (stats.moves() % NORMALIZE_PERIOD == 0){
            this->normalize();
        }
    }
	explicit Chain(void){}
    explicit Chain(bool circular, int r_length);
	explicit Chain(char const *filename, bool circular, int r_length);
	int dispChainCord();
	void clearAccepts();
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
	void snapshot(char *filename);
};

class CircularChain: public Chain{
public:
	int Lk;
protected:
	virtual int updateAllBangleKinkNum();
	virtual int updateBangleKinkNum(int i);
	void driftProof();
public:
	CircularChain(){}
	CircularChain(char const *filename,int length);
	virtual double calG_bSum();
	virtual int crankshaft(int m, int n, double a,bool trialMoveFlag=false);
	virtual double deltaE_TrialCrankshaft(int m, int n, double a);
};

class LinearChain: public Chain{
public:
	virtual int updateAllBangleKinkNum();
protected:
	virtual int updateBangleKinkNum(int i);

public:
    
    LinearChain(){}
    LinearChain(int length);
	LinearChain(char const *filename,int length);
	virtual double calG_bSum();
	virtual int crankshaft(int m, int n, double a,bool trialMoveFlag=false);
	virtual double deltaE_TrialCrankshaft(int m, int n, double a);
	virtual int halfChain(int m, double rv[3],double a,bool trialMoveFlag=false);
	virtual double deltaE_TrialHalfChain(int m,double rv[3],double a);
    bool trialLigateAfterHalfChainOK(int m, double rv[3],double a,
       double endToEndDistanceThreshold,
       double endToEndAngleThreshold);
    bool trialLigateAfterCrankshaftOK(int m, int n, double a,
       double endToEndDistanceThreshold,
       double endToEndAngleThreshold);
};

#endif /* CHAIN_H */