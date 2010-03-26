#ifndef MCBOX_H
#define MCBOX_H
#include <string> //used in map<string, string> and config reading.
#include <sstream>//string to number convertion

#include <iostream> 
#include <fstream>
#include <cstring>
#include "math\randplus.h"
#include "string\strplus.h"
#include "chain.h"


using std::ofstream;
using std::cout;
using std::endl;

class MCbox_circular{
protected:
	static const int strBufSize=180;
	char buf[strBufSize];
	char filePrefix[80];
	unsigned long anglenums[180];
	unsigned long N_kink[MAXKINKNUM+1];
	unsigned long seeding;
	ofstream fp_log;


public:
    CircularChain dnaChain;
    MCbox_circular(
        char const *r_filePrefix, 
        int const length, 
        double const r_theta_k, 
        double const r_h, 
        unsigned long r_seeding = 23UL);
    virtual ~MCbox_circular(){
        fp_log.close();
    }
    void logParameters(void);
    void logAngleDist(char *suffix = "");
    void clearAngleStats(void);
    void appendAngleStats(void);
    void logAccepts(void);
    void performMetropolisCircularCrankOnly(long monte_step);
};

class MCbox_linear{
protected:
	static const int strBufSize=180;
	char buf[strBufSize];
	char filePrefix[80];
	unsigned long anglenums[180];
	unsigned long N_kink[MAXKINKNUM+1];
	unsigned long seeding;
	ofstream fp_log;
    double crankratio,rotpercentage;


public:
    LinearChain dnaChain;
    MCbox_linear(
        char const *r_filePrefix, 
        int const length, 
        double const r_theta_k, 
        double const r_h, 
        unsigned long r_seeding = 23UL,
        double r_crankratio=0.5,
        double r_rotpercentage=const_rotpercentage);
    virtual ~MCbox_linear(){
        fp_log.close();
    }
    void logParameters(void);
    void logAngleDist(char *suffix = "");
    void getAngleDist(unsigned long *r_dist){
        for (int i=0;i<MAXKINKNUM;i++) r_dist[i]=this->N_kink[i];
    }
    void clearAngleStats(void);
    void appendAngleStats(void);
    void logAccepts(void);
    void logAcceptsEndToEndDistance(void);
    void logAcceptsEndToEndDistanceAngle(void);
    void performMetropolisLinearCrankOnly(long monte_step);
    void performMetropolisLinearCrankHalfChain(long monte_step);
    void performMetropolisLinearTryLigation(long monte_step,
                                        double endToEndDistanceThreshold,
                                        double r_endToEndAngleThresholdDeg);
    double performMetropolisLinearLigationCondP(long monte_step,
                                      double endToEndDistanceThreshold,
                                      double r_endToEndAngleThresholdDeg,
                                      double endToEndDistanceThreshold_ligate,
                                      double r_endToEndAngleThresholdDeg_ligate,
                                      bool conformationBinningOn=false);
    double performSimpleCondP(long monte_step,
                                  double endToEndDistanceThreshold,
                                  double r_endToEndAngleThresholdDeg,
                                  double endToEndDistanceThreshold_ligate,
                                  double r_endToEndAngleThresholdDeg_ligate,
                                  bool conformationBinningOn=false);

};

#endif /* MCBOX_H */