#ifndef MCBOX_H
#define MCBOX_H
#include <string> //used in map<string, string> and config reading.
#include <sstream>//string to number convertion

#include <iostream> 
#include <fstream>
#include <cstring>
#include "math/randplus.h"
#include "string/strplus.h"
#include "chain.h"
#include "file/configldr.h"
#include <map>
#include <sstream>


using std::ofstream;
using std::cout;
using std::endl;



class MCbox_circular{
protected:
	static const int strBufSize=180;
	char buf[strBufSize];
	char filePrefix[80];
	unsigned long anglenums[180];
	unsigned long seeding;
	ofstream * fp_log;

private:
	MCbox_circular ();

public:
    CircularChain * dnaChain;
	std::map <std::string,std::string> config;
    MCbox_circular(char const * configFile);
	virtual ~MCbox_circular(){
        (*fp_log).close();
    }
    void logParameters(void);
    void logAngleDist(char *suffix = "");
    void clearAngleStats(void);
    void pushAngleStats(void);
    void logAccepts(void);
	double calcGyration(void);
    void performMetropolisCircularCrankRept(long monte_step);
};

#endif /* MCBOX_H */