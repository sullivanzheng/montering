#ifndef RANDPLUS_H
#define	RANDPLUS_H
#include <cstdlib> //for rand() support,RAND_MAX=2 bytes.
#include <cmath>
#include "math\mtrand.h"
#if 0
static MTRand53 glmt(2351234UL);

inline double drand(double d_ceiling){//[0,ceiling)
	return double(floor(glmt()*(RAND_MAX+1)))/double(RAND_MAX+1)*d_ceiling;
}
inline double drand(double d_floor,double d_ceiling){//[floor,ceiling)
	return double(floor(glmt()*(RAND_MAX+1)))/double(RAND_MAX+1)*(d_ceiling-d_floor)+d_floor;
}
inline int irand(int i_ceiling){//integer random number [0, ceiling-1]
	return int(floor(drand(double(i_ceiling))));
}

inline int irand(int i_floor, int i_ceiling){//integer random number [floor, ceiling-1]
    return irand(i_ceiling-i_floor)+i_floor;
}
#endif

inline double drand(double d_ceiling){//[0,ceiling)
	return double(rand())/double(RAND_MAX+1)*d_ceiling;
}
inline double drand(double d_floor,double d_ceiling){//[floor,ceiling)
	return double(rand())/double(RAND_MAX+1)*(d_ceiling-d_floor)+d_floor;
}
inline int irand(int i_ceiling){//integer random number [0, ceiling-1]
	return int(floor(drand(double(i_ceiling))));
}

inline int irand(int i_floor, int i_ceiling){//integer random number [floor, ceiling-1]
    return irand(i_ceiling-i_floor)+i_floor;
}


#endif /* RANDPLUS_H */