#ifndef RANDPLUS_H
#define	RANDPLUS_H
#include <cstdlib> //for rand() support,long(RAND_MAX)=2 bytes.
#include <cmath>

#if 0
#include "math\mtrand.h"
static MTRand53 glmt(2351234UL);

inline double drand(double d_ceiling){//[0,ceiling)
	return double(floor(glmt()*(long(RAND_MAX)+1)))/double(long(RAND_MAX)+1)*d_ceiling;
}
inline double drand(double d_floor,double d_ceiling){//[floor,ceiling)
	return double(floor(glmt()*(long(RAND_MAX)+1)))/double(long(RAND_MAX)+1)*(d_ceiling-d_floor)+d_floor;
}
inline long irand(long i_ceiling){//integer random number [0, ceiling-1]
	return long(floor(drand(double(i_ceiling))));
}

inline long irand(long i_floor, long i_ceiling){//integer random number [floor, ceiling-1]
    return irand(i_ceiling-i_floor)+i_floor;
}
#endif
#define RANDRANGE (double(RAND_MAX)+1.0)

inline double drand(double d_ceiling){//[0,ceiling)
	return double(rand())/RANDRANGE*d_ceiling;
}
inline double drand(double d_floor,double d_ceiling){//[floor,ceiling)
	return double(rand())/RANDRANGE*(d_ceiling-d_floor)+d_floor;
}
inline long irand(long i_ceiling){//integer random number [0, ceiling-1]
	return long(floor(drand(double(i_ceiling))));
}

inline long irand(long i_floor, long i_ceiling){//integer random number [floor, ceiling-1]
    return irand(i_ceiling-i_floor)+i_floor;
}

#endif /* RANDPLUS_H */