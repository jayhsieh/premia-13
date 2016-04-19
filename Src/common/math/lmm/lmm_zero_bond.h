
#ifndef LMM_ZERO_BOND_H
#define LMM_ZERO_BOND_H


#include"lmm_header.h"


extern int mallocZeroBond( ZeroBond **ptZb ,int numOfMat, double tenorVal );
extern int freeZeroBond(ZeroBond **ptZb);
extern int printZeroBond(ZeroBond *ptZb);


#endif

