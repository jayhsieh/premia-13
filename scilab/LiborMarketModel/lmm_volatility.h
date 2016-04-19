
#ifndef LMM_VOLATILITY_H 
#define LMM_VOLATILITY_H

#include"lmm_header.h"

extern int mallocVolatility(Volatility** ptVol , int numOfFac );
extern double evalVolatility( Volatility* ptVol ,int factorNumber, double t, double T);
extern int freeParametre(Volatility **ptVol);               
extern int initVolatility(Volatility *ptVol);
extern int computeEigenVector(Volatility *ptVol);
extern int shiftVolatility(Volatility *ptVol);
extern int computeCovariance(Volatility *ptVol);
extern int copyVolatility(Volatility *ptSrc, Volatility *ptDest);
extern int mallocVolatilityInteger(Volatility **ptVol , int numOfFac, float  tenor );

#endif

