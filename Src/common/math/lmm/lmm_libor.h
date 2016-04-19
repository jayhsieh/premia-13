
#ifndef LMM_LIBOR_H
#define LMM_LIBOR_H


#include"lmm_header.h"


extern int mallocLibor(Libor **ptLib,int numOfMat, double tenorVal,double l0);
extern void freeLibor(Libor **ptLib);
extern int initLibor(Libor *ptLib,double l0);
extern int printLibor(Libor *ptLib);
extern int putLiborToZero(Libor *ptLib, int index);
extern int copyLibor(Libor *ptLibSrc , Libor *ptLibDest );
extern double computeZeroCouponSum(Libor* ptLib, int o, int s,int m );
extern double computeSwapRate(Libor* ptLib, int o, int s,int m );
extern double computeSwapPrice(Libor* ptLib, Swaption* ptSwp,int o, int s,int m );
extern int mallocHistLibor(HistLibor **ptLib ,int numOfMat, double tenorVal, int numberOfTraj );
extern double computeZeroCoupon(Libor* ptLib, int o,int s);


void Libor_To_ZeroCoupon(Libor* ptLib, PnlVect* zc); // Compute P(0, Ti) i=0:N
double Sum_ZC(Libor* ptLib, int i_first, int i_last); // Compute "sum P(0, T_i)" for "i" from "i_first" to "i_last".

#endif

