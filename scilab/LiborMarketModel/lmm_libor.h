
#ifndef LMM_LIBOR_H 
#define LMM_LIBOR_H


#include"lmm_header.h"

extern int mallocLibor(Libor **ptLib ,int numOfMat, double tenorVal );
extern int freeLibor(Libor **ptLib);
extern int initLibor(Libor *ptLib);
extern int printLibor(Libor *ptLib);
extern int putLiborToZero(Libor *ptLib, int index);
extern double computeZeroCouponSum(Libor* ptLib, int o, int s,int m );
extern double computeSwapRate(Libor* ptLib, int o, int s,int m );
extern double computeSwapPrice(Libor* ptLib, Swaption* ptSwp,int o, int s,int m );
extern int mallocHistLibor(HistLibor **ptLib ,int numOfMat, double tenorVal, int numberOfTraj );
extern double computeZeroCoupon(Libor* ptLib, int o,int s);

#endif  





