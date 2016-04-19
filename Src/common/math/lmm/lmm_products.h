
#ifndef LMM_PRODUCTS_H
#define LMM_PRODUCTS_H


#include"lmm_header.h"

extern int mallocSwaption(Swaption **ptSwpt ,double swptMat , double swpMat ,double vola , double K, double tenor);
extern void freeSwaption(Swaption **ptSwpt);
extern int printSwaption(Swaption *pt);
extern int mallocCaplet(Caplet **ptCplt ,double Mat , double priceVal , double K);
extern void freeCaplet(Caplet **ptCplt);
extern int printCaplet(Caplet *pt);


#endif

