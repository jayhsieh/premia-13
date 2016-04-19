#include <malloc.h>
#include <stdio.h>
#include <math.h>
#ifndef _DATAMARKET_
#define _DATAMARKET_

typedef struct DataMarket 
{
  int nbdata;
  int *TypeOpt;
  double *St0;
  double *K;
  double *T;
  double *r;
  double *d;
  double *wi;
  double *PrixObs;
  double *SigmaImpObs;
  double *PrixMod;
  double *SigmaImp;
}DataMarket;

DataMarket * CreateDataMarket(int N);
void AfficheDataMarket(DataMarket *DM);
void FreeDataMarket(DataMarket *DM);


#endif
