#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdio.h>

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
void PrintDataMarket(DataMarket *DM);
void PrintDataMarketFile(DataMarket *DM);

void affiche_sol(int dimx,double *sol,double *x,double fx,double *grad,double *grad0);
void faffiche_sol(FILE *ftest,int dimx,double *x,double fx,double *grad);
void init_bornes(int TypeModel,int dimx,int *nbd,double *xmin,double *xmax);

void readDataMarket(DataMarket *DM);

void read_params(int *TypeModel,int *typen,int *nbdata, int *nvsdata, double *xinit);

#endif
