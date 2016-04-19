#ifndef _VARSWAPS_H_
#define _VARSWAPS_H_

#include "utils1.h"

typedef struct VSMarket 
{
  int nbdata;
  double *T;
  double *wi;
  double *PriceObs;
  double *PriceMod;
}VSMarket;

VSMarket * CreateVSMarket(int N);
void PrintVSMarket(VSMarket *VS);
void FreeVSMarket(VSMarket *VS);
void readVSMarket(VSMarket *VS);

double VS_Heston(double T, double V0, double kappa, double theta);

double costFunctionVS(int dimx, double *x);
void  gradFunctionVS(int dimx, double *x, double *grad);
void initParamsVS(int nbdata, VSMarket *_VS, int _TypeNorme, int _TypeModel);
int fitEquationVS(int dimx, double *x);
void CreateSyntheticVSMarket(VSMarket *VS, int nt, double *TS, double *x);

void FreeCostFunctionVS();
void PrintVSPrices();
void PrintVSMarketFile();
//
#endif
