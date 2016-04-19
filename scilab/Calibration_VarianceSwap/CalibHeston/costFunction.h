#ifndef _COSTFUNCTION_H_
#define _COSTFUNCTION_H_

#include "utils1.h"

int FT_Call_Heston(double s, double strike, double t, double ri, double dividi, double sigma0,double ka,double theta,double sigma2,double rhow,double *ptprice, double *ptdelta);
int FT_Put_Heston(double s, double strike, double t, double ri, double dividi, double sigma0,double ka,double theta,double sigma2,double rhow,double *ptprice, double *ptdelta);

double costFunction(int dimx, double *x);

void initParams(int nbdata, DataMarket *_DM,int _TypeNorme,int _TypeModel,int _LogNorme,int _MelangeNormes);

int CreateSyntheticDM(DataMarket *DMS, int nt, double *_T, double *ri, double *di, double St0, int nk, double *K, int *TypeOpt, int dimx, double *x);

double NormeMarket(int TypeNorme,int TypeOpt,double St0,double T,double K,double r,double divid,double prix,double prixobs,double sigma_imp_obs,double *sigmaimp);
void checksigmaimp();
void FreeCostFunction();
void AffichePrixDataMarket();
void AfficheDataMarketFile();
//
#endif
