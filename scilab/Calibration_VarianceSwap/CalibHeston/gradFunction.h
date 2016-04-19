#include "utils1.h"
#include <math.h>
#include <stdio.h>

void  gradFunction(int dimx, double *x, double *grad);

void initParamsGrad(int nbdata,DataMarket *_DM,int _TypeNorme,int _TypeModel,int _LogNorme,int _MelangeNormes);

double dNormeMarket(int TypeNorme,int TypeOpt,double St0,double T,double K,double r,double divid,double prix,double prixobs,double sigma_imp_obs,int dimx,double *gradi);

void FreeGradFunction();

