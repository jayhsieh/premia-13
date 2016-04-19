#include "optype.h"
#include "./utils/datamarket.h"
#include "parameter.h"
#include "./utils/bs.h"
#include "./svj_std/svj.h"
//
double costFunction(int dimx, double *x);
void initParams(int nbdata,DataMarket *_DM,int _TypeNorme,int _TypeModel,int _LogNorme,int _MelangeNormes);
double NormeMarket(int TypeNorme,int TypeOpt,double St0,double T,double K,double r,double divid,double prix,double prixobs,double sigma_imp_obs,double *sigmaimp);
void checksigmaimp();
void FreeCostFunction();
void AffichePrixDataMarket();
//
