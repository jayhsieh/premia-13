#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "./utils/datamarket.h"
#include "optype.h"
#include "../svj_std/svj.h"
#include "costFunction.h"
#include "gradFunction.h"
#include "utils/plot.h"
#include "bs.h"
//
double uniforme(void);
void affiche_sol(int dimx,double *sol,double *x,double fx,double *grad,double *grad0);
void faffiche_sol(FILE *ftest,int dimx,double *sol,double *x,double fx,double *grad,double *grad0);
double norme(int dimx,double *x);
double dist(int dimx,double *x,double *y);

void init_sol(int TypeModel,int *dimx,double *x,int dimx1,int dimx2,int dimx3,double V0,double kappa,double theta,double sigmav,double rho,double lambda, double m0,double v);

void init_bornes(int TypeModel,int dimx,int *nbd,double *xmin,double *xmax);

void init_PrixObs(int nbdata1,int nbdata2,int TypeModel,double Kmin,double Kmax,double Tmin,double Tmax,double St0,double r,double divid,double V0,double kappa,double theta,double sigmav,double rho,double lambda, double m0,double v);

void trace_1d(int iplot1,int iplot2,int nbplot,int dimx,double *x0,double *y0,double *y1,double *xmin,double *xmax);

void trace_3d(int iplot1,int iplot2,int nbplot,int nbdata1,int dimx,double *x0,double *y0,double *y1,double fy0y1[2*nbplot+1][2*nbplot+1],double *xmin,double *xmax);
void convert_date(int jj,int mm, int aa,char *date);
void readDataMarket(DataMarket *DM);
void read_params(int TypeModel,int *typen,int *nbetapes,int *nbdata,double *xinit);
void ReInitParams(int nbdata,DataMarket *DM,int TypeNorme,int TypeModel,int LogNorme,int MelangeNormes);
void BloqueVariables(int dimx,int TypeModel,int lequel,int *nbd,double *x,double *xmin,double *xmax);
void BloqueVariable_i(int dimx,int i,int *nbd,double *x,double *xmin,double *xmax);
void ReInitBornes(int dimx,int *nbd0,double *xmin0,double *xmax0,int *nbd,double *xmin,double *xmax);
