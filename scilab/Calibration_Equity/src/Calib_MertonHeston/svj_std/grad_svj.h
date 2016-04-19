#ifndef _GRADSVJ_H_
#define _GRADSVJ_H_
//
#include "optype.h"
#include "../utils/my_integral.h"
//
int calc_grad_svj(SVJPARAMS *svj,int dimx, double *grad);
//
// a enlever une fois integree dans premia
int Grad_FT_Call_Heston(double St0, NumFunc_1  *p, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,int dimx,double *grad);
int Grad_FT_Put_Heston(double St0, NumFunc_1  *p, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,int dimx,double *grad);
int Grad_FT_Call_Merton(double St0, NumFunc_1  *p, double T, double r, double divid, double V0,double lambda,double m0,double v,int dimx, double *grad);
int Grad_FT_Put_Merton(double St0, NumFunc_1  *p, double T, double r, double divid, double V0,double lambda,double m0,double v,int dimx, double *grad);
int Grad_FT_Call_HestMert(double St0, NumFunc_1  *p, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,double lambda, double m0,double v,int dimx,double *grad);
int Grad_FT_Put_HestMert(double St0, NumFunc_1  *p, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,double lambda, double m0,double v,int dimx,double *grad);
//
#endif
