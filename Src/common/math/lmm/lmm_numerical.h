
#ifndef LMM_NUMERICAL_H 
#define LMM_NUMERICAL_H

#include "pnl/pnl_vector.h"
#include "pnl/pnl_cdf.h"
#include"lmm_header.h"
#include"lmm_volatility.h"
#include"lmm_libor.h"
#include"lmm_zero_bond.h"

double maxi(double a,double b);
double ppos(double x);
int Set_to_Zero(double *x,int dim);
void computeNumeraire(char*nameOfthemeasure,Libor* ptLib,Swaption* ptSwpt,double *Numeraire,int j,int k,double auxspot);
void Name_To_Measure(char *ErrorMessage, char *name,
                     int (**computeEvol)(const PnlVect* ptRand,Libor* ptLibOld, Libor* ptLibNew,
                                         Volatility* ptVol,double dt,double t, double sigma_cost));


double ps_lmm(int n, double *u, double *v); /*scalar product of 2 vectors*/
double BSFormula(Swaption *ptSwpt, Libor* ptLib, double evalTime, double blackVol);



#endif

