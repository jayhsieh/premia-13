#ifndef __tool_box__
#define __tool_box__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"

  
typedef struct X_Grid
{
  double Bnd_Left;
  double Bnd_Right;
  double step;
  int N;
}X_Grid;

extern X_Grid * x_grid_create(double dBnd_Left,double dBnd_Right,int dN);
extern double x_grid_point(X_Grid *grid ,int i);


  
// Math tools  
  
extern void progonka(int N, const double low,  const double diag, const double up, 
                     const PnlVect * right_side, PnlVect * result);
extern void quadratic_interpolation(double Fm1, double F0,double Fp1,double Xm1,double X0,double Xp1,double X,double * FX,double * dFX);

// All these functions ara not yet tested 
extern void polint(double* xa,double* ya, int n, double x, double *y, double *dy);
extern double trapzd( PnlFunc *func,void * params, double a, double b, int n,double old_s);
extern double qromb(PnlFunc *func,void * params, double a, double b);


#endif 

