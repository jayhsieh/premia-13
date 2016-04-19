

#include"cumulfunc.h"

double Chi2N(double x, double df, double ncparam)
{
  double P,Q,r;
  int status;
  int i=1;

  cdfchn(&i,&P,&Q,&x,&df,&ncparam,&status,&r);

  return P;
}

