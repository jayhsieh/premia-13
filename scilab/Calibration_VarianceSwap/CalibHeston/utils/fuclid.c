#include "optim.h"

int optim_fuclid (int *n, double *x, double *y, double *ps, opt_simul_data *optim_data)
{
  int i;
  *ps = 0.;
  for (i = 0 ; i < *n ; ++i) *ps += x[i] * y[i];
  return 0;
}	

