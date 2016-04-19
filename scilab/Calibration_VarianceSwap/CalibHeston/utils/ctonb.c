#include "optim.h"

int optim_ctonb (int *n, double *u, double *v, opt_simul_data *optim_data)
{
  int i;
  for (i = 0; i < *n; ++i)
    {
      v[i] = u[i];
    }
  return 0;
}
