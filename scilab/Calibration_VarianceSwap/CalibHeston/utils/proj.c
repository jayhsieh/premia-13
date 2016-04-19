#include "optim.h"

int optim_proj (int *n, double *binf, double *bsup, double *x)
{
  int i;
  for (i = 0; i < *n; ++i)
    {
      x[i] = Max (binf[i],Min (x[i],bsup[i]));
    }
  return 0;
}	

