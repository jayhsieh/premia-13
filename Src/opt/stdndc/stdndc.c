#include "stdndc.h"
#include "error_msg.h"

/**
 * determine the number of parameters for the given type of recovery
 *
 * @param n number of parameter (set on output)
 * @param t array containing on output the default parameters
 * @param recovery_value an integer describing the type of recovery
 * @param with_init if set to 1 t is initializes.
 * @return OK or WRONG
 */
int n_param_recovery (int *n, double *t, int recovery_value, int with_init)
{
  if (with_init && t==NULL) return WRONG;
  switch (recovery_value)
    {
      /* Constant */
    case 1: *n = 1;
      if (with_init) { t[0] = 0.4; }
      break;
      /* Uniform */
    case 2: *n = 2;
      if (with_init) { t[0] = 0.3; t[1]=0.5; }
      break;
    default: *n = 0; return WRONG; break;
    }
  return OK;
}

extern Option OPT(CDO);
extern Option OPT(CDO_HEDGING);
extern Option OPT(CDO_COPULA);
extern Option OPT(CDO_HAWKES_INTENSITY);

Option* OPT(family)[]={
  &OPT(CDO),
  &OPT(CDO_HEDGING),
  &OPT(CDO_COPULA),
  &OPT(CDO_HAWKES_INTENSITY),
  NULL
};
