#include "intg.h"

#include "pnl/pnl_integration.h"

typedef double (realf)(double);
static double intg_f (double x, void *p)
{
  realf *f = (realf *) p;
  return (*f)(x);
}

/*
 * this function is only a wrapper for pnl_integration_qag for function of
 * a single variate
 */
void intg(double a, double b, double (*f)(double), double ea, double er, double *val, double *abserr)
{
  int neval;
  PnlFunc F;
  F.function = intg_f;
  F.params = (void *) f;
  pnl_integration_qag (&F, a, b, ea, 0, er, val, abserr, &neval);
}


