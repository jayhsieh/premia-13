#ifndef  _STDNDC_H
#define _STDNDC_H

#include  "optype.h"
#include  "var.h"
#include "enums.h"
#include "option.h"

#include  "chk.h"
#include  "numfunc.h"

#define TYPEOPT STDNDc


typedef struct TYPEOPT
{         
  VAR Ncomp; /* number of companies inherited from the model */
  VAR maturity; /* maturity time */
  VAR t_nominal; /* Homogeneous nominals ? */
  VAR tranch; /* vector of the tranches */
  VAR t_recovery; /* Type of recovery in case of defaut */
  VAR p_recovery; /* Parameters of the recovery in case of defaut */
  VAR NbPayment; /* Number of coupon payments per year */
  VAR date; /* current date */
  VAR n_defaults; /* number of defaults at current date */
} TYPEOPT;

int n_param_recovery (int *n, double *t, int recovery_value, int with_init);

#endif
