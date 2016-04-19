
#ifndef __PRICE_CDO_H
#define __PRICE_CDO_H

#include "copula_stdndc.h"

extern double* get_t_recovery_arg (const VAR *x);
extern int      price_cdo(const int    *n_comp, 
                          const double *nominal,
                          const int     n_dates,
                          const double *dates,
                          const int     n_tranches,
                          const double *tr,
                          const double *intensity,
                          const int     n_rates,
                          const double *x_rates,
                          const double *y_rates,
                          const int    *t_recovery,
                          const double *recovery,
                          const int    *t_copula,
                          const double *p_copula,
                          const int    *t_method,
                          const int    *p_method,
                          double       *price,
                          double       *def_leg,
                          double       *pay_leg);


extern int
premia_interf_price_cdo ( TYPEOPT *ptOpt, TYPEMOD *ptMod, PricingMethod *Met,
                          PnlVect **nominal, PnlVect **intensity,
                          int *n_rates, PnlVect **x_rates, PnlVect **y_rates,
                          int *n_dates, PnlVect **dates, int *n_tranches, 
                          int **p_method, int *is_homo);

#endif /* __PRICE_CDO_H */

