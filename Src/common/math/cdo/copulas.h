#ifndef _COPULAS_H
#define _COPULAS_H

#include "cdo_math.h"
#include "structs.h"
#include "pnl/pnl_cdf.h"

typedef struct copula 
{
  char        *name;
  int         nfactor; /* number of factors */
  int         size;
  double      *points; /* array of size nfactor x size */
  double      *weights; /* array of size nfactor x size */
  
  double      (*density)(const struct copula *cop,const double x);
  double      *(*compute_cond_prob)(const struct copula   *cop,const double f_t);
  double     stu;/* uniquement pour la student */

  void        (*generate)(struct copula *cop);
  int         (*compute_default_time)(const struct copula *cop,
                                      const step_fun      *H,
                                      double          *time);
  void        *parameters;
} copula;

copula      *init_gaussian_copula(const double  rho);
void         free_copula (copula **cop);
copula      *init_clayton_copula(const double   theta);
copula      *init_student_copula(const double   rho,const double   t1);
copula      *init_nig_copula(const double       rho,
                             const double       alpha,
                             const double       beta);
copula      *init_double_t_copula(const double  rho,const double t1,const double t2);
copula      *init_copula (int t_copula, const double *p_copula);

#endif
