#ifndef __COMPANY_H
#define __COMPANY_H

#include        "cdo_math.h"
#include        "structs.h"

typedef struct company
{
  double      nominal;    
  step_fun        *h;     
  step_fun        *H;     
  double      mean_delta;
  pfun_void_R     *generate_delta;
  pfun_R_complex  *phi_recov; 
  void        *p_recovery;    
} company;

typedef struct params_cov_unif 
{
  double  a;
  double  b;
} params_recov_unif;

typedef struct params_cov_gauss
{
  double  m;
  double  s;
} params_recov_gauss;


company         *init_company(const double      nominal,
                              const int         intensity_n_step,
                              const double      *intensity_x,
                              const double      *intensity_h_x,
                              const double      mean_delta,
                              pfun_void_R       *generate_delta,
                              pfun_R_complex        *phi_recov,
                              void          *p_recovery);
company         *init_company_cov_cst(const double      nominal,
                                      const int         intensity_n_step,
                                      const double      *intensity_x,
                                      const double      *intensity_h_x,
                                      const double      delta);
company         *init_company_cov_unif(const double     nominal,
                                       const int        intensity_n_step,
                                       const double     *intensity_x,
                                       const double     *intensity_h_x,
                                       const double     a,
                                       const double     b);
company         *init_company_cov_gauss(const double        nominal,
                                        const int        intensity_n_step,
                                        const double     *intensity_x,
                                        const double     *intensity_h_x,
                                        const double     m,
                                        const double     s);

company         *homogenize_company(const company   *comp,
                                    const double    new_nominal,
                                    const double    new_delta);

void            free_company(company            *co);

#endif  /* __COMPANY_H */
