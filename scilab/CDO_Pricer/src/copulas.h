#include "maths.h"
#include "structs.h"


typedef struct copula 
{
  char        *name;
  int         size;
  double      *points;
  double      *weights;
  
  double      (*density)(const struct copula *cop,const double x);

  double      *(*compute_cond_prob)(const struct copula   *cop,const double f_t);
  double     stu;/* uniquement pour la student */

 
  void        (*generate)(struct copula           *cop);
   
  int         (*compute_default_time)(const struct copula *cop,
                                      const step_fun      *H,
                                      double          *time);
  void        *parameters;
} copula;

copula      *init_gaussian_copula(const double  rho);
copula      *init_clayton_copula(const double   theta);
copula      *init_student_copula(const double   rho,const double   t1);
copula      *init_nig_copula(const double       rho,
                             const double       alpha,
                             const double       beta);
copula	    *init_double_t_copula(const double	rho,const double t1,const double t2);
