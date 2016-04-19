#include "copulas.h"


copula *init_copula (int t_copula, const double *p_copula)
{
  copula *cop;
  switch (t_copula)
    {
    case 1 :
      cop = init_gaussian_copula (p_copula[0]);
      break;
    case 2 :
      cop = init_clayton_copula (p_copula[0]);
      break;
    case 3 :
      cop = init_nig_copula (p_copula[0], p_copula[1], p_copula[2]);
      break;
    case 4:
      cop = init_student_copula ( p_copula[0], p_copula[1]);
      break;
    case 5:
      cop =  init_double_t_copula ( p_copula[0], p_copula[1], p_copula[2]);
      break;
    default:
      return NULL;
    }
  return cop;
}


void free_copula (copula **cop)
{
  free ((*cop)->points);
  free ((*cop)->weights);
  free ((*cop)->parameters);
  free (*cop); *cop = NULL;
}
