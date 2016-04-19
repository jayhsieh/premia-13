/***************************/
/* CEV scilab interface    */
/***************************/


#include "lmm_cev_pricer.h"

int lmm_cev_sci(double * price,size_t* type_model,size_t* type_product,size_t* type_pricing, 
		 size_t* type_scheme, double* alpha, double* fzero, double* H,
	       double* expiry,double* swp_lenght,int* int_pts)
{
  cev_price(*type_model,*type_product, *type_pricing,* type_scheme, * alpha, *fzero, * H,
	    *expiry,*swp_lenght,* int_pts, price);
  
  return(EXIT_SUCCESS);
}

