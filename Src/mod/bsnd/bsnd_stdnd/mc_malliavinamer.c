
/***********************************************************************/
/* Author Lokman A. Abbas-Turki   <lokman.abbas-turki@laposte.net>      */
/*                                                                      */
/* Pricing American options using Malliavin calculus and non-parametric */
/* variance reduction methods based on conditionning and a judicious  	*/
/* choice of the number of paths used for the approximation of the		*/
/* numerator and the denomenator that intervene in the conditional		*/ 
/* expectation (the continuation).										*/ 
/*																		*/	
/* This method does not use any other variance reduction method except  */
/* the ones described in:												*/
/*																		*/	
/* American Options Based on Malliavin Calculus and Nonparametric		*/
/* Variance Reduction Methods, Lokman Abbas-Turki and Bernard Lapeyre,	*/
/* preprint on arXiv.org.												*/
/*																		*/
/* In this version of the program, the author provides only the pricing */
/* method for the multidimensional non-correlated Black & Scholes model.*/
/* As mentioned in the paper above, we do not use controle variate for  */
/* the variance reduction.												*/
/************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "bsnd_stdnd.h"
#include "black.h"
#include "optype.h"
#include "enums.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_vector.h"

static int CHK_OPT(MC_MalliavinAmer)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(MC_MalliavinAmer)(void *Opt, void *Mod, PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
      Met->HelpFilenameHint = "mc_longstaffschwatrz_nd";
      Met->Par[0].Val.V_LONG=1000;
      Met->Par[1].Val.V_ENUM.value=0;
      Met->Par[1].Val.V_ENUM.members=&PremiaEnumMCRNGs;
      Met->Par[2].Val.V_INT=10;

    }
  return OK;
}

PricingMethod MET(MC_MalliavinAmer)=
{
  "MC_MalliavinAmer",
  {{"N iterations",LONG,{100},ALLOW},
      {"RandomGenerator",ENUM,{0},ALLOW},
      {"Number of Exercise Dates",INT,{100},ALLOW},
      {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(MC_MalliavinAmer),
  {{"Price",DOUBLE,{100},FORBID},
      {"Error",DOUBLE,{100},FORBID},
      {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(MC_MalliavinAmer),
  CHK_mc,
  MET(Init)
};



