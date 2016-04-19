#include "costFunction.h"
#include <stdlib.h>
#include "utils/intg.h"
#include "pnl/pnl_finance.h"

// ===================================================================
// Derivee de la fonction de repartition de la loi Gaussienne :
// e^{-x^2/2}/sqrt(2.*PI) tout simplement !
// ===================================================================
static inline  double dN(double x)
{
   return (1./sqrt(2.*M_PI)) * exp(-x*x/2.);
}

// ===================================================================
// d(formule de BS pour un Call Euro) / d Sigma
// ===================================================================
int dCall_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *dprice){
  double sigmasqrt,d1,d2;
  double sqrtt,A,dd1,dd2;
  //
  sqrtt     = sqrt(t);
  sigmasqrt = sigma*sqrt(t);
  A         = (log(s/k)+(r-divid)*t)/sqrtt;
  //
  d1  = A /sigma + sigmasqrt/2.;
  d2  = d1-sigmasqrt;
  dd1 = - A /(sigma*sigma) + 0.5*sqrtt;
  dd2 = dd1 - sqrtt;
  //
  *dprice = s*exp(-divid*t)*dN(d1)*dd1 -exp(-r*t)*k*dN(d2)*dd2;
  //
  return 0;
}
// ===================================================================
// d(formule de BS pour un Put Euro) / d Sigma
// ===================================================================
int dPut_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *dprice)
{
  double sigmasqrt,d1,d2;
  double sqrtt,A,dd1,dd2;
  //
  sqrtt     = sqrt(t);
  sigmasqrt = sigma*sqrt(t);
  A         = (log(s/k)+(r-divid)*t)/sqrtt;
  //
  d1  = A /sigma + sigmasqrt/2.;
  d2  = d1-sigmasqrt;
  dd1 = - A /(sigma*sigma) + 0.5*sqrtt;
  dd2 = dd1 - sqrtt;
  //
  *dprice =-exp(-r*t)*k*dN(-d2)*dd2 + s*exp(-divid*t)*dN(-d1)*dd1;
  //
  return 0;
}
//
// ===================================================================
// Fonction pour calculer le Sigma Implicite dans la Call ou le Put
//  Europeen avec BS et de prix cbs. Methode utlisee : dichotomie.
// ===================================================================
double SigmaImplicite(int TypeOpt,double St,double T,double K,
			double r,double divid,double cbs,int *nbiters)
{

  int error;
  return pnl_bs_implicit_vol (TypeOpt, cbs, St, K, T, r, divid, &error);
}
