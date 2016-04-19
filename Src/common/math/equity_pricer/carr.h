#ifndef __carr__
#define __carr__

#include "tool_box.h"
#include "finance_tool_box.h"
#include "levy_process.h"
#include "levy_diffusion.h"
#include "pnl/pnl_complex.h"

extern int CarrMethod_VectStrike(PnlVect *K,
                                 PnlVect * Price,
                                 double S0,
                                 double T,
                                 double Kstep,
                                 double CallPut,
                                 double r,
                                 double divid,
                                 double sigma,
                                 void * Model,
                                 dcomplex (*ln_phi)(dcomplex u,double t,void * model)
                                 );

extern int CarrMethod_onStrikeList(PnlVect *K,
                                   PnlVect * Price,
                                   double S0,
                                   double T,
                                   double CallPut,
                                   double r,
                                   double divid,
                                   double sigma,
                                   Levy_diffusion * Model);

extern int CarrMethod(double S0,
                      double T,
                      double K,
                      double CallPut,
                      double r,
                      double divid,
                      double sigma,
                      //Levy_diffusion * Model,
                      void * Model,
                      dcomplex (*ln_phi)(dcomplex u,double t,void * model),
                      double *ptprice,
                      double *ptdelta);

extern int CarrMethod_Vanilla_option(Option_Eqd *opt,
                                     double sigma,
                                     Levy_process * Model);

extern int CarrMethod_Vanilla_option_LD(Option_Eqd *opt,
                                        double sigma,
                                        Levy_diffusion * Model);

extern int AttariMethod_Vanilla_option(Option_Eqd *opt,
                                       double sigma,
                                       Levy_process * Model);

extern int AttariMethod_Vanilla_option_LD(Option_Eqd *opt,
                                          double sigma,
                                          Levy_diffusion * Model);

extern int Var_Swap_Price_option(Option_Eqd *opt,Levy_process * Model);
extern int Var_Swap_Price_option_LD(Option_Eqd *opt,Levy_diffusion * Model);

#endif 


