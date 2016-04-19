#ifndef  _IMPLIED_BS_H
#define _IMPLIED_BS_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "pnl/pnl_matrix_int.h"
#include "pnl/pnl_cdf.h"

/*  Compute delta forward because this quantity  is sold/bought to hedge option with forward*/
extern double pnl_forward_price(double Spot,double r, double divid, double Maturity);
extern double pnl_bs_impli_call (double Vol,double Bond, double Forward, double Strike, double Maturity);
extern double pnl_bs_impli_put  (double Vol,double Bond, double Forward, double Strike, double Maturity);
extern double pnl_bs_impli_call_delta_forward (double Vol,double Bond, double Forward,
                                         double Strike, double Maturity); 
extern double pnl_bs_impli_put_delta_forward  (double Vol,double Bond, double Forward,
                                         double Strike, double Maturity); 
extern double pnl_bs_impli_call_put (int Is_Call, double Vol,double Bond, double Forward,
                               double Strike, double Maturity);
extern double pnl_bs_impli_call_put_delta_forward (int Is_Call, double Vol,double Bond,
                                             double Forward, double Strike, double Maturity);
extern double pnl_bs_impli_vega(double Vol,double Bond, double Forward, double Strike, double Maturity);
extern double pnl_bs_impli_gamma(double Vol,double Bond, double Forward, double Strike, double Maturity);
extern double pnl_bs_impli_s_square_gamma (double Vol,double Bond, double Forward, double Strike,double Maturity);
extern double pnl_bs_impli_implicit_vol(int Is_Call, double Price,double Bond,
                                  double Forward, double Strike, double Maturity);  
extern int pnl_bs_impli_matrix_implicit_vol(const PnlMatInt * Is_Call, const PnlMat * Price,
                                      double spot,double rate, double divid,
                                      const PnlVect * Strike,const PnlVect * Maturity,PnlMat * Vol);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
