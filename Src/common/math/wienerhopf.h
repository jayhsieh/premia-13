#ifndef _WIENERHOPF_H
#define _WIENERHOPF_H

#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_fft.h"
#include "pnl/pnl_matvect.h"

/*/////////////////////////////////////////////////////*/

 int fastwienerhopf(int model, double mu, double qu, double om, int am, int upordown, int ifCall, double Spot, double lm1, double lp1,
			      double num,double nup, double cnum,double cnup,
			      double r, double divid,
			      double T, double h, double Strike1,
			      double bar,double rebate,
			      double er, long int step,
			      double *ptprice, double *ptdelta);

int fastwienerhopfamerican(int model, double mu, double qu, double om, 
int ifCall, double Spot, double lm1, double lp1,
    double num,double nup, double cnum,double cnup,
    double r, double divid,
    double T, double h, double Strike1,
    double er, long int step,
                           double *ptprice, double *ptdelta);
//Code for Backward Fourier
int bi_american(double mu, double qu, double om, 
int ifCall, double Spot, double lm1, double lp1,
    double num,double nup, double cnum,double cnup,
    double r, double divid,
    double T, double h, double Strike1,
    double er, long int step,
                double *ptprice, double *ptdelta);
  int bi_barr(double mu, double qu, double om, int upordown, int ifCall, double Spot, double lm1, double lp1,
			      double num,double nup, double cnum,double cnup,
			      double r, double divid,
			      double T, double h, double Strike1,
			      double bar,double rebate,
			      double er, long int step,
              double *ptprice, double *ptdelta);
  // Code for lookback floating strike
  int lookback_fls(int model, double mu, double qu, double om, int ifCall, double Spot, double minmax, double lm1, double lp1,
			      double num,double nup, double cnum,double cnup,
			      double r, double divid,
			      double T, double h, 
			      double er, double *ptprice, double *ptdelta);
    // Code for lookback fixed strike
  int lookback_fxs(int model, double mu, double qu, double om, int ifCall, double Spot, double minmax, double lm1, double lp1,
			      double num,double nup, double cnum,double cnup,
			      double r, double divid, double Strike,
			      double T, double h, 
			      double er, double *ptprice, double *ptdelta);
      // Code for swing
  int swing(int model, double mu, double qu, double om, 
int ifCall, double Spot, double lm1, double lp1,
    double num,double nup, double cnum,double cnup,
    double r, double divid,
    double T, double h, double Strike1, double del, int Nd,
    double er, long int step,
    double *ptprice, double *ptdelta);
  // Code for var
  int var_fft(int model, double mu,  
    double Spot, double lm1, double lp1,
    double num,double nup, double cnum,double cnup,
    double T, double h, double Strike1, double er, double alpha,
    double *ptprice, double *ptdelta);
#endif
