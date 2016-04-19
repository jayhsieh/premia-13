
#ifndef LMM_HEADER_H
#define LMM_HEADER_H

#include<math.h>
#include"stdio.h"
#include"stdlib.h"
#include "pnl/pnl_vector.h"



typedef double (*funcVol)(double t, double T,double sigma) ; /*// volatility function*/

typedef struct
{
  int numberOfFactors;
  double flat_sigma;
  funcVol *vol; /*// vol is an array of volatility functions */
}Volatility;

typedef struct
{
  int numberOfMaturities ;
  double tenor;
  PnlVect *maturity ;
  PnlVect *libor ;
}Libor;

typedef struct
{
  int numberOfMaturities ;
  int numberOfTrajectories;
  double tenor;
  double *maturity ;
  double ***libor ;
}HistLibor;

typedef struct
{
  int numberOfFactors ;
  double *val ;

}RandomGenerator;


typedef struct
{
  double swaptionMaturity;
  double swapMaturity;
  double price;
  double strike;
  double tenor;
  int numberOfDates;

}Swaption;


typedef struct
{
  double maturity;
  double price;
  double strike;

}Caplet;


typedef struct
{
  int numberOfMaturities ;
  double tenor;
  double * value ;

}ZeroBond;






#endif

