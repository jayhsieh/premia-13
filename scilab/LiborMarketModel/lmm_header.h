#ifndef LMM_HEADER_H 
#define LMM_HEADER_H

#include<math.h>
#include"stdio.h"
#include"stdlib.h"



typedef double (*funcVol)(double t, double T) ; /*// volatility function*/

typedef struct 
{
  int numberOfFactors ;
  funcVol *vol; /*// vol is an array of volatility functions */
}Volatility;

typedef struct 
{
  int numberOfMaturities ; 
  double tenor;
  double *maturity ;
  double *libor ;
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

