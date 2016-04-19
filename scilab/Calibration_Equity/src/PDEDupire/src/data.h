#include <stdlib.h>
#ifndef DATA
#define DATA 


/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

/**************type definitions**********************/




/*marketData = price of the option for the couple (strike,maturity)*/
struct marketData{
  double strike;
  double maturity;
  double price;
  struct marketData *next;
};



#endif


