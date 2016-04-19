#ifndef _MRC30D_H
#define _MRC30D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD MRC30D

/* MRC30D World */ 
typedef struct TYPEMOD {
  VAR Size;
  VAR T;
  VAR R; 
  VAR kappa;
  VAR eta;
  VAR gama;
  VAR a;
  VAR InitialStocksWeights;
  VAR LocalVolatilities;
  VAR Basket_Correlation;
  VAR BasketLocalVolatility;
} TYPEMOD;

#endif

