#ifndef _LMM_HESTON1D_H
#define _LMM_HESTON1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD LMM_HESTON1D

/*1D Libor Market Model Stochastic Volatility World: with 1  or 2 factor*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR NbFactors;
  VAR l0;
  VAR Sigma;
  VAR Sigma0;
  VAR MeanReversion;
  VAR LongRunVariance;
  VAR Sigma2;
  VAR Rho1;  
  VAR Rho2; 
} TYPEMOD;

#endif
