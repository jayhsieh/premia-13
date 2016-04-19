#ifndef _INFLATION_LMM_HESTON1D_H
#define _INFLATION_LMM_HESTON1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD INFLATION_LMM_HESTON1D

/*1D INFLATION Libor Market Model Stochastic Volatility World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR I0;
  VAR SigmaI;
  VAR F0;
  VAR SigmaF;
  VAR Sigma0;
  VAR SpeedMeanReversion;
  VAR LongRunVariance;
  VAR Sigma2;
  VAR RhoFI;
  VAR RhoFV;
  VAR RhoIV;
  VAR RhoI;
} TYPEMOD;

#endif


 
