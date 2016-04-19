#ifndef _STEIN1D_H
#define _STEIN1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD STEIN1D

/*1D STEIN World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR Sigma0;
  VAR MeanReversion;
  VAR LongRunVariance;
  VAR Sigma;
  VAR Rho;  
} TYPEMOD;

#endif
