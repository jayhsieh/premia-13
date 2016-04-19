#ifndef _SCOTT1D_H
#define _SCOTT1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD SCOTT1D

/*1D SCOTT World*/     
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
