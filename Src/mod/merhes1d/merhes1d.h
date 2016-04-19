#ifndef _MERHES1D_H
#define _MERHES1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD MERHES1D

/*1D MERTON-HESTON World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR Sigma0;
  VAR MeanReversion;
  VAR LongRunVariance;
  VAR Sigma;
  VAR Lambda;  
  VAR Mean;
  VAR Variance;
  VAR Rho;  
} TYPEMOD;


#endif
