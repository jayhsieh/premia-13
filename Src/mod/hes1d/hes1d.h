#ifndef _HES1D_H
#define _HES1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD HES1D

/*1D HESTON World*/    
 
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
