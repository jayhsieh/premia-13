#ifndef _DPS_H
#define _DPS_H

#include "optype.h"
#include "var.h"

#define TYPEMOD DPS

/* DPS World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 

  VAR Rho;
  VAR Sigma0;
  VAR Kappa;
  VAR Eta;
  VAR Theta;

  VAR LambdaS;
  VAR MeanS;
  VAR SigmaS;
  VAR LambdaV;
  VAR MeanV;
  VAR LambdaSV;
  VAR MeanSV;
  VAR SigmaSV;
  VAR MeanVS;
  VAR RhoSV;
  
 
} TYPEMOD;

#endif

