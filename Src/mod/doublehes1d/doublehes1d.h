#ifndef _DOUBLEHES1D_H
#define _DOUBLEHES1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD DOUBLEHES1D

/* DOUBLEHES1D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR Sigma0;
  VAR MeanReversion;
  VAR LongRunVariance;
  VAR Sigma;
  VAR Sigma0V;
  VAR MeanReversionV;
  VAR LongRunVarianceV;
  VAR SigmaV;
   VAR Rho;
  VAR RhoSV2;
  //VAR RhoVV; 
} TYPEMOD;

#endif

