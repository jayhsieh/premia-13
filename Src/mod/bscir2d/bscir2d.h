#ifndef _BSCIR2D_H
#define _BSCIR2D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD BSCIR2D

/* BSCIR2D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0; 
  VAR Sigma;
  VAR r0;
  VAR k;
  VAR SigmaR;
  VAR theta;
  VAR Rho;
  VAR Mortality;
} TYPEMOD;

#endif

