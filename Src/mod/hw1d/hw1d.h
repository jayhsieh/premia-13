#ifndef _HW1D_H
#define _HW1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD HW1D

/*1D HULL-WHITE World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR Sigma0;
  VAR Mean;
  VAR Sigma;
  VAR Rho;  
} TYPEMOD;

#endif
