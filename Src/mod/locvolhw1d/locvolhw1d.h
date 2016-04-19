#ifndef _LOCVOLHW1D_H
#define _LOCVOLHW1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD LOCVOLHW1D

/* LOCVOLHW1D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR csi;
  VAR kappa; 
  VAR v;
  VAR beta; 
  VAR rho;
  VAR f0t; 

} TYPEMOD;

#endif

