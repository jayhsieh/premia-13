#ifndef _HESCIR1D_H
#define _HESCIR1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD HESCIR1D

/* HESCIR1D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR r0;
  VAR kr;
  VAR thetar;
  VAR Sigmar;
  VAR V0;
  VAR kV;
  VAR thetaV;
  VAR SigmaV;
  VAR RhoSr;
  VAR RhoSV;
  VAR RhorV;
  
} TYPEMOD;

#endif

