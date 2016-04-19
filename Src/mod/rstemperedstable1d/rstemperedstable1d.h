#ifndef _RSTEMPEREDSTABLE1D_H
#define _RSTEMPEREDSTABLE1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD RSTEMPEREDSTABLE1D

/*1D Regime Switching Temperedstable World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Transition_probabilities;
 
} TYPEMOD;

#endif
