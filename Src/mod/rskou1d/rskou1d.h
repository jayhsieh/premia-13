#ifndef _RSKOU1D_H
#define _RSKOU1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD RSKOU1D

/*1D Regime Switching Kou World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Transition_probabilities;
 
} TYPEMOD;

#endif
