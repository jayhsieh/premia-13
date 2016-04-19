#ifndef _NIG1FACT1D_H
#define _NIG1FACT1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD NIG1FACT1D

/* NIG1FACT1D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR R;
  VAR alpha;
  VAR beta;
  VAR delta;
  VAR mu;
  VAR Sigma;
  VAR lambda;
} TYPEMOD;

#endif

