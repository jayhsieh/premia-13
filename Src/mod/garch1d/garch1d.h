#ifndef _GARCH1D_H
#define _GARCH1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD GARCH1D

/* GARCH1D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR R; 
  VAR alpha0;
  VAR alpha1;
  VAR lambda;
  VAR beta1;
} TYPEMOD;

#endif

