#ifndef _OU1D_H
#define _OU1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD OU1D

/* OU1D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR R;
  VAR Speed;
  VAR Sigma;
  VAR a1;
  VAR a2;
  VAR a3;
} TYPEMOD;

#endif

