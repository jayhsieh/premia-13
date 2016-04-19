#ifndef _UVM1D_H
#define _UVM1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD UVM1D

/* UVM1D World */ 
typedef struct TYPEMOD {
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR sigmamin;
  VAR sigmamax;
} TYPEMOD;

#endif

