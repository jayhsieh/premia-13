#ifndef _CEV1D_H
#define _CEV1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD CEV1D

/* CEV1D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR R;
  VAR Divid;
  VAR v;
  VAR beta; 
} TYPEMOD;

#endif

