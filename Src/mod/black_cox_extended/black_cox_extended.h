#ifndef _BLACK_COX_EXTENDED_H
#define _BLACK_COX_EXTENDED_H

#include "optype.h"
#include "var.h"

#define TYPEMOD BLACK_COX_EXTENDED

/* BLACK_COX_EXTENDED World */ 
typedef struct TYPEMOD {
  VAR S0;
  VAR R; 
  VAR Sigma;
  VAR L;
  VAR alpha;
  VAR mu;
} TYPEMOD;

#endif

