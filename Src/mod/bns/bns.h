#ifndef _BNS_H
#define _BNS_H

#include "optype.h"
#include "var.h"

#define TYPEMOD BNS


typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR Sigma0;
  VAR Lambda;
  VAR Rho;
  VAR Beta;
  VAR Alpha;
} TYPEMOD;

#endif

