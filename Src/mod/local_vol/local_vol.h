#ifndef _LOCAL_VOL_H
#define _LOCAL_VOL_H

#include "optype.h"
#include "var.h"

#define TYPEMOD LOCAL_VOL

/* LOCAL_VOL World */ 
typedef struct TYPEMOD {
  VAR S0;
  VAR Interest;
  VAR Divid;
  VAR Eta;
  VAR Sigma;
} TYPEMOD;

#endif

