#ifndef _SCHWARTZ_H
#define _SCHWARTZ_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD SCHWARTZ

typedef struct TYPEMOD{ 
  VAR T;
  VAR R;
  VAR Divid;
  VAR sigmad;
  VAR sigmas;
  VAR alpha;
  VAR Rho;
} TYPEMOD;


#endif
