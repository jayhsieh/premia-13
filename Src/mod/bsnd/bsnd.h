#ifndef _BSND_H
#define _BSND_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD BSND

typedef struct TYPEMOD{ 
  VAR Size;
  VAR T;
  VAR S0;
  VAR Sigma;
  VAR Divid;
  VAR Rho;
  VAR R;       
} TYPEMOD;


#endif
