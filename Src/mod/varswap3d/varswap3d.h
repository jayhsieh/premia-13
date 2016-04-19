#ifndef _VARSWAP3D_H
#define _VARSWAP3D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD VARSWAP3D

typedef struct TYPEMOD{ 
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R;
  VAR V0;
  VAR Beta;
  VAR MeanReversion;
  VAR Rho;
} TYPEMOD;



#endif
