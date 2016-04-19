#ifndef _HullWhite2D_H
#define _HullWhite2D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD HullWhite2D

/*2D HullWhite World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR flat_flag;
  VAR InitialYieldsu;
  VAR aR;
  VAR SigmaR;
  VAR bu;
  VAR Sigmau;
  VAR Rho;
} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char* MOD(GetCurve)(TYPEMOD *pt);

#endif
