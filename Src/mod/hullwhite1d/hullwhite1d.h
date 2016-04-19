#ifndef _HullWhite1D_H
#define _HullWhite1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD HullWhite1D

/*1D HULL-WHITE World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR flat_flag;
  VAR a;
  VAR Sigma;
} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char* MOD(GetCurve)(TYPEMOD *pt);

#endif
