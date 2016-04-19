#ifndef _CirPlus1D_H
#define _CirPlus1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD CirPP1D


/*1D Cir++ World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR flat_flag;
  VAR a;
  VAR b;
  VAR Sigma;
  VAR InitialYields;
} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char* MOD(GetCurve)(TYPEMOD *pt);

#endif
