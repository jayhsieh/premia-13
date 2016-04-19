#ifndef _BlackKarasinski1D_H
#define _BlackKarasinski1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD BlackKarasinski1D



/*1D  BlackKarasinski World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR flat_flag;
  VAR r0;
  VAR a;
  VAR Sigma;
} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char* MOD(GetCurve)(TYPEMOD *pt);

#endif
