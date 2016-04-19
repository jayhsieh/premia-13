#ifndef _HullWhite1DGeneralized_H
#define _HullWhite1DGeneralized_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD HullWhite1DGeneralized

/*1D HULL-WHITE World*/
typedef struct TYPEMOD{
  VAR T;
  VAR flat_flag;
  VAR CapletCurve;
  VAR a;
} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char* MOD(GetCurve)(TYPEMOD *pt);

#endif
