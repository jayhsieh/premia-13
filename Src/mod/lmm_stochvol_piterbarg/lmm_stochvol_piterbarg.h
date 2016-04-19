#ifndef _LMM_STOCHVOL_PITERBARG_H_
#define _LMM_STOCHVOL_PITERBARG_H_

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD LMM_STOCHVOL_PITERBARG

typedef struct TYPEMOD{
  VAR T;
  VAR Flag_InitialYieldCurve;
  VAR Var_SpeedMeanReversion;
  VAR Var_Volatility;
  VAR SkewsParams_a;
  VAR SkewsParams_b;
  VAR SkewsParams_c;
  VAR SkewsParams_d;
  VAR VolsParams_a;
  VAR VolsParams_b;
  VAR VolsParams_c;
  VAR VolsParams_d;
} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char* MOD(GetCurve)(TYPEMOD *pt);

#endif
