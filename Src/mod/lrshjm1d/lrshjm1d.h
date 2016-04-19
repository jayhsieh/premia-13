#ifndef _LiRitchkenSankarasubramanian1D_H
#define _LiRitchkenSankarasubramanian1D_H

#include "optype.h"
#include "var.h"
#include "enums.h"

#define TYPEMOD LRSHJM1D

/*1D Li Ritchken Sankarasubramanian World*/
typedef struct TYPEMOD{
  VAR T;
  VAR flat_flag;
  VAR Sigma;
  VAR Kappa;
  VAR Rho;
  VAR Lambda;
} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char * MOD(GetCurve)(TYPEMOD *pt);

#endif
