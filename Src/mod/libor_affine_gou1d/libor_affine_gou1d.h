#ifndef _LIBOR_AFFINE_GOU1D_H
#define _LIBOR_AFFINE_GOU1D_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD LIBOR_AFFINE_GOU1D

/*1D Libor Market Model World*/
typedef struct TYPEMOD
{
    VAR T;
    VAR flat_flag;
    VAR x0;
    VAR lambda;
    VAR alpha;
    VAR beta;

} TYPEMOD;

extern double MOD(GetYield)(TYPEMOD *pt);
extern char* MOD(GetCurve)(TYPEMOD *pt);
#endif
