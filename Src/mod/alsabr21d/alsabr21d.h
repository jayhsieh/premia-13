#ifndef _ALSABR21D_H
#define _ALSABR21D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD ALSABR21D

/* ALSABR21D World */
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR z0;
  VAR Divid;
  VAR R;
  VAR mu;
  VAR eta;
  VAR a1;
  VAR a2;
  VAR c1;
  VAR c2;
} TYPEMOD;

#endif

