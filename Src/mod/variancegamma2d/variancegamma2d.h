#ifndef _VARIANCEGAMMA2D_H
#define _VARIANCEGAMMA2D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD VARIANCEGAMMA2D

/* VARIANCEGAMMA2D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S01;
  VAR S02;
  VAR R; 
  VAR ap;
  VAR am;
  VAR lambda;
  VAR alpha;
} TYPEMOD;

#endif

