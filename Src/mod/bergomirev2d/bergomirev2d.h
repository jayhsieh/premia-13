#ifndef _BERGOMIREV2D_H
#define _BERGOMIREV2D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD BERGOMIREV2D

/* BERGOMIREV2D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R;
  VAR theta;
  VAR k1;
  VAR k2; 
  //VAR rhoxy;
  VAR rhoSx;
  VAR rhoSy;
  VAR ForwardVarianceData;
  
} TYPEMOD;

#endif

