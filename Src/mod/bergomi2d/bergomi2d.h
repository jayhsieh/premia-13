#ifndef _BERGOMI2D_H
#define _BERGOMI2D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD BERGOMI2D

/* BERGOMI2D World */ 
typedef struct TYPEMOD {
  VAR T;
  VAR S0;
  VAR Divid;
  VAR R; 
  VAR csi0;
  VAR omega;
  VAR theta;
  VAR k1;
  VAR k2; 
  //VAR rhoxy;
  VAR rhoSx;
  VAR rhoSy;
  
} TYPEMOD;

#endif

