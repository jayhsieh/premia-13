#ifndef _HAWKES_INTENSITY_H
#define _HAWKES_INTENSITY_H

#include "optype.h"
#include "var.h"

#define TYPEMOD HAWKES_INTENSITY

/* HAWKES_INTENSITY World */ 
typedef struct TYPEMOD {
  VAR Ncomp; /* must be the first variable */
  VAR lambda0;
  VAR kappa;
  VAR c;
  VAR delta;
  VAR r; 
} TYPEMOD;

#endif

