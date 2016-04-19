#ifndef _SCHWARTZTROLLE_H
#define _SCHWARTZTROLLE_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD SCHWARTZTROLLE

/*SCHWARTZTROLLE World*/     
typedef struct TYPEMOD{ 
  VAR T;
  VAR P0T0;
  VAR f0T;
  VAR v0;
  VAR eta;
  VAR kappa;
  VAR sigma; 
  VAR rho;
  VAR alpha;
  VAR gammac;
} TYPEMOD;


#endif
