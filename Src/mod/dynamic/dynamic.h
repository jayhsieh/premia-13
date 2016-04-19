#ifndef _DYNAMIC
#define _DYNAMIC

#include "optype.h"
#include "var.h"
#include "error_msg.h"

#define TYPEMOD DYNAMIC

/* Copula model for CDO */     

typedef struct TYPEMOD{ 
  VAR Ncomp; /* number of companies */
  VAR r; /* interest rate */
} TYPEMOD;



#endif
