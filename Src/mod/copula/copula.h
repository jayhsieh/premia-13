#ifndef _COPULA_H
#define _COPULA_H

#include "optype.h"
#include "var.h"
#include "error_msg.h"
#include "enums.h"

#define TYPEMOD COPULA

/* Copula model for CDO */     

typedef struct TYPEMOD{ 
  VAR Ncomp; /* number of companies */
  VAR r; /* interest rate */
  VAR t_copula; /* type of copula */
  VAR t_intensity; /* type of intensity */
} TYPEMOD;



#endif
