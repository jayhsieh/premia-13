#ifndef  _DOUBLIM_H
#define _DOUBLIM_H

#include "optype.h"
#include "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h"

#define TYPEOPT DOUBLIM

/*Limit Option// General double barrier*/

typedef struct TYPEOPT{         
  /* var_setable */
  VAR    PayOff;
  VAR    Rebate;		
  VAR    LowerLimit;
  VAR    UpperLimit;
  VAR    Maturity;                                  
  /* var_fixed */
  VAR	 OutOrIn;
  VAR    Parisian;
  VAR    RebOrNo;
  VAR    EuOrAm;	         

} TYPEOPT;

int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);      
int OPT(Check)(int user,Planning *pt_plan,Option *opt);     

#endif
