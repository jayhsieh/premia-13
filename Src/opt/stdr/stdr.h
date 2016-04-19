#ifndef  _STDr_H
#define _STDr_H

#include  "optype.h"
#include  "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h" 

#define TYPEOPT STDr


typedef struct TYPEOPT
{         
  VAR	   Maturity;     
  VAR	   Strike;
  VAR      Alpha;
  VAR      NumberofCreditors;
  VAR      ConstantProbabilityofDefault;
} TYPEOPT;

int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);      
int OPT(Check)(int user,Planning *pt_plan,Option *opt);    

#endif
