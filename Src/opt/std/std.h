#ifndef  _STD_H
#define _STD_H

#include  "optype.h"
#include  "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h"

#define TYPEOPT STD

typedef struct TYPEOPT
{         
  VAR      PayOff;
  VAR	   Maturity;      
  VAR	   EuOrAm;   
} TYPEOPT;

int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);      
int OPT(Check)(int user,Planning *pt_plan,Option *opt);    

#endif
