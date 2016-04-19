#ifndef  _STDg_H
#define _STDg_H

#include  "optype.h"
#include  "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h"

#define TYPEOPT STDg

typedef struct TYPEOPT
{         
  VAR      PayOff;
  VAR	   Maturity; 
   VAR	   FutureMaturity; 
  VAR	   EuOrAm;   
  VAR      NbExerciseDate;
  VAR      RefractingPeriod;
  
} TYPEOPT;

int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);      
int OPT(Check)(int user,Planning *pt_plan,Option *opt);    

#endif
