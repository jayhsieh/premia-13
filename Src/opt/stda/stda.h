#ifndef  _STDa_H
#define _STDa_H

#include  "optype.h"
#include  "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h"

#define TYPEOPT STDa

typedef struct TYPEOPT
{
  VAR      PayOff;
  VAR	   EuOrAm; 
  VAR	   Maturity;      
  VAR      DeemedContribution;
  VAR      InitialAge;
  VAR      Premium;
  VAR      MinimumGuaranteedInterestRate;
   VAR     NumberofMonitoringDates;
  VAR Alpha;
} TYPEOPT;

int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);      
int OPT(Check)(int user,Planning *pt_plan,Option *opt);    

#endif
