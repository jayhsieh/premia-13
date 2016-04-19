#ifndef  _STDf_H
#define _STDf_H

#include  "optype.h"
#include  "var.h"
#include  "chk.h"
#include  "numfunc.h"
#include "option.h"
  
#define TYPEOPT STDf

typedef struct TYPEOPT
{
  VAR PayOff;
  VAR EuOrAm;   
  VAR OMaturity;
  VAR BMaturity;
  VAR Nominal;
  VAR FixedRate;
  VAR ResetPeriod;
  VAR FirstResetDate;
  VAR NbResetDate;
  VAR FloatingRate;

} TYPEOPT;

int OPT(Get)(int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt, Model *mod);    
int OPT(Show)(int user,Planning *pt_plan,Option *opt, Model *mod);      
int OPT(Check)(int user,Planning *pt_plan,Option *opt);    

#endif
