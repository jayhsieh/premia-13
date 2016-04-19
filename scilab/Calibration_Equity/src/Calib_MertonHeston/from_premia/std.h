#ifndef  _STD_H
#define _STD_H

#include  "optype.h"
#include  "var.h"

#include  "chk.h"
#include  "numfunc.h"

#define TYPEOPT STD

typedef struct TYPEOPT
	{         
	VAR      PayOff;
	VAR					  EuOrAm;   
	VAR								Maturity;      
	} TYPEOPT;

int OPT(Get)(int user,Planning *pt_plan,Option *opt);    
int OPT(FGet)(char **InputFile,int user,Planning *pt_plan,Option *opt);    
int OPT(Show)(int user,Planning *pt_plan,Option *opt);      
int OPT(Check)(int user,Planning *pt_plan,Option *opt);    

#endif
