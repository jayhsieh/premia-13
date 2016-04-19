#ifndef _BS1D_H
#define _BS1D_H

#include "optype.h"
#include "var.h"

#define TYPEMOD BS1D

/*1D BlackScholes World*/     
typedef struct  TYPEMOD{ 
	VAR				 T;
	VAR              S0;
	VAR				Mu;
	VAR             Sigma;
	VAR              Divid;
	VAR              R;       
	} TYPEMOD;

int MOD(Get)(int user,Planning* pt_plan,Model *model);     
int MOD(FGet)(char **InputFile,int user,Planning* pt_plan,Model *model);     
int MOD(Show)(int user,Planning *pt_plan,Model *model);    
int MOD(Check)(int user,Planning *pt_plan,Model *model);     
int MOD(Init)(Model *model) ;

#endif
