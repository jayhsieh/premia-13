/////////////////VAR and CTE////////////////////////////
#include  "stdr.h"

static TYPEOPT VaRisk=
  {
    {"Maturity",DATE,{0},ALLOW,SETABLE},
	{"Strike",PDOUBLE,{100},ALLOW,SETABLE},
	{"Alpha",PDOUBLE,{100},ALLOW,SETABLE},
	{"Number of Creditors",PINT,{100},IRRELEVANT,UNSETABLE} ,
    {"Constant Probability of Default",RGDOUBLE,{10},IRRELEVANT,UNSETABLE},
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if(opt->init == 0 )
    {
      opt->init = 1;
      opt->nvar = 5;
      opt->nvar_setable = 3;


      (pt->Maturity).Val.V_DATE=1.0;
	  (pt->Strike).Val.V_PDOUBLE=100.0;
	  (pt->NumberofCreditors).Val.V_PINT=100;
    }

  return OK;
}

MAKEOPT(VaRisk);
