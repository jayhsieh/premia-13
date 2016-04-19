/////////////////CREDIT VAR and CTE////////////////////////////
#include  "stdr.h"

static TYPEOPT CreditVaRisk=
  {
	{"Maturity",DATE,{0},IRRELEVANT,UNSETABLE},
	{"Strike",PDOUBLE,{100},IRRELEVANT,UNSETABLE},
	{"Confidence Level",RGDOUBLE,{100},ALLOW,SETABLE},
	{"Number of Borrowers",PINT,{100},ALLOW,SETABLE},
    {"Constant Probability of Default",RGDOUBLE,{10},ALLOW,SETABLE},
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(TYPEOPT*)(opt->TypeOpt);

  if(opt->init == 0)
    {
      opt->init = 1;
      opt->nvar =5;
      opt->nvar_setable =3;

	  opt->HelpFilenameHint = "credit_varisk";
	  (pt->NumberofCreditors).Val.V_PINT=100;
	  (pt->ConstantProbabilityofDefault).Val.V_RGDOUBLE=0.0021;
	  (pt->Alpha).Val.V_RGDOUBLE=0.999;
    }

  return OK;
}

MAKEOPT(CreditVaRisk);
