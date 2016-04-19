#include  "stda.h"

static NumFunc_1 put=
  {
    Put,
    {{"Strike",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static TYPEOPT Gap=
  {
    /*PayOff*/      {"Payoff",NUMFUNC_1,{0},ALLOW,SETABLE},
	 /*EuOrAm*/       {"Euro",BOOL,{AMER},IRRELEVANT,UNSETABLE},
    /*Maturity*/     {"Maturity",DATE,{0},ALLOW,SETABLE},
	/*DeemedContribution;*/ {"Deemed Contribution",PDOUBLE,{0},IRRELEVANT,UNSETABLE},
    /*InitialAge*/ {"Initial Age",PDOUBLE,{0},IRRELEVANT,UNSETABLE},
    /*Premium;*/ {"Premium",PDOUBLE,{0},IRRELEVANT,UNSETABLE},
     /*MinimumGuaranteedInterestRate*/ {"MinimumGuaranteedInterestRate",PDOUBLE,{0},IRRELEVANT,UNSETABLE},
	 {"Number of Monitoring Dates",PINT,{0},ALLOW,SETABLE},
	/*Alpha*/ {"Alpha",RGDOUBLE,{0},ALLOW,SETABLE},
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if (opt->init == 0 )
    {
      opt->init = 1;
      opt->nvar = 9;
      opt->nvar_setable = 4;

      pt->PayOff.Val.V_NUMFUNC_1=&put;

      (pt->Maturity).Val.V_DATE=1.0;
      (pt->PayOff.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=90.0;
	  (pt->Alpha).Val.V_RGDOUBLE=0.9;
	  (pt->NumberofMonitoringDates).Val.V_PINT=252;
    }

  return OK;
}

MAKEOPT(Gap);
