#include  "std2dg.h"

static NumFunc_2 callspread2d=
  {
    CallSpread2d,
    {{"Strike",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static TYPEOPT CallSpread2dEuro=
  {
    /*Maturity*/        {"Maturity",DATE,{0},ALLOW,SETABLE},
    /*PayOff*/      {"Payoff",NUMFUNC_2,{0},FORBID,SETABLE},
    /*EuOrAm*/       {"Euro",BOOL,{EURO},FORBID,UNSETABLE},

  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if ( opt->init == 0)
    {
      opt->init = 1;
      opt->nvar = 3;
      opt->nvar_setable=2;


      pt->PayOff.Val.V_NUMFUNC_2=&callspread2d;

      (pt->EuOrAm).Val.V_BOOL=EURO;
      (pt->Maturity).Val.V_DATE=1.0;
      (pt->PayOff.Val.V_NUMFUNC_2)->Par[0].Val.V_PDOUBLE=2.;

    }

  return OK;
}

MAKEOPT(CallSpread2dEuro);
