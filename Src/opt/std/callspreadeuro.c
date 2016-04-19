#include  "std.h"

static NumFunc_1 callspread=
  {
    CallSpread,
    {{"Strike 1",PDOUBLE,{100},ALLOW,SETABLE},
     {"Strike 2",PDOUBLE,{110},ALLOW,SETABLE},
     {" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}
    },
    CHK_callspread
  };

static TYPEOPT CallSpreadEuro=
  {
    /*PayOff*/      {"Payoff",NUMFUNC_1,{0},FORBID,SETABLE},
    /*Maturity*/    {"Maturity",DATE,{0},ALLOW,SETABLE},
    /*EuOrAm*/      {"Euro",BOOL,{EURO},FORBID,UNSETABLE}
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if (opt->init == 0 )
    {
      opt->init = 1;
      opt->nvar = 3;
      opt->nvar_setable = 2;

      pt->PayOff.Val.V_NUMFUNC_1=&callspread;

      (pt->EuOrAm).Val.V_BOOL=EURO;
      (pt->Maturity).Val.V_DATE=1.0;
      (pt->PayOff.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=100.0;
      (pt->PayOff.Val.V_NUMFUNC_1)->Par[1].Val.V_PDOUBLE=110.0;

    }

  return OK;
}

MAKEOPT(CallSpreadEuro);
