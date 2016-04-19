#include  "std2d.h"

static NumFunc_2 exchange=
  {
    Exchange,
    {
      {"Ratio",PDOUBLE,{100},ALLOW,SETABLE}
    },
    CHK_call
  };

static TYPEOPT ExchangeAmer=
  {
    /*Maturity*/        {"Maturity",DATE,{0},ALLOW,SETABLE},
    /*PayOff*/      {"Payoff",NUMFUNC_2,{0},FORBID,SETABLE},
    /*EuOrAm*/       {"Amer",BOOL,{AMER},FORBID,UNSETABLE},

  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if ( opt->init == 0)
    {
      opt->init = 1;
      opt->nvar = 3;
      opt->nvar_setable=2;

      pt->PayOff.Val.V_NUMFUNC_2=&exchange;
      
      (pt->EuOrAm).Val.V_BOOL=AMER;
      (pt->Maturity).Val.V_DATE=1.0;
      (pt->PayOff.Val.V_NUMFUNC_2)->Par[0].Val.V_PDOUBLE=0.75;

    }

  return OK;
}

MAKEOPT(ExchangeAmer);
