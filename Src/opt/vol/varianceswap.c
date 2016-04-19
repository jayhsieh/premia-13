#include  "vol.h"

static NumFunc_1 varianceswap=
  {
    Call,
    {{"Strike in volatilty terms",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static TYPEOPT VarianceSwap=
  {
    /*PayOff*/      {"Payoff",NUMFUNC_1,{0},FORBID,SETABLE},
    /*Maturity*/    {"Maturity",DATE,{0},ALLOW,SETABLE},
        /*Variance Budget*/    {"Variance Budget",PDOUBLE,{0},FORBID,UNSETABLE},

    /*EuOrAm*/      {"Euro",BOOL,{EURO},FORBID,UNSETABLE}
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if (opt->init == 0 )
    {
      opt->init = 1;
      opt->nvar = 4;
      opt->nvar_setable = 2;

      pt->PayOff.Val.V_NUMFUNC_1=&varianceswap;

      (pt->EuOrAm).Val.V_BOOL=EURO;
      (pt->Maturity).Val.V_DATE=1.0;
            (pt->VarianceBudget).Val.V_PDOUBLE=0.0625;

      (pt->PayOff.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=10.0;
    }

  return OK;
}

MAKEOPT(VarianceSwap);
