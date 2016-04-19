#include  "vol.h"

static NumFunc_1 timer=
  {
    Call,
    {{"Strike",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static TYPEOPT Timer=
  {
    /*PayOff*/      {"Payoff",NUMFUNC_1,{0},FORBID,SETABLE},
     /*Maturity*/    {"Maturity",DATE,{0},FORBID,UNSETABLE},
    /*Variance Budget*/    {"Variance Budget",PDOUBLE,{0},ALLOW,SETABLE},
    /*EuOrAm*/      {"",BOOL,{EURO},FORBID,UNSETABLE}
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if (opt->init == 0 )
    {
      opt->init = 1;
      opt->nvar = 4;
      opt->nvar_setable = 2;
      opt->HelpFilenameHint="Timer";

      pt->PayOff.Val.V_NUMFUNC_1=&timer;
      (pt->VarianceBudget).Val.V_PDOUBLE=0.0265;
      (pt->Maturity).Val.V_DATE=1.0;
      (pt->EuOrAm).Val.V_BOOL=EURO;
      (pt->PayOff.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=100.;
    }

  return OK;
}

MAKEOPT(Timer);
