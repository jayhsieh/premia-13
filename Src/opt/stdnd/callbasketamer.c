#include  "stdnd.h"

static NumFunc_nd callbasketamer_nd=
{
  CallBasket_nd,
  {{"Strike",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
  CHK_call
};

static TYPEOPT CallBasketAmer_nd=
{
  /*Size*/        {"Size", PINT, {1}, FORBID, UNSETABLE},
  /*Maturity*/    {"Maturity",DATE,{0},ALLOW,SETABLE},
  /*PayOff*/      {"Payoff",NUMFUNC_ND,{0},FORBID,SETABLE},
  /*EurOrAmer*/   {"Amer",BOOL,{AMER},FORBID,UNSETABLE},
};

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);
  VAR* ptMod=(VAR*)(mod->TypeModel);
  
  if ( opt->init == 0)
    {
      opt->init = 1;
       opt->HelpFilenameHint = "callbasketamer";
      opt->nvar = 4;
      opt->nvar_setable=2;

      pt->PayOff.Val.V_NUMFUNC_ND=&callbasketamer_nd;
      (pt->Maturity).Val.V_DATE=1.0;
      pt->EuOrAm.Val.V_BOOL=AMER;
      (pt->PayOff.Val.V_NUMFUNC_ND)->Par[0].Val.V_PDOUBLE=100.;
    }
  pt->Size.Val.V_PINT=ptMod[0].Val.V_INT;
  return OK;
}

MAKEOPT(CallBasketAmer_nd);
