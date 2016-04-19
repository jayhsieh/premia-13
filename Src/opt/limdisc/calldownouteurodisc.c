#include "limdisc.h"

static NumFunc_1 call=
  {
    Call,
    {{"Strike",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static NumFunc_1 rebate=
  {
    Const,
    {{"Rebate",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_digit
  };

static NumFunc_1 limit=
  {
    ConstLim,
    {
      {"StartingDate",DATE,{0},ALLOW,SETABLE},
      {"FinalDate",DATE,{0},IRRELEVANT,UNSETABLE},
      {"Frequency",INT2,{0},ALLOW,SETABLE},
      {"Limit",PDOUBLE,{90},ALLOW,SETABLE},
      {" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}
    },
    CHK_digit
  };

static TYPEOPT CallDownOutDiscEuro=
  {
    /*Maturity*/    {"Maturity",DATE,{0},ALLOW,SETABLE},
    /*Limit*/       {"Limit",NUMFUNC_1,{0},FORBID,SETABLE},
    /*PayOff*/      {"Payoff",NUMFUNC_1,{0},FORBID,SETABLE},
    /*Rebate*/      {"Rebate",NUMFUNC_1,{0},FORBID,SETABLE},

    /*OutOrIn*/     {"Out",BOOL,{OUT},FORBID,UNSETABLE},
    /*DownOrUp*/    {"Down",BOOL,{DOWN},FORBID,UNSETABLE},
    /*RebNo*/       {"Rebate",BOOL,{NOREBATE},FORBID,UNSETABLE},
    /*EuOrAm*/      {"Euro",BOOL,{EURO},FORBID,UNSETABLE},
    /*PartOrTot*/   {"Total",BOOL,{TOTAL},FORBID,UNSETABLE},
    /*ContOrDisc*/  {"Disc",BOOL,{DISC},FORBID,UNSETABLE},
    /*ConstLim*/    {"ConstLim",BOOL,{CONSTLIM},ALLOW,UNSETABLE},
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if (opt->init == 0 )
    {
      opt->init = 1;
      opt->HelpFilenameHint = "calldownouteurodisc";
      opt->nvar = 11;
      opt->nvar_setable = 4;

      pt->PayOff.Val.V_NUMFUNC_1=&call;
      pt->Rebate.Val.V_NUMFUNC_1=&rebate;
      pt->Limit.Val.V_NUMFUNC_1=&limit;

      (pt->PayOff.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=100.0;
      (pt->Rebate.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=10.0;

      (pt->OutOrIn).Val.V_BOOL=OUT;
      (pt->DownOrUp).Val.V_BOOL=DOWN;
      (pt->RebOrNo).Val.V_BOOL=NOREBATE;
      (pt->EuOrAm).Val.V_BOOL=EURO;
      (pt->PartOrTot).Val.V_BOOL=TOTAL;
      (pt->ContOrDisc).Val.V_BOOL=DISC;
      (pt->ConstLim).Val.V_BOOL=CONSTLIM;

      (pt->Limit.Val.V_NUMFUNC_1)->Par[0].Val.V_DATE=0.;
      (pt->Limit.Val.V_NUMFUNC_1)->Par[2].Val.V_INT2=100;
      (pt->Limit.Val.V_NUMFUNC_1)->Par[3].Val.V_PDOUBLE=90.0;

      (pt->Maturity).Val.V_DATE=1.0;

      /* test for setability */
      if ((pt->RebOrNo).Val.V_BOOL==REBATE) 
	pt->Rebate.Vsetable=SETABLE;
      else 
	pt->Rebate.Vsetable=UNSETABLE;

    }

  return OK;
}

MAKEOPT(CallDownOutDiscEuro);
