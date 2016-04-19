#include  "pad.h"

static NumFunc_2 put=
  {
    Put_OverSpot2,  /*(K-Minimum)+*/
    {{"Strike",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static NumFunc_2 minimum=
  {
    Minimum,
    {
      {"StartingDate",DATE,{0},IRRELEVANT,UNSETABLE},
      {"FinalDate",DATE,{0},IRRELEVANT,UNSETABLE},
      {"Frequency",PDOUBLE,{0},IRRELEVANT,UNSETABLE},
      {"InitialValue",PDOUBLE,{100},IRRELEVANT,UNSETABLE},
      {"Minimum",PDOUBLE,{100},ALLOW,SETABLE},
      {" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}
    },
    CHK_call
  };


TYPEOPT LookBackPutFixedEuro=
  {
    /*Maturity*/        {"Maturity",DATE,{0},ALLOW,SETABLE},
    /*PayOff*/      {"Payoff",NUMFUNC_2,{0},FORBID,SETABLE},
    /*PathDep*/     {"PathDep",NUMFUNC_2,{0},FORBID,SETABLE},

    /*MinOrElse*/    {"Minimum",PADE,{MINIMUM},ALLOW,UNSETABLE},
    /*EuOrAm*/       {"Euro",BOOL,{EURO},FORBID,UNSETABLE},
    /*PartOrTot*/    {"Total",BOOL,{TOTAL},FORBID,UNSETABLE},
    /*ContOrDisc*/   {"Continuous",BOOL,{CONT},FORBID,UNSETABLE},
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if ( opt->init == 0)
    {
      opt->init = 1;
      opt->nvar = 7;
      opt->nvar_setable=3;


      pt->PayOff.Val.V_NUMFUNC_2=&put;
      pt->PathDep.Val.V_NUMFUNC_2=&minimum;

      (pt->PayOff.Val.V_NUMFUNC_2)->Par[0].Val.V_PDOUBLE=100.0;

      (pt->MinOrElse).Val.V_PADE=MINIMUM;
      (pt->EuOrAm).Val.V_BOOL=EURO;
      (pt->PartOrTot).Val.V_BOOL=TOTAL;
      (pt->ContOrDisc).Val.V_BOOL=CONT;

      (pt->PathDep.Val.V_NUMFUNC_2)->Par[0].Val.V_DATE=0.0;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[1].Val.V_DATE=0.0;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[2].Val.V_PDOUBLE=0.0;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[3].Val.V_PDOUBLE=100.0;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[4].Val.V_PDOUBLE=100.0;

      (pt->Maturity).Val.V_DATE=1.0;

    }

  return OK;
}

MAKEOPT(LookBackPutFixedEuro);
