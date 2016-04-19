#include  "pad.h"

static NumFunc_2 put=
  {
    Put_OverSpot2,
    {{"Strike",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static NumFunc_2 MovingAverage=
  {
    Asian,
    {
      {"Nb Dates",PINT,{0},ALLOW,SETABLE},
      {"Window",PINT,{0},ALLOW,SETABLE},
      {"Delay",INT,{0},ALLOW,SETABLE},
      {" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}
    },
    CHK_call
  };

static TYPEOPT MovingAveragePutFixedAmer=
  {
    /*Maturity*/    {"Maturity",DATE,{0},ALLOW,SETABLE},
    /*PayOff*/      {"Payoff",NUMFUNC_2,{0},FORBID,SETABLE},
    /*PathDep*/     {"PathDep",NUMFUNC_2,{0},FORBID,SETABLE},

    /*MinOrElse*/   {"Average",PADE,{AVERAGE},ALLOW,UNSETABLE},
    /*EuOrAm*/      {"Amer",BOOL,{AMER},FORBID,UNSETABLE},
    /*PartOrTot*/   {"Total",BOOL,{TOTAL},FORBID,UNSETABLE},
    /*ContOrDisc*/  {"Continuous",BOOL,{CONT},FORBID,UNSETABLE},

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
      pt->PathDep.Val.V_NUMFUNC_2=&MovingAverage;

      (pt->MinOrElse).Val.V_PADE=AVERAGE;
      (pt->EuOrAm).Val.V_BOOL=AMER;
      (pt->PartOrTot).Val.V_BOOL=TOTAL;
      (pt->ContOrDisc).Val.V_BOOL=CONT;

      (pt->PayOff.Val.V_NUMFUNC_2)->Par[0].Val.V_DOUBLE=100;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[0].Val.V_PINT=50;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[1].Val.V_PINT=5;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[2].Val.V_INT=0;

      (pt->Maturity).Val.V_DATE=1.0;

    }

  return OK;
}

MAKEOPT(MovingAveragePutFixedAmer);
