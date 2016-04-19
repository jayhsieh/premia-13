#include  "pad.h"

static NumFunc_2 call=
  {
    Call_OverSpot2,
    {{"Strike",PDOUBLE,{100},FORBID,UNSETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static NumFunc_2 cliquet=
  {
    Asian,
    {
      {"Fg",PDOUBLE,{0},ALLOW,SETABLE},
      {"Cg",PDOUBLE,{0},ALLOW,SETABLE},
      {"Fl",PDOUBLE,{0},ALLOW,SETABLE},
      {"Cl",PDOUBLE,{100},ALLOW,SETABLE},
      {" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}
    },
    CHK_call
  };


static TYPEOPT Cliquet=
{
  {"Maturity(in years)",PINT,{0},ALLOW,SETABLE},     /*Maturity*/
  
   /*PayOff*/      {"Payoff",NUMFUNC_2,{0},FORBID,UNSETABLE},
  /*PathDep*/     {"PathDep",NUMFUNC_2,{0},FORBID,SETABLE},
  
  /*MinOrElse*/    {"Average",PADE,{AVERAGE},ALLOW,UNSETABLE},
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
      opt->nvar_setable=2;

      pt->PayOff.Val.V_NUMFUNC_2=&call;
      pt->PathDep.Val.V_NUMFUNC_2=&cliquet;

      (pt->PayOff.Val.V_NUMFUNC_2)->Par[0].Val.V_PDOUBLE=0.;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[0].Val.V_DATE=0.16;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[1].Val.V_DATE=1000.;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[2].Val.V_PDOUBLE=0.0;
      (pt->PathDep.Val.V_NUMFUNC_2)->Par[3].Val.V_PDOUBLE=0.08;
      

      (pt->Maturity).Val.V_PINT=5;
    
      (pt->MinOrElse).Val.V_PADE=AVERAGE;
      (pt->EuOrAm).Val.V_BOOL=EURO;
      (pt->PartOrTot).Val.V_BOOL=TOTAL;
      (pt->ContOrDisc).Val.V_BOOL=CONT;
      
      
      opt->HelpFilenameHint = "Cliquet";

    }
  return OK;
}

MAKEOPT(Cliquet);
