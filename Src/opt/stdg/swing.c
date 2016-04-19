#include  "stdg.h"

static NumFunc_1 put=
  {
    Put,
    {{"Strike",PDOUBLE,{100},ALLOW,SETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static TYPEOPT Swing=
  {
    /*PayOff*/      {"Payoff",NUMFUNC_1,{0},FORBID,SETABLE},
    /*Maturity*/        {"Maturity",DATE,{0},ALLOW,SETABLE},
        /*Future Maturity*/    {"Future Maturity",DATE,{0},ALLOW,SETABLE},
    /*EuOrAm*/       {"Amer",BOOL,{AMER},FORBID,UNSETABLE},
    /* NbExerciseDate;*/   {"Nb of Put Exercise",PINT,{0},ALLOW,SETABLE},   
    /* RefractingPeriod;*/ {"Refracting Period",SPDOUBLE,{0},ALLOW,SETABLE},   
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);
  
  if (opt->init == 0 )
    {
      opt->init = 1;
      opt->nvar = 6;
      opt->nvar_setable = 5;

      pt->PayOff.Val.V_NUMFUNC_1=&put;

      (pt->EuOrAm).Val.V_BOOL=AMER;
      (pt->Maturity).Val.V_DATE=1.0;
      (pt->FutureMaturity).Val.V_DATE=0.5;
      (pt->PayOff.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=100.0;
      (pt->NbExerciseDate).Val.V_PINT=2;
      (pt->RefractingPeriod).Val.V_SPDOUBLE=0.1;
      
      /* the following variables are not set interactively */            
      pt->PayOff.Vsetable=UNSETABLE;
      pt->EuOrAm.Vsetable=UNSETABLE;
      pt->Maturity.Vsetable=SETABLE;
      pt->FutureMaturity.Vsetable=UNSETABLE;
      pt->RefractingPeriod.Vsetable=SETABLE;
      pt->NbExerciseDate.Vsetable=SETABLE;

    }

  return OK;
}

MAKEOPT(Swing);
