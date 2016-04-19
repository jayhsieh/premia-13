#include  "stdi.h"
static NumFunc_1 call=
  {
    Call,
    {{"Strike",PDOUBLE,{100},FORBID,UNSETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static TYPEOPT Floor=
  {
    {"Payoff",NUMFUNC_1,{0},FORBID,SETABLE},     /* PayOff; */
    {"Euro",BOOL,{EURO},FORBID,SETABLE},         /* EuOrAm */
    {"Option Maturity",DATE,{0},FORBID,SETABLE}, /* OMaturity;*/  
    {"Contract Maturity",DATE,{0},ALLOW,SETABLE},/* BMaturity;*/  
    {"Nominal Value",PDOUBLE,{0},ALLOW,SETABLE}, /* Nominal;*/  
    {"Strike",PDOUBLE,{0},ALLOW,SETABLE},      /* FixedRate;*/  
    {"Reset Period",PDOUBLE,{0},ALLOW,SETABLE},  /* ResetPeriod;*/  
    {"First Reset Date",DATE,{0},ALLOW,SETABLE}, /* FirstResetDate;*/  
    {"Nb of Reset",PINT,{0},FORBID,SETABLE},     /* NbResetDate;*/  
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if ( opt->init == 0)
    {
      opt->init = 1;
      opt->nvar = 9;
      opt->nvar_setable=9;

      pt->PayOff.Val.V_NUMFUNC_1=&call;

      (pt->EuOrAm).Val.V_BOOL=AMER;
      (pt->OMaturity).Val.V_DATE=7.0;
      (pt->BMaturity).Val.V_DATE=7.0;
      (pt->Nominal).Val.V_PDOUBLE=1.0;
      (pt->FixedRate).Val.V_PDOUBLE=0.07;
      (pt->ResetPeriod).Val.V_PDOUBLE=0.5;
      (pt->FirstResetDate).Val.V_DATE=0.5;
      (pt->NbResetDate).Val.V_PINT=10;
      (pt->PayOff.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=0.22313;

      /* the following variables are set interactively or not */
            
      pt->PayOff.Vsetable=UNSETABLE; 
      pt->EuOrAm.Vsetable=UNSETABLE; 
      pt->OMaturity.Vsetable=UNSETABLE;
      pt->BMaturity.Vsetable=SETABLE;
      pt->Nominal.Vsetable=SETABLE;
      pt->FixedRate.Vsetable=SETABLE;
      pt->ResetPeriod.Vsetable=SETABLE;
      pt->FirstResetDate.Vsetable=SETABLE;
      pt->NbResetDate.Vsetable=UNSETABLE;


    }

  return OK;
}

MAKEOPT(Floor);
