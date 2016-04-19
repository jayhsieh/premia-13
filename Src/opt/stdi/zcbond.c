#include  "stdi.h"
static NumFunc_1 call=
  {
    Call,
    {{"Strike",PDOUBLE,{100},FORBID,UNSETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
  };

static TYPEOPT ZeroCouponBond=
  {
    {"Payoff",NUMFUNC_1,{0},FORBID,UNSETABLE},     /* PayOff; */
    {"Euro",BOOL,{EURO},FORBID,UNSETABLE},         /* EuOrAm */
    {"Option Maturity",DATE,{0},FORBID,UNSETABLE}, /* OMaturity;*/  
    {"Contract Maturity",DATE,{0},ALLOW,SETABLE},/* BMaturity;*/  
    {"Nominal Value",PDOUBLE,{0},ALLOW,UNSETABLE}, /* Nominal;*/  
    {"Strike",PDOUBLE,{0},ALLOW,UNSETABLE},      /* FixedRate;*/  
    {"Reset Period",PDOUBLE,{0},ALLOW,UNSETABLE},  /* ResetPeriod;*/  
    {"First Reset Date",DATE,{0},ALLOW,UNSETABLE}, /* FirstResetDate;*/  
    {"Nb of Reset",PINT,{0},FORBID,UNSETABLE},     /* NbResetDate;*/  
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

  if ( opt->init == 0)
    {
      opt->init = 1;
      opt->nvar = 9;
      opt->nvar_setable=1;

      pt->PayOff.Val.V_NUMFUNC_1=&call;
 
      (pt->EuOrAm).Val.V_BOOL=EURO;
      (pt->OMaturity).Val.V_DATE=1.0;
      (pt->BMaturity).Val.V_DATE=1.0;
      (pt->Nominal).Val.V_PDOUBLE=1.0;
      (pt->FixedRate).Val.V_PDOUBLE=1.0;
      (pt->ResetPeriod).Val.V_PDOUBLE=1.0;
      (pt->FirstResetDate).Val.V_DATE=1.0;
      (pt->NbResetDate).Val.V_INT=1;
      (pt->PayOff.Val.V_NUMFUNC_1)->Par[0].Val.V_PDOUBLE=0.22313;

	  opt->HelpFilenameHint = "zcbond";
  }

  return OK;
}

MAKEOPT(ZeroCouponBond);
