#include  "callable.h"

static TYPEOPT IntermittentCall =
{
    {"Strike",PDOUBLE,{0},ALLOW,SETABLE}, 
    {"Maturity (in days)",INT,{0},ALLOW,SETABLE}, 
    {"Nominal Coupon Rate",PDOUBLE,{0},ALLOW,SETABLE}, 
    {"Nominal Recovery",PDOUBLE,{0},ALLOW,SETABLE}, 
    {"PutStrike",PDOUBLE,{0},ALLOW,SETABLE}, 
    {"Calltrike",PDOUBLE,{0},ALLOW,SETABLE}, 
    {"LowerBarrier",PDOUBLE,{0},ALLOW,SETABLE}, 
    {"UpperBarrier",PDOUBLE,{0},ALLOW,SETABLE}, 
    {"Window (in days)",INT,{0},IRRELEVANT,UNSETABLE}, 
    {"Period (in days)",INT,{0},IRRELEVANT,UNSETABLE}, 
};

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(TYPEOPT*)(opt->TypeOpt);

  if (opt->init == 0 )
    {
      opt->init = 1;
      opt->nvar = 10;
      opt->nvar_setable = 8;

      pt->Strike.Val.V_PDOUBLE=100.; 
      pt->Maturity.Val.V_INT=180; 
      pt->Coupon.Val.V_PDOUBLE=14.4; 
      pt->Recovery.Val.V_PDOUBLE=0.; 
      pt->PutStrike.Val.V_PDOUBLE=0.; 
      pt->CallStrike.Val.V_PDOUBLE=103.; 
      pt->LowerBarrier.Val.V_PDOUBLE=97.; 
      pt->UpperBarrier.Val.V_PDOUBLE=103.; 
      pt->Window.Val.V_INT=0; 
      pt->Period.Val.V_INT=0; 
    }

  return OK;
}

MAKEOPT(IntermittentCall);
