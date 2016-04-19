#include  "stdc.h"


static TYPEOPT CreditDefaultSwap=
{
    /*Maturity*/        {"Maturity",DATE,{0},ALLOW,SETABLE},
    /* Payement Period in Month */   {"Nb of Payements",PINT,{0},ALLOW,SETABLE},   
    /* Recovery */ {"Recovery",PDOUBLE,{0},ALLOW,SETABLE},   
  };

static int OPT(Init)(Option *opt,Model *mod)
{
  TYPEOPT* pt=(TYPEOPT*)(opt->TypeOpt);
  
  if (opt->init == 0 )
    {
      
      opt->init = 1;
       opt->HelpFilenameHint = "cds";
      opt->nvar = 3;
      opt->nvar_setable = 3;

      (pt->Maturity).Val.V_DATE=5.0;
      (pt->NbPayement).Val.V_PINT=4;
      (pt->Recovery).Val.V_PDOUBLE=0.4;
      
      /* the following variables are set interactively */            
      pt->Maturity.Vsetable=SETABLE;
      pt->Recovery.Vsetable=SETABLE;
      pt->NbPayement.Vsetable=SETABLE;

    }

  return OK;
}

MAKEOPT(CreditDefaultSwap);
