#include "bharchiarella1d.h"
#include "chk.h"
#include "error_msg.h"
#include "model.h"

extern char* path_sep;


static int MOD(Init)(Model *model)       
{
  TYPEMOD* pt=(TYPEMOD*)(model->TypeModel);
  
  if (model->init == 0 )
    {
      model->init = 1;
      model->nvar=0;
      pt->T.Vname = "Current Date";
      pt->T.Vtype=DATE;
      pt->T.Val.V_DATE=0.0;
      pt->T.Viter=ALLOW;
      model->nvar++;
 
      pt->alpha0.Vname = "Volatility Parameter Alpha0";
      pt->alpha0.Vtype=PDOUBLE;
      pt->alpha0.Val.V_PDOUBLE=0.001;
      pt->alpha0.Viter=ALLOW;
      model->nvar++;
      
      pt->alphar.Vname = "Volatility Parameter Alphar";
      pt->alphar.Vtype=PDOUBLE;
      pt->alphar.Val.V_PDOUBLE=0.04;
      pt->alphar.Viter=ALLOW;
      model->nvar++;
      
      pt->alphaf.Vname = "Volatility Parameter Alphaf";
      pt->alphaf.Vtype=PDOUBLE;
      pt->alphaf.Val.V_PDOUBLE=0.0;
      pt->alphaf.Viter=ALLOW;
      model->nvar++;

      pt->gamm.Vname = "Volatility Parameter Gamma";
      pt->gamm.Vtype=PDOUBLE;
      pt->gamm.Val.V_PDOUBLE=0.5;
      pt->gamm.Viter=ALLOW;
      model->nvar++;

      
      pt->lambda.Vname = "Volatility Parameter Lambda";
      pt->lambda.Vtype=PDOUBLE;
      pt->lambda.Val.V_PDOUBLE=0.2;
      pt->lambda.Viter=ALLOW;
      model->nvar++;

      pt->tau.Vname = "Volatility Parameter  Tau";
      pt->tau.Vtype=PDOUBLE;
      pt->tau.Val.V_PDOUBLE=1.;
      pt->tau.Viter=ALLOW;
      model->nvar++;
      
      pt->beta0.Vname = "Forward Rate Parameter Beta0";
      pt->beta0.Vtype=PDOUBLE;
      pt->beta0.Val.V_PDOUBLE=0.04;
      pt->beta0.Viter=ALLOW;
      model->nvar++;

      pt->beta1.Vname = "Forward RateParameter Beta1";
      pt->beta1.Vtype=PDOUBLE;
      pt->beta1.Val.V_PDOUBLE=0.04;
      pt->beta1.Viter=ALLOW;
      model->nvar++;

      pt->eta.Vname = "Forward Rate Parameter Eta";
      pt->eta.Vtype=PDOUBLE;
      pt->eta.Val.V_PDOUBLE=0.05;
      pt->eta.Viter=ALLOW;
      model->nvar++;
      
    }
  return OK;
}

TYPEMOD BharChiarella1d;

MAKEMOD(BharChiarella1d);


