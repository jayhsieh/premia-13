#include "schwartz.h"
#include "chk.h"
#include "error_msg.h"
#include "model.h"

static int MOD(Init)(Model *model)       
{
  TYPEMOD* pt=(TYPEMOD*)(model->TypeModel);
    
  if (model->init == 0 )
    {
      model->init = 1;
      model->nvar=0;
      pt->T.Vname = "Current Date";
      pt->T.Vtype=DATE;
      pt->T.Val.V_DATE=0.;
      pt->T.Viter=ALLOW;
      model->nvar++;
     
      pt->Divid.Vname = "Annual Dividend Rate";
      pt->Divid.Vtype=DOUBLE;
      pt->Divid.Val.V_DOUBLE=0;
      pt->Divid.Viter=ALLOW;
      model->nvar++;
      
      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=4.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->sigmad.Vname = "Sigma Delta";
      pt->sigmad.Vtype=PDOUBLE;
      pt->sigmad.Val.V_PDOUBLE=0.0556;
      pt->sigmad.Viter=ALLOW;
      model->nvar++;

      pt->sigmas.Vname = "Sigma S";
      pt->sigmas.Vtype=PDOUBLE;
      pt->sigmas.Val.V_PDOUBLE=0.374;
      pt->sigmas.Viter=ALLOW;
      model->nvar++;

      pt->alpha.Vname = "Alpha";
      pt->alpha.Vtype=PDOUBLE;
      pt->alpha.Val.V_PDOUBLE=1.829;
      pt->alpha.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Correlation";
      pt->Rho.Vtype=RGDOUBLEM11;
      pt->Rho.Val.V_RGDOUBLEM11=0.882;
      pt->Rho.Viter=ALLOW;
      model->nvar++;

    }
  return OK;
}

TYPEMOD schwartz;
MAKEMOD(schwartz);
         

