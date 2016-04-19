#include "hescir1d.h"
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

      pt->S0.Vname = "Spot";
      pt->S0.Vtype=PDOUBLE;
      pt->S0.Val.V_PDOUBLE=100.;
      pt->S0.Viter=ALLOW;
      model->nvar++;
 
      pt->r0.Vname = "Current Interest Rate";
      pt->r0.Vtype=PDOUBLE;
      pt->r0.Val.V_PDOUBLE=0.04;
      pt->r0.Viter=ALLOW;
      model->nvar++;

      pt->kr.Vname = "Interest Rate Speed of Mean Reversion";
      pt->kr.Vtype=PDOUBLE;
      pt->kr.Val.V_PDOUBLE=0.3;
      pt->kr.Viter=ALLOW;
      model->nvar++;

      pt->thetar.Vname = "Interest Rate Long Term Mean";
      pt->thetar.Vtype=PDOUBLE;
      pt->thetar.Val.V_PDOUBLE=0.04;
      pt->thetar.Viter=ALLOW;
      model->nvar++;

      pt->Sigmar.Vname = "Interest Rate Volatility";
      pt->Sigmar.Vtype=PDOUBLE;
      pt->Sigmar.Val.V_PDOUBLE=0.1;
      pt->Sigmar.Viter=ALLOW;
      model->nvar++;
      
      pt->V0.Vname = "Current Variance";
      pt->V0.Vtype=DOUBLE;
      pt->V0.Val.V_DOUBLE=0.1;
      pt->V0.Viter=ALLOW;
      model->nvar++;

      pt->kV.Vname = "Mean Reversion Variance";
      pt->kV.Vtype=DOUBLE;
      pt->kV.Val.V_DOUBLE=1.5;
      pt->kV.Viter=ALLOW;
      model->nvar++;

      pt->thetaV.Vname = "Long-Run Variance";
      pt->thetaV.Vtype=DOUBLE;
      pt->thetaV.Val.V_DOUBLE=0.02;
      pt->thetaV.Viter=ALLOW;
      model->nvar++;

      pt->SigmaV.Vname = "Volatility of Variance";
      pt->SigmaV.Vtype=DOUBLE;
      pt->SigmaV.Val.V_DOUBLE=0.2;
      pt->SigmaV.Viter=ALLOW;
      model->nvar++;

      pt->RhoSr.Vname = "Rho S r";
      pt->RhoSr.Vtype=RGDOUBLEM11;
      pt->RhoSr.Val.V_RGDOUBLEM11=0.1;
      pt->RhoSr.Viter=ALLOW;
      model->nvar++;

      pt->RhoSV.Vname = "Rho S V";
      pt->RhoSV.Vtype=RGDOUBLEM11;
      pt->RhoSV.Val.V_RGDOUBLEM11=0.;
      pt->RhoSV.Viter=ALLOW;
      model->nvar++;

      pt->RhorV.Vname = "Rho r V";
      pt->RhorV.Vtype=RGDOUBLEM11;
      pt->RhorV.Val.V_RGDOUBLEM11=0.0;
      pt->RhorV.Viter=ALLOW;
      model->nvar++;
          }

  return OK;
}

TYPEMOD hescir1d;
MAKEMOD(hescir1d);
