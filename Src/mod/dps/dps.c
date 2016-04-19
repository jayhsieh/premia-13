#include "dps.h"
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

      pt->Divid.Vname = "Annual Dividend Rate";
      pt->Divid.Vtype=DOUBLE;
      pt->Divid.Val.V_DOUBLE=0.;
      pt->Divid.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=10.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Rho";
      pt->Rho.Vtype=DOUBLE;
      pt->Rho.Val.V_DOUBLE=0.5;
      pt->Rho.Viter=ALLOW;
      model->nvar++;
      
      pt->Sigma0.Vname = "Current Variance";
      pt->Sigma0.Vtype=DOUBLE;
      pt->Sigma0.Val.V_DOUBLE=0.01;
      pt->Sigma0.Viter=ALLOW;
      model->nvar++;

      pt->Kappa.Vname = "Mean Reversion";
      pt->Kappa.Vtype=DOUBLE;
      pt->Kappa.Val.V_DOUBLE=2.;
      pt->Kappa.Viter=ALLOW;
      model->nvar++;

      pt->Eta.Vname = "Long-Run Variance";
      pt->Eta.Vtype=DOUBLE;
      pt->Eta.Val.V_DOUBLE=0.01;
      pt->Eta.Viter=ALLOW;
      model->nvar++;

      pt->Theta.Vname = "Volatility of Volatility";
      pt->Theta.Vtype=DOUBLE;
      pt->Theta.Val.V_DOUBLE=0.2;
      pt->Theta.Viter=ALLOW;
      model->nvar++;

      pt->LambdaS.Vname = "Lambda Spot Jump";
      pt->LambdaS.Vtype=DOUBLE;
      pt->LambdaS.Val.V_DOUBLE=0.1382;
      pt->LambdaS.Viter=ALLOW;
      model->nvar++;

      pt->MeanS.Vname = "Mean Spot Jump";
      pt->MeanS.Vtype=DOUBLE;
      pt->MeanS.Val.V_DOUBLE=0.1791;
      pt->MeanS.Viter=ALLOW;
      model->nvar++;

      pt->SigmaS.Vname = "Variance Spot Jump";
      pt->SigmaS.Vtype=DOUBLE;
      pt->SigmaS.Val.V_DOUBLE=0.1346;
      pt->SigmaS.Viter=ALLOW;
      model->nvar++;


      pt->LambdaV.Vname = "Lambda Variance Jump";
      pt->LambdaV.Vtype=DOUBLE;
      pt->LambdaV.Val.V_DOUBLE=0.;
      pt->LambdaV.Viter=ALLOW;
      model->nvar++;

      pt->MeanV.Vname = "Mean Variance Jump";
      pt->MeanV.Vtype=DOUBLE;
      pt->MeanV.Val.V_DOUBLE=1.;
      pt->MeanV.Viter=ALLOW;
      model->nvar++;

      pt->LambdaSV.Vname = "Lambda Spot-Variance Jump correlated ";
      pt->LambdaSV.Vtype=DOUBLE;
      pt->LambdaSV.Val.V_DOUBLE=0.;
      pt->LambdaSV.Viter=ALLOW;
      model->nvar++;

      pt->MeanSV.Vname = "Mean Spot Jump correlated ";
      pt->MeanSV.Vtype=DOUBLE;
      pt->MeanSV.Val.V_DOUBLE=0.1;
      pt->MeanSV.Viter=ALLOW;
      model->nvar++;

      pt->SigmaSV.Vname = "Variance Spot Jump correlated ";
      pt->SigmaSV.Vtype=DOUBLE;
      pt->SigmaSV.Val.V_DOUBLE=0.16;
      pt->SigmaSV.Viter=ALLOW;
      model->nvar++;

      pt->MeanVS.Vname = "Mean Variance Jump correlated ";
      pt->MeanVS.Vtype=DOUBLE;
      pt->MeanVS.Val.V_DOUBLE=0.1;
      pt->MeanVS.Viter=ALLOW;
      model->nvar++;


      pt->RhoSV.Vname = "Rho_Jump";
      pt->RhoSV.Vtype=DOUBLE;
      pt->RhoSV.Val.V_DOUBLE=0.5;
      pt->RhoSV.Viter=ALLOW;
      model->nvar++;

     }

  return OK;
}

TYPEMOD dps;
MAKEMOD(dps);
