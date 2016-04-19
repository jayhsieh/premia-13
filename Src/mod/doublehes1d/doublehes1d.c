#include "doublehes1d.h"
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
      pt->R.Val.V_DOUBLE=3.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->Sigma0.Vname = "Current Variance 1";
      pt->Sigma0.Vtype=DOUBLE;
      pt->Sigma0.Val.V_DOUBLE=0.36;
      pt->Sigma0.Viter=ALLOW;
      model->nvar++;

      pt->MeanReversion.Vname = "Mean Reversion of Variance 1";
      pt->MeanReversion.Vtype=DOUBLE;
      pt->MeanReversion.Val.V_DOUBLE=0.9;
      pt->MeanReversion.Viter=ALLOW;
      model->nvar++;

      pt->LongRunVariance.Vname = "Long-Run of Variance 1";
      pt->LongRunVariance.Vtype=DOUBLE;
      pt->LongRunVariance.Val.V_DOUBLE=0.1;
      pt->LongRunVariance.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Volatility of Variance 1";
      pt->Sigma.Vtype=DOUBLE;
      pt->Sigma.Val.V_DOUBLE=0.1;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Rho Spot-Variance 1";
      pt->Rho.Vtype=DOUBLE;
      pt->Rho.Val.V_DOUBLE=-0.5;
      pt->Rho.Viter=ALLOW;
      model->nvar++;

      pt->Sigma0V.Vname = "Current Variance 2";
      pt->Sigma0V.Vtype=DOUBLE;
      pt->Sigma0V.Val.V_DOUBLE=0.49;
      pt->Sigma0V.Viter=ALLOW;
      model->nvar++;

      pt->MeanReversionV.Vname = "Mean Reversion  of Variance 2";
      pt->MeanReversionV.Vtype=DOUBLE;
      pt->MeanReversionV.Val.V_DOUBLE=1.2;
      pt->MeanReversionV.Viter=ALLOW;
      model->nvar++;

      pt->LongRunVarianceV.Vname = "Long-Run of Variance 2";
      pt->LongRunVarianceV.Vtype=DOUBLE;
      pt->LongRunVarianceV.Val.V_DOUBLE=0.15;
      pt->LongRunVarianceV.Viter=ALLOW;
      model->nvar++;

      pt->SigmaV.Vname = "Volatility of Variance 2";
      pt->SigmaV.Vtype=DOUBLE;
      pt->SigmaV.Val.V_DOUBLE=0.2;
      pt->SigmaV.Viter=ALLOW;
      model->nvar++;

      pt->RhoSV2.Vname = "Rho Spot-Variance 2";
      pt->RhoSV2.Vtype=DOUBLE;
      pt->RhoSV2.Val.V_DOUBLE=-0.5;
      pt->RhoSV2.Viter=ALLOW;
      model->nvar++;
      
      /* pt->RhoVV.Vname = "Rho Variance-Variance of Variance"; */
      /* pt->RhoVV.Vtype=DOUBLE; */
      /* pt->RhoVV.Val.V_DOUBLE=0.5; */
      /* pt->RhoVV.Viter=ALLOW; */
      //model->nvar++;
          }

  return OK;
}

TYPEMOD DoubleHeston1d;
MAKEMOD(DoubleHeston1d);
