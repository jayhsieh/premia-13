#include "inflation_lmm_heston1d.h"
#include "chk.h"
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

      pt->I0.Vname = "Flat Initial Libor Forward CPI";
      pt->I0.Vtype=PDOUBLE;
      pt->I0.Val.V_PDOUBLE=0.08;
      pt->I0.Viter=ALLOW;
      model->nvar++;

      pt->SigmaI.Vname = "Volatility Libor Forward CPI";
      pt->SigmaI.Vtype=PDOUBLE;
      pt->SigmaI.Val.V_PDOUBLE=0.01;
      pt->SigmaI.Viter=ALLOW;
      model->nvar++;

      pt->F0.Vname = "Flat Initial Libor Nominal Rates";
      pt->F0.Vtype=PDOUBLE;
      pt->F0.Val.V_PDOUBLE=0.04;
      pt->F0.Viter=ALLOW;
      model->nvar++;

      pt->SigmaF.Vname = "Volatility Libor Nominal Rates";
      pt->SigmaF.Vtype=PDOUBLE;
      pt->SigmaF.Val.V_PDOUBLE=0.75;
      pt->SigmaF.Viter=ALLOW;
      model->nvar++;
      
      pt->Sigma0.Vname = "Current Variance";
      pt->Sigma0.Vtype=DOUBLE;
      pt->Sigma0.Val.V_DOUBLE=0.7;
      pt->Sigma0.Viter=ALLOW;
      model->nvar++;
      
      pt->SpeedMeanReversion.Vname = "Speed of Mean Reversion";
      pt->SpeedMeanReversion.Vtype=DOUBLE;
      pt->SpeedMeanReversion.Val.V_DOUBLE=0.2;
      pt->SpeedMeanReversion.Viter=ALLOW;
      model->nvar++;

      pt->LongRunVariance.Vname = "Long-Run Variance";
      pt->LongRunVariance.Vtype=DOUBLE;
      pt->LongRunVariance.Val.V_DOUBLE=1.4;
      pt->LongRunVariance.Viter=ALLOW;
      model->nvar++;

      pt->Sigma2.Vname = "Volatility of Volatility";
      pt->Sigma2.Vtype=DOUBLE;
      pt->Sigma2.Val.V_DOUBLE=0.001;
      pt->Sigma2.Viter=ALLOW;
      model->nvar++;

      pt->RhoFI.Vname = "RhoFI";
      pt->RhoFI.Vtype=DOUBLE;
      pt->RhoFI.Val.V_DOUBLE=0.0;
      pt->RhoFI.Viter=ALLOW;
      model->nvar++;
      
      pt->RhoFV.Vname = "RhoFV";
      pt->RhoFV.Vtype=DOUBLE;
      pt->RhoFV.Val.V_DOUBLE=0.0;
      pt->RhoFV.Viter=ALLOW;
      model->nvar++;

      pt->RhoIV.Vname = "RhoIV";
      pt->RhoIV.Vtype=DOUBLE;
      pt->RhoIV.Val.V_DOUBLE=0.4;
      pt->RhoIV.Viter=ALLOW;
      model->nvar++;

      pt->RhoI.Vname = "RhoI";
      pt->RhoI.Vtype=DOUBLE;
      pt->RhoI.Val.V_DOUBLE=0.8;
      pt->RhoI.Viter=ALLOW;
      model->nvar++;

    }
  return OK;
}
TYPEMOD INFLATION_LMM_HESTON1d;
MAKEMOD(INFLATION_LMM_HESTON1d);


