#include "fps1d.h"
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

      pt->Sigma0.Vname = "Current Y";
      pt->Sigma0.Vtype=DOUBLE;
      pt->Sigma0.Val.V_DOUBLE=0.0;
      pt->Sigma0.Viter=ALLOW;
      model->nvar++;

      pt->MeanReversion.Vname = "Speed of Mean Reversion";
      pt->MeanReversion.Vtype=DOUBLE;
      pt->MeanReversion.Val.V_DOUBLE=1.;
      pt->MeanReversion.Viter=ALLOW;
      model->nvar++;

      pt->LongRunVariance.Vname = "Long-Run Standard Deviation";
      pt->LongRunVariance.Vtype=DOUBLE;
      pt->LongRunVariance.Val.V_DOUBLE=0.5;
      pt->LongRunVariance.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Rho";
      pt->Rho.Vtype=DOUBLE;
      pt->Rho.Val.V_DOUBLE=-0.15;
      pt->Rho.Viter=ALLOW;
      model->nvar++;

      pt->SigmaF.Vname = "Sigma F";
      pt->SigmaF.Vtype=DOUBLE;
      pt->SigmaF.Val.V_DOUBLE=0.2;
      pt->SigmaF.Viter=ALLOW;
      model->nvar++;

     model->HelpFilenameHint = "FPS1D";

    }

  return OK;
}


TYPEMOD FPS1dim;

MAKEMOD(FPS1dim);


