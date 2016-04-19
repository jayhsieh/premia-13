#include "fps2d.h"
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

      pt->InitialSlow.Vname = "Current Y";
      pt->InitialSlow.Vtype=DOUBLE;
      pt->InitialSlow.Val.V_DOUBLE=-1.;
      pt->InitialSlow.Viter=ALLOW;
      model->nvar++;
      
      pt->InitialFast.Vname = "Current Z";
      pt->InitialFast.Vtype=DOUBLE;
      pt->InitialFast.Val.V_DOUBLE=-1.;
      pt->InitialFast.Viter=ALLOW;
      model->nvar++;

      
      pt->SigmaSlow.Vname = "Sigma Y";
      pt->SigmaSlow.Vtype=DOUBLE;
      pt->SigmaSlow.Val.V_DOUBLE=0.8;
      pt->SigmaSlow.Viter=ALLOW;
      model->nvar++;

      pt->SigmaFast.Vname = "Sigma Z";
      pt->SigmaFast.Vtype=DOUBLE;
      pt->SigmaFast.Val.V_DOUBLE=0.5;
      pt->SigmaFast.Viter=ALLOW;
      model->nvar++;

      pt->MeanReversionSlow.Vname = "Speed of Mean Reversion Slow";
      pt->MeanReversionSlow.Vtype=DOUBLE;
      pt->MeanReversionSlow.Val.V_DOUBLE=0.01;
      pt->MeanReversionSlow.Viter=ALLOW;
      model->nvar++;

      pt->MeanReversionFast.Vname = "Speed of Mean Reversion Fast";
      pt->MeanReversionFast.Vtype=DOUBLE;
      pt->MeanReversionFast.Val.V_DOUBLE=100.;
      pt->MeanReversionFast.Viter=ALLOW;
      model->nvar++;

      pt->LongRunVarianceSlow.Vname = "Long-Run Term Slow";
      pt->LongRunVarianceSlow.Vtype=DOUBLE;
      pt->LongRunVarianceSlow.Val.V_DOUBLE=-0.8;
      pt->LongRunVarianceSlow.Viter=ALLOW;
      model->nvar++;

      pt->LongRunVarianceFast.Vname = "Long-Run Term Fast";
      pt->LongRunVarianceFast.Vtype=DOUBLE;
      pt->LongRunVarianceFast.Val.V_DOUBLE=-0.8;
      pt->LongRunVarianceFast.Viter=ALLOW;
      model->nvar++;

      pt->Rho1.Vname = "Rho 1";
      pt->Rho1.Vtype=RGDOUBLEM11;
      pt->Rho1.Val.V_RGDOUBLEM11=-0.2;
      pt->Rho1.Viter=ALLOW;
      model->nvar++;

      pt->Rho2.Vname = "Rho 2";
      pt->Rho2.Vtype=RGDOUBLEM11;
      pt->Rho2.Val.V_RGDOUBLEM11=-0.2;
      pt->Rho2.Viter=ALLOW;
      model->nvar++;

      pt->Rho12.Vname = "Rho 12";
      pt->Rho12.Vtype=RGDOUBLEM11;
      pt->Rho12.Val.V_RGDOUBLEM11=0.;
      pt->Rho12.Viter=ALLOW;
      model->nvar++;

      model->HelpFilenameHint = "FPS2D";
    }

  return OK;
}


TYPEMOD FPS2dim;

MAKEMOD(FPS2dim);


