#include "bs2d.h"
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

      pt->S01.Vname = "Spot 1";
      pt->S01.Vtype=PDOUBLE;
      pt->S01.Val.V_PDOUBLE=100.;
      pt->S01.Viter=ALLOW;
      model->nvar++;

      pt->Mu1.Vname = "Trend 1";
      pt->Mu1.Vtype=DOUBLE;
      pt->Mu1.Val.V_DOUBLE=0.;
      pt->Mu1.Viter=ALLOW;
      model->nvar++;

      pt->Sigma1.Vname = "Volatility 1";
      pt->Sigma1.Vtype=PDOUBLE;
      pt->Sigma1.Val.V_PDOUBLE=0.2;
      pt->Sigma1.Viter=ALLOW;
      model->nvar++;

      pt->Divid1.Vname = "Annual Dividend Rate 1";
      pt->Divid1.Vtype=DOUBLE;
      pt->Divid1.Val.V_DOUBLE=0.;
      pt->Divid1.Viter=ALLOW;
      model->nvar++;

      pt->S02.Vname = "Spot 2";
      pt->S02.Vtype=PDOUBLE;
      pt->S02.Val.V_PDOUBLE=100.;
      pt->S02.Viter=ALLOW;
      model->nvar++;

      pt->Mu2.Vname = "Trend 2";
      pt->Mu2.Vtype=DOUBLE;
      pt->Mu2.Val.V_DOUBLE=0.;
      pt->Mu2.Viter=ALLOW;
      model->nvar++;

      pt->Sigma2.Vname = "Volatility 2";
      pt->Sigma2.Vtype=PDOUBLE;
      pt->Sigma2.Val.V_PDOUBLE=0.2;
      pt->Sigma2.Viter=ALLOW;
      model->nvar++;

      pt->Divid2.Vname = "Annual Dividend Rate 2";
      pt->Divid2.Vtype=DOUBLE;
      pt->Divid2.Val.V_DOUBLE=0.;
      pt->Divid2.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Correlation";
      pt->Rho.Vtype=RGDOUBLEM11;
      pt->Rho.Val.V_RGDOUBLEM11=0.;
      pt->Rho.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=5.0;
      pt->R.Viter=ALLOW;
      model->nvar++;
    }

  return OK;
}

TYPEMOD BlackScholes2dim;

MAKEMOD(BlackScholes2dim);


