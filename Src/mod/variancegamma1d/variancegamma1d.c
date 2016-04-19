#include "variancegamma1d.h"
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
      pt->T.Val.V_DATE=0.;
      pt->T.Viter=ALLOW;
      model->nvar++;

      pt->S0.Vname = "Spot";
      pt->S0.Vtype=PDOUBLE;
      pt->S0.Val.V_PDOUBLE=100.;
      pt->S0.Viter=ALLOW;
      model->nvar++;

      pt->Mu.Vname = "Trend";
      pt->Mu.Vtype=DOUBLE;
      pt->Mu.Val.V_DOUBLE=0.;
      pt->Mu.Viter=ALLOW;
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

      pt->Sigma.Vname = "Sigma";
      pt->Sigma.Vtype=SPDOUBLE;
      pt->Sigma.Val.V_SPDOUBLE=0.12;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->Theta.Vname = "Theta";
      pt->Theta.Vtype=DOUBLE;
      pt->Theta.Val.V_DOUBLE=-0.33;
      pt->Theta.Viter=ALLOW;
      model->nvar++;

      pt->Kappa.Vname = "Kappa";
      pt->Kappa.Vtype=SPDOUBLE;
      pt->Kappa.Val.V_SPDOUBLE=0.16;
      pt->Kappa.Viter=ALLOW;
      model->nvar++;
    }

  return OK;
}


TYPEMOD VarianceGamma1dim;

MAKEMOD(VarianceGamma1dim);


