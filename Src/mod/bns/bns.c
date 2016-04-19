#include "bns.h"
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

      pt->Sigma0.Vname = "Current Variance";
      pt->Sigma0.Vtype=DOUBLE;
      pt->Sigma0.Val.V_DOUBLE=0.0433;
      pt->Sigma0.Viter=ALLOW;
      model->nvar++;

      pt->Lambda.Vname = "Lambda";
      pt->Lambda.Vtype=DOUBLE;
      pt->Lambda.Val.V_DOUBLE=0.5474;
      pt->Lambda.Viter=ALLOW;
      model->nvar++;

      pt->Beta.Vname = "Beta";
      pt->Beta.Vtype=DOUBLE;
      pt->Beta.Val.V_DOUBLE=18.6075;
      pt->Beta.Viter=ALLOW;
      model->nvar++;

      pt->Alpha.Vname = "Alpha";
      pt->Alpha.Vtype=DOUBLE;
      pt->Alpha.Val.V_DOUBLE=0.6069;
      pt->Alpha.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Rho";
      pt->Rho.Vtype=DOUBLE;
      pt->Rho.Val.V_DOUBLE=-4.6750;
      pt->Rho.Viter=ALLOW;
      model->nvar++;
    }

  return OK;
}

TYPEMOD bns;
MAKEMOD(bns);
