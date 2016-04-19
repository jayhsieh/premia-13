#include "locvolhw1d.h"
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

      pt->csi.Vname = "csi";
      pt->csi.Vtype=PDOUBLE;
      pt->csi.Val.V_PDOUBLE=0.007;
      pt->csi.Viter=ALLOW;
      model->nvar++;

      pt->kappa.Vname = "kappa";
      pt->kappa.Vtype=PDOUBLE;
      pt->kappa.Val.V_PDOUBLE=0.01;
      pt->kappa.Viter=ALLOW;
      model->nvar++;
     
      pt->v.Vname = "v";
      pt->v.Vtype=PDOUBLE;
      pt->v.Val.V_PDOUBLE=0.2;
      pt->v.Viter=ALLOW;
      model->nvar++;

      pt->beta.Vname = "beta";
      pt->beta.Vtype=PDOUBLE;
      pt->beta.Val.V_PDOUBLE=0.8;
      pt->beta.Viter=ALLOW;
      model->nvar++;

      pt->rho.Vname = "rho";
      pt->rho.Vtype=RGDOUBLEM11;
      pt->rho.Val.V_RGDOUBLEM11=0.15;
      pt->rho.Viter=ALLOW;
      model->nvar++;

      pt->f0t.Vname = "Forward rate";
      pt->f0t.Vtype=PDOUBLE;
      pt->f0t.Val.V_PDOUBLE=0.02;
      pt->f0t.Viter=ALLOW;
      model->nvar++;

    }

  return OK;
}

TYPEMOD locvolhw1d;
MAKEMOD(locvolhw1d);
