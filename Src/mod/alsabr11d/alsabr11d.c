#include "alsabr11d.h"
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
      model->nvar++;;

      pt->z0.Vname = "z0";
      pt->z0.Vtype=SPDOUBLE;
      pt->z0.Val.V_SPDOUBLE=10.;
      pt->z0.Viter=ALLOW;
      model->nvar++;;

      pt->gam.Vname = "Gamma Exponent";
      pt->gam.Vtype=RGDOUBLE02;
      pt->gam.Val.V_RGDOUBLE02=1.90;
      pt->gam.Viter=ALLOW;
      model->nvar++;

      pt->eta.Vname = "Eta";
      pt->eta.Vtype=SPDOUBLE;
      pt->eta.Val.V_SPDOUBLE=0.01;
      pt->eta.Viter=ALLOW;
      model->nvar++;

      pt->Divid.Vname = "Annual Dividend Rate";
      pt->Divid.Vtype=DOUBLE;
      pt->Divid.Val.V_DOUBLE=0.0;
      pt->Divid.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=10.;
      pt->R.Viter=ALLOW;
      model->nvar++;
    }

  return OK;
}

TYPEMOD alsabr11d;
MAKEMOD(alsabr11d);
