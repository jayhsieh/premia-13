#include "local_vol.h"
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

      pt->S0.Vname = "Spot";
      pt->S0.Vtype=PDOUBLE;
      pt->S0.Val.V_PDOUBLE=100.55;
      pt->S0.Viter=ALLOW;
      model->nvar++;

      pt->Interest.Vname = "Annual Interest Rate";
      pt->Interest.Vtype=PDOUBLE;
      pt->Interest.Val.V_PDOUBLE=5;
      pt->Interest.Viter=ALLOW;
      model->nvar++;

      pt->Eta.Vname = "Fractional Loss";
      pt->Eta.Vtype=PDOUBLE;
      pt->Eta.Val.V_PDOUBLE=1;
      pt->Eta.Viter=ALLOW;
      model->nvar++;

      pt->Divid.Vname = "Annual Dividend Rate";
      pt->Divid.Vtype=PDOUBLE;
      pt->Divid.Val.V_PDOUBLE=0.;
      pt->Divid.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Volatility";
      pt->Sigma.Vtype=PDOUBLE;
      pt->Sigma.Val.V_PDOUBLE=0.2;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;
    }

  return OK;
}

TYPEMOD local_vol;
MAKEMOD(local_vol);
