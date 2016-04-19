#include "uvm1d.h"
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

      pt->sigmamin.Vname = "Sigma Min";
      pt->sigmamin.Vtype=DOUBLE;
      pt->sigmamin.Val.V_DOUBLE=0.2;
      pt->sigmamin.Viter=ALLOW;
      model->nvar++;

      
      pt->sigmamax.Vname = "Sigma Max";
      pt->sigmamax.Vtype=DOUBLE;
      pt->sigmamax.Val.V_DOUBLE=0.3;
      pt->sigmamax.Viter=ALLOW;
      model->nvar++;
      
    }

  return OK;
}

TYPEMOD uvm1d;
MAKEMOD(uvm1d);
