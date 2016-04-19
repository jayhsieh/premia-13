#include "jump1d.h"
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

      pt->Mu.Vname = "Trend";
      pt->Mu.Vtype=DOUBLE;
      pt->Mu.Val.V_DOUBLE=0.;
      pt->Mu.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Volatility";
      pt->Sigma.Vtype=PDOUBLE;
      pt->Sigma.Val.V_PDOUBLE=0.3;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->Divid.Vname = "Annual Dividend Rate";
      pt->Divid.Vtype=DOUBLE;
      pt->Divid.Val.V_DOUBLE=0.;
      pt->Divid.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=5.;
      pt->R.Viter=ALLOW;
      model->nvar++;
      
      pt->Lambda.Vname = "Lambda";
      pt->Lambda.Vtype=DOUBLE;
      pt->Lambda.Val.V_DOUBLE=1;
      pt->Lambda.Viter=ALLOW;
      model->nvar++;

      pt->Mean.Vname = "Jump Size";
      pt->Mean.Vtype=DOUBLE;
      pt->Mean.Val.V_DOUBLE=0.1;
      pt->Mean.Viter=ALLOW;
      model->nvar++;

    }

  return OK;
}


TYPEMOD Jump1dim;

MAKEMOD(Jump1dim);


