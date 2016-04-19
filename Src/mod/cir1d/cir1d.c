#include "cir1d.h"
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
      pt->T.Val.V_DATE=0.0;
      pt->T.Viter=ALLOW;
      model->nvar++;

      pt->r0.Vname = "Current Rate";
      pt->r0.Vtype=PDOUBLE;
      pt->r0.Val.V_PDOUBLE=0.05;
      pt->r0.Viter=ALLOW;
      model->nvar++;

      pt->k.Vname = "Speed of Mean Reversion";
      pt->k.Vtype=PDOUBLE;
      pt->k.Val.V_PDOUBLE=0.15;
      pt->k.Viter=ALLOW;
      model->nvar++;

      pt->theta.Vname = "Long Term Mean";
      pt->theta.Vtype=PDOUBLE;
      pt->theta.Val.V_PDOUBLE=0.05;
      pt->theta.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Volatility";
      pt->Sigma.Vtype=PDOUBLE;
      pt->Sigma.Val.V_PDOUBLE=0.1;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;
    }
  return OK;
}
TYPEMOD Cir1d;
MAKEMOD(Cir1d);


