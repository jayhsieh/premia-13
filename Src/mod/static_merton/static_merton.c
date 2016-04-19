#include "static_merton.h"
#include "chk.h"
#include "error_msg.h"
#include "model.h"
#include "premia_obj.h"

static int MOD(Init)(Model *model)       
{
  TYPEMOD* pt=(TYPEMOD*)(model->TypeModel);
	
  if (model->init == 0 )
    {
      model->init = 1;
      model->nvar=0;

      pt->rho.Vname = "Positive Correlation between Borrowers";
      pt->rho.Vtype=RGDOUBLE;
      pt->rho.Val.V_RGDOUBLE=0.15;
      pt->rho.Viter=ALLOW;
      model->nvar++;

    }
  return OK;
}
TYPEMOD Static_Merton1d;
MAKEMOD(Static_Merton1d);
