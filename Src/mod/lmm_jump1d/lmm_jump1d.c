#include "lmm_jump1d.h"
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
      pt->T.Val.V_DATE=0.0;
      pt->T.Viter=ALLOW;
      model->nvar++;

      pt->l0.Vname = "Flat Initial Libor Rates";
      pt->l0.Vtype=PDOUBLE;
      pt->l0.Val.V_PDOUBLE=0.05;
      pt->l0.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Flat Volatility Libor Rates";
      pt->Sigma.Vtype=PDOUBLE;
      pt->Sigma.Val.V_PDOUBLE=0.2;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

    }
  return OK;
}
TYPEMOD LMM_JUMP1d;
MAKEMOD(LMM_JUMP1d);


