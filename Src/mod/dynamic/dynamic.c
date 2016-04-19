#include "dynamic.h"
#include "chk.h"
#include "model.h"
#include "premia_obj.h"
 
/**
 * Initialization of the Dynamic CDO model
 * @param model 
 */
static int MOD(Init)(Model *model)       
{
  TYPEMOD* pt=(TYPEMOD*)(model->TypeModel);
  if (model->init == 0 )
    {
      model->init = 1;
      model->nvar=0;
      pt->Ncomp.Vname = "Number of Companies";
      pt->Ncomp.Vtype=PINT;
      pt->Ncomp.Val.V_PINT=125;
      pt->Ncomp.Viter=ALLOW;
      model->nvar++;

      pt->r.Vname = "Interest rate";
      pt->r.Vtype=PDOUBLE;
      pt->r.Val.V_PDOUBLE=0.04;
      pt->r.Viter=ALLOW;
      model->nvar++;
    }
  return OK;
}

TYPEMOD Dynamic;

MAKEMOD(Dynamic);


