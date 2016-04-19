#include "nonpar1d.h"
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

      pt->S0.Vname = "Spot";
      pt->S0.Vtype=PDOUBLE;
      pt->S0.Val.V_PDOUBLE=100.;
      pt->S0.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=10.;
      pt->R.Viter=ALLOW;
      model->nvar++;
      
      pt->implied_volatility.Vname = "Implied Volatility";
      pt->implied_volatility.Vtype=FILENAME;
      pt->implied_volatility.Val.V_INT=0;
      pt->implied_volatility.Viter=FORBID;
      pt->implied_volatility.Vsetable = UNSETABLE;

      model->nvar++;
      if (( pt->implied_volatility.Val.V_FILENAME=malloc(sizeof(char)*MAX_PATH_LEN))==NULL)
        return MEMORY_ALLOCATION_FAILURE;
      sprintf( pt->implied_volatility.Val.V_FILENAME, "%s%simplied_volatility.dat", premia_data_dir,path_sep);

    }
  return OK;
}
TYPEMOD NonPar1d;
MAKEMOD(NonPar1d);
