#include "rskou1d.h"
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

      pt->Transition_probabilities.Vname = "Transition probabilities";
      pt->Transition_probabilities.Vtype=FILENAME;
      pt->Transition_probabilities.Val.V_INT=0;
      pt->Transition_probabilities.Viter=FORBID;
      pt->Transition_probabilities.Vsetable = SETABLE;

      model->nvar++;
      if (( pt->Transition_probabilities.Val.V_FILENAME=malloc(sizeof(char)*MAX_PATH_LEN))==NULL)
        return MEMORY_ALLOCATION_FAILURE;
      sprintf( pt->Transition_probabilities.Val.V_FILENAME, "%s%sKOU_Transition_probabilities.dat", premia_data_dir,path_sep);
      
    }

  return OK;
}

TYPEMOD RSKou1dim;

MAKEMOD(RSKou1dim);


