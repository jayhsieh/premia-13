#include "cirpp1d.h"
#include "chk.h"
#include "model.h"
#include "premia_obj.h"
#include "enums.h"


double MOD(GetYield)(TYPEMOD *pt)
{
  VAR *Par;
  Par = lookup_premia_enum_par (&(pt->flat_flag), 0);
  return Par[0].Val.V_PDOUBLE;
}

char* MOD(GetCurve)(TYPEMOD *pt)
{
  VAR *Par;
  Par = lookup_premia_enum_par (&(pt->flat_flag), 1);
  return Par[0].Val.V_FILENAME;
}

static int MOD(Init)(Model *model)       
{
  VAR *Par;
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

      pt->flat_flag.Vname = "Initial Yields Curve";
      pt->flat_flag.Vtype=ENUM;
      pt->flat_flag.Val.V_ENUM.value=0;
      pt->flat_flag.Val.V_ENUM.members=&PremiaEnumFlat;
      pt->flat_flag.Viter=ALLOW;
      model->nvar++;
      Par = lookup_premia_enum_par (&(pt->flat_flag), 0);
      Par[0].Vname = "Yield Value";
      Par[0].Vtype=PDOUBLE;
      Par[0].Val.V_PDOUBLE=0.05;
      Par[0].Viter=ALLOW;
      Par = lookup_premia_enum_par (&(pt->flat_flag), 1);
      Par[0].Vname = "Yield Curve";
      Par[0].Vtype=FILENAME;
      Par[0].Val.V_FILENAME=NULL;
      Par[0].Viter=FORBID;


      pt->a.Vname = "Speed of Mean Reversion";
      pt->a.Vtype=PDOUBLE;
      pt->a.Val.V_PDOUBLE=0.15;
      pt->a.Viter=ALLOW;
      model->nvar++;
      
      pt->b.Vname = "Long Term Mean";
      pt->b.Vtype=PDOUBLE;
      pt->b.Val.V_PDOUBLE=0.05;
      pt->b.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Volatility";
      pt->Sigma.Vtype=PDOUBLE;
      pt->Sigma.Val.V_PDOUBLE=0.1;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      model->HelpFilenameHint = "CIRPP1D";
    }

  Par = lookup_premia_enum_par (&(pt->flat_flag), 1);
  if (Par[0].Val.V_FILENAME==NULL)
    {
      if ((Par[0].Val.V_FILENAME=malloc(sizeof(char)*MAX_PATH_LEN))==NULL)
        return MEMORY_ALLOCATION_FAILURE;
      sprintf(Par[0].Val.V_FILENAME, "%s%sinitialyield.dat", premia_data_dir, path_sep);
    }

  return OK;
}
TYPEMOD CirPlus1d;
MAKEMOD(CirPlus1d);


