#include "jarrowyildirim1d.h"
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

      pt->I0.Vname = "Current CPI";
      pt->I0.Vtype=PDOUBLE;
      pt->I0.Val.V_PDOUBLE=0.04;
      pt->I0.Viter=ALLOW;
      model->nvar++;

      
      pt->an.Vname = "Volatility Parameter a Nominal";
      pt->an.Vtype=PDOUBLE;
      pt->an.Val.V_PDOUBLE=0.07676;
      pt->an.Viter=ALLOW;
      model->nvar++;
      
      pt->ar.Vname = "Volatility Parameter a Creal";
      pt->ar.Vtype=PDOUBLE;
      pt->ar.Val.V_PDOUBLE=1.057;
      pt->ar.Viter=ALLOW;
      model->nvar++;

       
      pt->sigman.Vname = "Volatility Parameter sigma Nominal";
      pt->sigman.Vtype=PDOUBLE;
      pt->sigman.Val.V_PDOUBLE=0.007;
      pt->sigman.Viter=ALLOW;
      model->nvar++;
      
      pt->sigmar.Vname = "Volatility Parameter sigma Creal";
      pt->sigmar.Vtype=PDOUBLE;
      pt->sigmar.Val.V_PDOUBLE=0.0111;
      pt->sigmar.Viter=ALLOW;
      model->nvar++;

      pt->sigma_cpi.Vname = "Volatility of CPI";
      pt->sigma_cpi.Vtype=PDOUBLE;
      pt->sigma_cpi.Val.V_PDOUBLE=0.038;
      pt->sigma_cpi.Viter=ALLOW;
      model->nvar++;
  
      pt->Rhonr.Vname = "Rho Nominal Creal";
      pt->Rhonr.Vtype=DOUBLE;
      pt->Rhonr.Val.V_DOUBLE=-0.11;
      pt->Rhonr.Viter=ALLOW;
      model->nvar++;

      
      pt->Rhoncpi.Vname = "Rho Nominal CPI";
      pt->Rhoncpi.Vtype=DOUBLE;
      pt->Rhoncpi.Val.V_DOUBLE=0.;
      pt->Rhoncpi.Viter=ALLOW;
      model->nvar++;

      pt->Rhorcpi.Vname = "Rho Creal CPI";
      pt->Rhorcpi.Vtype=DOUBLE;
      pt->Rhorcpi.Val.V_DOUBLE=0.0932;
      pt->Rhorcpi.Viter=ALLOW;
      model->nvar++;
      
    }
  return OK;
}

TYPEMOD JarrowYildirim1d;

MAKEMOD(JarrowYildirim1d);


