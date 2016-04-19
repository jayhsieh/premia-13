#include "qtsm2d.h"
#include "chk.h"
#include "error_msg.h"
#include "model.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

extern char* path_sep;

static int MOD(Init)(Model *model)       
{
  TYPEMOD* pt=(TYPEMOD*)(model->TypeModel);
  double tmp1[3],tmp2[4];
  double tmp[2];
  
  if (model->init == 0 )
    {
      model->init = 1;
      model->nvar=0;
      
      pt->T.Vname = "Current Date";
      pt->T.Vtype=DATE;
      pt->T.Val.V_DATE=0.0;
      pt->T.Viter=ALLOW;
      model->nvar++;

      pt->x.Vname = "Initial x";
      pt->x.Vtype=PNLVECT;
     
      tmp[0]=0.735;
      tmp[1]=-0.525;
      pt->x.Val.V_PNLVECT=pnl_vect_create_from_ptr(2,tmp);
      pt->x.Viter=FORBID;
      model->nvar++;

      pt->d0.Vname = "d0";
      pt->d0.Vtype=PDOUBLE;
      pt->d0.Val.V_PDOUBLE=0.0088;
      pt->d0.Viter=ALLOW;
      model->nvar++;

      pt->d.Vname = "Initial d";
      pt->d.Vtype=PNLVECT;
      tmp[0]=0.0066;
      tmp[1]=-0.022;
      pt->d.Val.V_PNLVECT=pnl_vect_create_from_ptr(2,tmp);
      pt->d.Viter=FORBID;
      model->nvar++;
     

      pt->theta.Vname = "Theta";
      pt->theta.Vtype=PNLVECT;
      tmp[0]=0.;
      tmp[1]=0.;
      pt->theta.Val.V_PNLVECT=pnl_vect_create_from_ptr(2,tmp);
      pt->theta.Viter=FORBID;
      model->nvar++;

      pt->GammaV.Vname = "Gamma11 Gamma12 Gamma22";
      pt->GammaV.Vtype=PNLVECT;
      tmp1[0]=0.0176;
      tmp1[1]=-0.0132;
      tmp1[2]=0.1100;
      pt->GammaV.Val.V_PNLVECT=pnl_vect_create_from_ptr(3,tmp1);
      pt->GammaV.Viter=FORBID;
      model->nvar++;
      
      pt->SigmaV.Vname = "Sigma11 Sigma12 Sigma22";
      pt->SigmaV.Vtype=PNLVECT;
      tmp1[0]=1.;
      tmp1[1]=0.;
      tmp1[2]=1.;
      pt->SigmaV.Val.V_PNLVECT=pnl_vect_create_from_ptr(3,tmp1);
      pt->SigmaV.Viter=FORBID;
      model->nvar++;

     
     
      pt->KappaV.Vname = "Kappa11 Kappa12 Kappa21 Kappa22";
      pt->KappaV.Vtype=PNLVECT;
      tmp2[0]=0.264;
      tmp2[1]=0.;
      tmp2[2]=0.1;
      tmp2[3]=0.66;
      pt->KappaV.Val.V_PNLVECT=pnl_vect_create_from_ptr(4,tmp2);
      pt->KappaV.Viter=FORBID;
      model->nvar++;
     
    }
  return OK;
       
}
TYPEMOD QTSM2d;
MAKEMOD(QTSM2d);


