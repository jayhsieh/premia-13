#include "nsp/interf.h"
#define __PARAMS(x)     x
#include "../demo/price_cdo.h"

int nsp_premia_cdo(Stack stack, int rhs, int opt, int lhs)
{
  NspMatrix       *ncomp = GetRealMatCopy(stack, 1);
  NspMatrix       *nominal = GetRealMatCopy(stack, 2);
  NspMatrix       *dates = GetRealMatCopy(stack, 3);
  NspMatrix       *tranches = GetRealMatCopy(stack, 4);
  NspMatrix       *intensity = GetRealMatCopy(stack, 5);
  NspMatrix       *xrates = GetRealMatCopy(stack, 6);
  NspMatrix       *yrates = GetRealMatCopy(stack, 7);
  NspMatrix       *t_recovery = GetRealMatCopy(stack, 8);
  NspMatrix       *recovery = GetRealMatCopy(stack, 9);
  NspMatrix       *t_copula = GetRealMatCopy(stack, 10);
  NspMatrix       *copula = GetRealMatCopy(stack, 11);
  NspMatrix       *t_method = GetRealMatCopy(stack, 12);
  NspMatrix       *method = GetRealMatCopy(stack, 13);
  NspMatrix       *price;
  NspMatrix       *dl;
  NspMatrix       *pl;
  int         m, n;
    
  CheckRhs(13, 13);
  CheckLhs(1, 3);
  m = Max(tranches->m, tranches->n)-1;
  n = 1;
  if ((price = nsp_matrix_create(NVOID,'r',m,1)) == NULLMAT) return RET_BUG;
  if ((dl = nsp_matrix_create(NVOID,'r',m,1)) == NULLMAT) return RET_BUG;
  if ((pl = nsp_matrix_create(NVOID,'r',m,1)) == NULLMAT) return RET_BUG;
  ncomp = Mat2int(ncomp);
  t_recovery = Mat2int(t_recovery);
  t_copula = Mat2int(t_copula);
  t_method = Mat2int(t_method);
  method = Mat2int(method);
  price_cdo(ncomp->I, 
            nominal->R, 
            Max(dates->m, dates->n), 
            dates->R, 
            Max(tranches->m, tranches->n), 
            tranches->R, 
            intensity->R, 
            Max(xrates->m, xrates->n), 
            xrates->R, 
            yrates->R, 
            t_recovery->I, 
            recovery->R, 
            t_copula->I, 
            copula->R, 
            t_method->I, 
            method->I, 
            price->R, 
            dl->R, 
            pl->R);
  NthObj(rhs+1) = NSP_OBJECT(price);
  NthObj(rhs+2) = NSP_OBJECT(dl);
  NthObj(rhs+3) = NSP_OBJECT(pl);

  NSP_OBJECT(dl)->ret_pos = 3;
  NSP_OBJECT(pl)->ret_pos = 2;
  NSP_OBJECT(price)->ret_pos = 1;

  return Max(lhs, 1);
}
