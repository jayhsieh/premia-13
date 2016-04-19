#include        "cdo.h"

double          **hw_numdef(const CDO       *cdo, 
                            const copula    *cop,
                            const grid      *t,
                            const cond_prob *cp)
{
  double      **nd;
  double      *U;
  double      *V;
  double      p0;
  double      w_jn;
  int         jv;
  int         jV;
  int         jt;
  int         jn;
  int         jk;
    
  nd = malloc(t->size * sizeof(double*));
  U = malloc((cdo->n_comp+1) * sizeof(double));
  V = malloc((cdo->n_comp+1) * sizeof(double));
  for (jt = 0; jt < t->size; jt++) {
    nd[jt] = malloc((cdo->n_comp+1) * sizeof(double));
    for (jV = 0; jV < (cdo->n_comp+1); jV++) 
      nd[jt][jV] = 0.;
    for (jv = 0; jv < cop->size; jv++) {
      p0 = 1.;
      for (jn = 0; jn < cdo->n_comp; jn++)
        p0 = p0 * (1. - cp->p[jn][jt][jv]);
      nd[jt][0] += p0 * cop->weights[jv];
      U[0] = 1.;
      for (jV = 1; jV < (cdo->n_comp+1); jV++) {
        V[jV] = 0;
        for (jn = 0; jn < cdo->n_comp; jn++) {
          w_jn = cp->p[jn][jt][jv] / (1. - cp->p[jn][jt][jv]);
          V[jV] += pow(w_jn, jV); /* explosion! */
        }
        U[jV] = 0;
        for (jk = 1; jk <= jV; jk++) {
          U[jV] += (pow(-1., jk+1) * V[jk] * U[jV-jk]);
        }
        U[jV] = U[jV] / jV;
        nd[jt][jV] += p0 * U[jV] * cop->weights[jv];
      }
    }
  }
  free(U);
  free(V);

  return (nd);
}


double          **hw_numdef1(const CDO      *cdo, 
                             const copula    *cop,
                             const grid      *t,
                             const cond_prob *cp)
{
  double      ***nd;
  double              **nd1;
  int                  jw;
  double      *U;
  double      *V;
  double      p0;
  double      w_jn;
  int         jv;
  int         jV;
  int         jt;
  int         jn;
  int         jk;
    
  nd1 = malloc(t->size * sizeof(double*));
  nd = malloc(t->size * sizeof(double**));
  U = malloc((cdo->n_comp+1) * sizeof(double));
  V = malloc((cdo->n_comp+1) * sizeof(double));
  for (jt = 0; jt < t->size; jt++) {
    nd[jt] = malloc((cdo->n_comp+1) * sizeof(double*));
    nd1[jt] = malloc((cdo->n_comp+1) * sizeof(double));
    for (jV = 0; jV < (cdo->n_comp+1); jV++) {
      nd[jt][jV] = malloc((cop->size) * sizeof(double));  
      nd1[jt][jV] = 0.;
    }

    for (jv = 0; jv < cop->size; jv++) {
      for(jV=0;jV<cdo->n_comp+1;jV++){
        nd[jt][jV][jv]=0.0;
      }
          
      for(jw=0;jw<cop->size;jw++){
        p0 = 1.;
        for (jn = 0; jn < cdo->n_comp; jn++){
          p0 = p0 * (1. - cp->p[jn][jt][jv+jw*cop->size]);
        }
            
        nd[jt][0][jv] += p0 * cop->weights[jw+cop->size];
        U[0] = 1.;
        for (jV = 1; jV < (cdo->n_comp+1); jV++) {
          V[jV] = 0;
          for (jn = 0; jn < cdo->n_comp; jn++) {
            w_jn = cp->p[jn][jt][jv+jw*cop->size] / (1. - cp->p[jn][jt][jv+jw*cop->size]);
            V[jV] += pow(w_jn, jV);
          }
          U[jV] = 0;
          for (jk = 1; jk <= jV; jk++) {
            U[jV] += (pow(-1., jk+1) * V[jk] * U[jV-jk]);
          }
          U[jV] = U[jV] / jV;
          nd[jt][jV][jv] += p0 * U[jV] * cop->weights[jw+cop->size];
        }
      }
      for (jV = 1; jV < (cdo->n_comp+1); jV++) {
        nd1[jt][jV]=nd1[jt][jV]+nd[jt][jV][jv]*(cop->weights[jv]);
      }    

    }
  }
  free(U);
  free(V);
  /*    free(nd);*/
  return (nd1);
}


double          **hw_losses_h(const CDO     *cdo,
                              const copula  *cop, 
                              const grid    *t,
                              const grid    *x,  
                              const cond_prob   *cp)
{
  double      **cond_losses;
  double      **losses;
  double      *delta;
  double      p_default;
  double      sum;
  int         jt;
  int         jx;
  int         jv;
  int         jn;
    
  cond_losses = malloc(x->size * sizeof(double*));
  for (jx = 0; jx < x->size; jx++) 
    cond_losses[jx] = malloc(cop->size * sizeof(double));
  delta = malloc(x->size * sizeof(double));
  losses = malloc(t->size * sizeof(double*));
  for (jt = 0; jt < t->size; jt++) {
    for (jv = 0; jv < cop->size; jv++)
      cond_losses[0][jv] = 1.;
    for (jx = 1; jx < x->size; jx++) {
      for (jv = 0; jv < cop->size; jv++)
        cond_losses[jx][jv] = 0.;
    }
    for (jv = 0; jv < cop->size; jv++) {
      for (jn = 0; jn < cdo->n_comp; jn++) {
        p_default = cp->p[jn][jt][jv];
        sum = 0;
        for (jx = 1; jx < x->size; jx++) {
          delta[jx] = p_default * (cond_losses[jx-1][jv] 
                                   - cond_losses[jx][jv]);
          sum += delta[jx];
        }
        cond_losses[0][jv] -= sum;
        for (jx = 1; jx < x->size; jx++) {
          cond_losses[jx][jv] += delta[jx];
        }
      }
    }
    losses[jt] = malloc(x->size * sizeof(double));
    for (jx = 0; jx < x->size; jx++) {
      losses[jt][jx] = 0;
      for (jv = 0; jv < cop->size; jv++) {
        losses[jt][jx] += cond_losses[jx][jv] * cop->weights[jv];   
      }
    }
  }
  for (jx = 0; jx < x->size; jx++) 
    free(cond_losses[jx]);
  free(cond_losses);
  free(delta);

  return (losses);
}


double          **hw_losses_h1(const CDO    *cdo,
                               const copula  *cop, 
                               const grid    *t,
                               const grid    *x,  
                               const cond_prob   *cp)
{
  double      ***cond_losses;
  double      **losses;
  double              ***losses1;
  double      *delta;
  double      p_default;
  double      sum;
  int         jt;
  int                 jw;
  int         jx;
  int         jv;
  int         jn;
    
  cond_losses = malloc(x->size * sizeof(double**));
  for (jx = 0; jx < x->size; jx++) 
    cond_losses[jx] = malloc(cop->size * sizeof(double*));
  delta = malloc(x->size * sizeof(double));
  losses = malloc(t->size * sizeof(double*));
  losses1 = malloc(t->size * sizeof(double**));

  for(jx=0;jx< x->size;jx++){
    for(jv=0;jv<cop->size;jv++){
      cond_losses[jx][jv] = malloc(cop->size * sizeof(double));
    }
  }
  for(jt=0;jt<t->size;jt++){
    losses1[jt] = malloc(x->size * sizeof(double*));
    losses[jt] = malloc(x->size * sizeof(double)); 
  }
  for(jt=0;jt<t->size;jt++){
    for(jx=0;jx<x->size;jx++){
      losses1[jt][jx] = malloc(cop->size * sizeof(double));
      losses[jt][jx]=0;
    }
  }  
    
  for(jt=0;jt<t->size;jt++){
    for(jx=0;jx<x->size;jx++){
      for(jv=0;jv<cop->size;jv++){
        losses1[jt][jx][jv] = 0.0;
     
      }
    }  
  }  
     
  for (jt = 0; jt < t->size; jt++) {
    for (jv = 0; jv < cop->size; jv++) {
      for(jw=0;jw<cop->size;jw++){   
        cond_losses[0][jv][jw] = 1.;
        
      }
      for (jx = 1; jx < x->size; jx++) {
        for(jw=0;jw<cop->size;jw++){
          cond_losses[jx][jv][jw] = 0.;
        }
      }  
       
      for(jw=0;jw<cop->size;jw++){  
            
        for (jn = 0; jn < cdo->n_comp; jn++) {
          p_default = cp->p[jn][jt][jv+jw*cop->size];
          sum = 0;
               
          for (jx = 1; jx < x->size; jx++) {
            delta[jx] = p_default * (cond_losses[jx-1][jv][jw]- cond_losses[jx][jv][jw]);
            sum = sum+delta[jx];
          }

          cond_losses[0][jv][jw] =cond_losses[0][jv][jw]- sum;
        
          for (jx = 1; jx < x->size; jx++) {
            cond_losses[jx][jv][jw] =cond_losses[jx][jv][jw]+ delta[jx];
          }
        }
      }
     
       
      for (jx = 0; jx < x->size; jx++) {
        for(jw=0;jw<cop->size;jw++){
          losses1[jt][jx][jv] =losses1[jt][jx][jv]+ cond_losses[jx][jv][jw] * cop->weights[jw+cop->size];   
        }
        losses[jt][jx]=losses[jt][jx]+losses1[jt][jx][jv]*cop->weights[jv]; 
      }
          
       
    }
  }
 
    
       
 
    


  return (losses);
}




double          **hw_losses_nh(const CDO        *cdo,
                               const copula        *cop, 
                               const grid      *t,
                               const grid      *x,  
                               const cond_prob     *cp)
{
  double      **cond_losses;
  double      **mean_cond_losses;
  double      **losses;
  double      *add_cond;
  double      *add_mean;
  double      L_j;
  double      A_kpL_j;
  double      p_default;
  int         jt;
  int         jx;
  int         ujx;
  int         jv;
  int         jn;
    
  cond_losses = malloc(x->size * sizeof(double*));
  mean_cond_losses = malloc(x->size * sizeof(double*));
  for (jx = 0; jx < x->size; jx++) {
    cond_losses[jx] = malloc(cop->size * sizeof(double));
    mean_cond_losses[jx] = malloc(cop->size * sizeof(double));
  }
  losses = malloc(t->size * sizeof(double*));
  add_cond = malloc(x->size * sizeof(double));
  add_mean = malloc(x->size * sizeof(double));
  for (jt = 0; jt < t->size; jt++) {
    for (jv = 0; jv < cop->size; jv++) {
      cond_losses[0][jv] = 1.;
      mean_cond_losses[0][jv] = 0.;
    }
    for (jx = 1; jx < x->size; jx++) {
      for (jv = 0; jv < cop->size; jv++) {
        cond_losses[jx][jv] = 0.;
        mean_cond_losses[jx][jv] = 0.;
      }
    }
    for (jn = 0; jn < cdo->n_comp; jn++) {
      L_j = cdo->C[jn]->nominal * (1 - RECOVERY(jn));
      for (jv = 0; jv < cop->size; jv++) {
        p_default = cp->p[jn][jt][jv];
        for (jx = 0; jx < x->size; jx++) {
          add_cond[jx] = 0.;
          add_mean[jx] = 0.;
        }
        for (jx = 0; jx < x->size; jx++) {
          A_kpL_j = mean_cond_losses[jx][jv] + L_j;
          ujx = jx;
          while ((ujx+1 < x->size) && (A_kpL_j >= x->data[ujx+1])) ujx++;
          if (ujx > jx) {
            add_cond[jx] -= cond_losses[jx][jv] * p_default;
            add_cond[ujx] += cond_losses[jx][jv] * p_default;
            if (cond_losses[ujx][jv] + cond_losses[jx][jv] * p_default == 0) 
              add_mean[ujx] = 0; 
            else
              add_mean[ujx] += (cond_losses[jx][jv] * p_default * (A_kpL_j - mean_cond_losses[ujx][jv])) / (cond_losses[ujx][jv] + cond_losses[jx][jv] * p_default);  
          }
          else {
            add_mean[jx] += p_default * L_j;
          }
        }
        for (jx = 0; jx < x->size; jx++) {
          cond_losses[jx][jv] += add_cond[jx];
          mean_cond_losses[jx][jv] += add_mean[jx];
        }
      }
    }
    losses[jt] = malloc(x->size * sizeof(double));
    for (jx = 0; jx < x->size; jx++) {
      losses[jt][jx] = 0;
      for (jv = 0; jv < cop->size; jv++) {
        losses[jt][jx] += cond_losses[jx][jv] * cop->weights[jv];
      }
    }
  }
  free(add_cond);
  free(add_mean);

  return (losses);
}


double          **hw_losses_nh1(const CDO       *cdo,
                                const copula        *cop, 
                                const grid      *t,
                                const grid      *x,  
                                const cond_prob     *cp)
{
  double      ***cond_losses;
  double      ***mean_cond_losses;
  double      **losses;
  double      *add_cond;
  double      ***losses1;
  double      *add_mean;
  double      L_j;
  double      A_kpL_j;
  double      p_default;
  int         jt;
  int         jx;
  int         ujx;
  int         jv;
  int         jn;
  int                 jw;
  losses1 = malloc(t->size * sizeof(double**)); 
  cond_losses = malloc(x->size * sizeof(double**));
  mean_cond_losses = malloc(x->size * sizeof(double**));
      
  for(jt=0;jt<t->size;jt++){
    losses1[jt] = malloc(x->size * sizeof(double*));
  }
  for(jt=0;jt<t->size;jt++){
    for(jx=0;jx<x->size;jx++){
      losses1[jt][jx] = malloc(cop->size * sizeof(double));
    }
  }
    
  for (jx = 0; jx < x->size; jx++) {
        
    cond_losses[jx] = malloc(cop->size * sizeof(double*));
    mean_cond_losses[jx] = malloc(cop->size * sizeof(double*));
  }

  for (jx = 0; jx < x->size; jx++) {
    for(jv=0;jv<cop->size;jv++){
      cond_losses[jx][jv] = malloc(cop->size * sizeof(double));
      mean_cond_losses[jx][jv] = malloc(cop->size * sizeof(double));
    }
  }
  losses = malloc(t->size * sizeof(double*));
  add_cond = malloc(x->size * sizeof(double));
  add_mean = malloc(x->size * sizeof(double));
  for (jt = 0; jt < t->size; jt++) {
    for (jv = 0; jv < cop->size; jv++) {
      for (jw = 0; jw < cop->size; jw++) {
        cond_losses[0][jv][jw] = 1.;
        mean_cond_losses[0][jv][jw] = 0.;
      }
    }

    for (jx = 1; jx < x->size; jx++) {
      for (jv = 0; jv < cop->size; jv++) {
        for (jw = 0; jw < cop->size; jw++) {
          cond_losses[jx][jv][jw] = 0.;
          mean_cond_losses[jx][jv][jw] = 0.;
        }
      }
    }
    for (jn = 0; jn < cdo->n_comp; jn++) {
      L_j = cdo->C[jn]->nominal * (1 - RECOVERY(jn));
      for (jv = 0; jv < cop->size; jv++) {
        for (jw = 0; jw < cop->size; jw++) {
          p_default = cp->p[jn][jt][jv+jw*cop->size];
          for (jx = 0; jx < x->size; jx++) {
            add_cond[jx] = 0.;
            add_mean[jx] = 0.;
          }
          for (jx = 0; jx < x->size; jx++) {
            A_kpL_j = mean_cond_losses[jx][jv][jw] + L_j;
            ujx = jx;
            while ((ujx+1 < x->size) && (A_kpL_j >= x->data[ujx+1])) ujx++;
            if (ujx > jx) {
              add_cond[jx] -= cond_losses[jx][jv][jw] * p_default;
              add_cond[ujx] += cond_losses[jx][jv][jw] * p_default;
              if (cond_losses[ujx][jv][jw] + cond_losses[jx][jv][jw] * p_default == 0) 
                add_mean[ujx] = 0; 
              else
                add_mean[ujx] += (cond_losses[jx][jv][jw] * p_default * (A_kpL_j - mean_cond_losses[ujx][jv][jw])) / (cond_losses[ujx][jv][jw] + cond_losses[jx][jv][jw] * p_default);  
            }
            else {
              add_mean[jx] += p_default * L_j;
            }
          }
          for (jx = 0; jx < x->size; jx++) {
            cond_losses[jx][jv][jw] += add_cond[jx];
            mean_cond_losses[jx][jv][jw] += add_mean[jx];
          }
        }
      }
    }
    losses[jt] = malloc(x->size * sizeof(double));
    
    for (jx = 0; jx < x->size; jx++) {
      for(jv=0;jv<cop->size;jv++){
        losses1[jt][jx][jv] = 0;
              
        for(jw=0;jw<cop->size;jw++){
          losses1[jt][jx][jv] =losses1[jt][jx][jv]+ cond_losses[jx][jv][jw] * cop->weights[jw+cop->size];
        }
      }
    }
    for (jx = 0; jx < x->size; jx++) {
      losses[jt][jx]=0;
      for(jv=0;jv<cop->size;jv++){
        
        losses[jt][jx]=losses[jt][jx]+losses1[jt][jx][jv]*cop->weights[jv];
      }
    }
  }
  free(add_cond);
  free(add_mean);

  return (losses);
}




double          **hw_losses_nh2(const CDO       *cdo,
                                const copula        *cop, 
                                const grid      *t,
                                const grid      *x,  
                                const cond_prob     *cp)
{
  double      **cond_losses;
  double      **mean_cond_losses;
  double      **losses;
  double      *add_cond;
  double      *add_mean;
  double      L_j;
  double      A_kpL_j;
  double      p_default;
  int         jt;
  int         jx;
  int         ujx;
  int         jv;
  int         jn;
    
  cond_losses = malloc(x->size * sizeof(double*));
  mean_cond_losses = malloc(x->size * sizeof(double*));
  for (jx = 0; jx < x->size; jx++) {
    cond_losses[jx] = malloc(cop->size * sizeof(double));
    mean_cond_losses[jx] = malloc(cop->size * sizeof(double));
  }
  losses = malloc(t->size * sizeof(double*));
  add_cond = malloc(x->size * sizeof(double));
  add_mean = malloc(x->size * sizeof(double));
  for (jt = 0; jt < t->size; jt++) {
    for (jv = 0; jv < cop->size; jv++) {
      cond_losses[0][jv] = 1.;
      mean_cond_losses[0][jv] = 0.;
    }
    for (jx = 1; jx < x->size; jx++) {
      for (jv = 0; jv < cop->size; jv++) {
        cond_losses[jx][jv] = 0.;
        mean_cond_losses[jx][jv] = 0.;
      }
    }
    for (jn = 0; jn < cdo->n_comp; jn++) {
      L_j = cdo->C[jn]->nominal * (1 - RECOVERY(jn));
      for (jv = 0; jv < cop->size; jv++) {
        p_default = cp->p[jn][jt][jv];
        for (jx = 0; jx < x->size; jx++) {
          add_cond[jx] = 0.;
          add_mean[jx] = 0.;
        }
        for (jx = 0; jx < x->size; jx++) {
          A_kpL_j = mean_cond_losses[jx][jv] + L_j;
          ujx = jx;
          while ((ujx+1 < x->size) && (A_kpL_j >= x->data[ujx+1])) ujx++;
          add_cond[jx] -= cond_losses[jx][jv] * p_default;
          add_cond[ujx] += cond_losses[jx][jv] * p_default;
          add_mean[ujx] += add_mean[ujx] + (mean_cond_losses[jx][jv] + L_j - mean_cond_losses[ujx][jv]) * p_default;
        }
        for (jx = 0; jx < x->size; jx++) {
          cond_losses[jx][jv] += add_cond[jx];
          mean_cond_losses[jx][jv] += add_mean[jx];
        }
      }
    }
    losses[jt] = malloc(x->size * sizeof(double));
    for (jx = 0; jx < x->size; jx++) {
      losses[jt][jx] = 0;
      for (jv = 0; jv < cop->size; jv++) {
        losses[jt][jx] += cond_losses[jx][jv] * cop->weights[jv];
      }
    }
  }
  free(add_cond);
  free(add_mean);

  return (losses);
}





