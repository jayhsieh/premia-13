#include 		"cdo.h"
#include                "maths.h"

double			**lg_numdef(const CDO		*cdo, 
				    const copula	*cop,
				    const grid		*t,
				    const cond_prob	*cp)
{
    double		**nd;
    double		*r;
    double		*r_cpy;
    int			jt;
    int			jr;
    int			jn;
    int			jv;
    
    nd = malloc(t->size * sizeof(double*));
    r = malloc((cdo->n_comp+1) * sizeof(double));
    r_cpy = malloc((cdo->n_comp+1) * sizeof(double));
    for (jt = 0; jt < t->size; jt++) {
	nd[jt] = malloc((cdo->n_comp+1) * sizeof(double));
	for (jr = 0; jr < cdo->n_comp+1; jr++) {
	    nd[jt][jr] = 0;
	}
	for (jv = 0; jv < cop->size; jv++) {
	    r[0] = 1.;
	    for (jr = 1; jr < cdo->n_comp+1; jr++) 
		r[jr] = 0;
	    for (jn = 0; jn < cdo->n_comp; jn++) {
		for (jr = 0; jr < jn+1; jr++) { 
		    r_cpy[jr] = r[jr];
		}
		r[0] = (1. - cp->p[jn][jt][jv]) * r_cpy[0];
		for (jr = 1; jr < jn+2; jr++) { 
		    r[jr] += cp->p[jn][jt][jv] * (r_cpy[jr-1] - r[jr]);
		}
	    }
	    for (jr = 0; jr < cdo->n_comp+1; jr++) {
		nd[jt][jr] += r[jr] * cop->weights[jv];
	    }
	}
    }
    free(r);
    

    return (nd);
}


double			**lg_losses(const CDO		*cdo, 
				    const copula	*cop,
				    const grid		*t,
				    const grid 		*x,
				    const cond_prob	*cp)
{
    double		**losses;
    grid		*u;
    complex		*phi_cov;
    complex		prod;
    int			jt;
    int			ju;
    int			jv;
    int			jn;
    fftw_complex	*fft_in;
    fftw_complex	*fft_out;
    fftw_plan		fft_plan;
    double		C;
    double		F;

    fft_in = fftw_malloc(x->size * sizeof(fftw_complex));
    fft_out = fftw_malloc(x->size * sizeof(fftw_complex));
    fft_plan = fftw_plan_dft_1d(x->size, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    u = create_grid(x->size);
    C = x->data[x->size-1];
    F = (C_2PI * pow(x->size-1.,2.)) / (2. * C * (double) x->size);
    for (ju = 0; ju < x->size; ju++)
	u->data[ju] = - F + 2.*F*(ju / ((double) x->size-1.));
    phi_cov = malloc(cdo->n_comp * sizeof(complex));
    losses = malloc(t->size * sizeof(double*)); 
    for (jt = 0; jt < t->size; jt++) {
	for (ju = 0; ju < x->size; ju++) {
	    for (jn = 0; jn < cdo->n_comp; jn++) {
		phi_cov[jn] = PHI_COV(jn, u->data[ju] * cdo->C[jn]->nominal);
	    }
	    fft_in[ju] = 0.;
	    for (jv = 0; jv < cop->size; jv++) {
		prod = 1.;
		for (jn = 0; jn < cdo->n_comp; jn++) {
		    prod = prod * 
			( 1. + cp->p[jn][jt][jv] * (phi_cov[jn] - 1.) );
		}
		fft_in[ju] += prod * cop->weights[jv];
	    }
	}
	fftw_execute(fft_plan);
	losses[jt] = malloc(x->size * sizeof(double));
	for (ju = 0; ju < x->size; ju++) {
	    losses[jt][ju] = (C / ((double) x->size-1.)) * creal(cexp(I*F*x->data[ju]) * fft_out[ju]) * 2.*F/(C_2PI * ((double) x->size-1.));
	}
    }

    return (losses);
}



double			**lg_numdef1(const CDO		*cdo, 
				    const copula	*cop,
				    const grid		*t,
				    const cond_prob	*cp)
{    
    double		**nd;
 ;
     double		**r;
    double		**r_cpy;
    double              *z ;
   
    int			jt;
    int			jr;
    int			jn;
    int			jv,jw;
  
   
    nd = malloc(t->size * sizeof(double*));
   
   r = malloc(((cdo->n_comp+1)) * sizeof(double*));
    r_cpy = malloc(((cdo->n_comp+1)) * sizeof(double*));
  
    for (jr=0;jr<cdo->n_comp+1;jr++){
     r[jr]=malloc((cop->size)*sizeof(double));
     r_cpy[jr]=malloc((cop->size)*sizeof(double));
    }
     z = malloc((cdo->n_comp+1) * sizeof(double));
    
    for (jt = 0; jt < t->size; jt++) {
	nd[jt] = malloc((cdo->n_comp+1) * sizeof(double));
       
	for (jr = 0; jr < cdo->n_comp+1; jr++) {
	    nd[jt][jr] = 0;
      
            z[jr]=0;
          
	   
        }
      
       for (jv = 0; jv < cop->size; jv++) { 
               
         for( jw=0;jw<cop->size;jw++){
               r[0][jw] = 1.;
              for (jr = 1; jr < cdo->n_comp+1; jr++) { 
               r[jr][jw] = 0;
              
               }
             
	      for (jn = 0; jn < cdo->n_comp; jn++) {
		   for (jr = 0; jr < jn+1; jr++) { 
		    r_cpy[jr][jw] = r[jr][jw];
		   }
		   r[0][jw] = (1. - cp->p[jn][jt][jv+jw*cop->size]) * r_cpy[0][jw];
		    for (jr = 1; jr < jn+2; jr++) { 
		    r[jr][jw] += cp->p[jn][jt][jv+jw*cop->size] * (r_cpy[jr-1][jw] - r[jr][jw]);
		    }
	     
             }
             for (jr = 0; jr < cdo->n_comp+1; jr++) { 
             z[jr]+=r[jr][jw]*(cop->weights[jw+cop->size]);
           
           }
         }

          for (jr = 0; jr < cdo->n_comp+1; jr++) { 
          nd[jt][jr] +=z[jr]* cop->weights[jv];
	  z[jr] =0;
         

          }
      }
    
   }

 

    free(r);
    
    

    return (nd);
}

double			**lg_losses1(const CDO		*cdo, 
				    const copula	*cop,
				    const grid		*t,
				    const grid 		*x,
				    const cond_prob	*cp)
{
    double		**losses;
    grid		*u;
    complex		*phi_cov;
    complex		*prod;
    complex             prod1=0;
    int			jt;
    int			ju;
    int			jv;
     int                iv;
    int			jn;
    fftw_complex	*fft_in;
    fftw_complex	*fft_out;
    fftw_plan		fft_plan;
    double		C;
    double		F;

    fft_in = fftw_malloc(x->size * sizeof(fftw_complex));
    fft_out = fftw_malloc(x->size * sizeof(fftw_complex));
    fft_plan = fftw_plan_dft_1d(x->size, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    prod=malloc(cop->size*sizeof(complex));
    u = create_grid(x->size);
    C = x->data[x->size-1];
    F = (C_2PI * pow(x->size-1.,2.)) / (2. * C * (double) x->size);
    for (ju = 0; ju < x->size; ju++)
	u->data[ju] = - F + 2.*F*(ju / ((double) x->size-1.));
    phi_cov = malloc(cdo->n_comp * sizeof(complex));
    losses = malloc(t->size * sizeof(double*)); 
    for (jt = 0; jt < t->size; jt++) {
	for (ju = 0; ju < x->size; ju++) {
	    for (jn = 0; jn < cdo->n_comp; jn++) {
		phi_cov[jn] = PHI_COV(jn, u->data[ju] * cdo->C[jn]->nominal);
	    }
	    fft_in[ju] = 0.;
	    for (jv = 0; jv < cop->size; jv++) {
               for(iv=0;iv<cop->size;iv++){
		prod[iv] = 1.0;
		for (jn = 0; jn < cdo->n_comp; jn++) {
		    prod[iv] = prod[iv] * 
			( 1. + cp->p[jn][jt][jv+iv*cop->size] * (phi_cov[jn] - 1.) );
		}
                prod1=0.0;
               }
                for(iv=0;iv<cop->size;iv++){
               prod1=prod1+prod[iv]*(cop->weights[iv+cop->size]);
                }
		fft_in[ju] += prod1 * cop->weights[jv];
	    }
	}
	fftw_execute(fft_plan);
	losses[jt] = malloc(x->size * sizeof(double));
	for (ju = 0; ju < x->size; ju++) {
	    losses[jt][ju] = (C / ((double) x->size-1.)) * creal(cexp(I*F*x->data[ju]) * fft_out[ju]) * 2.*F/(C_2PI * ((double) x->size-1.));
	}
    }

    return (losses);
}









