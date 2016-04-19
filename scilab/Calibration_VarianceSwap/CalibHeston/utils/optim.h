#ifndef OPTIMLIB_H
#define OPTIMLIB_H

#include <string.h> 
#include <stdio.h> 
#include "nsp-math.h"

typedef struct { double r, i; } doubleC;

#include "blas.h"

#define Sciprintf printf 

/* 
 * data used when nsp is used to evaluate an  
 * objective function  
 */

typedef struct _opt_simul_data opt_simul_data;

typedef void (*opt_simul)(int *ind, int *n, double *x, double *f, double *g,  opt_simul_data *optim_data);
typedef int (*opt_prosca)(int *n, double *x, double *y, double *ps,  opt_simul_data *optim_data);
typedef int (*opt_ct)(int *n, double *f, double *g, opt_simul_data *optim_data);

typedef void (*U_fp)();

/* for icse */
extern int C2F(dgefa) (double *, int *, int *, int *, int *);
extern int C2F(dgesl) (double *, int *, int *, int *, double *, int *);
extern int dmmul_scicos (double *,   int *,   double *,   int *,   double *, 
			 int *,   int *,   int *,   int *);

extern double optim_dlamch (char *, long int);
extern double optim_dpmpar (int *i__);
extern double optim_enorm (int *n, double *x);

extern int optim_ctcab (int *n, double *u, double *v,opt_simul_data *optim_data);
extern int optim_ctonb (int *n, double *u, double *v,opt_simul_data *optim_data);

extern int optim_dcube (double *t, double *f, double *fp, double *ta,double *fa, 
			double *fpa, double *tlower,double *tupper);

extern int optim_fcube (double *t, double *f, double *fp, double *ta,double *fa, 
			double *fpa, double *tlower,double *tupper);

extern int optim_fpq2 (int *inout, double *x, double *cx, double *fx,    double *gx, 
		       double *d__, double *sthalf,    double *penlty, int *iyflag, double *y, 
		       double *cy,    double *fy, double *gy, double *z__, double *cz, 
		       double *fz, double *gz, double *gg, double *hh,    double *s);

extern int optim_fuclid (int *n, double *x, double *y, double *ps, 
			 opt_simul_data *optim_data);

extern int optim_gcbd (int *indgc, opt_simul simul, char *nomf, int *n, double *x,    
		       double *f, double *g, int *imp, int *io, double *zero,   
		       int *napmax, int *itmax, double *epsf, double *epsg,    
		       double *epsx, double *df0, double *binf, double *bsup,    
		       int *nfac, double *vect, int *nvect, int *ivect,    
		       int *nivect, opt_simul_data *optim_data,  long int nomf_len);

extern int optim_icscof (int *ico, int *ntob, int *nex, int *nob, double *yob, double *ob, double *cof);

extern int optim_icse (int *ind, int *nu, double *u, double *co, double *g,    
		       int *itv, double *rtv, double *dtv, U_fp icsef,    
		       U_fp icsec2, U_fp icsei);

extern int optim_icse0 (int *nu, double *t0, double *tf, double *dti,double *dtf, 
			double *ermx, int *iu, int *nuc,int *nuv, int *ilin, int *nti, 
			int *ntf, int *ny,int *nea, int *itmx, int *nex, int *nob, 
			int *ntob,int *ntobi, int *nitu, int *ndtu, int *nitv,int *nrtv, int *ndtv);

extern int optim_icse1 (int *ind, int *nu, double *u, U_fp icsef, double *y0,double *ytot, 
			double *f, double *b, double *fy,double *fu, int *ipv1, double *dm,
			double *yold,double *smold, double *yint, double *yerr,double *dif1,
			double *dif2, double *dif3, int *itu,double *dtu, double *t0, double *tf,
			double *dti,double *dtf, double *ermx, int *iu, int *nuc,int *nuv, int *ilin,
			int *nti, int *ntf, int *ny,int *nea, int *itmx, int *nex, int *nob, 
			int *ntob,int *ntobi, int *nitu, int *ndtu);

extern int optim_icse2 (int *ind, int *nu, double *u, double *co, double *g,U_fp icsef, U_fp icsec2, U_fp icsei, double *y0,double *tob, double *obs, double *ob, double *ytot,double *f, double *b, double *fy, double *fu,int *ipv2, int *itob, double *cof, double *ytob,double *c2y, double *y0u, double *dmy, double *smy,double *oldmu, double *y, double *oldp, double *p,double *p0, double *gt, double *yob, double *d__,int *itu, double *dtu, double *t0, double *tf,double *dti, double *dtf, double *ermx, int *iu,int *nuc, int *nuv, int *ilin, int *nti, int *ntf,int *ny, int *nea, int *itmx, int *nex, int *nob,int *ntob, int *ntobi, int *nitu, int *ndtu,char *nomf, char *nomc, char *nomi);

extern int optim_icsec2 (int *indc, int *nu, double *tob, double *obs, double *cof, double *ytob, double *ob, double *u, double *c__, double *cy, double *g, double *yob, double *d__, int *itu, double *dtu, double *t0, double *tf, double *dti, double *dtf, double *ermx, int *iu, int *nuc, int *nuv, int *ilin, int *nti, int *ntf, int *ny, int *nea, int *itmx, int *nex, int *nob, int *ntob, int *ntobi, int *nitu, int *ndtu);
extern int optim_icsei (int *indi, int *nui, double *u, double *y0,double *y0u, int *itu, double *dtu, double *t0,double *tf, double *dti, double *dtf, double *ermx,int *iu, int *nuc, int *nuv, int *ilin, int *nti,int *ntf, int *ny, int *nea, int *itmx, int *nex,int *nob, int *ntob, int *ntobi, int *nitu,int *ndtu);

extern int optim_majour (double *hm, double *hd, double *dd, int *n, double *hno, int *ir,
			 int *indic, double *eps);
extern int optim_majysa (int *n, int *nt, int *np, double *y, double *s, double *ys, int *lb, 
			 double *g, double *x, double *g1, double *x1, int *index, int *ialg, int *nb);

extern int optim_n1fc1 (opt_simul simul, opt_prosca prosca, int *n, double *xn,double *fn,
			double *g, double *dxmin, double *df1,double *epsf, double *zero,
			int *imp, int *io,int *mode, int *iter, int *nsim, int *memax,
			int *iz,double *rz, double *dz, opt_simul_data *optim_data);


extern int optim_n1gc2 (opt_simul simul, opt_prosca prosca, int *n, double *x, 
			double *f,double *g, double *dxmin, double *df1, 
			double *epsrel,int *imp, int *io, int *mode, 
			int *niter, int *nsim,double *rz, int *nrz,
			opt_simul_data *optim_data);

extern int optim_n1qn1 (opt_simul simul, int *n, double *x, double *f, double *g,double *var, double *eps, int *mode, int *niter,int *nsim, int *imp, int *lp, double *zm,
			opt_simul_data *optim_data);

extern int optim_n1qn2 (opt_simul simul, opt_prosca prosca, int *n, double *x, double *f,double *g, double *dxmin, double *df1, double *epsg,int *impres, int *io, int *mode, int *niter,int *nsim, double *dz, int *ndz,opt_simul_data *optim_data); 

extern int optim_n1qn3(opt_simul simul,opt_prosca prosca, opt_ct ctonb, opt_ct ctcab,int *n, double *x, double *f, double *g,double *dxmin, double *df1, double *epsg, int *impres,int *io, int *mode, int *niter, int *nsim, double *dz,int *ndz,opt_simul_data *optim_data); 

extern int optim_nlis0 (int *n, opt_simul simul, opt_prosca prosca, double *xn, double *fn,
			double *fpn, double *t, double *tmin, double *tmax, double *d__,
			double *g, double *amd, double *amf, int *imp, int *io,
			int *logic, int *nap, int *napmax, double *x,opt_simul_data *optim_data); 

extern int optim_proj (int *n, double *binf, double *bsup, double *x);

extern int optim_qnbd (int *indqn, opt_simul simul, int *n, double *x, double *f,    
		       double *g, int *imp, int *io, double *zero,   
		       int *napmax, int *itmax, double *epsf, double *epsg,
		       double *epsx, double *df0, double *binf, double *bsup, 
		       int *nfac, double *trav, int *ntrav, int *itrav, 
		       int *nitrav,opt_simul_data *optim_data);

extern int optim_satur (int *n, double *x, double *binf, double *bsup,double *d__, 
			double *ttmin, double *ttsup,double *topt, double *tg, 
			double *td, double *tmi,int *icoi, int *icos, int *iproj);

int optim_rlbd (int *indrl, int *n, opt_simul simul, double *x, double *binf,
		double *bsup, double *f, double *hp, double *t, double *tmax,
		double *d__, double *gn, double *tproj, double *amd, double *amf,
		int *imp, int *io, double *zero, int *nap, int *napmax,
		double *xn,opt_simul_data *optim_data);

#endif /* OPTIMLIB_H */
