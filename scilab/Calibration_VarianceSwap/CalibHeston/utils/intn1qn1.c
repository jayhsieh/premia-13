/* Nsp
 * Copyright (C) 2006 Jean-Philippe Chancelier Enpc/Cermics
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * n1qn1 : interface 
 *--------------------------------------------------------------------------*/

#include "nsp/interf.h"
#include "optim.h"
#include "nsp/gtk/gobject.h" /* FIXME: nsp_gtk_eval_function */

static void n1qn1_clean(opt_simul_data *obj);
static int n1qn1_prepare(int n,NspObject *fcn,NspObject *args,opt_simul_data *obj);
static void n1qn1_lfcn(int *ind, int *n, double *x, double *f, double *g,void *n1qn1_obj_d);

static NspObject *get_function(Stack stack, int pos,opt_simul_data *data);

typedef enum { algo_qn=0, algo_gc=1, algo_nd=2,algo_unknown=3 } optim_algo;
static optim_algo check_optim_algo(Stack stack,const char *fname,const char *algo); 

int int_optim_n1qn1 (Stack stack, int rhs, int opt, int lhs)
{ 
  int constraints = FALSE;/* TRUE if bound constraints are given */
  optim_algo algo ;
  NspObject *args=NULL; /* extra arguments to be passed to the simulator */
  NspMatrix *var=NULL,*var1=NULL,*g=NULL,*epsx=NULL,*xinf=NULL,*xsup=NULL;
  NspObject *fcn;
  int imp=0,itermax=100,simmax=100,mode=1,retf=0,ret = RET_BUG,lp=0,mem=10,indic=4,i,io=0;
  double epsf=0.0,epsg,df0=1.0,epsm,f,zng,dxmin;
  char *alg=NULL;
  nsp_option opts[] ={{ "args",obj,NULLOBJ,-1},
		      { "imp",s_int,NULLOBJ,-1},
		      { "itermax",s_int,NULLOBJ,-1},
		      { "simmax",s_int,NULLOBJ,-1},
		      { "var",realmat,NULLOBJ,-1},
		      { "alg",string,NULLOBJ,-1},
		      { "df0",s_double,NULLOBJ,-1},
		      { "mem",s_int,NULLOBJ,-1},
		      { "epsf",s_double,NULLOBJ,-1},
		      { "epsg",s_double,NULLOBJ,-1},
		      { "epsx",realmat,NULLOBJ,-1},
		      { "xinf",realmat,NULLOBJ,-1},
		      { "xsup",realmat,NULLOBJ,-1},
		      { NULL,t_end,NULLOBJ,-1}};

  NspMatrix *X, *work1=NULLMAT, *work2=NULLMAT, *work3=NULLMAT;
  opt_simul_data optim_data;
  
  CheckStdRhs(2,2);
  CheckLhs(1,3);

  if ((fcn = get_function(stack,1,&optim_data))==NULL)  return RET_BUG;
  if ((X = GetRealMatCopy(stack,2)) == NULLMAT) return RET_BUG;

  epsg = epsm = nsp_dlamch("eps");

  if ( get_optional_args(stack,rhs,opt,opts,&args,&imp,&itermax,&simmax,&var,&alg,&df0,&mem,
			 &epsf,&epsg,&epsx,&xinf,&xsup) == FAIL) return RET_BUG;

  algo = check_optim_algo(stack,NspFname(stack),alg);
  
  if ( n1qn1_prepare(X->mn,fcn,args,&optim_data)==FAIL) return RET_BUG;

  if ((g =nsp_matrix_create(NVOID,'r',X->mn,1)) == NULLMAT) goto bug;

  /* check constraints */
  
  if ( xinf != NULLMAT)
    {
      if ( xsup == NULLMAT )
	{
	  Scierror("Error: you must provide xsup if xinf is given\n");
	  goto bug;
	}
      if ( xinf->mn != X->mn )
	{
	  Scierror("Error: xinf should be of lenght %d\n",X->mn);
	  goto bug;
	}
      constraints = TRUE;
    }

  if ( xsup != NULLMAT)
    {
      if ( xinf == NULLMAT )
	{
	  Scierror("Error: you must provide xinf if xsup is given\n");
	  goto bug;
	}
      if ( xsup->mn != X->mn )
	{
	  Scierror("Error: xsup should be of lenght %d\n",X->mn);
	  goto bug;
	}
      constraints = TRUE;
    }
  

  switch ( algo )
    {
    case algo_qn : 
      if ( constraints == FALSE ) 
	{
	  /* qn without constraints 
	   * n1qn1 with (epsg and var as stop criteria)
	   */
	  if ((work1 =nsp_matrix_create(NVOID,'r',X->mn*(X->mn+13)/2,1)) == NULLMAT) goto bug;
	  if ( var == NULLMAT) 
	    {
	      if ((var1 =nsp_matrix_create(NVOID,'r',X->mn,1)) == NULLMAT) goto bug;
	      nsp_mat_set_rval(var1,0.1);
	    }
	  else 
	    {
	      if ( var->mn != X->mn) 
		{
		  Scierror("Error: optional argument var should be of size %d\n",X->mn);
		  goto bug;
		}
	    }
	  retf = optim_n1qn1(optim_data.f_fcn,&X->mn,X->R,&f,g->R,
			     (var == NULL) ? var1->R : var->R,
			     &epsg,&mode,&itermax,&simmax,&imp,&lp,
			     work1->R,&optim_data);
	  n1qn1_clean(&optim_data);
	  if ( retf == FAIL ) 
	    {
	      Scierror("Error: execution of %s aborted\n",NspFname(stack));
	      goto bug;
	    }
	  break;
	}
      else 
	{
	  int nfac = 0,indopt;
	  /* qn with constraints  */
	  if ( epsx == NULLMAT ) 
	    {
	      if ((epsx =nsp_matrix_create(NVOID,'r',X->mn,1)) == NULLMAT) goto bug;
	      nsp_mat_set_rval(epsx,epsm);
	    }
	  else
	    {
	      if ( epsx->mn != X->mn ) 
		{
		  Scierror("Error: optional argument epsx should be of size %d\n",X->mn);
		  goto bug;
		}
	    }

	  if ((work1 =nsp_matrix_create(NVOID,'r',X->mn*(X->mn+1)/2 + 4*X->mn + 1,1)) == NULLMAT)
	    goto bug;
	  if ((work2 =nsp_matrix_create(NVOID,'r',2*X->mn,1)) == NULLMAT)
	    goto bug;
	  indic = 4;
	  optim_data.f_fcn(&indic,&X->mn,X->R,&f,g->R,&optim_data);
	  if ( indic <= 0) 
	    {
	      Scierror("Error: optim cost evaluation failed \n");
	      goto bug;
	    }
	  indopt=1;
	  optim_qnbd(&indopt,optim_data.f_fcn,&X->mn,X->R,&f,g->R,&imp,&io,&epsm,&itermax,
		     &simmax,&epsf,&epsg,epsx->R,&df0,xinf->R,xsup->R,&nfac,
		     work1->R,&work1->mn, work2->I,&work2->mn,&optim_data);
	  break;
	}
    case algo_gc : 
      if ( constraints == FALSE ) 
	{
	  /*  n1qn3: conjugate gradient without constraints */
	  indic = 4;
	  optim_data.f_fcn(&indic,&X->mn,X->R,&f,g->R,&optim_data);
	  if ( indic <= 0) 
	    {
	      Scierror("Error: optim cost evaluation failed \n");
	      goto bug;
	    }
	  /* compute epsrel */
	  zng=0.0;
	  for ( i = 0 ; i < X->mn ; i++) zng += g->R[i]*g->R[i];
	  zng=sqrt(zng);
	  if (zng > 0.0) epsg /=zng;
	  /*  dxmin */
	  dxmin= epsm;
	  if ( epsx != NULL) 
	    {
	      if ( epsx->mn != X->mn ) 
		{
		  Scierror("Error: optional argument epsx should be of size %d\n",X->mn);
		  goto bug;
		}
	      dxmin=epsx->R[0];
	      for ( i=1; i < X->mn ; i++) dxmin=Min(dxmin,epsx->R[i]);
	    }
	  /* working arrays  */
	  if ((work1 =nsp_matrix_create(NVOID,'r',4*X->mn + mem*(2*X->mn + 1),1)) == NULLMAT) 
	    goto bug;
	  optim_n1qn3(optim_data.f_fcn,optim_fuclid,optim_ctonb,optim_ctcab,
		      &X->mn,X->R,&f,g->R,&dxmin,&df0,&epsg,&imp,&io,&mode,&itermax,&simmax,
		      work1->R,&work1->mn,&optim_data);
	  if ( mode <= 0) 
	    {
	      Scierror("Error: optim cost evaluation failed \n");
	      goto bug;
	    }
	}
      else
	{
	  int nt,ntv,nitv,indopt=1,nfac=0;
	  indic = 4;
	  optim_data.f_fcn(&indic,&X->mn,X->R,&f,g->R,&optim_data);
	  if ( indic <= 0) 
	    {
	      Scierror("Error: optim cost evaluation failed \n");
	      goto bug;
	    }
	  if ( epsx == NULLMAT ) 
	    {
	      if ((epsx =nsp_matrix_create(NVOID,'r',X->mn,1)) == NULLMAT) goto bug;
	      nsp_mat_set_rval(epsx,epsm);
	    }
	  else
	    {
	      if ( epsx->mn != X->mn ) 
		{
		  Scierror("Error: optional argument epsx should be of size %d\n",X->mn);
		  goto bug;
		}
	    }
	  nt= Max(1,mem/3);
	  ntv=X->mn*(5 + 3*nt) + 2*nt +1;
	  nitv=X->mn + nt + 1;
	  if ((work1 =nsp_matrix_create(NVOID,'r',ntv,1)) == NULLMAT)
	    goto bug;
	  if ((work2 =nsp_matrix_create(NVOID,'r',nitv,1)) == NULLMAT)
	    goto bug;
	  /* "" should be the name of the simulator */
	  optim_gcbd(&indopt,optim_data.f_fcn,"",&X->mn,X->R,
		     &f,g->R,&imp,&io,&epsm,&simmax,&itermax,&epsf,&epsg,epsx->R,
		     &df0,xinf->R,xsup->R,&nfac,work1->R,&work1->mn,work2->I,&work2->mn,
		     &optim_data,strlen(""));
	  break;
	}
      break;
    case algo_nd : 
      /*  n1fc1:  'nd' without constraints */
      indic = 4;
      optim_data.f_fcn(&indic,&X->mn,X->R,&f,g->R,&optim_data);
      if ( epsx == NULLMAT ) 
	{
	  if ((epsx =nsp_matrix_create(NVOID,'r',X->mn,1)) == NULLMAT) goto bug;
	  nsp_mat_set_rval(epsx,epsm);
	  
	}
      else
	{
	  if ( epsx->mn != X->mn ) 
	    {
	      Scierror("Error: optional argument epsx should be of size %d\n",X->mn);
	      goto bug;
	    }
	}
      if ((work1 =nsp_matrix_create(NVOID,'r',2*mem + 2,1)) == NULLMAT) goto bug;
      if ((work2 =nsp_matrix_create(NVOID,'r',5*X->mn + (X->mn+4)*mem,1)) == NULLMAT) goto bug;
      if ((work3 =nsp_matrix_create(NVOID,'r',(mem+9)*mem + 8,1)) == NULLMAT) goto bug;
      optim_n1fc1(optim_data.f_fcn,optim_fuclid,&X->mn,X->R,&f,g->R,epsx->R,&df0,&epsf,
		  &epsm,&imp,&io,&mode,&itermax,&simmax,&mem,
		  work1->I,work2->R,work3->R,&optim_data);
      if ( mode <= 0 )
	{
	  Scierror("Error: stop in optim cost function\n");
	  goto bug;
	}
      switch ( mode ) 
	{
	case 1: if ( imp != 0 ) 
	    Sciprintf("normal stop, at last iteration f decreases by less than %g\n",epsf);break;
	case 2: Scierror("Error: incoherent optim call\n"); goto bug;
	case 3: if ( imp != 0 ) Sciprintf("reduce the x scale\n"); break;
	case 4: if ( imp != 0 ) Sciprintf("maximum number of iterations is reached\n"); break;
	case 5: if ( imp != 0 ) Sciprintf("maximum number of calls to f is reached\n"); break;
	case 6: if ( imp != 0 ) 
	    Sciprintf("optimization stops because too small variations for x\n");break;
	case 7: if ( imp != 0 ) Sciprintf("fprf2 failed\n");break;
	case 8: if ( imp != 0 ) Sciprintf("we are starting to loop\n");break;
	}
      break;
    case algo_unknown: 
      Scierror("Error: unknown algorithm\n");
      goto bug;
    }
  NSP_OBJECT(X)->ret_pos=1;
  ret = Max(lhs,1);
 bug:
  if ( work1 != NULL) nsp_matrix_destroy(work1);
  if ( g != NULL) nsp_matrix_destroy(g);
  if ( var1 != NULL) nsp_matrix_destroy(var1);
  return ret;
}


static optim_algo check_optim_algo(Stack stack,const char *fname,const char *algo)
{
  static char *Table[] = {"qn",  "gc",  "nd", NULL};
  int rep ; 
  if ( algo == NULL ) return algo_qn;
  rep = is_string_in_array(algo,Table,1);
  if ( rep < 0 ) 
    {
      Scierror("Error:\toptional argument alg of function optim has a wrong value %s\n",algo);
      Scierror("\texpected values are 'qn',or 'gc,  or 'nd' \n");
      return algo_unknown;
    }
  return rep;
}



/*--------------------------------------------------------------
 * objective 
 *-------------------------------------------------------------*/

/*
 * n1qn1_prepare:
 * if m < 0 prepare for n1qn1d or n1qn1j 
 * if m >= 0 then for lmdiff or lmderr
 **/

static int n1qn1_prepare(int n,NspObject *fcn,NspObject *args,opt_simul_data *obj)
{
  if (( obj->fcn =nsp_object_copy(fcn)) == NULL) return FAIL;
  if (( nsp_object_set_name(obj->fcn,"n1qn1_fcn")== FAIL)) return FAIL;
  if ( args != NULL ) 
    {
      if (( obj->args = nsp_object_copy(args)) == NULL ) return FAIL;
      if (( nsp_object_set_name((NspObject *) obj->args,"arg")== FAIL)) return FAIL;
    }
  else 
    {
      obj->args = NULL;
    }
  if ((obj->x = nsp_matrix_create("x",'r',n,1))== NULL) return FAIL;
  if ((obj->ind = nsp_matrix_create("ind",'r',1,1))== NULL) return FAIL;
  return OK;
}

/*
 * n1qn1_clean:
 **/

static void n1qn1_clean(opt_simul_data *obj)
{
  if ( obj->args != NULL) nsp_object_destroy(&obj->args);
  nsp_object_destroy(&obj->fcn);
  nsp_matrix_destroy(obj->x);
  nsp_matrix_destroy(obj->ind);
}

/*
 * soft evaluation of fcn 
 * [f,g,ind]=cost(x,ind [,args]) 
 */

static void n1qn1_lfcn(int *ind, int *n, double *x, double *f, double *g,void *n1qn1_obj_d)
{
  opt_simul_data *n1qn1_obj = n1qn1_obj_d;
  NspObject *targs[3];/* arguments to be transmited to n1qn1_obj->objective */
  NspObject *nsp_ret[3];
  int nret = 3,nargs = 2,ret=FAIL,i;
  
  targs[0]= NSP_OBJECT(n1qn1_obj->x); 
  memcpy(n1qn1_obj->x->R,x,n1qn1_obj->x->mn*sizeof(double));
  targs[1]= NSP_OBJECT(n1qn1_obj->ind); 
  n1qn1_obj->ind->R[0]= (int) *ind;
  if (n1qn1_obj->args != NULL ) 
    {
      targs[2]= NSP_OBJECT(n1qn1_obj->args);
      nargs= 3;
    }
  /* FIXME : a changer pour metre une fonction eval standard */
  if ( nsp_gtk_eval_function((NspPList *)n1qn1_obj->fcn ,targs,nargs,nsp_ret,&nret)== FAIL)
    {
      goto stop;
    }
  if ( nret != 3)
    {
      Scierror("Error: optim cost function should return three values\n");
      goto stop;
    }
  if ( IsMat(nsp_ret[0]) && ((NspMatrix *) nsp_ret[0])->rc_type == 'r' 
       &&((NspMatrix *) nsp_ret[0])->mn == 1 ) 
    {
      *f =((NspMatrix *)  nsp_ret[0])->R[0];
    }
  else 
    {
      Scierror("Error: first returned argument of optim cost function should be scalar\n");
      goto stop;
    }
  if ( IsMat(nsp_ret[1]) && ((NspMatrix *) nsp_ret[1])->rc_type == 'r' 
       &&((NspMatrix *) nsp_ret[1])->mn == n1qn1_obj->x->mn ) 
    {
      memcpy(g,((NspMatrix *) nsp_ret[1])->R,n1qn1_obj->x->mn*sizeof(double));
    }
  else 
    {
      Scierror("Error: second returned  argument of optim cost function should be a vector of size %d\n",n1qn1_obj->x->mn);
      goto stop;
    }
  if ( IsMat(nsp_ret[2]) && ((NspMatrix *) nsp_ret[2])->rc_type == 'r' 
       &&((NspMatrix *) nsp_ret[2])->mn == 1 ) 
    {
      *ind = (int) ((NspMatrix *) nsp_ret[2])->R[0];
    }
  else 
    {
      Scierror("Error: third returned argument of optim cost function should be scalar\n");
      goto stop;
    }
  ret=OK;
 stop: 
  for( i= 0 ; i < nret ; i++) nsp_void_object_destroy(&nsp_ret[i]);
  if (ret == FAIL) 
    {
      Scierror("Error: optim cost function evaluation failed\n");
      Scierror("\tsetting ind to zero to force stop\n");
      *ind= 0 ;
    }
}

/* FIXME: should be in a .h */

extern int SearchInDynLinks (char *op, int (**realop)());

static void optim_genros(int *ind, int *n, double *x, double *f, double *g, opt_simul_data *optim_data);


static NspObject *get_function(Stack stack, int pos,opt_simul_data *data)
{
  NspObject *obj;
  if ( (obj=nsp_get_object(stack,pos)) == NULLOBJ ) return NULLOBJ;
  if ( IsNspPList(obj) )
    {
      data->f_fcn = (opt_simul) n1qn1_lfcn;
      return obj;
    }
  else if ( IsString(obj) )
    {
      char *str = ((NspSMatrix *)obj)->S[0];
      int (*func) (void);
      if ( strcmp(str,"genros")==0 )
	{
	  data->f_fcn= (opt_simul) optim_genros;
	  return obj;
	}
      /* search string in the dynamically linked functions */
      if ( SearchInDynLinks(str, &func) == -1 )
	{
	  Scierror("Error: function %s is not dynamically linked in nsp\n",str);
	  return NULL;
	}
      data->f_fcn= (opt_simul) func;
      return obj;
    }
  else
    {
      Scierror("Error: argument %d should be a function or a string\n",pos,NspFname(stack));
      return NULL;
    }
  return NULL;
}


/* interface used to call hard coded examples (fcn and grad)
 * at nsp level
 * [f,g,ind]=optim_test(x,ind,'genros');
 */

static char *test_names[]={ "genros" , NULL };
static opt_simul test_functions[]={ optim_genros , NULL};

int int_minpack_n1qn1d_test (Stack stack, int rhs, int opt, int lhs)
{ 
  NspMatrix *X,*f,*g,*indic;
  int rep,ind;
  CheckStdRhs(2,2);
  CheckLhs(1,3);

  if ((X = GetRealMatCopy(stack,1)) == NULLMAT) return RET_BUG;
  if (GetScalarInt(stack,2,&ind) == FAIL) return RET_BUG;
  if ((rep= GetStringInArray(stack,3,test_names,1)) == -1) return RET_BUG; 
  if ((g =nsp_matrix_create(NVOID,'r',X->mn,1)) == NULLMAT) return RET_BUG;
  if ((f =nsp_matrix_create(NVOID,'r',1,1)) == NULLMAT) return RET_BUG;
  if ((indic =nsp_matrix_create(NVOID,'r',1,1)) == NULLMAT) return RET_BUG;
  (test_functions[rep])(&ind,&X->mn,X->R,f->R,g->R,NULL);
  if ( ind < 0 ) 
    {
      Scierror("Error: execution of %s stopped by ind <= 0\n",NspFname(stack));
      return RET_BUG;
    }
  indic->R[0]= ind;
  MoveObj(stack,1,NSP_OBJECT(f));
  if ( lhs >= 2) 
    MoveObj(stack,2,NSP_OBJECT(g));
  if ( lhs >= 3)
    MoveObj(stack,3,NSP_OBJECT(indic));
  return Max(lhs,1);
}

/*
 * Example of hard code cost function for optim 
 *  if n<=2 returns ind=0 
 */

static void optim_genros(int *ind, int *n, double *x, double *f, double *g, opt_simul_data *optim_data)
{
  double d1,ddzs=100;
  int i;
  if (*n < 3)
    {
      *ind = 0;
      return;
    }

  if (*ind == 2 || *ind == 4 )
    {
      /* cost */
      *f = 1.;
      for (i = 1; i < *n ; ++i)
	{
	  d1 = x[i] - x[i - 1]*x[i - 1];
	  *f +=  ddzs*(d1*d1) + (1 - x[i])*(1 - x[i]);
	}
    }
  if (*ind == 3 || *ind == 4 )
    {
      /* gradient */
      d1 = x[0];
      g[0] = ddzs * (-4.0) * (x[1] - d1 * d1) * x[0];
      for (i = 1; i < *n -1 ; ++i)
	{
	  g[i] = ddzs*(2.0) * (x[i] - x[i - 1]*x[i - 1]);
	  g[i] +=  - ddzs*(4.0) * (x[i + 1] - x[i]*x[i])*x[i] - (1 - x[i])*2;
	}
      g[*n-1] = ddzs * 2. * (x[*n-1] - x[*n-2]*x[*n-2]) - (1 - x[*n-1])* 2;
    }
}

