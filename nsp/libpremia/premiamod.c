/* Nsp
 * Copyright (C) 2007 Jean-Philippe Chancelier Enpc/Cermics
 * Jerome Lelong Enpc Inria
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
 * Interface for Premia
 *
 */

#include <nsp/object.h> 
#include <nsp/matrix.h> 
#include <nsp/smatrix.h> 
#include <nsp/list.h> 
#include <nsp/hash.h> 
#include <nsp/hobj.h> 
#include <nsp/type.h> 
#include <nsp/file.h> 
#include "nsp/serial.h"
#include "nsp/interf.h"
#include "nsp/list.h"
#include "nsp/nsptcl.h"
#include "nsp/dstring.h"

#define PremiaModel_Private
#include "premiamod.h"


#include "premia_vars.h"
#include "enums.h"
#include "var.h"
/* 
 * NspPremiaModel inherits from NspObject 
 */

int nsp_type_premiamodel_id=0;
NspTypePremiaModel *nsp_type_premiamodel=NULL;

/*
 * Type object for PremiaModel 
 * all the instance of NspTypePremiaModel share the same id. 
 * nsp_type_premiamodel: is an instance of NspTypePremiaModel 
 *    used for objects of NspPremiaModel type (i.e built with new_premiamodel) 
 * other instances are used for derived classes 
 */
NspTypePremiaModel *new_type_premiamodel(type_mode mode)
{
  NspTypePremiaModel *type= NULL;
  NspTypeObject *top;
  if (  nsp_type_premiamodel != 0 && mode == T_BASE ) 
    {
      /* initialization performed and T_BASE requested */
      return nsp_type_premiamodel;
    }
  if ((type =  malloc(sizeof(NspTypePremiaModel))) == NULL) return NULL;
  type->interface = NULL;
  type->surtype = (NspTypeBase *) new_type_object(T_DERIVED);
  if ( type->surtype == NULL) return NULL;
  type->attrs = premiamodel_attrs ; 
  type->get_attrs = (attrs_func *) int_get_attribute;
  type->set_attrs = (attrs_func *) int_set_attribute; 
  type->methods = premiamodel_get_methods; 
  type->new = (new_func *) new_premiamodel;

  
  top = NSP_TYPE_OBJECT(type->surtype);
  while ( top->surtype != NULL ) top= NSP_TYPE_OBJECT(top->surtype);
  
  /* object methods redefined for premiamodel */ 

  top->pr = (print_func *) nsp_premiamodel_print;                  
  top->dealloc = (dealloc_func *) nsp_premiamodel_destroy;
  top->copy  =  (copy_func *) nsp_premiamodel_copy;                 
  top->size  = (size_func *) nsp_premiamodel_size;                
  top->s_type =  (s_type_func *) nsp_premiamodel_type_as_string;  
  top->sh_type = (sh_type_func *) nsp_premiamodel_type_short_string;
  top->info = (info_func *) nsp_premiamodel_info ;                  
  /* top->is_true = (is_true_func  *) nsp_premiamodel_is_true; */
  /* top->loop =(loop_func *) nsp_premiamodel_loop;*/
  top->path_extract = (path_func *)  object_path_extract; 
  top->get_from_obj = (get_from_obj_func *) nsp_premiamodel_object;
  top->eq  = (eq_func *) nsp_premiamodel_eq;
  top->neq  = (eq_func *) nsp_premiamodel_neq;
  top->save  = (save_func *) nsp_premiamodel_xdr_save;
  top->load  = (load_func *) nsp_premiamodel_xdr_load;
  top->create = (create_func*) int_premiamodel_create;
  
  /* specific methods for premiamodel */
      
  type->init = (init_func *) init_premiamodel;

  /* 
   * PremiaModel interfaces can be added here 
   * type->interface = (NspTypeBase *) new_type_b();
   * type->interface->interface = (NspTypeBase *) new_type_C()
   * ....
   */
  if ( nsp_type_premiamodel_id == 0 ) 
    {
      /* 
       * the first time we get here we initialize the type id and
       * an instance of NspTypePremiaModel called nsp_type_premiamodel
       */
      type->id =  nsp_type_premiamodel_id = nsp_new_type_id();
      nsp_type_premiamodel = type;
      if ( nsp_register_type(nsp_type_premiamodel) == FALSE) return NULL;
      return ( mode == T_BASE ) ? type : new_type_premiamodel(mode);
    }
  else 
    {
      type->id = nsp_type_premiamodel_id;
      return type;
    }
}

/*
 * initialize PremiaModel instances 
 * locally and by calling initializer on parent class 
 */

static int init_premiamodel(NspPremiaModel *o,NspTypePremiaModel *type)
{
  /* jump the first surtype */ 
  if ( type->surtype->init(&o->father,type->surtype) == FAIL) return FAIL;
  o->type = type; 
  NSP_OBJECT(o)->basetype = (NspTypeBase *)type;
  /* specific */  
  premia_set_global_vars();
  return OK;
}

/*
 * new instance of PremiaModel 
 */

NspPremiaModel *new_premiamodel() 
{
  NspPremiaModel *loc; 
  /* type must exists */
  nsp_type_premiamodel = new_type_premiamodel(T_BASE);
  if ( (loc = malloc(sizeof(NspPremiaModel)))== NULLPREMIAMODEL) return loc;
  /* initialize object */
  if ( init_premiamodel(loc,nsp_type_premiamodel) == FAIL) return NULLPREMIAMODEL;
  return loc;
}

/*----------------------------------------------
 * Object method redefined for PremiaModel 
 *-----------------------------------------------*/
/*
 * size 
 */

static int nsp_premiamodel_size(NspPremiaModel *Mat, int flag)
{
  return 0;
}

/*
 * type as string 
 */
static char premiamodel_type_name[]="PremiaModel";
static char premiamodel_short_type_name[]="premiamodel";

static char *nsp_premiamodel_type_as_string(void)
{
  return(premiamodel_type_name);
}

static char *nsp_premiamodel_type_short_string(NspObject *v)
{
  return(premiamodel_short_type_name);
}

/*
 * A == B 
 */
static int nsp_premiamodel_eq(NspPremiaModel *A, NspObject *B)
{
  NspPremiaModel *loc = (NspPremiaModel *) B;
  if ( check_cast(B,nsp_type_premiamodel_id) == FALSE) return FALSE ;
  if ( A->obj == loc->obj ) return TRUE;
  return FALSE;
}

/*
 * A != B 
 */
static int nsp_premiamodel_neq(NspPremiaModel *A, NspObject *B)
{
  return ( nsp_premiamodel_eq(A,B) == TRUE ) ? FALSE : TRUE;
}

/*
 * delete 
 */
void nsp_premiamodel_destroy(NspPremiaModel *H)
{
  nsp_object_destroy_name(NSP_OBJECT(H));
  H->obj->ref_count--;
  if ( H->obj->ref_count == 0 )
    {      
      /* FREE(H->obj->mod.TypeModel); */
      if (H->obj->mod.TypeModel != NULL) 
        nsp_premia_free_vars(H->obj->mod.TypeModel,TRUE,H->obj->mod.nvar);
      if (H->obj->opt.TypeOpt != NULL) 
        nsp_premia_free_vars(H->obj->opt.TypeOpt,TRUE,H->obj->opt.nvar);
      if ( H->obj->meth.Name != NULL) 
        {
          nsp_premia_free_vars(H->obj->meth.Res,FALSE,
                               nsp_premia_get_nvar(H->obj->meth.Res));
          nsp_premia_free_vars(H->obj->meth.Par,FALSE,
                               nsp_premia_get_nvar(H->obj->meth.Par));
        }
      if (H->obj->it != NULL)
        {
          if (H->obj->it->range1 != NULL) nsp_matrix_destroy(H->obj->it->range1);
          if (H->obj->it->range2 != NULL) nsp_matrix_destroy(H->obj->it->range2);
          nsp_matrix_destroy(H->obj->it->prix);
          FREE(H->obj->it);
        }
      FREE(H->obj);
    }
  FREE(H);
}

/*
 * info 
 */
int nsp_premiamodel_info(NspPremiaModel *M, int indent,const char *name, int rec_level)
{
  const char *pname = (name != NULL) ? name : NSP_OBJECT(M)->name;
  Sciprintf1(indent,"%s\t= ...\t\tpremia\n",pname);
  return TRUE;
}

/*-----------------------------------
 * print an nsp_premia object
 *---------------------------------- */
static void nsp_premiamodel_print_vars(NspPremiaModel *self,const VAR *vars, int nvar, int indent, int depth)
{
  NspList *L=NULL;
  if (depth <= 0) return;
  if ((L = nsp_premia_get_var_names(self,vars, nvar)) == NULL) return;
  nsp_object_print(NSP_OBJECT(L), indent, NULL, 2);
  if (L!=NULL) nsp_list_destroy(L);
}

static void nsp_premiamodel_print_gen(NspPremiaModel *M,int indent,const char *name, int depth)
{
  if (depth<=0)
    {
      nsp_premiamodel_info(M,indent,name,0);
      return;
    }
  else
    {
      Sciprintf1(indent, "%s = \n", name);
      Sciprintf1(indent+2, "Asset type : %s\n",
                 (M->obj->asset.name != NULL) ? M->obj->asset.name : "undef");
      /* Model */
      Sciprintf1(indent+2, "Model : %s\n", (M->obj->mod.Name != NULL) ? M->obj->mod.Name : "undef");
      if (M->obj->mod.Name != NULL) 
        nsp_premiamodel_print_vars(M,M->obj->mod.TypeModel, M->obj->mod.nvar, indent+1, depth-1);
      /* Option */
      Sciprintf1(indent+2, "Option : %s\n", (M->obj->opt.Name != NULL) ? M->obj->opt.Name : "undef");
      if (M->obj->opt.Name != NULL)
        nsp_premiamodel_print_vars(M,M->obj->opt.TypeOpt, M->obj->opt.nvar, indent+1, depth-1);
      /* Method */
      Sciprintf1(indent+2, "Method : %s\n", (M->obj->meth.Name != NULL) ? M->obj->meth.Name : "undef");
      if (M->obj->meth.Name != NULL)
        nsp_premiamodel_print_vars(M,M->obj->meth.Par, nsp_premia_get_nvar(M->obj->meth.Par),
                                   indent+1, depth-1);
    }
}

/*
 * print 
 */
int nsp_premiamodel_print(NspPremiaModel *M,int indent,const char *name, int rec_level)
{
  const char *pname = (name != NULL) ? name : NSP_OBJECT(M)->name;
  if (user_pref.pr_as_read_syntax)
    {
      if ( strcmp(pname,NVOID) != 0) 
        Sciprintf1(indent,"%s=premiamodel_create();",pname);
      else 
        Sciprintf1(indent,"premiamodel_create();");
    }
  else 
    {
      nsp_premiamodel_print_gen(M, indent, pname, user_pref.pr_depth);
    }
  return TRUE;
}

static int nsp_premiamodel_save_vars (XDR *xdrs, NspPremiaModel *self, const VAR *vars, int nvar)
{
  NspList *L;
  if ((L = nsp_premia_get_var_names_forsave(self,vars, nvar)) == NULL) return FAIL;
  nsp_object_xdr_save (xdrs, NSP_OBJECT(L));
  if (L!=NULL) nsp_list_destroy(L);
  return OK;
}

/*
 * save method
 */
int nsp_premiamodel_xdr_save(XDR *xdrs, NspPremiaModel *M)
{
  /* ID */
  if (nsp_xdr_save_id(xdrs,(NspTypeBase *) nsp_type_premiamodel) == FAIL) return FAIL;
  if (nsp_xdr_save_string(xdrs, NSP_OBJECT(M)->name) == FAIL) return FAIL;
  /* Asset */
  if (nsp_xdr_save_string (xdrs, M->obj->asset.name) == FAIL) return FAIL;
  /* Model */
  if (nsp_xdr_save_string (xdrs, M->obj->mod.Name) == FAIL) return FAIL;
  if (nsp_premiamodel_save_vars(xdrs, M, M->obj->mod.TypeModel, M->obj->mod.nvar) == FAIL) return FAIL;
  /* Option */
  if (nsp_xdr_save_string (xdrs, M->obj->opt.Name) == FAIL) return FAIL;
  if (nsp_premiamodel_save_vars(xdrs, M, M->obj->opt.TypeOpt, M->obj->opt.nvar) == FAIL) return FAIL;
  /* Method */
  if (nsp_xdr_save_string (xdrs, M->obj->meth.Name) == FAIL) return FAIL;
  if (nsp_premiamodel_save_vars(xdrs, M, M->obj->meth.Par, nsp_premia_get_nvar(M->obj->meth.Par)) == FAIL) return FAIL;
  /* force init to make sure all dynamically allocated vars exist */
  M->obj->meth.Init(&(M->obj->meth), &(M->obj->opt));
  /* compute_err flag */
  if (nsp_xdr_save_i(xdrs,M->obj->compute_err) == FAIL) return FAIL;
  /* Results */
  if (nsp_premiamodel_save_vars(xdrs, M, M->obj->meth.Res, nsp_premia_get_nvar(M->obj->meth.Res)) == FAIL) return FAIL;
  return OK;
}

/*
 * Load method
 */
NspPremiaModel* nsp_premiamodel_xdr_load (XDR *xdrs)
{
  NspPremiaModel *M = NULL;
  NspList *L = NULL;
  char *buf = NULL;
  int len = 256, nvar;
  VAR *vars;
  static char name[NAME_MAXL];

  if (nsp_xdr_load_string(xdrs,name,NAME_MAXL) == FAIL) return NULL;

  if ((buf=malloc(len*sizeof(char)))==NULL) return NULL;
  nsp_type_premiamodel = new_type_premiamodel(T_BASE);
  if(( M = premiamodel_create(name,(NspTypeBase *) nsp_type_premiamodel)) == NULLPREMIAMODEL) 
    return NULL;
  InitErrorMsg();
  InitVar();

  /* Asset */
  nsp_xdr_load_string (xdrs, buf, len);
  nsp_premiamodel_set_asset (M, buf);
  
  /* Model */
  nsp_xdr_load_string (xdrs, buf, len);

  if (nsp_premiamodel_set_model_with_str (M, buf) != OK) return NULL;
  L = (NspList *) nsp_object_xdr_load (xdrs);
  if (nsp_premia_set_var_names_forload(M, M->obj->mod.TypeModel, M->obj->mod.nvar, L, 0,M->obj->mod.TypeModel) == FAIL) 
    {
      Scierror("Error: while trying to set model values\n");
      return NULL;
    }
  
  /* Option */
  nsp_xdr_load_string (xdrs, buf, len);
  if (nsp_premiamodel_set_option_with_str (M, buf) != OK) return NULL;
  L = (NspList *) nsp_object_xdr_load (xdrs);
  if (nsp_premia_set_var_names_forload(M, M->obj->opt.TypeOpt, M->obj->opt.nvar, L, 0,M->obj->opt.TypeOpt) == FAIL) 
    {
      Scierror("Error: while trying to set option values\n");
      return NULL;
    }
  M->obj->opt.Init(&(M->obj->opt), &(M->obj->mod));
  
  /* Method */
  nsp_xdr_load_string (xdrs, buf, len);
  if (nsp_premiamodel_set_method_with_str (M, buf) != OK) return NULL;
  L = (NspList *) nsp_object_xdr_load (xdrs);
  vars = M->obj->meth.Par;
  nvar = nsp_premia_get_nvar(vars);
  if (nsp_premia_set_var_names_forload(M, vars, nvar, L, 0, NULL) == FAIL) 
    {
      Scierror("Error: while trying to set method values\n");
      return NULL;
    }

  /* Results */
  nsp_xdr_load_i (xdrs, &(M->obj->compute_err));
  L = (NspList *) nsp_object_xdr_load (xdrs);
  vars = M->obj->meth.Res;
  nvar = nsp_premia_get_nvar(vars);
  if (nsp_premia_set_var_names_forload(M, vars, nvar, L, 0, NULL) == FAIL) 
    {
      Scierror("Error: while trying to set result values\n");
      return NULL;
    }
  free(buf); buf=NULL;
  return M;
}

/*-----------------------------------------------------
 * a set of functions used when writing interfaces 
 * for PremiaModel objects 
 * Note that some of these functions could become MACROS XXXXX 
 *-----------------------------------------------------*/
NspPremiaModel   *nsp_premiamodel_object(NspObject *O)
{
  /* Follow pointer */
  if ( check_cast(O,nsp_type_hobj_id) == TRUE)  O = ((NspHobj *) O)->O ;
  /* Check type */
  if ( check_cast (O,nsp_type_premiamodel_id) == TRUE ) return ((NspPremiaModel *) O);
  else 
    Scierror("Error:    Argument should be a %s\n",type_get_name(nsp_type_premiamodel));
  return NULL;
}

int IsPremiaModelObj(Stack stack, int i)
{
  return nsp_object_type(NthObj(i) , nsp_type_premiamodel_id);
}

int IsPremiaModel(NspObject *O)
{
  return nsp_object_type(O,nsp_type_premiamodel_id);
}

NspPremiaModel  *GetPremiaModelCopy(Stack stack, int i)
{
  if (  GetPremiaModel(stack,i) == NULL ) return NULL;
  return MaybeObjCopy(&NthObj(i));
}

NspPremiaModel  *GetPremiaModel(Stack stack, int i)
{
  NspPremiaModel *M;
  if (( M = nsp_premiamodel_object(NthObj(i))) == NULLPREMIAMODEL)
    ArgMessage(stack,i);
  return M;
}

/*-----------------------------------------------------
 * constructor 
 * if type is non NULL it is a subtype which can be used to 
 * create a NspClassB instance 
 *-----------------------------------------------------*/

static NspPremiaModel *premiamodel_create_void(char *name,NspTypeBase *type)
{
  NspPremiaModel *H  = (type == NULL) ? new_premiamodel() : type->new();
  if ( H ==  NULLPREMIAMODEL)
    {
      Sciprintf("No more memory\n");
      return NULLPREMIAMODEL;
    }
  if ( nsp_object_set_initial_name(NSP_OBJECT(H),name) == NULLSTRING) return NULLPREMIAMODEL;
  NSP_OBJECT(H)->ret_pos = -1 ;
  H->obj = NULL;
  return H;
}

NspPremiaModel *premiamodel_create(char *name,NspTypeBase *type)
{
  NspPremiaModel *H  = premiamodel_create_void(name,type);
  if ( H ==  NULLPREMIAMODEL) return NULLPREMIAMODEL;
  if ((H->obj = malloc(sizeof(nsp_premiamodel))) == NULL) return NULL;
  H->obj->ref_count=1;
  /* can be used to check that objects are non-initialized */
  H->obj->asset.name = NULL;
  H->obj->mod.TypeModel = NULL;
  H->obj->mod.Name = NULL;
  H->obj->opt.TypeOpt = NULL;
  H->obj->opt.Name= NULL;
  H->obj->meth.Name = NULL;
  H->obj->compat = FALSE;
  H->obj->it = NULL;
  H->obj->compute_err=-1;
  return H;
}

/*
 * copy for gobject derived class  
 */
NspPremiaModel *nsp_premiamodel_copy(NspPremiaModel *self)
{
  NspPremiaModel *H  =premiamodel_create_void(NVOID,(NspTypeBase *) nsp_type_premiamodel);
  if ( H ==  NULLPREMIAMODEL) return NULLPREMIAMODEL;
  H->obj = self->obj;
  self->obj->ref_count++;
  return H;
}

/*
 * Create a first object 
 */
int int_premiamodel_init (Stack stack, int rhs, int opt, int lhs)
{
  CheckStdRhs(0,0);
  /* want to be sure that type premiamodel is initialized */
  nsp_type_premiamodel = new_type_premiamodel(T_BASE);
  return 0;
}

/*
 * create a premiamodel object
 */
int int_premiamodel_create(Stack stack, int rhs, int opt, int lhs)
{
  NspPremiaModel *H;
  CheckStdRhs(0,0);
  /* want to be sure that type premiamodel is initialized */
  nsp_type_premiamodel = new_type_premiamodel(T_BASE);
  if(( H = premiamodel_create(NVOID,(NspTypeBase *) nsp_type_premiamodel)) == NULLPREMIAMODEL) 
    return RET_BUG;
  InitErrorMsg();
  InitVar();
  MoveObj(stack,1,(NspObject  *) H);
  return 1;
} 

int int_serial_test_unserialize(Stack stack, int rhs, int opt, int lhs) 
{
  NspObject *Obj;
  NspSerial *a;
  NspBMatrix *B;
  CheckRhs(1,1);
  CheckLhs(0,1);

  B = nsp_bmatrix_create (NVOID, 1, 1);
  if (( a= GetSerial(stack,1))== NULLSERIAL ) return RET_BUG;
  Obj = nsp_object_unserialize(a);
  if (Obj == NULLOBJ)
    {
      B->B[0] = FALSE;
    }
  else
    {
      if (nsp_object_set_name(Obj,NVOID) == FAIL) return RET_BUG;
      B->B[0] = TRUE;
    }
  MoveObj(stack,1,(NspObject *)B); 
  return 1;
}

/* -------------------------------- */
/* End of PremiaModel definition    */
/* -------------------------------- */


static int nsp_premia_set_var_names_forload(NspPremiaModel *self,VAR *vars,int n,NspList *L, int depth, void *obj)
{
  return nsp_premia_set_var_names_util(self,vars,n,L,depth,TRUE,obj);
}

static NspList* nsp_premia_get_var_names_forsave(NspPremiaModel *self, const VAR *vars,int n)
{
  return nsp_premia_get_var_names_util(self,vars,n, TRUE);
}
