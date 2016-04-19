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

/*
 * this file defines all the method implemented for a Premia
 * object
 */
#include "nsp/object.h"
#include "premiamod.h"
#include "nsp/interf.h"
#include "premia_vars.h"

#define CheckIsAssetSet(self)                                   \
  if ( self->obj->asset.name == NULL)                           \
    {                                                           \
      Scierror("Error: you must first set an asset type\n");    \
      return RET_BUG;                                           \
    }

#define CheckIsModelSet(self)                           \
  if ( self->obj->mod.TypeModel == NULL)                \
    {                                                   \
      Scierror("Error: you must first set a model\n");  \
      return RET_BUG;                                   \
    }

#define CheckIsOptionSet(self)                              \
  if ( self->obj->opt.TypeOpt == NULL)                      \
    {                                                       \
      Scierror("Error: you must first set an option\n");    \
      return RET_BUG;                                       \
    }

#define CheckIsMethodSet(self)                          \
  if ( self->obj->meth.Name == NULL)                    \
    {                                                   \
      Scierror("Error: you must first set a method\n"); \
      return RET_BUG;                                   \
    }

typedef enum { p_model, p_option, p_method_in, p_method_out } p_objs;

/* Get a list of (var-name,value,tag) 
 * where tag is true for variables which can be interactively changed.
 */
static int _wrap_premia_pb_get_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  NspList *L;
  int nvar;
  VAR *vars;
  CheckLhs(0,1);
  CheckRhs(0,0);
  switch (type ) 
    {
    case p_model : 
      vars = self->obj->mod.TypeModel;
      nvar = self->obj->mod.nvar;
      break;
    case p_option :
      vars = self->obj->opt.TypeOpt;
      nvar = self->obj->opt.nvar;
      break;
    case p_method_in: 
      vars=self->obj->meth.Par;
      /* nvar is detected with PREMIA_NULLTYPE tag */
      nvar = nsp_premia_get_nvar(vars);
      break;
    case p_method_out:
      if (self->obj->compute_err != OK)
        {
          if ((L=nsp_list_create(NVOID))==NULLLIST) return RET_BUG;
          MoveObj(stack,1,NSP_OBJECT(L));
          return Max(1,lhs);
        }
      vars=self->obj->meth.Res;
      /* nvar is detected with PREMIA_NULLTYPE tag */
      nvar = nsp_premia_get_nvar(vars);
      break;
    default:
      Scierror("Warning: to be done\n");
      return RET_BUG;
    }
  if ((L= nsp_premia_get_var_names(self,vars,nvar)) == NULL)
    return RET_BUG;
  MoveObj(stack,1,NSP_OBJECT(L));
  return Max(1,lhs);
}

static int _wrap_premiamodel_get_model_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  return _wrap_premia_pb_get_values(self,stack,rhs,opt,lhs,p_model);
}

static int _wrap_premiamodel_get_option_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  return _wrap_premia_pb_get_values(self,stack,rhs,opt,lhs,p_option);
}

static int _wrap_premiamodel_get_method_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  return _wrap_premia_pb_get_values(self,stack,rhs,opt,lhs,p_method_in);
}

static int _wrap_premiamodel_get_method_results(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  return _wrap_premia_pb_get_values(self,stack,rhs,opt,lhs,p_method_out);
}


static int _wrap_premiamodel_get_method_results_iter(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  NspObject *X, *Y, *Z;
  NspList *L;
  NspSMatrix *vars;
  CheckLhs(0,2);
  CheckRhs(0,0);
  if (self->obj->it == NULL)
    {
      Scierror("Error: while trying to get iteration results. No iteration structures\n");
      return RET_BUG;
    }
  if (self->obj->it->prix == NULLMAT)
    {
      Scierror("Error: while trying to get iteration results. Run compute first.\n");
      return RET_BUG;
    }
  if ((X=nsp_object_copy_and_name("X",(NspObject*)self->obj->it->range1))==NULLOBJ) return RET_BUG;
  if ((L=nsp_list_create("X")) == NULLLIST ) return RET_BUG;
  if (nsp_list_end_insert(L,X) == FAIL) return RET_BUG;
  if (self->obj->it->range2!=NULLMAT)
    {
      if ((Y=nsp_object_copy_and_name("Y",(NspObject*)self->obj->it->range2))==NULLOBJ) return RET_BUG;
      if (nsp_list_end_insert(L,Y) == FAIL) return RET_BUG;
      if ((vars=nsp_smatrix_create_with_length(NVOID,2,1,-1))==NULLSMAT) return RET_BUG;
      if ((vars->S[1] =nsp_string_copy(self->obj->it->location2->Vname)) == (nsp_string) 0) return RET_BUG; 
    }
  else
    {
      if ((vars=nsp_smatrix_create_with_length(NVOID,1,1,-1))==NULLSMAT) return RET_BUG;
    }
  if ((Z=nsp_object_copy_and_name("Z",(NspObject*)self->obj->it->prix))==NULLOBJ) return RET_BUG;
  if (nsp_list_end_insert(L,Z) == FAIL) return RET_BUG;
  if ((vars->S[0] =nsp_string_copy(self->obj->it->location1->Vname)) == (nsp_string) 0) return RET_BUG; 
  MoveObj(stack,1,NSP_OBJECT(L));
  MoveObj(stack,2,NSP_OBJECT(vars));  
  return Max(2, lhs);
  
}


static int _wrap_premia_pb_set_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  NspList *A, *B;
  VAR *vars;
  int nvar;
  char *msg = NULL;
  void *obj = NULL;
  CheckLhs(0,0);
  CheckRhs(1,1);
  if (( A=GetList(stack,1)) == NULLLIST) return RET_BUG;

  /* need to copy, because the list is deleted step by step
     in  nsp_premia_set_var_names. */
  B = nsp_list_copy(A);
  switch (type ) 
    {
    case p_model : 
      vars = self->obj->mod.TypeModel;
      nvar = self->obj->mod.nvar;
      obj = self->obj->mod.TypeModel;
      msg = "Error: while trying to set model values in a premia problem\n";
      break;
    case p_option :
      vars = self->obj->opt.TypeOpt;
      nvar = self->obj->opt.nvar;
      self->obj->opt.Init(&(self->obj->opt), &(self->obj->mod));
      obj = self->obj->opt.TypeOpt;
      msg = "Error: while trying to set option values in a premia problem\n";
      break;
    case p_method_in: 
      vars=self->obj->meth.Par;
      nvar = nsp_premia_get_nvar(vars);
      msg = "Error: while trying to set method parameters in a premia problem\n";
      break;
    case p_method_out: 
      vars=self->obj->meth.Res;
      nvar = nsp_premia_get_nvar(vars);
      msg = "Error: while trying to set method results in a premia problem\n";
      break;
    default:
      Scierror("Warning: to be done\n");
      return RET_BUG;
    }
  if(  nsp_premia_set_var_names(self,vars,nvar,B,0,obj) == FAIL) 
    {
      Scierror(msg);
      return RET_BUG;
    }

  return 0;
}

static int _wrap_premiamodel_set_option_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs) 
{
  return _wrap_premia_pb_set_values(self,stack,rhs,opt,lhs,p_option);
}

static int _wrap_premiamodel_set_model_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs) 
{
  if (self->obj->it!=NULL)
    {
      if (self->obj->it->range1!=NULLMAT) nsp_matrix_destroy(self->obj->it->range1);
      if (self->obj->it->range2!=NULLMAT) nsp_matrix_destroy(self->obj->it->range2);
      if (self->obj->it->prix != NULLMAT) nsp_matrix_destroy(self->obj->it->prix);
      free(self->obj->it);
      self->obj->it = NULL;
    }
  
  return _wrap_premia_pb_set_values(self,stack,rhs,opt,lhs,p_model);
}

static int _wrap_premiamodel_set_method_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs) 
{
  return _wrap_premia_pb_set_values(self,stack,rhs,opt,lhs,p_method_in);
}


/*
 * Check that parameters are correctly set. When a parameter is not a first
 * level parameter and is wrong, the sub variables are explored to discover
 * which parameter is wrong.
 */
static int _wrap_premia_pb_check_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,
                                        p_objs type)
{
  NspSMatrix *S;
  char *error,*format;
  int i,i_type,nvar;
  VAR *vars;
  CheckLhs(0,1);
  CheckRhs(0,0);

  switch (type ) 
    {
    case p_model : 
      vars = self->obj->mod.TypeModel;
      nvar = self->obj->mod.nvar;
      break;
    case p_option :
      vars = self->obj->opt.TypeOpt;
      nvar = self->obj->opt.nvar;
      break;
    case p_method_in: 
      vars=self->obj->meth.Par;
      nvar = nsp_premia_get_nvar(vars);
      break;
    case p_method_out: 
      vars=self->obj->meth.Res;
      nvar = nsp_premia_get_nvar(vars);
      break;
    default:
      Scierror("Warning: to be done\n");
      return RET_BUG;
    }
  for ( i = 0 ; i < nvar ; i++ )
    {
      VAR *vars1=NULL;
      int status =  (ChkVar1(NULL,&vars[i],WRONG) == 0);
      if ( status == FALSE ) 
        {
          /* If a variable is not a first level one, we try to 
           * get more precise informations 
           */
          switch( vars[i].Vtype)
            {
            case NUMFUNC_1:
              vars1 = (vars[i].Val.V_NUMFUNC_1)->Par;
              break;
            case NUMFUNC_2:
              vars1 = (vars[i].Val.V_NUMFUNC_2)->Par;
              break;
            case PTVAR:
              vars1 = (vars[i].Val.V_PTVAR)->Par;
              break;
            case PNLVECT: 
            default:
              break;
            }
          if ( vars1 != NULL) 
            {
              while (vars1->Vtype!=PREMIA_NULLTYPE)
                { 
                  int status1 =  (ChkVar1(NULL,vars1,WRONG) == 0);
                  if ( status1 == FALSE ) break;
                  vars1++;
                } 
            }
          else 
            {
              vars1 = &vars[i];
            }
          premia_Vtype_info(vars1,&format,&error,&i_type);
          /* Scierror("Error: %s is wrong, %s\n",vars1->Vname,error); */
          if ((S=nsp_smatrix_create_with_length(NVOID,1,2,-1))== NULLSMAT) return RET_BUG;
          if ((S->S[0] =nsp_string_copy(vars1->Vname)) == (nsp_string) 0) return RET_BUG;
          if ((S->S[1] =nsp_string_copy(error)) == (nsp_string) 0) return RET_BUG;
          MoveObj(stack,1,NSP_OBJECT(S));
          return 1;
        }
    }
  if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
    return RET_BUG;
  MoveObj(stack,1,NSP_OBJECT(S));
  return 1;
}

static int _wrap_premiamodel_check_model_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  return _wrap_premia_pb_check_values(self,stack,rhs,opt,lhs,p_model);
}

static int _wrap_premiamodel_check_option_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  return _wrap_premia_pb_check_values(self,stack,rhs,opt,lhs,p_option);
}

static int _wrap_premiamodel_check_method_values(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  return _wrap_premia_pb_check_values(self,stack,rhs,opt,lhs,p_method_in);
}


/*
 * set the type of asset using the syntax asset=asset_name.
 * asset_name can be as defined in premia_asset_names variable
 */
static int _wrap_premiamodel_set_asset(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  char *type=NULL;
  int ret;
  nsp_option set_model_opts[] = {{ "str",string,NULLOBJ,-1},
                                 { NULL,t_end,NULLOBJ,-1}};
  
  CheckStdRhs(0,0);
  CheckOptRhs(1,1);
  /* asset type */
  if ( get_optional_args(stack,rhs,opt,set_model_opts, &type) == FAIL)
    return RET_BUG;

  if ((ret=nsp_premiamodel_set_asset (self, type))!=OK) return ret;

  premia_free_attr (self, p_mod);
  premia_free_attr (self, p_opt);
  premia_free_attr (self, p_meth);
  return 0;
}

/*
 * get the asset type if already set
 */
static int _wrap_premiamodel_get_asset(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  NspSMatrix *S;
  CheckStdRhs(0,0);
  if ( self->obj->asset.name == NULL ) 
    {
      if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
        return RET_BUG;
      MoveObj(stack,1,(NspObject  *) S);
      return 1;
    }
  if ( nsp_move_string(stack,1,self->obj->asset.name,-1) == FAIL)
    return RET_BUG;
  return 1;
}

/*
 * model can be chosen using a number or a string in this
 * latter case the syntax is str=model_name where model_name
 * is given in the models variable
 */
static int _wrap_premiamodel_set_model(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  int model,ret;
  char *str=NULL;
  nsp_option set_model_opts[] = {{ "str",string,NULLOBJ,-1},
                                 { NULL,t_end,NULLOBJ,-1}};
  
  CheckStdRhs(0,1);

  CheckIsAssetSet (self);

  /* Free previously set model, option, method */
  premia_free_attr (self, p_mod);
  premia_free_attr (self, p_opt);
  premia_free_attr (self, p_meth);
  
  /* set model called with a named argument */
  if (opt==1)
    {
      if ( get_optional_args(stack,rhs,opt,set_model_opts, &str) == FAIL)
        return RET_BUG;
      if (str==NULL) return RET_BUG; /* no parameters found */
      if ((ret=nsp_premiamodel_set_model_with_str (self, str))!=OK) return ret;
      return 0;
    }

  /* set model called with an integer argument : the model number */
  if (GetScalarInt(stack,1,&model) == FAIL) return RET_BUG;
  model--;
  if ((ret=nsp_premiamodel_set_model (self, model))!=OK) return ret;
  return 0;
} 

static int _wrap_premiamodel_get_model(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  NspSMatrix *S;
  CheckStdRhs(0,0);
  if ( self->obj->mod.TypeModel == NULL ) 
    {
      if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
        return RET_BUG;
      MoveObj(stack,1,(NspObject  *) S);
      return 1;
    }
  if ( nsp_move_string(stack,1,self->obj->mod.Name,-1) == FAIL) return RET_BUG;
  return 1;
} 

/** returns the available models for the asset type already
 *   set.
 *  
 *  @returns RET_BUG if asset_type not set
 */
static int _wrap_premiamodel_get_models(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  PremiaAsset asset;
  NspSMatrix *S;
  CheckStdRhs(0,0);
  asset = self->obj->asset;
  if ( asset.name == NULL ) 
    {
      if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
        return RET_BUG;
      MoveObj(stack,1,(NspObject  *) S);
      return 1;
    }

  if ((S = premia_get_models_aux (asset)) == NULL) return RET_BUG;
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
} 

/*
 * Sets the option either using a Name (syntax str=Name) or using two
 * numbers. Note that these numbers are the indices of the familiy and option
 * within the restricted lists given the model choice
 */
static int _wrap_premiamodel_set_option(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  Option **loc;
  Family **families;
  Pricing **pricings;
  int option,family,i_family=-1,i_option=-1,ret,nf=0,m;
  int count_opt,count_fam;
  char *str=NULL;
  nsp_option set_option_opts[] = {{ "str",string,NULLOBJ,-1},
                                  { NULL,t_end,NULLOBJ,-1}};
  CheckStdRhs(0,2);

  CheckIsAssetSet (self);
  CheckIsModelSet (self);

  premia_free_attr (self, p_opt);
  premia_free_attr (self, p_meth);
  
  /* set_option is called using a named argument */
  if (opt==1)
    {
      if ( get_optional_args(stack,rhs,opt,set_option_opts, &str) == FAIL) return RET_BUG;
      if (str==NULL) return RET_BUG; /* no parameters found */
      if ((ret=nsp_premiamodel_set_option_with_str (self, str))!=OK) return ret;
      return 0;
    }

  /* set_option is called using two integer arguments : the family and
     option numbers */
  if (GetScalarInt(stack,1,&family) == FAIL) return RET_BUG;
  if (GetScalarInt(stack,2,&option) == FAIL) return RET_BUG;
  family--;
  option--;

  families = self->obj->asset.families;
  pricings = self->obj->asset.pricings;

  while ( families[nf] != NULL) nf++;
  count_fam=-1;
  for ( m=0; m<nf; m++ )
    {
      int i, fsize = 0;
      loc = (*families[m]);
      while (loc[fsize] != NULL) fsize++;

      count_opt=0;
      /* if n is given we check if the family is compatible with 
       * the model. 
       */
      for ( i=0 ; i < fsize ; i++) 
        {
          if (Premia_match_model_option(&(self->obj->mod), loc[i], pricings)==0)
            {
              if (count_opt == 0) count_fam++;
              if (count_fam == family && count_opt == option)
                {
                  i_option = i; i_family = m; break;
                }
              count_opt++;
            }
        }
    }
  
  if ((ret=nsp_premiamodel_set_option(self, i_family, i_option))!=OK) return ret;
  return 0;
} 

static int _wrap_premiamodel_get_option(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  NspSMatrix *S;
  CheckStdRhs(0,0);
  if ( self->obj->opt.TypeOpt == NULL ) 
    {
      if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
        return RET_BUG;
      MoveObj(stack,1,(NspObject  *) S);
      return 1;
    }
  if ( nsp_move_string(stack,1,self->obj->opt.Name,-1) == FAIL) return RET_BUG;
  return 1;
} 

static int _wrap_premiamodel_get_family(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  NspSMatrix *S;
  int index;
  CheckStdRhs(1,1);

 
  if ( self->obj->asset.name == NULL || self->obj->mod.TypeModel == NULL )
    {
      if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
        return RET_BUG;
      MoveObj(stack,1,(NspObject  *) S);
      return 1;
    }

  if (GetScalarInt(stack,1,&index) == FAIL) return RET_BUG;
  index--;

  if ((S=premia_get_family_aux(self->obj->asset, &(self->obj->mod), index))==NULL) return RET_BUG;
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
  
}

static int _wrap_premiamodel_set_method(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  int method,ret;
  char *str=NULL;
  nsp_option set_method_opts[] = {{ "str",string,NULLOBJ,-1},
                                  { NULL,t_end,NULLOBJ,-1}};
  CheckStdRhs(0,1);
  
  CheckIsAssetSet (self);
  CheckIsModelSet (self);
  CheckIsOptionSet (self);

  /* set_method called using the named argument str*/
  if (opt==1)
    {
      if ( get_optional_args(stack,rhs,opt,set_method_opts, &str) == FAIL)
        return RET_BUG;
      if (str==NULL) return RET_BUG; /* no parameters found */
      if ((ret=nsp_premiamodel_set_method_with_str (self, str))!=OK) return ret;
      return 0;
    }
  /* set_method is called using an integer argument : the method number */
  if (GetScalarInt(stack,1,&method) == FAIL) return RET_BUG;
  method--;
  
  if ((ret=nsp_premiamodel_set_method (self, method))!=OK) return ret;
  return 0;
} 

static int _wrap_premiamodel_get_method(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  NspSMatrix *S;
  CheckStdRhs(0,0);
  if ( self->obj->meth.Name == NULL ) 
    {
      if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
        return RET_BUG;
      MoveObj(stack,1,(NspObject  *) S);
      return 1;
      Scierror("Error: method is not set\n");
      return RET_BUG;
    }
  if ( nsp_move_string(stack,1,self->obj->meth.Name,-1) == FAIL) return RET_BUG;
  return 1;
} 

static int _wrap_premiamodel_get_methods(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  NspSMatrix *S;
  CheckStdRhs(0,0);

  CheckIsAssetSet (self);
  CheckIsModelSet (self);
  CheckIsOptionSet (self);

  if ((S = premia_get_methods_aux (self->obj->asset,&(self->obj->mod),&(self->obj->opt))) == NULL) return RET_BUG;
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
  
}

static int _wrap_premiamodel_compute(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  int user=NO_PAR, i, j;
  double *ptr1=NULL, *ptr2=NULL;
  Planning pt_plan;
  Pricing **pricings, *res;
  CheckStdRhs(0,0);

  CheckIsAssetSet (self);
  CheckIsModelSet (self);
  CheckIsOptionSet (self);
  CheckIsMethodSet (self);


  if ((self->obj->mod.Check)(user,&pt_plan,&self->obj->mod) != OK ) 
    {
      Scierror("Error: model is not correct\n");
      return RET_BUG;
    }
  if ((self->obj->opt.Check)(user,&pt_plan,&self->obj->opt) !=OK)
    {
      Scierror("Error: option is not correct\n");
      return RET_BUG;
    }
  /* Il faut ici un pricing */
  pricings = self->obj->asset.pricings;
  if (SelectPricing(0,&self->obj->mod,&self->obj->opt,pricings,&res) != OK) 
    {
      Scierror("Error: option and model have incompatibilities\n");
      return RET_BUG;
    }

  /* only one computation */
  if (self->obj->it == NULL)
    {
      
      /* force initialisation of method */
      (self->obj->mod.Init)(&self->obj->mod);
      (self->obj->opt.Init)(&self->obj->opt, &self->obj->mod);
      (self->obj->meth.Init)(&self->obj->meth, &self->obj->opt);
      if ((self->obj->meth.Check)(user,&pt_plan,&self->obj->meth)!=OK)
        {
          Scierror("Error: method is not correct\n");
          return RET_BUG;
        }
      self->obj->compute_err=(self->obj->meth.Compute)(self->obj->opt.TypeOpt,self->obj->mod.TypeModel,&self->obj->meth);
      return 0;
    }

  /* iteration is required */
  if (self->obj->it->range1 == NULLMAT) return FAIL;
  ptr1 = self->obj->it->range1->R;
  if (self->obj->it->range2 == NULLMAT)
    {
      /* 1d iteration */
      if (self->obj->it->prix == NULLMAT)
        {
          if ((self->obj->it->prix = nsp_matrix_create("X", 'r', 1,self->obj->it->range1->mn))== NULLMAT) return FAIL;
        }
      else
        if (nsp_matrix_resize(self->obj->it->prix, self->obj->it->range1->mn,1) == FAIL) return FAIL;
      for (i=0; i<self->obj->it->range1->mn; i++, ptr1++)
        {
          _nsp_premia_set_value(self->obj->it->location1,*ptr1);
          (self->obj->mod.Init)(&self->obj->mod);
          (self->obj->opt.Init)(&self->obj->opt, &self->obj->mod);
          (self->obj->meth.Init)(&self->obj->meth, &self->obj->opt);
          if ((self->obj->meth.Check)(user,&pt_plan,&self->obj->meth)!=OK)
            {
              Scierror("Error: method is not correct\n");
              return RET_BUG;
            }
          self->obj->compute_err=(self->obj->meth.Compute)(self->obj->opt.TypeOpt,self->obj->mod.TypeModel,&self->obj->meth);
          if (self->obj->compute_err != OK) return 0;
          self->obj->it->prix->R[i] = self->obj->meth.Res[0].Val.V_DOUBLE;
        }
      return 0;
    }

  /* 2d iteration */
  if (self->obj->it->prix == NULLMAT)
    {
      if ((self->obj->it->prix = nsp_matrix_create("X", 'r', self->obj->it->range1->mn,self->obj->it->range2->mn))== NULLMAT) return FAIL;
    }
  else
    if (nsp_matrix_resize(self->obj->it->prix, self->obj->it->range1->mn,self->obj->it->range2->mn) == FAIL) return FAIL;
  for (i=0; i<self->obj->it->range1->mn; i++, ptr1++)
    {
      ptr2 = self->obj->it->range2->R;
      _nsp_premia_set_value(self->obj->it->location1,*ptr1);
      for (j=0; j<self->obj->it->range2->mn; j++, ptr2++)
        {
          _nsp_premia_set_value(self->obj->it->location2,*ptr2);
          (self->obj->mod.Init)(&self->obj->mod);
          (self->obj->opt.Init)(&self->obj->opt, &self->obj->mod);
          (self->obj->meth.Init)(&self->obj->meth, &self->obj->opt);
          if ((self->obj->meth.Check)(user,&pt_plan,&self->obj->meth)!=OK)
            {
              Scierror("Error: method is not correct\n");
              return RET_BUG;
            }
          self->obj->compute_err=(self->obj->meth.Compute)(self->obj->opt.TypeOpt,self->obj->mod.TypeModel,&self->obj->meth);
          if (self->obj->compute_err != OK) return 0;
          self->obj->it->prix->R[j*self->obj->it->range1->mn+i] = self->obj->meth.Res[0].Val.V_DOUBLE;
        }
    }
  return 0;
}

/*
 * If the computation did not end successfuly, this method returns the error
 * message. When everything was OK, an empty Smatrix is returned.
 */
static int _wrap_premiamodel_get_compute_err(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs)
{
  CheckStdRhs(0,1);
  NspSMatrix *S;
  if (self->obj->compute_err>=2)
    {
      if ((S =nsp_smatrix_create (NVOID, 1, 1, error_msg[self->obj->compute_err], 1))
          == NULLSMAT) return RET_BUG;
    }
  else
    {
      if ((S =nsp_smatrix_create (NVOID, 0, 0, "", -1)) == NULLSMAT) return RET_BUG;
    }
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
}

static int _wrap_premia_get_help_model(NspPremiaModel *self);
static int _wrap_premia_get_help_option(NspPremiaModel *self);
static int _wrap_premia_get_help_method(NspPremiaModel *self);

/** _wrap_premia_get_help
 * displays the help concerning the model, the option, the
 * pricing method depending on the string argument
 */
int _wrap_premia_get_help(NspPremiaModel *self, Stack stack, int rhs, int opt, int lhs)
{
  char *type;
  CheckRhs(1,1);
  CheckLhs(0,0);

  /* make sure all the paths are properly initialized */
  premia_set_global_vars();
  
  if ((type=GetString(stack,1)) == ((char *) 0)) return RET_BUG;
  if (strcasecmp(type, "model")==0)
    return _wrap_premia_get_help_model(self);
  else if (strcasecmp(type, "option")==0)
    return _wrap_premia_get_help_option(self);
  else if (strcasecmp(type, "method")==0)
    return _wrap_premia_get_help_method(self);
  else
    Scierror("only accepts \"model\", \"option\", or \"method\"\n");

  return 0;
}


/*
 * displays the help concerning the model
 */
static int _wrap_premia_get_help_model(NspPremiaModel *self)
{
  char fhelp_name[MAX_PATH_LEN];
  
  if ( self->obj->mod.TypeModel == NULL )
    {
      if ((strlen(premiamandir)+strlen(path_sep)+strlen("mod_doc.pdf"))>=MAX_PATH_LEN)
        Scierror("path too long\n");
      strcpy(fhelp_name,premiamandir);
      strcat(fhelp_name,path_sep);
      strcat(fhelp_name,"mod_doc.pdf");
    }
  else
    {
      get_model_helpfile (&(self->obj->mod), fhelp_name);
      if ( fhelp_name[0] == '\0' ) Scierror("path too long\n");
    }
  premia_spawnlp(fhelp_name);

  return 0;
}

/*
 * displays the help concerning the option
 */
static int _wrap_premia_get_help_option(NspPremiaModel *self)
{
  char fhelp_name[MAX_PATH_LEN];
  

  if ( self->obj->opt.TypeOpt == NULL )
    {
      if ((strlen(premiamandir)+strlen(path_sep)+strlen("opt_doc.pdf"))>=MAX_PATH_LEN)
        Scierror("path too long\n");
      strcpy(fhelp_name,premiamandir);
      strcat(fhelp_name,path_sep);
      strcat(fhelp_name,"opt_doc.pdf");
    }
  else
    {
      get_option_helpfile (&(self->obj->opt), fhelp_name);
      if ( fhelp_name[0] == '\0' ) Scierror("path too long\n");
    }
  premia_spawnlp(fhelp_name);

  return 0;
}

/*
 * displays the help concerning the method
 */
static int _wrap_premia_get_help_method(NspPremiaModel *self)
{
  char fhelp_name[MAX_PATH_LEN];
  
  if ( self->obj->meth.Name == NULL )
    {
      Scierror("No help available : choose a method first.\n");
    }
  else
    {
      get_method_helpfile_with_ids (&(self->obj->meth), self->obj->mod.ID,self->obj->opt.ID, fhelp_name);
      if ( fhelp_name[0] == '\0' ) Scierror("path too long\n");
    }
  premia_spawnlp(fhelp_name);

  return 0;
}


static int _wrap_premiamodel_is_with_iter(NspPremiaModel *self,Stack stack,int rhs,int opt,int lhs,p_objs type)
{
  NspBMatrix *M;
  CheckRhs(0,0);
  CheckLhs(1,1);
  if ((M=nsp_bmatrix_create("X", 1, 1))==NULLBMAT) return RET_BUG;
  if (self->obj->it == NULL) M->B[0] = 0;
  else M->B[0] = 1;
  MoveObj(stack,1,(NspObject  *) M);
  return 1;
}


static NspMethods premiamodel_methods[] = {
  {"set_asset",(nsp_method *) _wrap_premiamodel_set_asset},
  {"get_asset",(nsp_method *) _wrap_premiamodel_get_asset},  
  {"set_model",(nsp_method *) _wrap_premiamodel_set_model},
  {"get_model",(nsp_method *) _wrap_premiamodel_get_model},
  {"get_models",(nsp_method *) _wrap_premiamodel_get_models},
  {"set_option",(nsp_method *) _wrap_premiamodel_set_option},
  {"get_option",(nsp_method *) _wrap_premiamodel_get_option},
  {"get_family",(nsp_method *) _wrap_premiamodel_get_family},
  {"set_method",(nsp_method *) _wrap_premiamodel_set_method},
  {"get_method",(nsp_method *) _wrap_premiamodel_get_method},
  {"get_methods",(nsp_method *) _wrap_premiamodel_get_methods},

  {"get_model_values",(nsp_method *) _wrap_premiamodel_get_model_values},
  {"set_model_values",(nsp_method *) _wrap_premiamodel_set_model_values},
  {"model_check",(nsp_method *) _wrap_premiamodel_check_model_values},

  {"get_option_values",(nsp_method *) _wrap_premiamodel_get_option_values},
  {"set_option_values",(nsp_method *) _wrap_premiamodel_set_option_values},
  {"option_check",(nsp_method *) _wrap_premiamodel_check_option_values},

  {"get_method_values",(nsp_method *) _wrap_premiamodel_get_method_values},
  {"set_method_values",(nsp_method *) _wrap_premiamodel_set_method_values},
  {"method_check",(nsp_method *) _wrap_premiamodel_check_method_values},

  {"get_method_results",(nsp_method *) _wrap_premiamodel_get_method_results},
  {"get_method_results_iter",(nsp_method *) _wrap_premiamodel_get_method_results_iter},
  {"is_with_iter",(nsp_method *) _wrap_premiamodel_is_with_iter},

  {"compute",(nsp_method *) _wrap_premiamodel_compute},
  {"get_compute_err",(nsp_method *) _wrap_premiamodel_get_compute_err},

  {"get_help", (nsp_method *) _wrap_premia_get_help},
  { NULL, NULL}
};

NspMethods *premiamodel_get_methods(void) { return premiamodel_methods;};
/*-------------------------------------------
 * Attributes
 *-------------------------------------------*/

AttrTab premiamodel_attrs[] = {
  { NULL,NULL,NULL,NULL },
};


#undef CheckIsAssetSet
#undef CheckIsModelSet
#undef CheckIsOptionSet
#undef CheckIsMethodSet
