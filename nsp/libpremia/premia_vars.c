#include "nsp/interf.h"
#include "nsp/list.h"
#include "nsp/nsptcl.h"
#include "nsp/dstring.h"
#include "premia_vars.h"
#include "premiamod.h"


/*
 * returns a String Matrix containing the names of all the
 * asset types available
 */
int int_premia_get_assets(Stack stack, int rhs, int opt, int lhs)
{
  int i;
  NspSMatrix *S;
  int nassets=0;
  CheckStdRhs(0,0);
  CheckOptRhs(0,1);
  
  while (premia_assets[nassets].name != NULL) nassets++;
  if ((S=nsp_smatrix_create_with_length(NVOID,nassets,1,-1))== NULLSMAT) 
    return RET_BUG;
  for ( i=0 ; i < nassets ; i++) 
    if ((S->S[ i] =nsp_string_copy(premia_assets[i].name)) == (nsp_string) 0) 
      return RET_BUG;
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
}

/**
 * prints the available asset types according to
 * @premia_asset_names
 */
void premia_accepted_assets()
{
  PremiaAsset *dummy_asset=premia_assets;
  while (dummy_asset->name!=NULL)
    {
      Scierror("%s, ", dummy_asset->name); 
      dummy_asset++;
    }
}

/**
 * determine from the optionnal argument the asset type
 * amongst the ones defined by @premia_asset_names.
 *
 * this function is not available at Nsp level.
 */
int premia_get_asset_type(Stack stack, int rhs, int opt, int lhs, PremiaAsset *asset)
{
  char *asset_type = NULL;
  PremiaAsset *dummy_asset=premia_assets;
  
  nsp_option asset_type_opts[] = {{"asset",string,NULLOBJ,-1},
                                  { NULL,t_end,NULLOBJ,-1}};

  if ( get_optional_args(stack,rhs,opt, asset_type_opts,&asset_type) == FAIL)
    {
      Scierror("argument must be given in the form asset=str where str \n can only be ");
      premia_accepted_assets();
      Scierror("\n");
      return FAIL;
    }
  if ( asset_type == NULL)
    {
      *asset = *premia_assets; return OK;
    }

  while (dummy_asset->name!=NULL)
    {
      if (strcmp(asset_type, dummy_asset->name)==0)
        {
          *asset = *dummy_asset; return OK;
        }
      else
        dummy_asset++;
    }
  
  Scierror ("asset type unknown, accepted values are : ");
  premia_accepted_assets();
  Scierror("\n");
  return FAIL;
}

/*
 * Return the model which has index n in the restricted list of models having
 * some options availble for pricing
 */
Model* premia_get_model_from_index (PremiaAsset asset, int n)
{
  Model *model;
  Family **families;
  Model **models;
  Pricing **pricings;
  int nmodels=0, nmod=0;

  models = asset.models;
  families = asset.families;
  pricings = asset.pricings;
  model = NULL;

  
  while (models[nmodels] != NULL)
    {
      if (Premia_model_has_products (models[nmodels], families, pricings) == OK)
        {
          if (n == nmod) { model = models[nmodels]; break; }
          nmod++;
        }
      nmodels++;
    }
  return model;
}

/*
 * Returns the list of all the models available for the given asset as a SMat
 */
NspSMatrix* premia_get_models_aux (PremiaAsset asset)
{
  int i, nmodels, count_mod;
  NspSMatrix *S;
  Model **models;
  
  nmodels=0;
  count_mod=0;
    
  models = asset.models;
  
  while (models[nmodels] != NULL) nmodels++;
  if ((S=nsp_smatrix_create_with_length(NVOID,nmodels,1,-1))== NULLSMAT) return NULL;
  for (i=0 ; i<nmodels ; i++)
    {
      if (Premia_model_has_products(models[i], asset.families, asset.pricings) == OK)
        {
          if ((S->S[count_mod]=nsp_string_copy(models[i]->Name)) == (nsp_string) 0) 
            return NULL;
          count_mod++;
        }
    }
  if (nsp_smatrix_resize(S, count_mod, 1) == FAIL) return NULL;
  return S;
}

/*
 * Returns the list of available methods to price a given product
 */
NspSMatrix* premia_get_methods_aux (PremiaAsset asset, Model *mod, Option *opt)
{
  NspSMatrix *S;
  Pricing **pricings, *res;
  int npm=0, npm_ok=0, i, status;

  pricings = asset.pricings;
  status = SelectPricing(0, mod, opt, pricings, &res);
  if (status == PREMIA_NONE) goto empty;
  /* walk on the selected pricing and check correct methods */
  while ( res->Methods[npm] != NULL) npm++;
  npm_ok=0;
  for ( i=0 ; i < npm ; i++) 
    {
      if ( res->Methods[i]->CheckOpt(opt, mod) == OK) npm_ok++;
    }
  if ((S=nsp_smatrix_create_with_length(NVOID,npm_ok,1,-1))== NULLSMAT) 
    return NULL;
  npm_ok=0;
  for ( i=0 ; i < npm ; i++) 
    {
      if ( res->Methods[i]->CheckOpt(opt, mod) == OK) 
        {
          if ((S->S[npm_ok] =nsp_string_copy(res->Methods[i]->Name)) == (nsp_string) 0) 
            return NULL;
          npm_ok++;
        }
    }
  return S;
 empty:
  if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) return NULL;
  return S;
}

/*
 * Returns the list of available options in family n and in model mod. If
 * mod==NULL, all the options of family n are returned
 */
NspSMatrix* premia_get_family_aux(PremiaAsset asset, Model *mod, int n)
{
  NspSMatrix *S;
  Option **loc;
  int i,nf=0,fsize=0,count=0;
  Family **families;
  Pricing **pricings;

  families = asset.families;
  pricings = asset.pricings;

  /* count available families */
  while ( families[nf] != NULL) nf++;
  if ( n < 0 || n > nf -1 ) 
    {
      if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
        return NULL;
      return S;
    }

  /* now, we really find the available options for the model
     already set and the given family n*/
  loc = (*families[n]);
  while ( loc[fsize] != NULL) fsize++;
  if ((S=nsp_smatrix_create_with_length(NVOID,fsize,1,-1))== NULLSMAT) 
    return NULL;
  /* we check if the family is compatible with  the model. */
  for ( i=0 ; i < fsize ; i++) 
    {
      if (mod==NULL || Premia_match_model_option(mod, loc[i], pricings)==0)
        {
          if ((S->S[count] =nsp_string_copy(loc[i]->Name)) == (nsp_string) 0) 
            return NULL;
          count++;
        }
    }
  if (count < fsize) { nsp_smatrix_resize(S, count, 1); }
  return S;  
}

/*
 * Returns the list of all the models available (for the given asset if
 * specified with asset="poo") as a SMat
 */
int int_premia_get_models(Stack stack, int rhs, int opt, int lhs)
{
  NspSMatrix *S;
  PremiaAsset asset;
  
  CheckStdRhs(0,0);
  CheckOptRhs(0,1);
    
  if (premia_get_asset_type (stack, rhs, opt, lhs, &asset) == FAIL) return RET_BUG;

  if ((S = premia_get_models_aux (asset)) == NULL) return RET_BUG;
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
}

/*
 * Gets all the family names available in Premia. If called
 * with a second optional argument, the name of the family
 * whose index is given as optional argument is returned
 */
int int_premia_get_families(Stack stack, int rhs, int opt, int lhs)
{
  int nf=0;
  int index=0;
  PremiaAsset asset;
  Family **families;
  NspSMatrix *S;
  CheckStdRhs(0,1);
  CheckOptRhs(0,1);

  if (premia_get_asset_type (stack, rhs, opt, lhs, &asset) == FAIL) return RET_BUG;
  families =  asset.families;

  while ( families[nf] != NULL) nf++;

  if (rhs-opt==1)
    {
      /* a family index is given*/
      if (GetScalarInt(stack,1,&index) == FAIL) return RET_BUG;
      index--;
      if (index < 0 || index > nf-1)
        {
          Scierror("Error: family %d does not exist\n",index+1);
          return RET_BUG;
        }
      if ((S=nsp_smatrix_create_with_length(NVOID,1,1,-1))== NULLSMAT) return RET_BUG;
      if ((S->S[0] = nsp_string_copy((*families[index][0])->ID)) == (nsp_string) 0) 
        return RET_BUG;
    }
  else
    {
      if ((S=nsp_smatrix_create_with_length(NVOID,nf,1,-1))== NULLSMAT) return RET_BUG;
      for (nf=0; families[nf] != NULL; nf++)
        {
          if ((S->S[nf] =nsp_string_copy((*families[nf][0])->ID)) == (nsp_string) 0) 
            return RET_BUG;
        }
    }
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
}
  
/*
 * Gets the options in family m if a second optional argument is given the
 * subset of the family compatible with model n is returned. The asset type
 * considered can be specified using asset="poo". The index of the model if
 * given is to be considered in the restricted list of available models.
 */
int int_premia_get_family(Stack stack, int rhs, int opt, int lhs)
{
  Model *poo = NULL;
  NspSMatrix *S;
  int m,n=-1;
  PremiaAsset asset;
  CheckStdRhs(1,2);
  CheckOptRhs(0,1);
    

  /* determine the asset type */
  if (premia_get_asset_type (stack, rhs, opt, lhs, &asset) == FAIL) return RET_BUG;
  
  /* the family */
  if (GetScalarInt(stack,1,&m) == FAIL) return RET_BUG;
  m--; 
  if ( rhs-opt == 2 )
    {
      /* a model id is given too  */
      if (GetScalarInt(stack,2,&n) == FAIL) return RET_BUG;
      n--;
      if ((poo = premia_get_model_from_index (asset, n)) == NULL)
        {
          Scierror("Error: model %d does not exist\n",n+1);
          return RET_BUG;
        }
    }
  if ((S=premia_get_family_aux(asset, poo, m))==NULL) return RET_BUG;
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
}

/* 
 * get possible methods. Takes 3 arguments (integers)
 * - family : index of the family
 * - option : index of the option whithin the list of available options in the
 * given model and family
 * - model : index of the model.
 */
int int_premia_get_methods(Stack stack, int rhs, int opt, int lhs)
{
  Model *model;
  Option **loc, *option;
  NspSMatrix *S;
  PremiaAsset asset;
  Family **families;
  Pricing **pricings;
  int nf=0,i,m,n=-1,nopt=0;
  CheckOptRhs(0,1);
  CheckStdRhs(3,3);

  model = NULL;
  option = NULL;
  /* determine the asset type */
  if (premia_get_asset_type (stack, rhs, opt, lhs, &asset) == FAIL) return RET_BUG;
  families = asset.families;
  pricings = asset.pricings;

  
  /* the family */
  if (GetScalarInt(stack,1,&m) == FAIL) return RET_BUG;
  m--;
  /* the option */
  if (GetScalarInt(stack,2,&opt) == FAIL) return RET_BUG;
  opt--;
  /* the model */
  if (GetScalarInt(stack,3,&n) == FAIL) return RET_BUG;
  n--;
  if ((model = premia_get_model_from_index (asset, n)) == NULL)
    {
      Scierror("Error: model %d does not exist\n",n+1);
      return RET_BUG;
    }
 
  model->Init(model);
  while ( families[nf] != NULL) nf++;
  if ( m < 0 || m > nf -1 ) goto empty;
  loc = (*families[m]);

  /* find the option in the restricted list of available options */
  nopt = 0; i=0;
  while (loc[i] != NULL) 
    {
      if (Premia_match_model_option(model, loc[i], pricings)==0)
        {
          if (nopt==opt) { option = loc[i]; break; }
          nopt++;
        }
      i++;
    }
  if (option==NULL) goto empty;

  
  /* select the first matching pricing */
  InitErrorMsg();
  InitVar();
  option->Init(option, model);

  if ((S = premia_get_methods_aux (asset, model, option)) == NULL) return RET_BUG;
  MoveObj(stack,1,(NspObject  *) S);
  return 1;

 empty:
  if ((S=nsp_smatrix_create_with_length(NVOID,0,0,-1))== NULLSMAT) 
    return RET_BUG;
  MoveObj(stack,1,(NspObject  *) S);
  return 1;
}

/* 
 * This function is used to clone an array of vars of size n.
 */
int nsp_premia_clone_vars(VAR **res,int flag,const VAR *vars,int n)
{
  int i;
  VAR *loc,*loc2; 
  if ( flag == TRUE ) 
    {
      if ((loc = malloc(n*sizeof(VAR))) == NULL) 
        return FAIL;
      *res = loc;
    }
  else 
    {
      loc = *res;
    }
  /* now we have to check if recursive allocation is needed */
  for (i=0 ; i < n ; i++)
    {
      int count =0;
      loc[i]=vars[i];
      switch( vars[i].Vtype)
        {
        case NUMFUNC_1:
          loc2 = (vars[i].Val.V_NUMFUNC_1)->Par; 
          /* count how many vars are present in vars[i] */
          while (loc2->Vtype!=PREMIA_NULLTYPE) { count++; loc2++;} 
          /* allocate */
          if ((loc[i].Val.V_NUMFUNC_1 = malloc(sizeof(NumFunc_1)))==NULL) return FAIL;
          *(loc[i].Val.V_NUMFUNC_1) = *(vars[i].Val.V_NUMFUNC_1);
          /* recursive allocation */
          loc2 =loc[i].Val.V_NUMFUNC_1->Par;
          nsp_premia_clone_vars(&loc2,FALSE,vars[i].Val.V_NUMFUNC_1->Par,count);
          break;
        case NUMFUNC_2:
          loc2 = (vars[i].Val.V_NUMFUNC_2)->Par;
          while (loc2->Vtype!=PREMIA_NULLTYPE) { count++; loc2++;} 
          /* allocate */
          if ((loc[i].Val.V_NUMFUNC_2 = malloc(sizeof(NumFunc_2)))==NULL) return FAIL;
          *(loc[i].Val.V_NUMFUNC_2) = *(vars[i].Val.V_NUMFUNC_2);
          /* recursive allocation */
          loc2 =loc[i].Val.V_NUMFUNC_2->Par;
          nsp_premia_clone_vars(&loc2,FALSE, vars[i].Val.V_NUMFUNC_2->Par,count);
          break;
        case NUMFUNC_ND:
          loc2 = (vars[i].Val.V_NUMFUNC_ND)->Par;
          while (loc2->Vtype!=PREMIA_NULLTYPE) { count++; loc2++;} 
          /* allocate */
          if ((loc[i].Val.V_NUMFUNC_ND = malloc(sizeof(NumFunc_nd)))==NULL) return FAIL;
          *(loc[i].Val.V_NUMFUNC_ND) = *(vars[i].Val.V_NUMFUNC_ND);
          /* recursive allocation */
          loc2 =loc[i].Val.V_NUMFUNC_ND->Par;
          nsp_premia_clone_vars(&loc2,FALSE, vars[i].Val.V_NUMFUNC_ND->Par,count);
          break;
        case PTVAR:
          loc2 = (vars[i].Val.V_PTVAR)->Par;
          while (loc2->Vtype!=PREMIA_NULLTYPE) { count++; loc2++;} 
          if ((loc[i].Val.V_PTVAR = malloc(sizeof(PTVAR)))==NULL) return FAIL;
          *(loc[i].Val.V_PTVAR) = *(vars[i].Val.V_PTVAR);
          loc2 =loc[i].Val.V_PTVAR->Par;
          nsp_premia_clone_vars(&loc2,FALSE,vars[i].Val.V_PTVAR->Par,count);
          break;
        case PNLVECT:
          /* when cloning Met->Res, it may happen to have unintialized
             PNLVECT until the Compute function is called */
          if (vars[i].Val.V_PNLVECT != NULL)
            loc[i].Val.V_PNLVECT = pnl_vect_copy(vars[i].Val.V_PNLVECT);
          else
            loc[i].Val.V_PNLVECT=NULL;
          break;
        case PNLVECTCOMPACT:
          if (vars[i].Val.V_PNLVECTCOMPACT != NULL)
            loc[i].Val.V_PNLVECTCOMPACT = pnl_vect_compact_copy(vars[i].Val.V_PNLVECTCOMPACT);
          break;
        case FILENAME:
          if (vars[i].Val.V_FILENAME != NULL)
            {
              count=strlen(vars[i].Val.V_FILENAME);
              if ((loc[i].Val.V_FILENAME= malloc(count+1)) ==NULL) return FAIL;
              memcpy(loc[i].Val.V_FILENAME, vars[i].Val.V_FILENAME, count+1);
            }
          break;
        default: 
          break;
        }
    }
  return OK;
}

/* deallocate.
 * if flag == TRUE then vars is also deallocated
 */
void nsp_premia_free_vars(VAR *vars,int flag,int n)
{
  int i;
  VAR *loc2; 
  for (i=0 ; i < n ; i++)
    {
      int count =0;
      switch( vars[i].Vtype)
        {
        case NUMFUNC_1:
          loc2 = (vars[i].Val.V_NUMFUNC_1)->Par; 
          /* count how many vars are present in vars[i] */
          while (loc2->Vtype!=PREMIA_NULLTYPE) { count++; loc2++;} 
          nsp_premia_free_vars(loc2,FALSE,count);
          break;
        case NUMFUNC_2:
          loc2 = (vars[i].Val.V_NUMFUNC_2)->Par;
          while (loc2->Vtype!=PREMIA_NULLTYPE) { count++; loc2++;} 
          nsp_premia_free_vars(loc2,FALSE,count);
          break;
        case NUMFUNC_ND:
          loc2 = (vars[i].Val.V_NUMFUNC_ND)->Par;
          while (loc2->Vtype!=PREMIA_NULLTYPE) { count++; loc2++;} 
          nsp_premia_free_vars(loc2,FALSE,count);
          break;
        case PTVAR:
          loc2 = (vars[i].Val.V_PTVAR)->Par;
          while (loc2->Vtype!=PREMIA_NULLTYPE) { count++; loc2++;} 
          nsp_premia_free_vars(loc2,FALSE,count);
          break;
        case PNLVECT:
          pnl_vect_free (&(vars[i].Val.V_PNLVECT));
          break;
        case PNLVECTCOMPACT:
          pnl_vect_compact_free (&(vars[i].Val.V_PNLVECTCOMPACT));
          break;
        case FILENAME:
          FREE(vars[i].Val.V_FILENAME);
        default: 
          break;
        }
    }
  if ( flag == TRUE ) { FREE(vars); }
}


extern int *true_typeV;


static int _nsp_premia_get_value(const VAR *x,double *val)
{
  int vt;
  if ( x->Vtype >= FIRSTLEVEL)
    {
      *val =0;
      return OK;
    }
  vt=true_typeV[x->Vtype];
  switch(vt)
    {
    case DOUBLE: *val = x->Val.V_DOUBLE; break;
    case INT:    *val = (double) x->Val.V_INT; break;
    case LONG:   *val = (double) x->Val.V_LONG;break;
    default:
      Scierror("Warning: unknown truetype in the var system\n");
      return FAIL;
      break;
    }
  return OK;
}

int _nsp_premia_set_value(VAR *x,double val)
{
  int vt;
  if ( x->Vtype >= FIRSTLEVEL) return OK;
  vt=true_typeV[x->Vtype];
  switch(vt)
    {
    case DOUBLE: x->Val.V_DOUBLE = val ; break;
    case INT:    x->Val.V_INT = val ; break;
    case LONG:   x->Val.V_LONG = val ;break;
    default:
      Scierror("Warning: unknown truetype in the var system\n");
      return FAIL;
      break;
    }
  return OK;
}

/*
 * Creates a list L = (len, str) where
 * len is the length of str, which is a binary string containing
 * the content of file.
 */
static NspObject* premia_dump_file_forsave (char *file)
{
  FILE *DATA;
  NspObject *Obj;
  int i, len;
  char c, *str;
  str = NULL;
  if ((DATA = fopen(file, "rb")) == NULL) return NULL;
  len = 0;
  while ((c=fgetc(DATA)) != EOF) { len++; }
  if ((str=malloc(sizeof(char)*len)) == NULL) goto err;
  rewind (DATA);
  for ( i=0 ; i<len ; i++) { str[i] = fgetc(DATA); }
  fclose (DATA);
  Obj = (NspObject *) nsp_serial_create("X",str,len);
  if (str != NULL) free (str);
  return Obj;
 err:
  if (str != NULL) free (str);
  fclose (DATA);
  return NULL;
}

static int premia_dump_file_forload (char *file, NspSerial *S)
{
  FILE *DATA;
  int len, i, hl;

  len = S->nbytes;
  hl = strlen (nsp_serial_header);
  if ((DATA = fopen (file, "wb")) == NULL) return FAIL;
  for ( i=hl; i<len ; i++ ) { fputc ( S->val[i], DATA); }
  fclose (DATA);
  return OK;
}

NspList* nsp_premia_get_var_names(NspPremiaModel *self, const VAR *vars,int n)
{
  return nsp_premia_get_var_names_util(self,vars,n, FALSE);
}


/** 
 * Constructs a NspList from the content of VAR. The length of the list is
 * n and each entry of the list has the format
 *
 *    - type (ARRAY, DOUBLE, ENUM, FILENAME)
 *    - name of the variable
 *    - content of the variable (the type of this entry depends on the
 *    field type)
 *    - a boolean (%t if the entry is setable and %f otherwise)
 *
 *    Extra fields for enums only
 *
 *    - vector of labels
 *    - vector of associated keys
 *    - a list of extra args of the form
 *
 *        list( list ( "List of the current type "), ...
 *              list ( "List of the current type "), ...
 *              list ( "List of the current type "));
 *
 *      Each entry of the list is a NspList as produced by
 *      nsp_premia_get_var_names_util. Even if the extra arg is a single
 *      value, we must encapsulate it in a list in case the extra arg is a
 *      real array.
 *
 *    Extra filed for FILENAME only
 *
 *    - the content of the file as a binary string
 *
 * 
 * @param self a PremiaModel
 * @param vars an array of VAR
 * @param n the length of vars
 * @param forsave a flag TRUE or FALSE. If TRUE, when an entry of vars has
 * type FILENAME, the content of the file is saved in binary format in a
 * string
 * 
 * @return a NspList
 */
NspList* nsp_premia_get_var_names_util(NspPremiaModel *self, const VAR *vars,int n, int forsave)
{
  int i;
  NspList *Res=NULL,*loc1=NULL;
  NspObject *Ob;
  NspMatrix *M;
  NspSMatrix *S;
  VAR *vars1,*vars2;
  if ((Res= nsp_list_create("X")) == NULLLIST ) return NULLLIST;
  for (i=0 ; i < n ; i++)
    {
      double val;
      /* type, name, data, setable_flag */
      int_types Ret_default[]={string,string,s_double,s_bool,t_end};
      int_types Ret_array[]={string,string,realmat,s_bool,t_end};
      /* for ENUM : type, name, data, setable_flag, label_vector, key_vector */
      int_types Ret_enum[]={string,string,s_int,s_bool,smat,realmat,t_end};
      int_types Ret_filename[]={string,string,stringcopy,s_bool,t_end};
      int count =0;
      PremiaEnumMember *em=NULL;
      if (vars[i].Viter == IRRELEVANT) continue;
      switch( vars[i].Vtype)
        {
        case NUMFUNC_1:
          vars1=vars2 = (vars[i].Val.V_NUMFUNC_1)->Par;
          while (vars2->Vtype!=PREMIA_NULLTYPE) { count++; vars2++;} 
          if ((loc1=nsp_premia_get_var_names_util(self,vars1,count,forsave))== NULL) goto err;
          if ( nsp_list_concat(Res,loc1) == FAIL) goto err;
          break;
        case NUMFUNC_2:
          vars1=vars2 = (vars[i].Val.V_NUMFUNC_2)->Par;
          while (vars2->Vtype!=PREMIA_NULLTYPE) { count++; vars2++;} 
          if ((loc1=nsp_premia_get_var_names_util(self,vars1,count,forsave))== NULL) goto err;
          if ( nsp_list_concat(Res,loc1) == FAIL) goto err;
          break;
        case NUMFUNC_ND:
          vars1=vars2 = (vars[i].Val.V_NUMFUNC_ND)->Par;
          while (vars2->Vtype!=PREMIA_NULLTYPE) { count++; vars2++;} 
          if ((loc1=nsp_premia_get_var_names_util(self,vars1,count,forsave))== NULL) goto err;
          if ( nsp_list_concat(Res,loc1) == FAIL) goto err;
          break;
        case PTVAR:
          /* list(name,sublist_of-args) */
          vars1=vars2 = (vars[i].Val.V_PTVAR)->Par;
          while (vars2->Vtype!=PREMIA_NULLTYPE) { count++; vars2++;} 
          if ((loc1=nsp_premia_get_var_names_util(self,vars1,count,forsave))== NULL) goto err;
          if ( nsp_list_concat(Res,loc1) == FAIL) goto err;
          break;
        case PNLVECT:
          if (vars[i].Val.V_PNLVECT == NULL)
            {
              if (( M = nsp_matrix_create("M",'R',0,0))== NULLMAT) goto err;
            }
          else
            {
              if (( M = nsp_matrix_create_from_array
                    ("M",1,vars[i].Val.V_PNLVECT->size,
                     vars[i].Val.V_PNLVECT->array,NULL))== NULLMAT) goto err;
           }
          if (( Ob = (NspObject *)
                BuildListFromArgs("lel",Ret_array,"ARRAY", vars[i].Vname,
                                  M,vars[i].Vsetable==SETABLE))== NULLOBJ )
            goto err;
          if ( nsp_list_end_insert(Res,Ob) == FAIL) goto err;
          break;
        case PNLVECTCOMPACT:
          if (vars[i].Val.V_PNLVECTCOMPACT->convert == 'a')
            {
              if (( M = nsp_matrix_create_from_array
                    ("M",1,vars[i].Val.V_PNLVECTCOMPACT->size,
                     vars[i].Val.V_PNLVECTCOMPACT->array,NULL))== NULLMAT) goto err;
              if (( Ob = (NspObject *)
                    BuildListFromArgs("lel",Ret_array,"ARRAY", vars[i].Vname,
                                      M,vars[i].Vsetable==SETABLE))== NULLOBJ )
                goto err;
            }
          else
            {
              if (( Ob = (NspObject *)
                    BuildListFromArgs ("lel",Ret_default, "DOUBLE",vars[i].Vname,
                                       vars[i].Val.V_PNLVECTCOMPACT->val,
                                       vars[i].Vsetable==SETABLE))== NULLOBJ )
                goto err;
            }
          if ( nsp_list_end_insert(Res,Ob) == FAIL) goto err;

          break;
        case ENUM:
          /* first create the vector of labels */
          if ((S = nsp_smatrix_create_from_struct
               ("S", vars[i].Val.V_ENUM.members->members,
                vars[i].Val.V_ENUM.members->size ))==NULLSMAT) goto err;
          /* create the vector of keys */
          em=vars[i].Val.V_ENUM.members->members;
          for ( em=vars[i].Val.V_ENUM.members->members, count=0 ; em->label!=NULL ; em++, count++ );
          if (( M = nsp_matrix_create("M",'r',count,1))== NULLMAT) goto err;
          for ( em=vars[i].Val.V_ENUM.members->members, count=0 ; em->label!=NULL ; em++, count++ )
            { M->R[count]=em->key; }
          if (( Ob = (NspObject *)BuildListFromArgs
                ("lel",Ret_enum,"ENUM", vars[i].Vname,vars[i].Val.V_ENUM.value,
                 vars[i].Vsetable==SETABLE,S,M))== NULLOBJ ) goto err;
          if ( nsp_list_end_insert(Res,Ob) == FAIL) goto err;
          loc1 = nsp_list_create("X");
          for ( em=vars[i].Val.V_ENUM.members->members ; em->label!=NULL ; em++ )
            {
              NspList *loc2 = NULL;
              if ((loc2=nsp_premia_get_var_names_util(self,em->Par,em->nvar,forsave))== NULL) goto err;
              if (nsp_list_end_insert(loc1,(NspObject*)loc2) == FAIL) goto err;
            }
          if ( nsp_list_end_insert((NspList*) Ob,(NspObject*)loc1) == FAIL) goto err;

          break;
        case FILENAME:
          if (forsave == TRUE && vars[i].Vsetable == SETABLE)
            {
              nsp_string basename;
              basename = nsp_tail (vars[i].Val.V_FILENAME);
              if ((Ob = (NspObject *)BuildListFromArgs
                   ("lel",Ret_filename,"FILENAME",vars[i].Vname,basename,
                    vars[i].Vsetable==SETABLE))== NULLOBJ) goto err;
              if ( nsp_list_end_insert(Res,Ob) == FAIL) goto err;
              nsp_string_destroy (&basename);
              if ((Ob = premia_dump_file_forsave (vars[i].Val.V_FILENAME)) == NULL) goto err;
              if ( nsp_list_end_insert(Res,Ob) == FAIL) goto err;
            }
          else
            {
              if ((Ob = (NspObject *)BuildListFromArgs
                   ("lel",Ret_filename,"FILENAME", vars[i].Vname,vars[i].Val.V_FILENAME,
                    vars[i].Vsetable==SETABLE))== NULLOBJ) goto err;
              if ( nsp_list_end_insert(Res,Ob) == FAIL) goto err;
            }
          break;
        default:
          if ( _nsp_premia_get_value(&vars[i],&val) == FAIL)  goto err;
          M = NULL;
          if (self->obj->it != NULL && self->obj->it->range1 != NULL
              && self->obj->it->location1 == &(vars[i]))
            M = self->obj->it->range1; 
          else if (self->obj->it != NULL && self->obj->it->range2 != NULL
                   && self->obj->it->location2 == &(vars[i]))
            M = self->obj->it->range2;

          if (M != NULL)
            {
              if (( Ob = (NspObject *)BuildListFromArgs
                    ("lel",Ret_array, "DOUBLE",vars[i].Vname,M,
                     vars[i].Vsetable==SETABLE))== NULLOBJ )
                goto err;
            }
          else
            {             
              if (( Ob = (NspObject *)BuildListFromArgs
                    ("lel",Ret_default, "DOUBLE",vars[i].Vname,val,
                     vars[i].Vsetable==SETABLE))== NULLOBJ )
                goto err;
            }
          if ( nsp_list_end_insert(Res,Ob) == FAIL) goto err;
          break;
        }
    }
  return Res;
 err:
  if ( Res != NULL ) nsp_list_destroy(Res);
  return NULL;
}

/* 
 * get number of vars: end is detected by PREMIA_NULLTYPE key word 
 */
int nsp_premia_get_nvar(const VAR *vars)
{
  int count=0;
  while (1) 
    {
      if( vars[count].Vtype == PREMIA_NULLTYPE ) break;
      count++;
    }
  return count;
}

int nsp_premia_set_var_names(NspPremiaModel *self,VAR *vars,int n,NspList *L, int depth, void *obj)
{
  return nsp_premia_set_var_names_util(self,vars,n,L,depth,FALSE,obj);
}

/* 
 * get names and values in a list
 */
int nsp_premia_set_var_names_util(NspPremiaModel *self,VAR *vars,int n,NspList *L, int depth, int forload, void *obj)
{
  VAR *vars1,*vars2;
  NspObject *Obj, *Obj1, *Name;
  int i;
  for (i=0 ; i < n ; i++)
    {
      double val;
      int count =0;
      /* We skip IRRELEVANT variables because they are not in L*/
      if ( vars[i].Viter == IRRELEVANT )  continue;
      if (nsp_list_length(L) == 0) return FAIL;
      Obj = nsp_list_get_element(L,1);
      if ( Obj == NULLOBJ || ! IsList(Obj) ) return FAIL;
      if ( nsp_list_length((NspList *)Obj) < 3 ) return FAIL;
      /* we ignore variables which are supposed to be UNSETABLE */
      if ( vars[i].Vsetable == UNSETABLE ) 
        {
          nsp_list_remove_first (L);
          continue;
        }
      Name = nsp_list_get_element((NspList *)Obj,2);
      if ( ! IsString(Name) ) return FAIL;
      Obj1 = nsp_list_get_element((NspList *)Obj,3);
      switch( vars[i].Vtype)
        {
        case NUMFUNC_1:
          vars1=vars2 = (vars[i].Val.V_NUMFUNC_1)->Par;
          while (vars2->Vtype!=PREMIA_NULLTYPE) { count++; vars2++;} 
          if (nsp_premia_set_var_names_util(self,vars1,count,(NspList *) L,depth+1,FALSE,obj) == FAIL)
            return FAIL;
          break;
        case NUMFUNC_2:
          vars1=vars2 = (vars[i].Val.V_NUMFUNC_2)->Par;
          while (vars2->Vtype!=PREMIA_NULLTYPE) { count++; vars2++;} 
          if (nsp_premia_set_var_names_util(self,vars1,count,(NspList *) L,depth+1,FALSE,obj) == FAIL)
            return FAIL;
          break;
        case NUMFUNC_ND:
          vars1=vars2 = (vars[i].Val.V_NUMFUNC_2)->Par;
          while (vars2->Vtype!=PREMIA_NULLTYPE) { count++; vars2++;} 
          if (nsp_premia_set_var_names_util(self,vars1,count,(NspList *) L,depth+1, FALSE,obj) == FAIL)
            return FAIL;
          break;
        case PTVAR:
          /* list(name,sublist_of-args) */
          vars1=vars2 = (vars[i].Val.V_PTVAR)->Par;
          while (vars2->Vtype!=PREMIA_NULLTYPE) { count++; vars2++;} 
          if (nsp_premia_set_var_names_util(self,vars1,count,(NspList *) L,depth+1, FALSE,obj) == FAIL)
            return FAIL;
          break;
        case PNLVECT:
          if ( ! IsMat(Obj1) ) return FAIL;
          if (strcmp(((NspSMatrix *) Name)->S[0],vars[i].Vname) != 0)  return FAIL;
          if (pnl_vect_resize( vars[i].Val.V_PNLVECT, ((NspMatrix *) Obj1)->mn) == FAIL)
            return FAIL;
          memcpy(vars[i].Val.V_PNLVECT->array,((NspMatrix *) Obj1)->R,
                 ((NspMatrix *) Obj1)->mn*sizeof(double));
          nsp_list_remove_first(L);
          break;
        case PNLVECTCOMPACT:
          if ( ! IsMat(Obj1) ) return FAIL;
          if (strcmp(((NspSMatrix *) Name)->S[0],vars[i].Vname) != 0)  return FAIL;
          if (((NspMatrix *) Obj1)->mn==1)
            {
              pnl_vect_compact_set_double (vars[i].Val.V_PNLVECTCOMPACT,
                                           ((NspMatrix *) Obj1)->R[0]);
            }
          else
            {
              if (vars[i].Val.V_PNLVECTCOMPACT->size != ((NspMatrix *) Obj1)->mn) return FAIL;
              pnl_vect_compact_set_ptr (vars[i].Val.V_PNLVECTCOMPACT,
                                        ((NspMatrix *) Obj1)->R);
            }
          nsp_list_remove_first(L);
          break;
        case FILENAME:
          if ( ! IsString(Obj1) ) return FAIL;
          if (strcmp(((NspSMatrix *) Name)->S[0],vars[i].Vname) != 0)  return FAIL;
          FREE(vars[i].Val.V_FILENAME);
          if (((NspSMatrix *)Obj1)->mn != 1) return FAIL;
          if (forload == FALSE)
            {
              vars[i].Val.V_FILENAME = nsp_string_copy(((NspSMatrix *)Obj1)->S[0]);
              nsp_list_remove_first(L);
            }
          else
            {
              char *tmpdir;
              int len;
              tmpdir = getenv("NSP_TMPDIR");
              len = strlen(tmpdir) + strlen (((NspSMatrix *)Obj1)->S[0]) + 2;
              vars[i].Val.V_FILENAME = malloc (len * sizeof(char));
              strcpy (vars[i].Val.V_FILENAME, tmpdir);
              strcat (vars[i].Val.V_FILENAME, path_sep);
              strcat (vars[i].Val.V_FILENAME, ((NspSMatrix *)Obj1)->S[0]);
              nsp_list_remove_first(L);
              if ( (Obj1 = nsp_list_get_element(L,1)) == NULLOBJ ) return FAIL;
              if ( ! IsSerial(Obj1) ) return FAIL;
              premia_dump_file_forload (vars[i].Val.V_FILENAME, (NspSerial *)Obj1);
              nsp_list_remove_first(L);
            }
          break;
        case ENUM:
            {  
              PremiaEnumMember *em;
              NspObject *Obj2, *Obj3;
              int index;
              if ( ! IsMat(Obj1) ) return FAIL;
              if (strcmp(((NspSMatrix *) Name)->S[0],vars[i].Vname) != 0)  return FAIL;
              if ( ((NspMatrix *) Obj1)->mn != 1  ) return FAIL;
              val = ((NspMatrix *) Obj1)->R[0];
              vars[i].Val.V_ENUM.value = (int)val;
              /* Get the Par array */
              em = lookup_premia_enum_with_index (&(vars[i]), (int)val, &index);
              if ( em->nvar > 0 )
                {
                  vars1 = em->Par;
                  /* extract the list of Par args */
                  Obj2 = nsp_list_get_element ((NspList *)Obj, 7);
                  if ( Obj2 == NULLOBJ || ! IsList(Obj2) ) return FAIL;
                  /* Choose the correct Par arg */
                  Obj3 = nsp_list_get_element ((NspList *)Obj2, index + 1);
                  if ( Obj3 == NULLOBJ || ! IsList(Obj3) ) return FAIL;
                  if (nsp_premia_set_var_names_util(self,vars1,em->nvar,(NspList *) Obj3,depth+1, forload,obj) == FAIL) return FAIL;
                }
              nsp_list_remove_first(L);
            }  
          break;
        default: 
          if ( ! IsMat(Obj1) ) return FAIL;
          if (strcmp(((NspSMatrix *) Name)->S[0],vars[i].Vname) != 0)  return FAIL;
          if ( ((NspMatrix *) Obj1)->mn < 1 ) return FAIL;
          /* always set a value, otherwise the test functions will crash */
          /* make sure the Matrix is expanded if ever it has been given as
             [1:10] for instance, see impl member of NspMatrix */
          Mat2double((NspMatrix *) Obj1);
          val = ((NspMatrix *) Obj1)->R[0];
          _nsp_premia_set_value(&vars[i],val);

          /* check if some iteration is required */
          if ( (((NspMatrix *) Obj1)->m == 1 || ((NspMatrix *) Obj1)->n == 1) 
               && ((NspMatrix *) Obj1)->mn > 1)
            {
              NspMatrix **range;
              VAR **location;
              int empty_it = 0;
              if (self->obj->it==NULL)
                {
                  if ((self->obj->it=malloc(sizeof(nsp_premia_iterator)))==NULL) return FAIL;
                  self->obj->it->location1 = NULL;
                  self->obj->it->location2 = NULL;
                  self->obj->it->range1 = NULLMAT;
                  self->obj->it->range2 = NULLMAT;
                  range = &(self->obj->it->range1);
                  location = &(self->obj->it->location1);
                  empty_it = 1;
                }
              else
                {
                  if (self->obj->it->range1 == NULLMAT)
                    {
                      range = &(self->obj->it->range1);
                      location = &(self->obj->it->location1);
                      empty_it = 1;
                    }
                  else if (self->obj->it->range2 == NULLMAT)
                    {
                      range = &(self->obj->it->range2);
                      location = &(self->obj->it->location2);
                      empty_it = 1;
                    }
                  else empty_it = 0;
                }

              if (empty_it == 1)
                {
                  self->obj->it->prix = NULLMAT;
                  *location = &(vars[i]);
                  if ((*range = (NspMatrix *)
                       nsp_object_copy_and_name(vars[i].Vname, Obj1)) == NULLMAT) return FAIL;
                  /* convert range to column matrix */
                  (*range)->m = (*range)->mn; (*range)->n = 1;
                }
            }
          nsp_list_remove_first(L);
          break;
        }
      if ( vars[i].setter ) vars[i].setter(obj);
    }
  return OK;
}

void premia_free_attr (NspPremiaModel *self, p_attr type)
{
  switch (type)
    {
    case p_mod :
      if ( self->obj->mod.TypeModel == NULL ) return;
      self->obj->mod.init = 0;
      nsp_premia_free_vars(self->obj->mod.TypeModel,TRUE,self->obj->mod.nvar);
      self->obj->mod.TypeModel = NULL;
      self->obj->mod.Name = NULL;
      break;
    case p_opt :
      if (self->obj->opt.TypeOpt == NULL) return;
      self->obj->opt.init = 0;
      nsp_premia_free_vars(self->obj->opt.TypeOpt,TRUE,self->obj->opt.nvar);
      self->obj->opt.TypeOpt= NULL;
      self->obj->opt.Name= NULL;
      break;
    case p_meth:
      if ( self->obj->meth.Name == NULL) return;
      self->obj->meth.init = 0;
      nsp_premia_free_vars(self->obj->meth.Res,FALSE,
                           nsp_premia_get_nvar(self->obj->meth.Res));
      nsp_premia_free_vars(self->obj->meth.Par,FALSE,
                           nsp_premia_get_nvar(self->obj->meth.Par));
      self->obj->meth.Name = NULL;
      break;
    default: break;
    }
}

/*
 * sets the asset type of a NspPremiaModel using its name
 */
int nsp_premiamodel_set_asset (NspPremiaModel *self, const char *type)
{
  PremiaAsset *dummy_asset=premia_assets;

  while (dummy_asset->name!=NULL)
    {
      if (strcmp(type, dummy_asset->name)==0)
        {
          self->obj->asset = *dummy_asset;
          break;
        }
      else
        dummy_asset++;
    }

  if (dummy_asset->name==NULL)
    {
      Scierror ("asset type unknown, accepted values are : ");
      premia_accepted_assets();
      return RET_BUG;
    }
  return OK;
}

/*
 * sets the model of a NspPremiaModel using its name
 */
int nsp_premiamodel_set_model_with_str (NspPremiaModel *self, const char *str)
{
  VAR *var;
  Model *premia_model;
  Model **models;
  int model=0;

  premia_model = NULL;
  models = self->obj->asset.models;

  while (models[model] != NULL)
    {
      if (strcmp(str,models[model]->Name)==0)
        {
          premia_model = models[model]; break;
        }
      model++;
    }
  if (premia_model == NULL)   return RET_BUG;

  /* be sure that model is initialized to have correct nvar */
  premia_model->Init(premia_model);
  /* be sure that type premiamodel is initialized */
  self->obj->mod = *premia_model;
  /* clone vars recursively, even if it's not usefull for models */
  if ( nsp_premia_clone_vars(&var,TRUE,(VAR *) premia_model->TypeModel,premia_model->nvar)==FAIL) 
    return RET_BUG;
  self->obj->mod.TypeModel=var;
  return OK;  
} 

/*
 * sets the model of a NspPremiaModel using its number
 */
int nsp_premiamodel_set_model (NspPremiaModel *self, int n)
{
  VAR *var;
  Model *premia_model;

  if ( (premia_model = premia_get_model_from_index (self->obj->asset, n)) == NULL )
    {
      Scierror("Error: model %d does not exist\n",n+1);
      return RET_BUG;
    }

  /* be sure that model is initialized to have correct nvar */
  premia_model->Init(premia_model);
  /* be sure that type premiamodel is initialized */
  self->obj->mod = *premia_model;
  /* clone vars recursively, even if it's not usefull for models */
  if ( nsp_premia_clone_vars(&var,TRUE,(VAR *) premia_model->TypeModel,premia_model->nvar)==FAIL) 
    return RET_BUG;
  self->obj->mod.TypeModel=var;
  return OK;  
}

/*
 * sets the option of a NspPremiaModel using its name. Its relies on the
 * uniqueness of option names in Premia.
 */
int nsp_premiamodel_set_option_with_str (NspPremiaModel *self, const char *str)
{
  Option **loc;
  Family **families;
  int option=0,n_option=0,family=0,n_family=0,ret;


  families = self->obj->asset.families;
  
  while ( families[n_family] != NULL) n_family++;
  for (family=0; family<n_family; family++)
    {
      loc = (*families[family]);
      n_option = 0;
      while ( loc[n_option] != NULL) n_option++;
      for (option=0; option<n_option; option++)
        {
          if (strcmp(loc[option]->Name, str) == 0)
            {
              if ((ret=nsp_premiamodel_set_option (self, family, option))!=OK) return RET_BUG;
              return OK;
            }
        }
    }
  Scierror("Error: option %s does not exist in any family \n",str);
  return RET_BUG;
} 

/*
 * sets the option of a NspPremiaModel using the family and option
 * indices. These indices are computed in full lists of families and
 * options.
 */
int nsp_premiamodel_set_option (NspPremiaModel *self, int family, int option)
{
  Option *premia_option, **loc;
  VAR *var;
  Family **families;
  Pricing **pricings;
  int n_option=0,n_family=0;

  families = self->obj->asset.families;
  pricings = self->obj->asset.pricings;
  
  while ( families[n_family] != NULL) n_family++;

  if ( family  < 0 || family > n_family -1 ) 
    {
      Scierror("Error: family %d does not exist\n",family+1);
      return RET_BUG;
    }
  loc = (*families[family]);
  while ( loc[n_option] != NULL) n_option++;
  if ( option < 0 || option > n_option -1 ) 
    {
      Scierror("Error: option %d does not exist in family %d\n",option+1,family+1);
      return RET_BUG;
    }
  premia_option=(*families[family])[option];
  /* be sure that option is initialized according to the model.
     Forcing a new initialization. */
  premia_option->init = 0; 
  premia_option->Init(premia_option, &self->obj->mod);
  /* check that the option is compatible with the model */
  if ( Premia_match_model_option(&self->obj->mod,premia_option,pricings) != 0) 
    {
      Scierror("Error: option %d is not compatible with model\n",option+1);
      return RET_BUG;
    }
  self->obj->opt = *premia_option;
  /* clone vars recursively */
  if ( nsp_premia_clone_vars(&var,TRUE,(VAR *) premia_option->TypeOpt,premia_option->nvar)==FAIL) 
    return RET_BUG;
  self->obj->opt.TypeOpt=var;
  return OK;
}

/*
 * sets the method of a NspPremiaModel using its name
 */
int nsp_premiamodel_set_method_with_str (NspPremiaModel *self, const char* str)
{
  Pricing *res;
  Pricing **pricings;
  int n_method=0,i,ret,i_method;

  pricings = self->obj->asset.pricings;

  /*
   * we know that MatcingPricing is ok.
   * we don't check the status of SelectPricing because it may happen that
   * model and option parameters are not compatible at the current step. It is
   * checked right before computing the result
   */
  SelectPricing(0,&self->obj->mod,&self->obj->opt,pricings,&res);
  while (res->Methods[n_method] != NULL) n_method++;

  i_method=0;
  for ( i=0 ; i < n_method ; i++) 
    {
      if ( res->Methods[i]->CheckOpt(&self->obj->opt,&self->obj->mod)== OK )
        {
          if ( strcmp (res->Methods[i]->Name, str) == 0 )
            {
              if ((ret=nsp_premiamodel_set_method (self, i_method))!=OK) return ret;
              return OK;
            }
          i_method++;
        }
    }

  Scierror("Error: method %s does not exist\n", str);
  return RET_BUG;
} 

/*
 * sets the method of a NspPremiaModel using its number. Take care that the
 * number integer is the index of the method in the list of available pricings
 * for the model and option (not in the full list)
 */
int nsp_premiamodel_set_method (NspPremiaModel *self, int method)
{
  Pricing *res;
  VAR *var;
  Pricing **pricings;
  int n_method=0,npar,nres,method_id=-1,i,count=0;

  pricings = self->obj->asset.pricings;
  
  /*
   * we know that MatcingPricing is ok.
   * we don't check the status of SelectPricing because it may happen that
   * model and option parameters are not compatible at the current step. It is
   * checked right before computing the result
   */
  SelectPricing(0,&self->obj->mod,&self->obj->opt,pricings,&res);
  while (res->Methods[n_method] != NULL) n_method++;

  /* find the index of the method in the full list. This index will be
     method_id */
  for ( i=0 ; i < n_method ; i++) 
    {
      if ( res->Methods[i]->CheckOpt(&self->obj->opt,&self->obj->mod)== OK) 
        {
          if ( count == method ) 
            {
              method_id=i;break;
            }
          count++;
        }
    }

  if (method_id == -1) 
    {
      Scierror("Error: method %d does not exist\n", method);
      return RET_BUG;
    }

  /* free previous method */
  premia_free_attr (self, p_meth);
  res->Methods[method_id]->Init(res->Methods[method_id], &self->obj->opt);
  self->obj->meth = *(res->Methods[method_id]);

  /* clone vars recursively */
  var=self->obj->meth.Par;
  npar = nsp_premia_get_nvar(var);
  if ( nsp_premia_clone_vars(&var,FALSE,(VAR *) res->Methods[method_id]->Par,npar)==FAIL) 
    return RET_BUG;
  var = self->obj->meth.Res;
  nres = nsp_premia_get_nvar(var);
  if ( nsp_premia_clone_vars(&var,FALSE,(VAR *) res->Methods[method_id]->Res,nres)==FAIL) 
    return RET_BUG;
  self->obj->meth.Init(&self->obj->meth, &self->obj->opt);
  return OK;
}

