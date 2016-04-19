/* -*- Mode: C -*- */
#ifndef NSP_INC_PremiaModel
#define NSP_INC_PremiaModel

/*
 * This Software is GPL (Copyright ENPC 1998-2005) 
 * Jean-Philippe Chancelier Enpc/Cermics
 * Jerome Lelong Enpc Inria
 */

/* PremiaModel */

#include <nsp/object.h> 
#include <nsp/matrix.h> 
#include <nsp/bmatrix.h> 
#include <nsp/smatrix.h> 
#include <nsp/list.h> 
#include <nsp/hash.h> 
#include <nsp/type.h> 
#include "nsp/serial.h"

#include "optype.h" 
#include "var.h" 
#include "tools.h"
#include "error_msg.h"
#include "premia_obj.h"

/*
 * NspPremiaModel inherits from Nsp+Object
 */

typedef struct _NspPremiaModel NspPremiaModel ;
typedef struct _NspTypePremiaModel NspTypePremiaModel ;

struct _NspTypePremiaModel {
  /*< private >*/
  NSP_TYPE_OBJECT__
  /*< public >*/
};

/* iterator structure for a Premiamodel. Only two iterators are
   allowed. location is the address of the variable to modify. range is the
   sequence of values to successively use. */
typedef struct _nsp_premia_iterator nsp_premia_iterator;
struct _nsp_premia_iterator {
  VAR *location1;
  NspMatrix *range1;
  VAR *location2;
  NspMatrix *range2;
  NspMatrix *prix;  
};

typedef struct _nsp_premiamodel nsp_premiamodel;
struct _nsp_premiamodel {
  int ref_count;
  PremiaAsset asset;
  Model mod;
  Option opt;
  PricingMethod meth;
  int compat;
  nsp_premia_iterator *it;
  int compute_err; /* return value of the compute method */
};

struct _NspPremiaModel {
  /*< private >*/
  NspObject father;
  NspTypePremiaModel*type;
  /*< public >*/
  nsp_premiamodel *obj;
};

extern int nsp_type_premiamodel_id;
extern NspTypePremiaModel *nsp_type_premiamodel;

/* type instances for object */

NspTypePremiaModel *new_type_premiamodel(type_mode mode);

/* instance for PremiaModel */

NspPremiaModel *new_premiamodel();

/*
 * Object methods redefined for premiamodel 
 */


#define NULLPREMIAMODEL (NspPremiaModel*) 0

extern NspPremiaModel *premiamodel_create(char *name,NspTypeBase *type);

/* from PremiaModelObj.c */

extern NspPremiaModel *nsp_premiamodel_copy(NspPremiaModel *H);
extern void nsp_premiamodel_destroy(NspPremiaModel *H);
int nsp_premiamodel_print(NspPremiaModel *M,int indent,const char *name, int rec_level);
int nsp_premiamodel_info(NspPremiaModel *M, int indent,const char *name, int rec_level);
extern NspPremiaModel *nsp_premiamodel_object (NspObject *O); 
extern int IsPremiaModelObj (Stack stack, int i); 
extern int IsPremiaModel(NspObject *O);
extern NspPremiaModel *GetPremiaModelCopy (Stack stack, int i); 
extern NspPremiaModel *GetPremiaModel (Stack stack, int i); 
int int_premiamodel_create(Stack stack, int rhs, int opt, int lhs);

extern int _nsp_premia_set_value(VAR *x,double val);
extern int nsp_premia_set_var_names(NspPremiaModel *M,VAR *vars,int n,NspList *L,int depth,void *opt);
extern int nsp_premia_get_nvar(const VAR *vars);
extern NspList* nsp_premia_get_var_names(NspPremiaModel *M, const VAR *vars,int n);
extern  void nsp_premia_free_vars(VAR *vars,int flag,int n);
extern int nsp_premia_clone_vars(VAR **res,int flag,const VAR *vars,int n);
extern int nsp_premiamodel_set_asset (NspPremiaModel *self, const char *type);
extern int nsp_premiamodel_set_model_with_str (NspPremiaModel *self, const char *str);
extern int nsp_premiamodel_set_model (NspPremiaModel *self, int model); 
extern int nsp_premiamodel_set_option_with_str (NspPremiaModel *self, const char *str);
extern int nsp_premiamodel_set_option (NspPremiaModel *self, int family, int option);
extern int nsp_premiamodel_set_method_with_str (NspPremiaModel *self, const char *str);
extern int nsp_premiamodel_set_method (NspPremiaModel *self, int method);
typedef enum { p_asset, p_mod, p_opt, p_meth } p_attr;
extern void premia_free_attr (NspPremiaModel *self, p_attr type);
NspList* nsp_premia_get_var_names_util(NspPremiaModel *self, const VAR *vars,int n, int forsave);
int nsp_premia_set_var_names_util(NspPremiaModel *self,VAR *vars,int n,NspList *L, int depth, int forload, void *obj);

#ifdef PremiaModel_Private 
static int init_premiamodel(NspPremiaModel *o,NspTypePremiaModel *type);
static int nsp_premiamodel_size(NspPremiaModel *Mat, int flag);
static char *nsp_premiamodel_type_as_string(void);
static char *nsp_premiamodel_type_short_string(NspObject *v);
static int nsp_premiamodel_eq(NspPremiaModel *A, NspObject *B);
static int nsp_premiamodel_neq(NspPremiaModel *A, NspObject *B);
extern AttrTab premiamodel_attrs[];
extern NspMethods *premiamodel_get_methods(void);
static NspPremiaModel *premiamodel_create_void(char *name,NspTypeBase *type);
static int nsp_premiamodel_xdr_save(XDR *xdrs, NspPremiaModel *M);
static NspPremiaModel* nsp_premiamodel_xdr_load (XDR *xdrs);
static int nsp_premia_set_var_names_forload(NspPremiaModel *self,VAR *vars,int n,NspList *L, int depth, void *obj);
static NspList* nsp_premia_get_var_names_forsave(NspPremiaModel *self, const VAR *vars,int n);

#endif /* PremiaModel_Private */

#endif /* NSP_INC_PremiaModel */ 
