#ifndef _PREMIA_VARS
#define _PREMIA_VARS

#include "optype.h"
#include "tools.h"
#include "premia_obj.h"
#include "nsp/object.h"
#include "nsp/interf.h"



extern int premia_get_asset_type(Stack stack, int rhs, int opt, int lhs, PremiaAsset *asset);
extern void premia_accepted_assets();
extern Model* premia_get_model_from_index (PremiaAsset asset, int n);
extern NspSMatrix* premia_get_models_aux (PremiaAsset asset);
extern NspSMatrix* premia_get_methods_aux (PremiaAsset asset, Model *mod, Option *opt);
extern NspSMatrix* premia_get_family_aux(PremiaAsset asset, Model *mod, int n);
#endif /* _PREMIA_VARS */
