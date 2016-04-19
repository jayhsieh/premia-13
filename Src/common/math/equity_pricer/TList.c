#include <stdlib.h>

#include "config.h"
#include "gd_list.h"

#define MY_NEW_BASE /* to avoid an unknown type error in templates_on.h */

#define BASE 
#define KEY int
#define VALUE double
#define CONTAIN PremiaContains
#define NODE PremiaNode
#define NODE_SHORT node
#define CONTAIN_SHORT contains 
#define SHORT  
#define LSHORT sort_list
#define IN_FORMAT "%lf"
#define OUT_FORMAT "%f"
#define OUT_PUT_FORMAT(a) a

#include "pnl/pnl_templates_on.h"
#define FUNCTION_NODE(dir,name) CONCAT3(dir,NODE_SHORT,name)
#define FUNCTION_CONTAIN(dir,name) CONCAT3(dir,CONTAIN_SHORT,name)

#include "TList_source.c"
#undef CONTAIN
#undef NODE
#undef NODE_SHORT
#undef CONTAIN_SHORT
#undef FUNCTION_CONTAIN
#undef FUNCTION_NODE
#include "pnl/pnl_templates_off.h"



#define BASE 
#define KEY PnlVectInt *
#define VALUE int
#define CONTAIN PremiaSparsePoint
#define NODE PremiaNodeSparsePoint
#define NODE_SHORT node_sparse_point
#define CONTAIN_SHORT sparse_point
#define SHORT SparsePoint  
#define LSHORT sort_list_sparse_point
#define IN_FORMAT "%lf"
#define OUT_FORMAT "%f"
#define OUT_PUT_FORMAT(a) a
#define FUNCTION_NODE(dir,name) CONCAT3(dir,NODE_SHORT,name)
#define FUNCTION_CONTAIN(dir,name) CONCAT3(dir,CONTAIN_SHORT,name)

#include "pnl/pnl_templates_on.h"
#include "TList_source.c"
#undef CONTAIN
#undef NODE
#undef NODE_SHORT
#undef CONTAIN_SHORT
#undef FUNCTION_CONTAIN
#undef FUNCTION_NODE
#include "pnl/pnl_templates_off.h"


