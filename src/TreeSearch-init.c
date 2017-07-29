#define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>

#include "ape_reorder.h"
#include "phangorn_fitch.h"
#include "renumber_tree.h"

// Abandoned: attempts to link to ape / phangorn functions
// static R_NativePrimitiveArgType ape_node_depth_t[] = {
//   INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP
// };
// void ape_node_depth(int *ntip, int *nnode, int *e1, int *e2,
// 		int *nedge, double *xx, int *method) {
//   typedef int (*Fun)(int*, int*, int*, int*, int*, double*, int*);
//   Fun fun = (Fun) R_GetCCallable( "ape", "node_depth" );   
//   fun(ntip, nnode, e1, e2, nedge, xx, method);
// }



static const R_CMethodDef cMethods[] = {
  {"order_edges_number_nodes", (DL_FUNC) &order_edges_number_nodes, 3, order_edges_number_nodes_t},
  {"ape_neworder_phylo", (DL_FUNC) &ape_neworder_phylo, 3, ape_neworder_phylo_t},
  {"ape_node_depth", (DL_FUNC) &ape_node_depth, 7, ape_node_depth_t},
  {"ape_neworder_pruningwise", (DL_FUNC) &ape_neworder_pruningwise, 3, ape_neworder_pruningwise_t},
  {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
  {"FITCH", (DL_FUNC) &FITCH, 8},
  {"RENUMBER_TREE", (DL_FUNC) &RENUMBER_TREE, 3},
  {NULL, NULL, 0}
};

void R_init_TreeSearch(DllInfo *info) {
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}

