#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void initialize_delta(int *,double *);
extern void add_row_wrapper(int *, int *);
extern void boost_wrapper(int *, int *, int *, double *, double *, int *, double *);
extern void delete_node_wrapper(int *);
extern void delete_oob();
extern void free_daop();
extern void free_predict();
extern void get_data_options_info(int *, int *, int *);
extern void get_oob(int *);
extern void get_rows(int *, int *, int *, int *);
extern void get_tree(int *, double *, int *);
extern void get_tree_all(int *, double *, int *);
extern void get_weights(int *, int *);
extern void permute_wrapper(int *, int *);
extern void predict_trees_all(int *, double *);
extern void predict_trees_fast(int *, double *);
extern void print_node_members(int *);
extern void quantile_aux(double *, double *, int *, double *, double *);
extern void quantile_R(double *, int *, double *, int *, double *, int *, double *);
extern void read_data(double *, int *, int *, double *, int *, int *, double *, int *, int *, int *, int *, int *,int *);
extern void read_predict(double *, int *, int *, double *, int *, int *, double *, int *, int *, int *,int *,int *);
extern void sample_indx_wrapper(int *, int *);
extern void set_data(double *, int *, int *, double *, int *, int *);
extern void set_lambda(double *);
extern void set_mtry(int *);
extern void set_oob(int *);
extern void set_tnode_row();
extern void set_train(int *);
extern void set_variable_status(int *, int *, int *, int *);
extern void varimp(int *, double *);
extern void varimp_boost(int *, double *);
extern void check_meansum_summary(int *,double *,int *,double *,double *,int *);
extern void set_n_max(int *);
extern void get_tree_all_gini(int *,int *,double *,int *);
extern void set_ncat(int *);
extern void predict_trees_gini(int *,int *,double *);
extern void set_minnodesize(int *);
extern void set_method_family(int *);

static const R_CMethodDef CEntries[] = {
    {"initialize_delta",      (DL_FUNC) &initialize_delta,        2},
    {"add_row_wrapper",       (DL_FUNC) &add_row_wrapper,        2},
    {"boost_wrapper",         (DL_FUNC) &boost_wrapper,          7},
    {"delete_node_wrapper",   (DL_FUNC) &delete_node_wrapper,    1},
    {"delete_oob",            (DL_FUNC) &delete_oob,             0},
    {"free_daop",             (DL_FUNC) &free_daop,              0},
    {"free_predict",          (DL_FUNC) &free_predict,           0},
    {"get_data_options_info", (DL_FUNC) &get_data_options_info,  3},
    {"get_oob",               (DL_FUNC) &get_oob,                1},
    {"get_rows",              (DL_FUNC) &get_rows,               4},
    {"get_tree",              (DL_FUNC) &get_tree,               3},
    {"get_tree_all",          (DL_FUNC) &get_tree_all,           3},
    {"get_weights",           (DL_FUNC) &get_weights,            2},
    {"permute_wrapper",       (DL_FUNC) &permute_wrapper,        2},
    {"predict_trees_all",     (DL_FUNC) &predict_trees_all,      2},
    {"predict_trees_fast",    (DL_FUNC) &predict_trees_fast,     2},
    {"print_node_members",    (DL_FUNC) &print_node_members,     1},
    {"quantile_aux",          (DL_FUNC) &quantile_aux,           5},
    {"quantile_R",            (DL_FUNC) &quantile_R,             7},
    {"read_data",             (DL_FUNC) &read_data,             13},
    {"read_predict",          (DL_FUNC) &read_predict,          12},
    {"sample_indx_wrapper",   (DL_FUNC) &sample_indx_wrapper,    2},
    {"set_data",              (DL_FUNC) &set_data,               6},
    {"set_lambda",            (DL_FUNC) &set_lambda,             1},
    {"set_mtry",              (DL_FUNC) &set_mtry,               1},
    {"set_oob",               (DL_FUNC) &set_oob,                1},
    {"set_tnode_row",         (DL_FUNC) &set_tnode_row,          0},
    {"set_train",             (DL_FUNC) &set_train,              1},
    {"set_variable_status",   (DL_FUNC) &set_variable_status,    4},
    {"varimp",                (DL_FUNC) &varimp,                 2},
    {"varimp_boost",          (DL_FUNC) &varimp_boost,           2},
    {"check_meansum_summary", (DL_FUNC) &check_meansum_summary,   6},
    {"set_n_max",             (DL_FUNC) &set_n_max,   1},
    {"get_tree_all_gini",     (DL_FUNC) &get_tree_all_gini,   4},
    {"set_ncat",              (DL_FUNC) &set_ncat,   1},
    {"predict_trees_gini",    (DL_FUNC) &predict_trees_gini,   3},
    {"set_minnodesize",       (DL_FUNC) &set_minnodesize,   1},
    {"set_method_family",     (DL_FUNC) &set_method_family,   1},
    {NULL, NULL, 0}
};

void R_init_htree(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
