
#include "global_variables.h"

void generate_random_data(instance* inst);

void read_data_alpha(instance* inst);

void read_data_beta(instance* inst);

void read_data_theta(instance* inst);

void read_data_x(instance* inst);

void free_data(instance* inst);

void print_data(instance* inst);


//------------------------------------------------
// Variables positions
//------------------------------------------------

int position_eta(instance* inst);

int position_x(instance* inst, int i, int c);

int position_w(instance* inst, int i, int j, int c);


//------------------------------------------------
// Constraints
//------------------------------------------------

void build_c1(instance* inst);

void build_constr_tau(instance* inst);

void build_constr_w(instance* inst);

void build_constr_homogeneiry(instance* inst); 

void build_constr_opt_bound(instance* inst);

void build_constr_symmetry_breaking(instance* inst);

void build_constr_min_sum_Taylor_QA(instance* inst);

void build_constr_min_max_Taylor_QA(instance* inst);

void build_constr_min_sum_first_cut(instance* inst);

void build_constr_min_max_first_cut(instance* inst);

//------------------------------------------------
// Lagrangian multipliers
//------------------------------------------------

void get_zeta_MW(instance* inst);

void get_zeta_MW_multiple(instance* inst);

void get_zeta(instance* inst, double rr);

void get_zeta_first_cut(instance* inst, double rr);

void get_zeta_first_cut_multiple(instance* inst, double rr);


//------------------------------------------------
// Heuristics
//------------------------------------------------


void R_MN_SAP_ConstructuveHeur(instance* inst);

void R_MN_SAP_LocalSearch(instance* inst);

void R_MN_SAP_TabuSearch(instance* inst);

void R_MN_SAP_SimulatedAnnealing(instance* inst);

void objective(instance* inst, int method, int sol);

