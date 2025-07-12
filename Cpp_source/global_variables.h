#ifndef VARIABLE_local_HEADER
#define VARIABLE_local_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include <random>

#include <ilcplex/cplex.h>

using namespace std;


typedef struct {

    string output_file_data;
    string output_file_x;

    char alpha_file_data[256], beta_file_data[256], theta_file_data[256], x_file_data[256];

    //////////////////////////////////////////////////////////////////////////
    // NETWORK INFLUENCE PROBLEM PARAMETERS
    //////////////////////////////////////////////////////////////////////////

    int n_individuals;
    int n_groups;
    int n_features;
    int n_sol;
    int tau_max;
    int tau_min;
    int method;
    int problem;
    int Heuristic_only;
    int LB_heur;
    int UB_heur;
    int MagnantiW;

    int* n_neighbors_up;    // [n_individuals] 
    int* n_neighbors_down;  // [n_individuals] 

    int** neighbors_up;     // [n_individuals] x [n_neighbors_up], a list of neighbors for each individuals
    int** neighbors_down;   // [n_individuals] x [n_neighbors_down], a list of neighbors for each individuals
    int** group_features;   // [n_features] x [n_groups]

    int** alpha_data;       // [n_individuals] x [n_individuals]
    int** beta_data;        // [n_individuals] x [n_features]
    int* theta_data;        // [n_features]


    /////////////////////////////////////CPLEX////////////////////////////////

    CPXENVptr env_MILP;
    CPXLPptr  lp_MILP;

    ///////////////////////////////////////////////////////////////////////////////

    double tol_epsilon_lazy, tol_epsilon_user;
    int max_cuts_lazy, max_cuts_user, model, iteration_BEN, n_variables_BEN;
    int status, status_sub, ccnt, ccnt_sub, rcnt, rcnt_sub, nzcnt, nzcnt_sub, nzcnt_BEN;
    int* rmatbeg, * rmatind, * cmatbeg, * cmatind, * effortlevel;
    double* rmatval, * cmatval, * x, * pi, * obj, * lb, * ub, * rhs;
    double* sol_callback;
    char* c_type, * sense;
    char** colname, ** rmatname;
    double objval, bestobjval;
    int lpstat, nodecount;
    int number_of_CPU;
    int timelimit;
    double TOLL_OPTIMALITY;
    int option;
    double Theta_val;
    int opt_bound;
    double rhs_BEN;
    int* cut_BEN_rmatind;
    double* cut_BEN_rmatval;
    int* cut_BEN_Tay_rmatind;
    double* cut_BEN_Tay_rmatval;
    double** x_ic;
    int*** X_ic;
    int*** Heur_X_ic;
    int** Tabu_X_ic;
    int Tabu_not_improving_max;
    int SA_Iter_max;
    double**** zeta;
    double**** zeta_i;
    double**** zeta_j;
    double** Zeta0;
    double** Zeta1;
    double incumbent_old;
    double incumbent_new;
    int external_sol;
    int n_cuts;
    bool warm_start;
    bool get_feasible;
    int Heu_sol;

    int sol_val;
    int* sol_val_EXTERNAL;
    int* sol_val_UP;
    int* sol_val_DOWN;
    int* sol_val_HEUR_DOWN;

    double solution_time;
    double solution_time_HEUR;

} instance;




#endif
