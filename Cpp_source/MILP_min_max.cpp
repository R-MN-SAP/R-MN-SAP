
#define _CRT_SECURE_NO_WARNINGS

#include "global_functions.h"
#include "global_variables.h"
#include "MILP_min_sum.h"

#define write_prob
#include "global_functions.h"
#include "MILP_min_max.h"

//------------------------------------------------
// Constraints
//------------------------------------------------


void build_obj_min_max(instance* inst)
{

	inst->rcnt = 1;
	inst->nzcnt = (1 + (int)pow(inst->n_individuals, 2));

	for (int c = 0; c < inst->n_groups; c++)
	{
		// * creating the knapsack time 0 constraint *

		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

		inst->rhs[0] = 0.0;
		inst->sense[0] = 'L';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		int counter = 0;

		inst->rmatind[counter] = position_eta(inst);
		inst->rmatval[counter] = (double)1.0;

		counter++;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int j = 0; j < inst->n_individuals; j++)
			{
				if (i == j) {
					inst->rmatind[counter] = position_w(inst, i, j, c);
					inst->rmatval[counter] = (double)(0.0);
				}
				else {
					inst->rmatind[counter] = position_w(inst, i, j, c);
					inst->rmatval[counter] = (double)(-inst->alpha_data[i][j]);
				}

				//cout << position_w(inst, i, j, c) << endl;

				counter++;
			}
		}

		inst->rmatbeg[0] = 0;

		inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, 1, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);

		if (inst->status != 0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}

		free(inst->rmatbeg);
		free(inst->rmatval);
		free(inst->rmatind);
		free(inst->rhs);
		free(inst->sense);

	}
};



//------------------------------------------------
// Build MILP
//------------------------------------------------

void build_MILP_min_max(instance* inst) {

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * setting CPLEX environment and problem

	//opening the environment
	inst->env_MILP = CPXopenCPLEX(&(inst->status));
	if (inst->status != 0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	//opening the pointer to the problem
	inst->lp_MILP = CPXcreateprob(inst->env_MILP, &(inst->status), "LP");
	if (inst->status != 0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the variables and linear part of the objective function*
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	inst->ccnt = 1 + ((int)(pow(inst->n_individuals, 2) * inst->n_groups) + inst->n_individuals * inst->n_groups);

	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));

	inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(2000, sizeof(char)); }

	int counter = 0;

	for (int c = 0; c < inst->n_groups; c++)
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			inst->obj[counter] = (double)0.0;
			inst->lb[counter] = (double)0.0;
			inst->ub[counter] = 1.0;
			inst->c_type[counter] = 'B';

			sprintf(inst->colname[counter], "x_%d_%d", i + 1, c + 1);

			//cout << "(" << i + 1 << ", " << c + 1 << "): " << position_x(inst, i, c) << " counter : " << counter << endl;

			counter++;

		}
	}


	for (int c = 0; c < inst->n_groups; c++)
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int j = 0; j < inst->n_individuals; j++)
			{

				inst->obj[counter] = (double)0.0;
				inst->lb[counter] = (double)0.0;
				inst->ub[counter] = (double)1.0;
				inst->c_type[counter] = 'C';

				sprintf(inst->colname[counter], "w_%d_%d_%d", i + 1, j + 1, c + 1);

				//cout << "(" << j + 1 << ","  << i + 1 << ", " << c + 1 << "): " << position_w(inst, i, j, c) << " counter : " << counter << endl;

				counter++;
			}
		}
	}

	inst->obj[counter] = (double)1.0;

	if (inst->LB_heur >= 1)
	{
		inst->lb[counter] = inst->sol_val_DOWN[0];
	}
	else
	{
		inst->lb[counter] = -CPX_INFBOUND;
	}

	if (inst->UB_heur >= 1)
	{
		inst->ub[counter] = inst->sol_val_UP[0];
	}
	else
	{
		inst->ub[counter] = CPX_INFBOUND;
	}

	inst->c_type[counter] = 'C';

	sprintf(inst->colname[counter], "eta");

	counter++;

	inst->status = CPXnewcols(inst->env_MILP, inst->lp_MILP, inst->ccnt, inst->obj, inst->lb, inst->ub, inst->c_type, inst->colname);
	if (inst->status != 0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->c_type);

	for (int i = 0; i < inst->ccnt; i++) { free(inst->colname[i]); }
	free(inst->colname);

	// * setting the objective function in the maximization form
	CPXchgobjsen(inst->env_MILP, inst->lp_MILP, CPX_MAX);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	cout << "VARIABLE CREATED (" << inst->ccnt << ", " << counter << ")" << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	inst->zeta = new double*** [inst->n_sol];
	inst->zeta_i = new double*** [inst->n_sol];
	inst->zeta_j = new double*** [inst->n_sol];

	for (int ell = 0; ell < inst->n_sol; ell++) {

		inst->zeta[ell] = new double** [inst->n_individuals];
		inst->zeta_i[ell] = new double** [inst->n_individuals];
		inst->zeta_j[ell] = new double** [inst->n_individuals];

		for (int i = 0; i < inst->n_individuals; i++)
		{
			inst->zeta[ell][i] = new double* [inst->n_individuals];
			inst->zeta_i[ell][i] = new double* [inst->n_individuals];
			inst->zeta_j[ell][i] = new double* [inst->n_individuals];

			for (int j = 0; j < inst->n_individuals; j++)
			{
				inst->zeta[ell][i][j] = new double[inst->n_groups];
				inst->zeta_i[ell][i][j] = new double[inst->n_groups];
				inst->zeta_j[ell][i][j] = new double[inst->n_groups];
			}
		}
	}

	inst->x_ic = new double* [inst->n_individuals];

	for (int i = 0; i < inst->n_individuals; i++)
	{
		inst->x_ic[i] = new double[inst->n_groups];
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	build_obj_min_max(inst);

	cout << "Constraint 'min max' inserted \n";

	build_c1(inst);

	cout << "Constraint 'c1' inserted \n";

	build_constr_tau(inst);

	cout << "Constraint 'constr_tau' inserted \n";

	build_constr_w(inst);

	cout << "Constraint 'constr_w' inserted \n";

	build_constr_homogeneiry(inst);

	cout << "Constraint 'constr_homogeneiry' inserted \n";

	if (inst->MagnantiW == 1) {
		build_constr_min_max_first_cut(inst);
		cout << "Constraint 'constr_first_cut' inserted \n";
	}

	if (inst->external_sol <= 0)
	{
		//build_constr_BEN_min_sum_symmetry_breaking(inst);
		cout << "Constraint 'constr_symmetry_breaking' inserted \n";
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef write_prob
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * writing the created ILP model on a file *
	inst->status = CPXwriteprob(inst->env_MILP, inst->lp_MILP, "LP_MODEL_MINMAX.lp", NULL);

	if (inst->status != 0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif




};

void solve_MILP_min_max(instance* inst) {

	cout << endl;
	cout << "-------------------------------------------" << endl;
	cout << "Solving the MIN-MAX problem using B&C (Initial MG cut = " << inst->MagnantiW << ")" << endl;
	cout << "-------------------------------------------" << endl;

	CPXsetintparam(inst->env_MILP, CPX_PARAM_SCRIND, CPX_ON);

	//	// * Set relative tolerance *
	//	inst->status = CPXsetdblparam (inst->env_MILP, CPX_PARAM_EPAGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPAGAP\n");
	//	}
	//
	//	// * Set a tolerance *
	//	inst->status = CPXsetdblparam (inst->env_MILP, CPX_PARAM_EPGAP, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPGAP\n");
	//	}
	//
	//	// * Set mip tolerances integrality *
	//	inst->status = CPXsetdblparam (inst->env_MILP, CPX_PARAM_EPINT, 0.0);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPINTP\n");
	//	}
	//
	//	// * Set Feasibility tolerance *
	//	inst->status = CPXsetdblparam (inst->env_MILP, CPX_PARAM_EPRHS, 1e-9);
	//	if (inst->status)
	//	{
	//		printf ("error for CPX_PARAM_EPRHS\n");
	//	}g

	// * Set number of CPU*
	inst->status = CPXsetintparam(inst->env_MILP, CPX_PARAM_THREADS, inst->number_of_CPU);
	if (inst->status)
	{
		printf("error for CPX_PARAM_EPRHS\n");
	}

	// * Set time limit *
	inst->status = CPXsetdblparam(inst->env_MILP, CPX_PARAM_TILIM, inst->timelimit);
	if (inst->status)
	{
		printf("error for CPX_PARAM_EPRHS\n");
	}

	// * Set Memory limit *
	inst->status = CPXsetdblparam(inst->env_MILP, CPX_PARAM_TRELIM, 10000);
	if (inst->status)
	{
		printf("error for CPX_PARAM_TRELIM\n");
	}

	// * Set a tolerance *
	inst->status = CPXsetdblparam(inst->env_MILP, CPX_PARAM_EPGAP, inst->TOLL_OPTIMALITY);
	if (inst->status)
	{
		printf("error for CPX_PARAM_EPGAP\n");
	}

	if (inst->option == 2)
	{

		//this is the only one necessary to avoid the removal of all continuous variables
		CPXsetintparam(inst->env_MILP, CPX_PARAM_PREIND, CPX_OFF);
		//		CPXsetintparam(inst->env_MILP, CPX_PARAM_AGGIND, CPX_OFF);
		//		CPXsetintparam(inst->env_MILP, CPX_PARAM_BNDSTRENIND, CPX_OFF);
		//		CPXsetintparam(inst->env_MILP, CPX_PARAM_COEREDIND, CPX_OFF);
		//		CPXsetintparam(inst->env_MILP, CPX_PARAM_RELAXPREIND, CPX_OFF);
		//		CPXsetintparam(inst->env_MILP, CPX_PARAM_REDUCE, CPX_OFF);
		//		CPXsetintparam(inst->env_MILP, CPX_PARAM_PREPASS, CPX_OFF);
		//		CPXsetintparam(inst->env_MILP, CPX_PARAM_REPEATPRESOLVE, CPX_OFF);

		cout << "***********\n\n AUTOMATIC BENDER'S DECOMPOSITION\n\n";

		inst->status = CPXsetintparam(inst->env_MILP, CPXPARAM_Benders_Strategy, 3);
		if (inst->status)
		{
			printf("error for CPXPARAM_Benders_Strategy\n");
		}

		inst->status = CPXbendersopt(inst->env_MILP, inst->lp_MILP);

		//inst->status = CPXwritebendersannotation(inst->env_MILP, inst->lp_MILP,"TEST.ann");
		//if (inst->status)
		//{
		//	printf("error for CPXwritebendersannotation\n");
		//	exit(-1);
		//}
	}
	if (inst->option == 3)
	{
		//	int n_val;
		//	//n_val = inst->ccnt;
		//	n_val = 1;

		//	inst->rmatbeg = (int*)calloc(1, sizeof(int));
		//	inst->rmatind = (int*)calloc(n_val, sizeof(int));
		//	inst->rmatval = (double*)calloc(n_val, sizeof(double));
		//	inst->effortlevel = (int*)calloc(1, sizeof(int));

		//	inst->rmatbeg[0] = 0;
		//	inst->rmatind[0] = 0;
		//	inst->rmatval[0] = 0;
		//	inst->effortlevel[0] = 3;

		//	// https://www.ibm.com/docs/en/icos/22.1.0?topic=mip-starting-from-solution-starts
		//		
		//	cout << "Point 1" << endl;
		//	inst->status = CPXaddmipstarts(inst->env_MILP, inst->lp_MILP, 1, inst->ccnt, inst->rmatbeg, inst->rmatind, inst->rmatval, inst->effortlevel, NULL);
		//	if (inst->status)
		//	{
		//		printf("error for CPXaddmipstarts\n");
		//	}
		//	cout << "Point 2" << endl;
			//CPXsetintparam(inst->env_MILP, CPXPARAM_MIP_Strategy_FPHeur, 2);

		CPXsetintparam(inst->env_MILP, CPXPARAM_MIP_Strategy_HeuristicFreq, 2);
		cout << "***********\n\n FEASIBILITY PUMP\n\n";

	}

	//CPXsetintparam(inst->env_MILP, CPX_PARAM_LPMETHOD,4);

	///////////////////////////////////////////////////////////////////////////////////
	// * solving the MIP model
	clock_t time_start = clock();

	cout << "\nCPXmipopt:\n";
	inst->status = CPXmipopt(inst->env_MILP, inst->lp_MILP);
	if (inst->status != 0)
	{
		printf("error in CPXmipopt\n");
		//exit(-1);
	}
	clock_t time_end = clock();
	inst->solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;
	///////////////////////////////////////////////////////////////////////////////////


	bool sol_found = true;

	// * getting the solution
	inst->x = (double*)calloc(inst->ccnt, sizeof(double));


	inst->status = CPXgetmipx(inst->env_MILP, inst->lp_MILP, inst->x, 0, inst->ccnt - 1);
	if (inst->status != 0)
	{
		sol_found = false;
		printf("error in CPXgetmipx\n");
	}

	inst->objval = -1;
	inst->status = CPXgetmipobjval(inst->env_MILP, inst->lp_MILP, &(inst->objval));
	if (inst->status != 0)
	{
		sol_found = false;
		printf("error in CPXgetmipobjval\n");
	}

	//printf("\n\nMIP solution value ->\t\%f\n",inst->objval);

	//for (int c = 0; c < inst->n_groups; c++)
	//{
	//	for (int i = 0; i < inst->n_individuals; i++)
	//	{
	//		cout << "x(" << i + 1 << ", " << c + 1 << ") = " << inst->x[position_x(inst, i, c)] << endl;
	//	}
	//}

	///////////////////////////////////////////////////////////////////////////////////

	int cur_numcols = 0;
	int cur_numrows = 0;

	inst->bestobjval = -1;
	inst->status = CPXgetbestobjval(inst->env_MILP, inst->lp_MILP, &(inst->bestobjval));
	if (inst->status != 0)
	{
		printf("error in CPXgetbestobjval\n");
	}
	else {
		cur_numcols = CPXgetnumcols(inst->env_MILP, inst->lp_MILP);
		cur_numrows = CPXgetnumrows(inst->env_MILP, inst->lp_MILP);
	}

	inst->nodecount = 0;
	inst->lpstat = CPXgetstat(inst->env_MILP, inst->lp_MILP);
	inst->nodecount = CPXgetnodecnt(inst->env_MILP, inst->lp_MILP);

	cout << "\n\nlpstat\t" << inst->lpstat << endl;

	cout << "\n\nSTAT:\tobjval\t" << inst->objval << "\tbest-bound\t" << inst->bestobjval << "\tlpstat\t" << inst->lpstat << "\ttime\t" << inst->solution_time << endl << endl;

	printf("\nnumcols\t%d\n", cur_numcols);
	printf("\nnumrows\t%d\n", cur_numrows);

	//printf("\n\nMIP solution value ->\t\%f\n",inst->objval);

	//for (int c = 0; c < inst->n_groups; c++)
	//{
	//	for (int i = 0; i < inst->n_individuals; i++)
	//	{
	//		cout << "x(" << i + 1 << ", " << c + 1 << ") = " << inst->x[position_BEN_min_sum_x(inst, i, c)] << endl;
	//	}
	//}

	///////////////////////////////////////////////////////////////////////////////////

	//for (int c = 0; c < inst->n_groups; c++)
	//{
	//	for (int i = 0; i < inst->n_individuals; i++)
	//	{
	//		cout << "x(" << i + 1 << ", " << c + 1 << ") = " << inst->x[position_BEN_min_sum_x(inst, i, c)] <<  endl;
	//	}
	//}

	for (int s = 0; s < inst->n_sol; s++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			for (int i = 0; i < inst->n_individuals; i++)
			{
				inst->X_ic[s][i][c] = inst->x[position_x(inst, i, c)];
				inst->Heur_X_ic[s][i][c] = inst->x[position_x(inst, i, c)];
			}
		}
	}

};

void clean_MILP_min_max(instance* inst) {

	inst->status = CPXfreeprob(inst->env_MILP, &(inst->lp_MILP));
	if (inst->status != 0) { printf("error in CPXfreeprob\n"); exit(-1); }

	inst->status = CPXcloseCPLEX(&(inst->env_MILP));
	if (inst->status != 0) { printf("error in CPXcloseCPLEX\n"); exit(-1); }

};