
#define _CRT_SECURE_NO_WARNINGS

#define write_prob
#include "global_functions.h"
#include "global_variables.h"
#include "BEN_min_max.h"



//------------------------------------------------
// Variables positions
//------------------------------------------------


int position_BEN_min_max_eta(instance* inst)
{
	return (int)(inst->n_individuals * inst->n_groups);
}


int position_BEN_min_max_x(instance* inst, int i, int c)
{
	return  (c * (inst->n_individuals) + i);
}



//------------------------------------------------
// Constraints
//------------------------------------------------


void build_constr_BEN_min_max_c1(instance* inst)
{
	for (int i = 0; i < inst->n_individuals; i++)
	{
		// * creating the knapsack time 0 constraint *
		inst->rcnt = 1;
		inst->nzcnt = inst->n_groups;
		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

		inst->rhs[0] = 1.0;
		inst->sense[0] = 'E';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int c = 0; c < inst->n_groups; c++)
		{
			inst->rmatind[c] = position_BEN_min_max_x(inst, i, c);
			inst->rmatval[c] = (double)1.0;
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

void build_constr_BEN_min_max_tau(instance* inst)
{
	for (int c = 0; c < inst->n_groups; c++)
	{

		//--------------------------------------------------
		// Tau_min
		//--------------------------------------------------

		// * creating the knapsack time 0 constraint *
		inst->rcnt = 1;
		inst->nzcnt = inst->n_individuals;
		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

		inst->rhs[0] = inst->tau_min;
		inst->sense[0] = 'G';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[i] = position_BEN_min_max_x(inst, i, c);
			inst->rmatval[i] = 1;

		}

		inst->rmatbeg[0] = 0;

		inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
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

		//--------------------------------------------------
		// Tau_max
		//--------------------------------------------------

		// * creating the knapsack time 0 constraint *
		inst->rcnt = 1;
		inst->nzcnt = inst->n_individuals;
		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

		inst->rhs[0] = inst->tau_max;
		inst->sense[0] = 'L';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[i] = position_BEN_min_max_x(inst, i, c);
			inst->rmatval[i] = 1;

		}

		inst->rmatbeg[0] = 0;

		inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
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

void build_constr_BEN_min_max_homogeneiry(instance* inst)
{

	int count_var;
	inst->rcnt = 1;
	inst->nzcnt = 2 * (inst->n_individuals);

	for (int f = 0; f < inst->n_features; f++)
	{
		for (int c0 = 0; c0 < inst->n_groups; c0++)
		{
			for (int c1 = 0; c1 < inst->n_groups; c1++)
			{
				if (c1 == c0) {
					continue;
				}

				inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
				inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

				inst->rhs[0] = inst->theta_data[f];
				inst->sense[0] = 'L';

				inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
				inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
				inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

				count_var = 0;
				for (int i = 0; i < inst->n_individuals; i++) {

					inst->rmatind[count_var] = position_BEN_min_max_x(inst, i, c0);
					inst->rmatval[count_var] = (double)inst->beta_data[i][f];
					//inst->rmatval[count_var] = 5.0;
					//cout << "positive" << count_var << endl;

					count_var++;

				}
				for (int i = 0; i < inst->n_individuals; i++) {

					inst->rmatind[count_var] = position_BEN_min_max_x(inst, i, c1);
					inst->rmatval[count_var] = -1.0 * ((double)inst->beta_data[i][f]);

					//inst->rmatval[count_var] = -5.0;
					//cout << "negative" << count_var << endl;

					count_var++;

				}

				inst->rmatbeg[0] = 0;

				inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);
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
		}
	}
};

void build_constr_BEN_min_max_opt_bound(instance* inst)
{
	//int count_var;

	//inst->rcnt = 1;
	//inst->nzcnt = (inst->n_groups) * (inst->n_individuals) ;

	//inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	//inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

	//inst->rhs[0] = inst->opt_bound;
	//inst->sense[0] = 'G';

	//inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	//inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	//inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	//count_var = 0;

	//for (int c = 0; c < inst->n_groups; c++) {
	//	for (int i = 0; i < inst->n_individuals; i++) {

	//	}
	//}
	//inst->rmatbeg[0] = 0;

	//inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);

	//if (inst->status != 0)
	//{
	//	printf("error in CPXaddrows\n");
	//	exit(-1);
	//}

	//free(inst->rmatbeg);
	//free(inst->rmatval);
	//free(inst->rmatind);
	//free(inst->rhs);
	//free(inst->sense);

};

void build_constr_BEN_min_max_symmetry_breaking(instance* inst)
{
	int count_var;

	for (int c = 0; c < inst->n_groups - 1; c++)
	{
		count_var = 0;

		inst->rcnt = 1;
		inst->nzcnt = 2 * (inst->n_individuals);

		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

		inst->rhs[0] = 0;
		inst->sense[0] = 'L';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[count_var] = position_BEN_min_max_x(inst, i, c);
			inst->rmatval[count_var] = (double)i;

			count_var++;

		}
		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[count_var] = position_BEN_min_max_x(inst, i, c + 1);
			inst->rmatval[count_var] = (-1.0) * (double)i;

			count_var++;

		}

		inst->rmatbeg[0] = 0;

		inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);

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



void build_constr_BEN_min_max_Taylor_QA(instance* inst)
{
	inst->rcnt = 1;
	inst->nzcnt = 1 + (inst->n_individuals);

	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

	inst->sense[0] = 'L';

	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	double prim;
	double dual0;
	double dual1;
	int counter;

	for (int c = 0; c < inst->n_groups; c++)
	{
	//int c = 8;

		prim = 0.0;
		dual0 = 0.0;
		dual1 = 0.0;

		counter = 0;

		inst->rmatind[counter] = position_BEN_min_max_eta(inst);
		inst->rmatval[counter] = -(double)1.0;

		counter++;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			dual0 = 0.0;
			dual1 = 0.0;

			for (int j = 0; j < inst->n_individuals; j++)
			{
				if (j == i) {
					continue;
				}

				dual0 += (inst->alpha_data[i][j] * (double)inst->X_ic[0][i][c]);
				dual1 += (inst->alpha_data[j][i] * (double)inst->X_ic[0][j][c]);
				prim += (double)inst->alpha_data[j][i] * (double)inst->X_ic[0][i][c] * (double)inst->X_ic[0][j][c];
			}

			inst->rmatind[counter] = position_BEN_min_max_x(inst, i, c);
			inst->rmatval[counter] = (dual0 + dual1);

			counter++;

		}

		//inst->rhs[0] = 429.0;
		//inst->rhs[0] = 60.0;
		inst->rhs[0] = prim;
		//inst->rhs[0] = 0.0;

		inst->rmatbeg[0] = 0;

		inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);

		if (inst->status != 0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}


		
	}//for c

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);

};


void build_constr_BEN_min_max_first_cut(instance* inst)
{
	inst->rcnt = 1;
	inst->nzcnt = 1 + (inst->n_individuals);

	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

	inst->sense[0] = 'L';

	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	double dual0;
	double dual1;
	double dual2;

	int counter;

	get_zeta_first_cut(inst, 0.5);

	for (int c = 0; c < inst->n_groups; c++)
	{
		dual0 = 0.0;
		dual1 = 0.0;
		dual2 = 0.0;

		counter = 0;

		inst->rmatind[counter] = position_BEN_min_max_eta(inst);
		inst->rmatval[counter] = (double)1.0;

		counter++;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			dual0 = 0.0;
			dual1 = 0.0;

			for (int j = 0; j < inst->n_individuals; j++)
			{
				if (j == i) {
					continue;
				}

				dual0 = dual0 - (inst->zeta_i[0][i][j][c] - inst->zeta[0][i][j][c]);
				dual1 = dual1 - (inst->zeta_j[0][j][i][c] - inst->zeta[0][j][i][c]);
				dual2 = dual2 + (inst->zeta[0][i][j][c]);

			}

			inst->rmatind[counter] = position_BEN_min_max_x(inst, i, c);
			inst->rmatval[counter] = (dual0 + dual1);

			counter++;

		}

		inst->rhs[0] = dual2;

		inst->rmatbeg[0] = 0;

		inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);

		if (inst->status != 0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}

	}//for c



	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);



};

void build_constr_BEN_min_max_first_cut_multiple(instance* inst)
{

	get_zeta_first_cut_multiple(inst, 0.5);

	for (int ell = 0; ell < inst->n_sol; ell++)
	{
		inst->rcnt = 1;
		inst->nzcnt = 1 + (inst->n_individuals);

		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

		inst->sense[0] = 'L';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		double dual0;
		double dual1;
		double dual2;

		int counter;	

		for (int c = 0; c < inst->n_groups; c++)
		{
			dual0 = 0.0;
			dual1 = 0.0;
			dual2 = 0.0;

			counter = 0;

			inst->rmatind[counter] = position_BEN_min_max_eta(inst);
			inst->rmatval[counter] = (double)1.0;

			counter++;

			for (int i = 0; i < inst->n_individuals; i++)
			{
				dual0 = 0.0;
				dual1 = 0.0;

				for (int j = 0; j < inst->n_individuals; j++)
				{
					if (j == i) {
						continue;
					}

					dual0 = dual0 - (inst->zeta_i[ell][i][j][c] - inst->zeta[ell][i][j][c]);
					dual1 = dual1 - (inst->zeta_j[ell][j][i][c] - inst->zeta[ell][j][i][c]);
					dual2 = dual2 + (inst->zeta[ell][i][j][c]);

				}

				inst->rmatind[counter] = position_BEN_min_max_x(inst, i, c);
				inst->rmatval[counter] = (dual0 + dual1);

				counter++;

			}

			inst->rhs[0] = dual2;

			inst->rmatbeg[0] = 0;

			inst->status = CPXaddrows(inst->env_MILP, inst->lp_MILP, 0, inst->rcnt, inst->nzcnt, inst->rhs, inst->sense, inst->rmatbeg, inst->rmatind, inst->rmatval, NULL, NULL);

			if (inst->status != 0)
			{
				printf("error in CPXaddrows\n");
				exit(-1);
			}

		}//for c

		free(inst->rmatbeg);
		free(inst->rmatval);
		free(inst->rmatind);
		free(inst->rhs);
		free(inst->sense);

	}



};

/*****************************************************************/
int CPXPUBLIC mycutcallback_IMAST_BEN_min_max_single(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
/*****************************************************************/
{
	(*useraction_p) = CPX_CALLBACK_DEFAULT;

	instance* inst = (instance*)cbhandle;

	int status;

	double prim;

	double dual0;
	double dual1;
	double dual2;

	// get the array of master solution (variables v) at the current B&B node and store it in inst->sol_callback

	status = CPXgetcallbacknodex(env, cbdata, wherefrom, inst->sol_callback, 0, inst->n_variables_BEN - 1);  // can retrieve less than those
	if (status != 0) {
		printf("cannot get the v\n");
		exit(-1);
	}

	for (int i = 0; i < inst->n_individuals; i++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			inst->x_ic[i][c] = (double)inst->sol_callback[position_BEN_min_max_x(inst, i, c)];
		}
	}

	//------------------------------------------
	// Generate Lagrangian multipliers
	//------------------------------------------

	int counter;
	double eta_val = (double)inst->sol_callback[position_BEN_min_max_eta(inst)];


	//cout << "CALLBACK = " << endl;

	if (inst->external_sol >= 1) // no initial solution is used
	{
		get_zeta_MW(inst);
		//get_zeta_MW_multiple(inst);
	}
	else {
		get_zeta(inst, 0.5);
	}

	for (int c = 0; c < inst->n_groups; c++)
	{
		prim = 0.0;
		dual2 = 0.0;

		counter = 0;

		inst->cut_BEN_rmatind[counter] = position_BEN_min_max_eta(inst);
		inst->cut_BEN_rmatval[counter] = (double)1.0;

		counter++;

		//cut_val_primal = (double)0.0;
		//cut_val_dual = (double)0.0;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			dual0 = (double)0.0;
			dual1 = (double)0.0;

			for (int j = 0; j < inst->n_individuals; j++)
			{
				if (j == i) {
					continue;
				}
				prim += inst->x_ic[i][c] * inst->x_ic[j][c] * ((double)inst->alpha_data[i][j]);

				//dual0 = dual0 + inst->alpha_data[i][j] * inst->x_ic[i][c] * inst->zeta_i[i][j][c];
				//dual1 = dual1 + inst->alpha_data[j][i] * inst->x_ic[i][c] * inst->zeta_j[j][i][c];
				//dual2 = dual2 + inst->alpha_data[i][j] * ((1 - inst->x_ic[i][c] - inst->x_ic[j][c]) * inst->zeta[i][j][c]);

				//dual0 = dual0 + inst->alpha_data[i][j] * inst->x_ic[i][c] * (inst->zeta_i[i][j][c] - inst->zeta[i][j][c]);
				//dual1 = dual1 + inst->alpha_data[j][i] * inst->x_ic[i][c] * (inst->zeta_j[j][i][c] - inst->zeta[j][i][c]);
				//dual2 = dual2 + inst->alpha_data[i][j] * (inst->zeta[i][j][c]);

				dual0 -= (inst->zeta_i[0][i][j][c] - inst->zeta[0][i][j][c]);
				dual1 -= (inst->zeta_j[0][j][i][c] - inst->zeta[0][j][i][c]);
				dual2 += (inst->zeta[0][i][j][c]);

			}

			inst->cut_BEN_rmatind[counter] = position_BEN_min_max_x(inst, i, c);
			inst->cut_BEN_rmatval[counter] = (dual0 + dual1);

			counter++;

		}

		status = CPXcutcallbackadd(env, cbdata, wherefrom, inst->nzcnt_BEN, dual2, 'L', inst->cut_BEN_rmatind, inst->cut_BEN_rmatval, 1);
		if (status != 0)
		{
			printf("error in CPXcutcallbackadd\n");
			exit(-1);
		}


	}//for c



	inst->iteration_BEN++;

	////////////////////////////////////////////////////////

	//counter = 0;
	//prim = 0.0;

	//inst->cut_BEN_rmatind[counter] = position_BEN_min_sum_eta(inst);
	//inst->cut_BEN_rmatval[counter] = -(double)1.0;

	//counter++;

	//for (int c = 0; c < inst->n_groups; c++)
	//{
	//	for (int i = 0; i < inst->n_individuals; i++)
	//	{
	//		dual0 = 0.0;
	//		dual1 = 0.0;

	//		for (int j = 0; j < inst->n_individuals; j++)
	//		{
	//			if (j == i) {
	//				continue;
	//			}

	//			dual0 += (inst->alpha_data[i][j] * inst->x_ic[i][c]);
	//			dual1 += (inst->alpha_data[j][i] * inst->x_ic[j][c]);

	//			prim += inst->alpha_data[j][i] * inst->x_ic[i][c] * inst->x_ic[j][c];

	//		}

	//		inst->cut_BEN_rmatind[counter] = position_BEN_min_sum_x(inst, i, c);
	//		inst->cut_BEN_rmatval[counter] = (dual0 + dual1);

	//		counter++;

	//	}

	//}//for c

	//dual2 = (double)prim;

	//if (prim > eta_val) {
	//	status = CPXcutcallbackadd(env, cbdata, wherefrom, (1 + inst->n_groups * inst->n_individuals), dual2, 'L', inst->cut_BEN_Tay_rmatind, inst->cut_BEN_Tay_rmatval, 1);
	//	if (status != 0)
	//	{
	//		printf("error in CPXcutcallbackadd\n");
	//		exit(-1);
	//	}
	//}

	(*useraction_p) = CPX_CALLBACK_SET;

	if (inst->iteration_BEN > 10000000)
	{
		exit(-1);
	}

	return 0;

}


/*****************************************************************/
int CPXPUBLIC mycutcallback_IMAST_BEN_min_max(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
/*****************************************************************/
{
	(*useraction_p) = CPX_CALLBACK_DEFAULT;

	instance* inst = (instance*)cbhandle;

	int status;

	double prim;

	double dual0;
	double dual1;
	double dual2;

	// get the array of master solution (variables v) at the current B&B node and store it in inst->sol_callback

	status = CPXgetcallbacknodex(env, cbdata, wherefrom, inst->sol_callback, 0, inst->n_variables_BEN - 1);  // can retrieve less than those
	if (status != 0) {
		printf("cannot get the v\n");
		exit(-1);
	}

	for (int i = 0; i < inst->n_individuals; i++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			inst->x_ic[i][c] = (double)inst->sol_callback[position_BEN_min_max_x(inst, i, c)];
		}
	}

	//------------------------------------------
	// Generate Lagrangian multipliers
	//------------------------------------------

	int counter;
	double eta_val = (double)inst->sol_callback[position_BEN_min_max_eta(inst)];

	//cout << "CALLBACK = " << endl;
	if (inst->external_sol >= 1) // no initial solution is used
	{
		//get_zeta_MW(inst);
		get_zeta_MW_multiple(inst);
	}
	else {
		get_zeta(inst, 0.5);
	}
	for (int ell = 0; ell < inst->n_sol; ell++)
	{

		for (int c = 0; c < inst->n_groups; c++)
		{
			prim = 0.0;
			dual2 = 0.0;

			counter = 0;

			inst->cut_BEN_rmatind[counter] = position_BEN_min_max_eta(inst);
			inst->cut_BEN_rmatval[counter] = (double)1.0;

			counter++;

			//cut_val_primal = (double)0.0;
			//cut_val_dual = (double)0.0;

			for (int i = 0; i < inst->n_individuals; i++)
			{
				dual0 = (double)0.0;
				dual1 = (double)0.0;

				for (int j = 0; j < inst->n_individuals; j++)
				{
					if (j == i) {
						continue;
					}
					prim += inst->x_ic[i][c] * inst->x_ic[j][c] * ((double)inst->alpha_data[i][j]);

					//dual0 = dual0 + inst->alpha_data[i][j] * inst->x_ic[i][c] * inst->zeta_i[i][j][c];
					//dual1 = dual1 + inst->alpha_data[j][i] * inst->x_ic[i][c] * inst->zeta_j[j][i][c];
					//dual2 = dual2 + inst->alpha_data[i][j] * ((1 - inst->x_ic[i][c] - inst->x_ic[j][c]) * inst->zeta[i][j][c]);

					//dual0 = dual0 + inst->alpha_data[i][j] * inst->x_ic[i][c] * (inst->zeta_i[i][j][c] - inst->zeta[i][j][c]);
					//dual1 = dual1 + inst->alpha_data[j][i] * inst->x_ic[i][c] * (inst->zeta_j[j][i][c] - inst->zeta[j][i][c]);
					//dual2 = dual2 + inst->alpha_data[i][j] * (inst->zeta[i][j][c]);

					dual0 -= (inst->zeta_i[ell][i][j][c] - inst->zeta[ell][i][j][c]);
					dual1 -= (inst->zeta_j[ell][j][i][c] - inst->zeta[ell][j][i][c]);
					dual2 += (inst->zeta[ell][i][j][c]);

				}

				inst->cut_BEN_rmatind[counter] = position_BEN_min_max_x(inst, i, c);
				inst->cut_BEN_rmatval[counter] = (dual0 + dual1);

				counter++;

			}

			status = CPXcutcallbackadd(env, cbdata, wherefrom, inst->nzcnt_BEN, dual2, 'L', inst->cut_BEN_rmatind, inst->cut_BEN_rmatval, 1);
			if (status != 0)
			{
				printf("error in CPXcutcallbackadd\n");
				exit(-1);
			}
		}//for c
	}

	inst->iteration_BEN++;

	////////////////////////////////////////////////////////

	//counter = 0;
	//prim = 0.0;

	//inst->cut_BEN_rmatind[counter] = position_BEN_min_sum_eta(inst);
	//inst->cut_BEN_rmatval[counter] = -(double)1.0;

	//counter++;

	//for (int c = 0; c < inst->n_groups; c++)
	//{
	//	for (int i = 0; i < inst->n_individuals; i++)
	//	{
	//		dual0 = 0.0;
	//		dual1 = 0.0;

	//		for (int j = 0; j < inst->n_individuals; j++)
	//		{
	//			if (j == i) {
	//				continue;
	//			}

	//			dual0 += (inst->alpha_data[i][j] * inst->x_ic[i][c]);
	//			dual1 += (inst->alpha_data[j][i] * inst->x_ic[j][c]);

	//			prim += inst->alpha_data[j][i] * inst->x_ic[i][c] * inst->x_ic[j][c];

	//		}

	//		inst->cut_BEN_rmatind[counter] = position_BEN_min_sum_x(inst, i, c);
	//		inst->cut_BEN_rmatval[counter] = (dual0 + dual1);

	//		counter++;

	//	}

	//}//for c

	//dual2 = (double)prim;

	//if (prim > eta_val) {
	//	status = CPXcutcallbackadd(env, cbdata, wherefrom, (1 + inst->n_groups * inst->n_individuals), dual2, 'L', inst->cut_BEN_Tay_rmatind, inst->cut_BEN_Tay_rmatval, 1);
	//	if (status != 0)
	//	{
	//		printf("error in CPXcutcallbackadd\n");
	//		exit(-1);
	//	}
	//}

	(*useraction_p) = CPX_CALLBACK_SET;

	if (inst->iteration_BEN > 10000000)
	{
		exit(-1);
	}

	return 0;

}



//------------------------------------------------
// Build MILP
//------------------------------------------------

void build_BEN_min_max(instance* inst) {

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

	inst->ccnt = 1 + ((int)(inst->n_individuals * inst->n_groups));

	inst->obj = (double*)calloc(inst->ccnt, sizeof(double));
	inst->lb = (double*)calloc(inst->ccnt, sizeof(double));
	inst->ub = (double*)calloc(inst->ccnt, sizeof(double));
	inst->c_type = (char*)calloc(inst->ccnt, sizeof(char));

	inst->colname = (char**)calloc(inst->ccnt, sizeof(char*));
	for (int i = 0; i < inst->ccnt; i++) { inst->colname[i] = (char*)calloc(2000, sizeof(char)); }

	int counter = 0;
	bool test;

	for (int c = 0; c < inst->n_groups; c++)
	{
		test = true;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			inst->obj[counter] = (double)0.0;
			inst->lb[counter] = (double)0.0;
			inst->ub[counter] = (double)1.0;

			if (inst->external_sol > 0)
			{
				if ((inst->X_ic[0][i][c] >= 1) & test)
				{
					inst->lb[counter] = (double)inst->X_ic[0][i][c];
					inst->ub[counter] = (double)inst->X_ic[0][i][c];
					test = false;
				}
				else {
					inst->lb[counter] = (double)0.0;
					inst->ub[counter] = (double)1.0;
				}
			}

			inst->c_type[counter] = 'B';

			sprintf(inst->colname[counter], "x_%d_%d", i + 1, c + 1);

			//cout << "(" << i + 1 << ", " << c + 1 << "): " << position_BEN_min_sum_x(inst, i, c) << " counter : " << counter << endl;

			counter++;

		}
	}

	inst->obj[counter] = (double)1.0;
	if (inst->LB_heur >= 1)
	{
		inst->lb[counter] = inst->sol_val_DOWN[0];
	}
	else
	{
		inst->lb[counter] = -(inst->n_individuals*inst->n_groups);
	}

	if (inst->UB_heur >= 1)
	{
		inst->ub[counter] = inst->sol_val_UP[0];
	}
	else
	{
		inst->ub[counter] = (inst->n_individuals * inst->n_groups);
	}

	inst->c_type[counter] = 'C';

	//printf(inst->colname[counter], "eta");

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

	cout << "VARIABLE CREATED (" << inst->ccnt << ", " << counter << ")" << endl;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


	cout << "Constraint 'min max' inserted \n";

	build_constr_BEN_min_max_c1(inst);

	cout << "Constraint 'c1' inserted \n";

	build_constr_BEN_min_max_tau(inst);

	cout << "Constraint 'constr_tau' inserted \n";

	build_constr_BEN_min_max_homogeneiry(inst);

	cout << "Constraint 'constr_homogeneiry' inserted \n";

	if (inst->MagnantiW == 1) {
		build_constr_BEN_min_max_first_cut(inst);
		cout << "Constraint 'constr_first_cut' inserted \n";
	}
	if (inst->MagnantiW == 2) 
	{
		build_constr_BEN_min_max_first_cut_multiple(inst);
		cout << "Constraint 'constr_first_cut' inserted \n";
	}

	//if (inst->external_sol <= 0)
	//{
	//	//build_constr_BEN_min_sum_symmetry_breaking(inst);
	//	cout << "Constraint 'constr_symmetry_breaking' inserted \n";
	//}
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

void solve_BEN_min_max(instance* inst) {

	cout << endl;
	cout << "-------------------------------------------" << endl;
	cout << "Solving the MIN-MAX problem using Benders decomposition (Initial MG cut = "<< inst->MagnantiW <<")" << endl;
	cout << "-------------------------------------------" << endl;

	inst->iteration_BEN = 0;

	inst->n_variables_BEN = 1 + ((int)(inst->n_individuals * inst->n_groups));

	inst->sol_callback = new double[inst->n_variables_BEN];

	inst->nzcnt_BEN = 1 + (int)(inst->n_individuals);

	inst->cut_BEN_rmatind = new int[inst->nzcnt_BEN];
	inst->cut_BEN_rmatval = new double[inst->nzcnt_BEN];

	inst->cut_BEN_Tay_rmatind = new int[inst->nzcnt_BEN];
	inst->cut_BEN_Tay_rmatval = new double[inst->nzcnt_BEN];

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
	//	}

	// * Set number of CPU*
	inst->status = CPXsetintparam(inst->env_MILP, CPX_PARAM_THREADS, inst->number_of_CPU);
	if (inst->status)
	{
		printf("error for CPX_PARAM_THREADS\n");
	}

	// * Set time limit *
	inst->status = CPXsetdblparam(inst->env_MILP, CPX_PARAM_TILIM, inst->timelimit);
	if (inst->status)
	{
		printf("error for CPX_PARAM_TILIM\n");
	}

	// * Set Memory limit *
	//inst->status =  CPXsetdblparam(inst->env_MILP,CPX_PARAM_TRELIM,10000);
	//if (inst->status)
	//{
	//	printf ("error for CPX_PARAM_TRELIM\n");
	//}

	// * Set a tolerance *
	inst->status = CPXsetdblparam(inst->env_MILP, CPX_PARAM_EPGAP, inst->TOLL_OPTIMALITY);
	if (inst->status)
	{
		printf("error for CPX_PARAM_EPGAP\n");
	}

	//--------------------------------------------------------
	// Emphasize the achievement of a feasible solution
	//--------------------------------------------------------

	if (inst->get_feasible)
	{
		// * Set mip tolerances integrality *
		inst->status = CPXsetdblparam(inst->env_MILP, CPX_PARAM_EPINT, 1e-3);
		if (inst->status)
		{
			printf("error for CPX_PARAM_EPINTP\n");
		}

		// * Set Feasibility tolerance *
		//inst->status = CPXsetdblparam(inst->env_MILP_LB, CPX_PARAM_EPRHS, 0.01);
		//if (inst->status)
		//{
		//	printf("error for CPX_PARAM_EPRHS\n");
		//}

		//inst->status = CPXsetintparam(inst->env_MILP, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_FEASIBILITY);
		//if (inst->status) {
		//	printf("Error setting CPX_PARAM_MIPEMPHASIS\n");
		//}

		//inst->status = CPXsetintparam(inst->env_MILP_LB, CPX_PARAM_HEURFREQ, 2);
		//if (inst->status) {
		//	printf("Error setting CPX_PARAM_HEURFREQ\n");
		//}

	}
	else {
	
		// * Set mip tolerances integrality *
		inst->status = CPXsetdblparam(inst->env_MILP, CPX_PARAM_EPINT, 1e-7);
		if (inst->status)
		{
			printf("error for CPX_PARAM_EPINTP\n");
		}

		inst->status = CPXsetintparam(inst->env_MILP, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);
		if (inst->status) {
			printf("Error setting CPX_PARAM_MIPEMPHASIS\n");
		}

	}

	//--------------------------------------------------------
	// Warm start
	//--------------------------------------------------------

	if (inst->warm_start)
	{
		//// Define the warm start
		//int num_vars = ((int)(inst->n_individuals * inst->n_groups));		// Number of variables with non-zero values
		//int* indices = new int[num_vars];		// Array of variable indices
		//double* solution = new double[num_vars];	// Array of corresponding variable values

		//// Parameters for CPXaddmipstarts
		//int mcnt = 1;                   // One MIP start
		//int nzcnt = num_vars;           // Number of variables in this start
		//int beg[] = { 0 };                // Start of the first (and only) MIP start
		//int effortlevel[] = { CPX_MIPSTART_AUTO }; // Let CPLEX decide the effort level
		//char** mipstartname = NULL;     // No names for the MIP start

		//int counter;

		//for (int ell = 0; ell < inst->n_sol; ell++)
		//{

		//	counter = 0;

		//	for (int i = 0; i < inst->n_individuals; i++)
		//	{
		//		for (int c = 0; c < inst->n_groups; c++)
		//		{
		//			indices[counter] = position_x(inst, i, c);
		//			solution[counter] = (double)inst->X_ic[ell][i][c];
		//			counter++;
		//		}
		//	}

		//	// Add the MIP start to the model
		//	inst->status = CPXaddmipstarts(
		//		inst->env_MILP,
		//		inst->lp_MILP,
		//		mcnt,
		//		nzcnt,
		//		beg,
		//		indices,
		//		solution,
		//		effortlevel,
		//		mipstartname
		//	);

		//	if (inst->status) {
		//		printf("Error adding MIP start: %d\n", inst->status);
		//	}
		//}

		// Define the warm starts
		int num_vars = ((int)(inst->n_individuals * inst->n_groups)); // Number of variables with non-zero values per MIP start
		int n_sol = inst->n_sol; // Number of MIP starts
		int total_nzcnt = num_vars * n_sol; // Total number of variables across all MIP starts

		// Allocate memory for all MIP starts
		int* indices = new int[total_nzcnt]; // Indices for all MIP starts
		double* solution = new double[total_nzcnt]; // Values for all MIP starts
		int* beg = new int[n_sol]; // Starting points for each MIP start
		int* effortlevel = new int[n_sol]; // Effort levels for each MIP start

		// Populate the arrays
		int counter = 0;
		for (int ell = 0; ell < n_sol; ell++) {
			beg[ell] = counter; // Start of this MIP start
			effortlevel[ell] = CPX_MIPSTART_AUTO; // Use default effort level

			for (int i = 0; i < inst->n_individuals; i++) {
				for (int c = 0; c < inst->n_groups; c++) {
					indices[counter] = position_x(inst, i, c);
					solution[counter] = (double)inst->X_ic[ell][i][c];
					counter++;
				}
			}
		}

		// Add all MIP starts to the model in a single call
		inst->status = CPXaddmipstarts(
			inst->env_MILP,
			inst->lp_MILP,
			n_sol,       // Number of MIP starts
			total_nzcnt, // Total number of non-zero values
			beg,         // Start indices of each MIP start
			indices,     // Variable indices
			solution,    // Variable values
			effortlevel, // Effort levels
			NULL         // No names for MIP starts
		);

		if (inst->status) {
			printf("Error adding MIP starts: %d\n", inst->status);
		}

		// Clean up dynamically allocated memory
		delete[] indices;
		delete[] solution;
		delete[] beg;
		delete[] effortlevel;


	}

	///////////////////////////////////////////////////////////////////////////////////
	//
	// The following chunck of code is preventing CPLEX to pre-reduce the set of variables
	// right before returning the current solution. In fact, we need to preserve the lengh of 
	// of the vector of variables. By default CPLEX would drop the variables which are not needed.
	//
	///////////////////////////////////////////////////////////////////////////////////

	CPXsetintparam(inst->env_MILP, CPX_PARAM_MIPCBREDLP, CPX_OFF);			// let MIP callbacks work on the original model
	CPXsetintparam(inst->env_MILP, CPX_PARAM_PRELINEAR, CPX_OFF);			// assure linear mappings between the presolved and original models
	//CPXsetintparam(inst->env_LP_MODEL_IBENDERS, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);

	///////////////////////////////////////////////////////////////////////////////////
	//
	// CPLEX callback:
	//
	// The following chunck of code is asking CPLEX to suspend the procedure every time  
	// it gets an integer solution (for the case CPXsetlazyconstraintcallbackfunc) and
	// every time it gets a fractional solution (for the case CPXsetusercutcallbackfunc)
	///////////////////////////////////////////////////////////////////////////////////


	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_MILP, mycutcallback_IMAST_BEN_min_max, inst);
	//inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_MILP, mycutcallback_IMAST_BEN_min_max_single, inst);
	if (inst->status)
	{
		printf("error for CPXsetlazyconstraintcallbackfunc\n");
	}

	/*if (inst->option > 1)
	{

		cout << "******FRACTIONAL SEPARATION******\n\n";

		inst->status = CPXsetusercutcallbackfunc(inst->env_MILP, myusercutcallback_IBEN, inst);

		if (inst->status)
		{
			printf("error for CPXsetuserconstraintcallbackfunc\n");
		}

	}*/

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

//#ifdef write_prob
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//// * writing the created ILP model on a file *
//    inst->status = CPXwriteprob(inst->env_MILP, inst->lp_MILP, "IP_MODEL_IBENDERS.lp", NULL);
//    if (inst->status != 0)
//    {
//        printf("error in CPXwriteprob\n");
//        exit(-1);
//    }
//    //exit(-1);
//    /////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//#endif

	bool sol_found = true;

	// * getting the solution
	inst->x = (double*)calloc(inst->n_variables_BEN, sizeof(double));
	memset(inst->x, 0, inst->n_variables_BEN * sizeof(double));

	inst->status = CPXgetmipx(inst->env_MILP, inst->lp_MILP, inst->x, 0, inst->n_variables_BEN - 1);
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
	cout << "\n\nSTAT:\tobjval\t" << inst->objval << "\tbestobjval\t" << inst->bestobjval << "\tlpstat\t" << inst->lpstat << "\ttime\t" << inst->solution_time << endl << endl;

	printf("\nnumcols\t%d\n", cur_numcols);
	printf("\nnumrows\t%d\n", cur_numrows);

	//for (int c = 0; c < inst->n_groups; c++)
	//{
	//	for (int i = 0; i < inst->n_individuals; i++)
	//	{
	//		cout << "x(" << i + 1 << ", " << c + 1 << ") = " << inst->x[position_BEN_min_max_x(inst, i, c)] << endl;
	//	}
	//}

	for (int s = 0; s < inst->n_sol; s++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			for (int i = 0; i < inst->n_individuals; i++)
			{
				inst->X_ic[s][i][c] = inst->x[position_BEN_min_max_x(inst, i, c)];
				inst->Heur_X_ic[s][i][c] = inst->x[position_BEN_min_max_x(inst, i, c)];
			}
		}
	}

};

void clean_BEN_min_max(instance* inst) {

	inst->status = CPXfreeprob(inst->env_MILP, &(inst->lp_MILP));
	if (inst->status != 0) { printf("error in CPXfreeprob\n"); exit(-1); }

	inst->status = CPXcloseCPLEX(&(inst->env_MILP));
	if (inst->status != 0) { printf("error in CPXcloseCPLEX\n"); exit(-1); }

	delete[] inst->sol_callback;
	delete[] inst->cut_BEN_rmatind;
	delete[] inst->cut_BEN_rmatval;
	delete[] inst->cut_BEN_Tay_rmatind;
	delete[] inst->cut_BEN_Tay_rmatval;

	inst->sol_callback = nullptr;
	inst->cut_BEN_rmatind = nullptr;
	inst->cut_BEN_rmatval = nullptr;
	inst->cut_BEN_Tay_rmatind = nullptr;
	inst->cut_BEN_Tay_rmatval = nullptr;

};