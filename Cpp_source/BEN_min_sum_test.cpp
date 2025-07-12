
#define _CRT_SECURE_NO_WARNINGS

#define write_prob
#include <algorithm>
#include "global_functions.h"
#include "global_variables.h"
#include "BEN_min_sum_test.h"


//------------------------------------------------
// Variables positions
//------------------------------------------------


int position_BEN_min_sum_eta_test(instance* inst, int c)
{
	return (int)(c + (inst->n_individuals * inst->n_groups));
}


int position_BEN_min_sum_x_test(instance* inst, int i, int c)
{
	return  (c * (inst->n_individuals) + i);
}



//------------------------------------------------
// Constraints
//------------------------------------------------


void build_constr_BEN_min_sum_c1_test(instance* inst)
{
	for (int i = 0; i < inst->n_individuals; i++)
	{
		// * creating the knapsack time 0 constraint *
		inst->rcnt = 1;
		inst->nzcnt = inst->n_groups;
		inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
		inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

		inst->rhs[0] = (double)1.0;
		inst->sense[0] = 'E';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int c = 0; c < inst->n_groups; c++)
		{
			inst->rmatind[c] = position_BEN_min_sum_x_test(inst, i, c);
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

void build_constr_BEN_min_sum_tau_test(instance* inst)
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

		inst->rhs[0] = (double)inst->tau_min;
		inst->sense[0] = 'G';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[i] = position_BEN_min_sum_x_test(inst, i, c);
			inst->rmatval[i] = (double)1.0;

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

		inst->rhs[0] = (double)inst->tau_max;
		inst->sense[0] = 'L';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[i] = position_BEN_min_sum_x_test(inst, i, c);
			inst->rmatval[i] = (double)1.0;

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

void build_constr_BEN_min_sum_homogeneiry_test(instance* inst)
{

	int count_var;

	for (int f = 0; f < inst->n_features; f++)
	{
		for (int c0 = 0; c0 < inst->n_groups; c0++)
		{
			for (int c1 = 0; c1 < inst->n_groups; c1++)
			{
				if (c1 == c0) {
					continue;
				}

				inst->rcnt = 1;
				inst->nzcnt = 2 * (inst->n_individuals);
				inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
				inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

				inst->rhs[0] = inst->theta_data[f];
				inst->sense[0] = 'L';

				inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
				inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
				inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

				count_var = 0;
				for (int i = 0; i < inst->n_individuals; i++) {

					inst->rmatind[count_var] = position_BEN_min_sum_x_test(inst, i, c0);
					inst->rmatval[count_var] = (double)inst->beta_data[i][f];
					//inst->rmatval[count_var] = 5.0;
					//cout << "positive" << count_var << endl;

					count_var++;

				}
				for (int i = 0; i < inst->n_individuals; i++) {

					inst->rmatind[count_var] = position_BEN_min_sum_x_test(inst, i, c1);
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

void build_constr_BEN_min_sum_first_cut_test(instance* inst)
{
	inst->rcnt = 1;
	inst->nzcnt = 1 + inst->n_individuals;

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
		counter = 0;

		inst->rmatind[counter] = position_BEN_min_sum_eta_test(inst, c);
		inst->rmatval[counter] = (double)1.0;

		counter++;

		dual2 = 0.0;

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

			inst->rmatind[counter] = position_BEN_min_sum_x_test(inst, i, c);
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

void build_constr_BEN_min_sum_symmetry_breaking_test(instance* inst)
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

			inst->rmatind[count_var] = position_BEN_min_sum_x_test(inst, i, c);
			inst->rmatval[count_var] = (double)i;

			count_var++;

		}
		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[count_var] = position_BEN_min_sum_x_test(inst, i, c + 1);
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

void build_constr_BEN_min_sum_symmetry_breaking_QA_test(instance* inst)
{
	int count_var;
	double RHS;

	inst->rcnt = 1;
	inst->nzcnt = (inst->n_individuals);

	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

	inst->sense[0] = 'G';

	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	for (int c = 0; c < inst->n_groups; c++)
	{
		count_var = 0;
		RHS = (double)0.0;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			RHS += (double)i * (double)inst->X_ic[0][i][c];
		}

		inst->rhs[0] = RHS;

		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[count_var] = position_BEN_min_sum_x_test(inst, i, c);
			inst->rmatval[count_var] = (double)i;

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

void build_constr_BEN_min_sum_Taylor_QA_test(instance* inst)
{
	inst->rcnt = 1;
	inst->nzcnt = 1 + ( inst->n_individuals);

	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

	inst->sense[0] = 'L';

	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	double prim = 0.0;
	double dual0 = 0.0;
	double dual1 = 0.0;

	int counter;

	

	for (int c = 0; c < inst->n_groups; c++)
	{
		counter = 0;
		prim = 0.0;

		inst->rmatind[counter] = position_BEN_min_sum_eta_test(inst, c);
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

			inst->rmatind[counter] = position_BEN_min_sum_x_test(inst, i, c);
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

void build_constr_BEN_min_sum_tau_QA_test(instance* inst)
{

	double RHS;

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

		RHS = 0;
		for (int i = 0; i < inst->n_individuals; i++) {
			RHS = RHS + (double)inst->X_ic[0][i][c];
		}

		inst->rhs[0] = (double)RHS;
		inst->sense[0] = 'E';

		inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
		inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
		inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[i] = position_BEN_min_sum_x_test(inst, i, c);
			inst->rmatval[i] = (double)1.0;

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


/*****************************************************************/
int CPXPUBLIC mycutcallback_IMAST_BEN_min_sum_test(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
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
			inst->x_ic[i][c] = (double)inst->sol_callback[position_BEN_min_sum_x_test(inst, i, c)];
		}
	}

	//------------------------------------------
	// Generate Lagrangian multipliers
	//------------------------------------------

	int counter;
	double eta_val;

	//get_zeta(inst,0.5);
	get_zeta_MW(inst);

	for (int c = 0; c < inst->n_groups; c++)
	{
		counter = 0;
		
		eta_val = (double)inst->sol_callback[position_BEN_min_sum_eta_test(inst, c)];

		inst->cut_BEN_rmatind[counter] = position_BEN_min_sum_eta_test(inst, c);
		inst->cut_BEN_rmatval[counter] = (double)1.0;

		counter++;

		dual2 = 0.0;
		prim = 0.0;

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

			inst->cut_BEN_rmatind[counter] = position_BEN_min_sum_x_test(inst, i, c);
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

	//inst->cut_BEN_rmatind[counter] = position_BEN_min_sum_eta_test(inst);
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

	//		inst->cut_BEN_rmatind[counter] = position_BEN_min_sum_x_test(inst, i, c);
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

///*****************************************************************/
//int CPXPUBLIC mycutcallback_FMAST_BEN_min_sum(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
///*****************************************************************/
//{
//	//if (wherefrom == CPX_CALLBACK_MIP_CUT_FEAS) {
//	//	// generate cutting planes in integer nodes
//	//}
//	int status;
//	double prim;
//	double dual0;
//	double dual1;
//	double dual2;
//
//	(*useraction_p) = CPX_CALLBACK_DEFAULT;
//
//	instance* inst = (instance*)cbhandle;
//
//	// get the array of master solution (variables v) at the current B&B node and store it in inst->sol_callback
//
//	status = CPXgetcallbacknodex(env, cbdata, wherefrom, inst->sol_callback, 0, inst->n_variables_BEN - 1);  // can retrieve less than those
//	if (status != 0) {
//		printf("cannot get the v\n");
//		exit(-1);
//	}
//
//	for (int i = 0; i < inst->n_individuals; i++)
//	{
//		for (int c = 0; c < inst->n_groups; c++)
//		{
//			inst->x_ic[i][c] = (double)inst->sol_callback[position_BEN_min_sum_x_test(inst, i, c)];
//		}
//	}
//
//	int counter;
//	double eta_val = (double)inst->sol_callback[position_BEN_min_sum_eta_test(inst)];
//
//	counter = 0;
//
//	counter = 0;
//	prim = 0.0;
//
//	inst->cut_BEN_Tay_rmatind[counter] = position_BEN_min_sum_eta_test(inst);
//	inst->cut_BEN_Tay_rmatval[counter] = -(double)1.0;
//
//	counter++;
//
//	for (int c = 0; c < inst->n_groups; c++)
//	{
//		for (int i = 0; i < inst->n_individuals; i++)
//		{
//			dual0 = 0.0;
//			dual1 = 0.0;
//
//			for (int j = 0; j < inst->n_individuals; j++)
//			{
//				if (j == i) {
//					continue;
//				}
//
//				dual0 += (inst->alpha_data[i][j] * inst->x_ic[i][c]);
//				dual1 += (inst->alpha_data[j][i] * inst->x_ic[j][c]);
//
//				prim += inst->alpha_data[j][i] * inst->x_ic[i][c] * inst->x_ic[j][c];
//
//			}
//
//			inst->cut_BEN_Tay_rmatind[counter] = position_BEN_min_sum_x_test(inst, i, c);
//			inst->cut_BEN_Tay_rmatval[counter] = (dual0 + dual1);
//
//			counter++;
//
//		}
//
//	}//for c
//
//	dual2 = (double)prim;
//
//	if (prim >= eta_val)
//	{
//		status = CPXcutcallbackadd(env, cbdata, wherefrom, (1 + inst->n_groups * inst->n_individuals), dual2, 'L', inst->cut_BEN_Tay_rmatind, inst->cut_BEN_Tay_rmatval, 1);
//		if (status != 0)
//		{
//			printf("error in CPXcutcallbackadd\n");
//			exit(-1);
//		}
//	}
//
//	(*useraction_p) = CPX_CALLBACK_SET;
//
//	//if (inst->iteration_BEN > 10000000) {
//	//	exit(-1);
//	//}
//
//	return 0;
//
//}
//

//------------------------------------------------
// Build - Solve - Clean
//------------------------------------------------

void build_BEN_min_sum_test(instance* inst) {

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

	//inst->incumbent_old = (double)(-1000 * pow(inst->n_individuals * inst->n_groups, 2));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * creating the variables and linear part of the objective function*
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	inst->ccnt = inst->n_groups + ((int)(inst->n_individuals * inst->n_groups));

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

			}

			inst->c_type[counter] = 'B';

			sprintf(inst->colname[counter], "x_%d_%d", i + 1, c + 1);

			//cout << "(" << i + 1 << ", " << c + 1 << "): " << position_BEN_min_sum_x_test(inst, i, c) << " counter : " << counter << endl;

			counter++;

		}
	}

	for (int c = 0; c < inst->n_groups; c++)
	{
		inst->obj[counter] = (double)1.0;
		inst->lb[counter] = -CPX_INFBOUND;
		//inst->ub[counter] = (double)(1000*pow(inst->n_individuals*inst->n_groups,2));
		inst->ub[counter] = (double)50.0;
		inst->c_type[counter] = 'C';

		//sprintf(inst->colname[counter], "eta_%d", c + 1);
		sprintf(inst->colname[counter], "eta");

		counter++;

	}

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

	build_constr_BEN_min_sum_c1_test(inst);

	cout << "Constraint 'c1' inserted \n";
	
	build_constr_BEN_min_sum_tau_test(inst);

	cout << "Constraint 'constr_tau' inserted \n";

	build_constr_BEN_min_sum_homogeneiry_test(inst);

	cout << "Constraint 'constr_homogeneiry' inserted \n";

	if (inst->external_sol > 0) {

		build_constr_BEN_min_sum_first_cut_test(inst);

		//cout << "Constraint 'constr_first_cut' inserted \n";

		build_constr_BEN_min_sum_Taylor_QA_test(inst);

		//cout << "Constraint 'constr_Taylor_QA' inserted \n";

		//build_constr_BEN_min_sum_tau_QA(inst);

		//cout << "Constraint 'constr_tau' inserted \n";

		//build_constr_BEN_min_sum_symmetry_breaking_QA(inst);

		//cout << "Constraint 'constr_symmetry_breaking_QA' inserted \n";

	}
	else {

		//build_constr_BEN_min_sum_symmetry_breaking(inst);

		cout << "Constraint 'constr_symmetry_breaking' inserted \n";

	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef write_prob
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// * writing the created ILP model on a file *
	inst->status = CPXwriteprob(inst->env_MILP, inst->lp_MILP, "LP_MODEL_MINSUM.lp", NULL);

	if (inst->status != 0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif




};

void solve_BEN_min_sum_test(instance* inst) {

	inst->iteration_BEN = 0;

	inst->n_variables_BEN = inst->n_groups + ((int)(inst->n_individuals * inst->n_groups));

	inst->sol_callback = new double[inst->n_variables_BEN];

	inst->nzcnt_BEN = 1 + (inst->n_individuals);

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

	int status = CPXsetintparam(inst->env_MILP, CPX_PARAM_VARSEL, CPX_VARSEL_STRONG);
	if (status != 0) {
		// Handle error
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

	inst->status = CPXsetlazyconstraintcallbackfunc(inst->env_MILP, mycutcallback_IMAST_BEN_min_sum_test, inst);

	if (inst->status)
	{
		printf("error for CPXsetlazyconstraintcallbackfunc\n");
	}

	//inst->status = CPXsetusercutcallbackfunc(inst->env_MILP, mycutcallback_FMAST_BEN_min_sum, inst);

	if (inst->status)
	{
		printf("error for CPXsetusercutcallbackfunc\n");
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
	double solution_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;
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

	int cur_numcols = CPXgetnumcols(inst->env_MILP, inst->lp_MILP);
	int cur_numrows = CPXgetnumrows(inst->env_MILP, inst->lp_MILP);

	inst->bestobjval = -1;
	inst->status = CPXgetbestobjval(inst->env_MILP, inst->lp_MILP, &(inst->bestobjval));
	if (inst->status != 0)
	{
		printf("error in CPXgetbestobjval\n");
	}

	inst->lpstat = CPXgetstat(inst->env_MILP, inst->lp_MILP);
	inst->nodecount = CPXgetnodecnt(inst->env_MILP, inst->lp_MILP);

	cout << "\n\nlpstat\t" << inst->lpstat << endl;

	cout << "\n\nSTAT:\tobjval\t" << inst->objval << "\tbestobjval\t" << inst->bestobjval << "\tlpstat\t" << inst->lpstat << "\ttime\t" << solution_time << endl << endl;

	printf("\nnumcols\t%d\n", cur_numcols);
	printf("\nnumrows\t%d\n", cur_numrows);

	//printf("\n\nMIP solution value ->\t\%f\n",inst->objval);

	//for (int c = 0; c < inst->n_groups; c++)
	//{
	//	for (int i = 0; i < inst->n_individuals; i++)
	//	{
	//		cout << "x(" << i + 1 << ", " << c + 1 << ") = " << inst->x[position_BEN_min_sum_x_test(inst, i, c)] << endl;
	//	}
	//}

	///////////////////////////////////////////////////////////////////////////////////

	//for (int c = 0; c < inst->n_groups; c++)
	//{
	//	for (int i = 0; i < inst->n_individuals; i++)
	//	{
	//		cout << "x(" << i + 1 << ", " << c + 1 << ") = " << inst->x[position_BEN_min_sum_x_test(inst, i, c)] <<  endl;
	//	}
	//}

	//cout << "ETA = " << inst->x[position_BEN_min_sum_eta_test(inst)] << endl;


	//ofstream compact_file;
	//compact_file.open("info_Results.txt", ios::app);
	//compact_file << fixed
	//	<< inst->input_file_data << "\t"
	//	<< inst->n_scenarios << "\t"
	//	<< inst->n_category << "\t"
	//	<< inst->n_levels << "\t"
	//	<< inst->objval << "\t"
	//	<< inst->bestobjval << "\t"
	//	<< inst->lpstat << "\t"
	//	<< solution_time << "\t"
	//	<< inst->nodecount << "\t"
	//	<< inst->algorithm << "\t"
	//	<< inst->option << "\t"
	//	<< inst->TOLL_OPTIMALITY << "\t"
	//	<< cur_numcols << "\t"
	//	<< cur_numrows << "\t"

	//	<< endl;
	//compact_file.close();


	free(inst->x);


};

void clean_BEN_min_sum_test(instance* inst) {

	inst->status = CPXfreeprob(inst->env_MILP, &(inst->lp_MILP));
	if (inst->status != 0) { printf("error in CPXfreeprob\n"); exit(-1); }

	inst->status = CPXcloseCPLEX(&(inst->env_MILP));
	if (inst->status != 0) { printf("error in CPXcloseCPLEX\n"); exit(-1); }

};