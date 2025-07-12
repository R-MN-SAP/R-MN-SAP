
#include "global_variables.h"
#include "global_functions.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "math.h"
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string> 
#include <limits>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <omp.h>  // OpenMP header for parallel programming


using namespace std;

std::mt19937 gen(21051983); // Fixed seed


//------------------------------------------------
// Read or generate data functions
//------------------------------------------------

void generate_random_data(instance* inst)
{

	// Somewhere in your main or initialization function:
	srand(static_cast<unsigned int>(123)); // Seed the random number generator

	//-------------------------------------------
	// Alpha data
	//-------------------------------------------

	inst->alpha_data = new int* [inst->n_individuals];
	int random_value;

	for (int i = 0; i < inst->n_individuals; i++) {

		inst->alpha_data[i] = new int[inst->n_individuals];

		for (int j = 0; j < inst->n_individuals; j++) {

			if (i == j) {
				inst->alpha_data[i][j] = 0;
			}
			else {
				// Generate a random number between 0 and 2, then map it to {-1, 0, 1}
				random_value = rand() % 3 - 1;
				inst->alpha_data[i][j] = random_value;
			}


		}
	}

	// Write the array to a text file
	std::ofstream outfile_alpha("Random_alpha.txt");  // Open a file named "alpha_data.txt" for writing

	if (outfile_alpha.is_open())
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int j = 0; j < inst->n_individuals; j++)
			{
				int value = inst->alpha_data[i][j];
				if (value != 0)
				{
					outfile_alpha << (i + 1) << " " << (j + 1) << " " << value << std::endl;  // Output in the format (i, j, val)
				}
			}
		}
		outfile_alpha.close();  // Close the file
		std::cout << "Array written to Random_alpha.txt successfully." << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing." << std::endl;
	}


	//-------------------------------------------
	// Beta data
	//-------------------------------------------

	inst->beta_data = new int* [inst->n_individuals];

	for (int i = 0; i < inst->n_individuals; i++)
	{
		inst->beta_data[i] = new int[inst->n_features];

		for (int f = 0; f < inst->n_features; f++)
		{

			random_value = rand() % 2;
			inst->beta_data[i][f] = random_value; // Assigns a random value of 0 or 1
		}
	}

	// Write the array to a text file
	std::ofstream outfile_beta("Random_beta.txt");  // Open a file named "alpha_data.txt" for writing

	if (outfile_beta.is_open())
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int f = 0; f < inst->n_features; f++)
			{
				outfile_beta << inst->beta_data[i][f] << " ";
			}
			outfile_beta << endl;

		}
		outfile_beta.close();  // Close the file
		cout << "Array written to Random_beta.txt successfully." << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing." << std::endl;
	}


	//-------------------------------------------
	// Theta data
	//-------------------------------------------

	inst->tau_max = (int)round((double)(inst->n_individuals / inst->n_groups)) + 1;
	inst->tau_min = (int)round((double)(inst->n_individuals / inst->n_groups)) - 1;

	inst->theta_data = new int[inst->n_features];

	for (int f = 0; f < inst->n_features; f++) {

		inst->theta_data[f] = inst->tau_min - 1;

	}

	// Write the array to a text file
	std::ofstream outfile_theta("Random_theta.txt");  // Open a file named "alpha_data.txt" for writing

	if (outfile_theta.is_open())
	{
		for (int f = 0; f < inst->n_features; f++)
		{
			outfile_theta << inst->theta_data[f] << endl;
		}

		outfile_theta << inst->tau_min << endl;
		outfile_theta << inst->tau_max << endl;

		outfile_theta.close();  // Close the file
		cout << "Array written to Random_theta.txt successfully." << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing." << std::endl;
	}


	//-------------------------------------------
	// X data
	//-------------------------------------------

	// Initialize adjacency matrix and list of neighbors
	inst->X_ic = new int** [inst->n_sol];

	// Initialize adjacency matrix and list of neighbors
	inst->X_ic[0] = new int* [inst->n_individuals];

	for (int i = 0; i < inst->n_individuals; i++)
	{
		inst->X_ic[0][i] = new int[inst->n_groups];
	}

	for (int i = 0; i < inst->n_individuals; i++)
	{
		inst->X_ic[0][i] = new int[inst->n_groups];

		for (int c = 0; c < inst->n_groups; c++)
		{
			inst->X_ic[0][i][c] = 0;
		}
	}

	// Write the array to a text file
	std::ofstream outfile_x("Random_x.txt");  // Open a file named "alpha_data.txt" for writing

	for (int i = 0; i < inst->n_individuals; i++)
	{
		int random_group = rand() % (inst->n_groups);
		inst->X_ic[0][i][random_group] = 1;

		outfile_x << (i + 1) << " " << random_group + 1 << endl;

	}

	outfile_x.close();  // Close the file

}

void read_data_alpha(instance* inst) {

	// Initialize adjacency matrix and list of neighbors

	inst->alpha_data = new int* [inst->n_individuals];

	for (int i = 0; i < inst->n_individuals; i++) {

		inst->alpha_data[i] = new int[inst->n_individuals];

		for (int j = 0; j < inst->n_individuals; j++) {

			inst->alpha_data[i][j] = 0;
		}
	}
	// number of neighbors

	inst->n_neighbors_up = new int[inst->n_individuals];
	inst->n_neighbors_down = new int[inst->n_individuals];

	// open file for reading data

	char filename[256];
	//strcpy_s(filename, sizeof(filename), "alpha_2.txt");
	strcpy_s(filename, sizeof(filename), inst->alpha_file_data);

	//inst->alpha_file_data

	ifstream graphfile(filename);

	//ifstream graphfile("alpha_2.txt");
	if (!graphfile) {
		cerr << "Can not open the file " << inst->alpha_file_data << endl;
		exit(-1);
	}

	string line;

	// Read the non-zero entries from the file and store them in the dense matrix

	int ii, jj, value;

	while (getline(graphfile, line)) {

		stringstream ss(line);

		ss >> ii >> jj >> value;

		inst->alpha_data[ii - 1][jj - 1] = value;

		//cout << "(" << ii << ", " << jj << ") = "<< inst->alpha_data[ii - 1][jj - 1] << endl;
	}

	graphfile.close();

	int accumulate_up;
	int accumulate_down;

	for (int i = 0; i < inst->n_individuals; i++) {

		accumulate_up = 0;
		accumulate_down = 0;

		for (int j = 0; j < inst->n_individuals; j++) {

			if (inst->alpha_data[i][j] > 0)
			{
				accumulate_up++;
			}
			if (inst->alpha_data[i][j] <= 0)
			{
				accumulate_down++;
			}
		}
		inst->n_neighbors_up[i] = (int)accumulate_up;
		inst->n_neighbors_down[i] = (int)accumulate_down;

		//cout << "Number of friends of " << i + 1 << " is " << inst->n_neighbors_up[i] << endl;
		//cout << "Number of enemies of " << i + 1 << " is " << inst->n_neighbors_down[i] << endl;

	}

	// Initialize the list of neighbors

	inst->neighbors_up = new int* [inst->n_individuals];
	inst->neighbors_down = new int* [inst->n_individuals];

	int neighbor_index_up;
	int neighbor_index_down;

	// Fill in the list of neighbors
	for (int i = 0; i < inst->n_individuals; i++) {

		inst->neighbors_up[i] = new int[inst->n_neighbors_up[i]];
		inst->neighbors_down[i] = new int[inst->n_neighbors_down[i]];

		neighbor_index_up = 0;
		neighbor_index_down = 0;

		for (int j = 0; j < inst->n_individuals; j++) {

			if (inst->alpha_data[i][j] > 0)
			{
				if (neighbor_index_up >= inst->n_neighbors_up[i]) {
					cerr << "Friend index is larger than the number of friends for node " << i << endl;
					exit(-1);
				}
				else {
					inst->neighbors_up[i][neighbor_index_up] = j;
					neighbor_index_up++;
				}
			}
			if (inst->alpha_data[i][j] <= 0)
			{
				if (neighbor_index_down >= inst->n_neighbors_down[i]) {
					cerr << "Enemy index is larger than the number of enemies for node " << i << endl;
					exit(-1);
				}
				else {
					inst->neighbors_down[i][neighbor_index_down] = j;
					neighbor_index_down++;
				}
			}
		}
	}

}

void read_data_x(instance* inst) {

	// open file for reading data

	char filename[256];
	strcpy_s(filename, sizeof(filename), inst->x_file_data);

	ifstream graphfile(filename);

	if (!graphfile) {
		cerr << "Can not open the file " << inst->x_file_data << endl;
		exit(-1);
	}

	string line;

	// Read the non-zero entries from the file and store them in the dense matrix

	getline(graphfile, line);

	stringstream ss(line);

	ss >> inst->n_sol;

	// Initialize adjacency matrix and list of neighbors
	inst->X_ic = new int** [inst->n_sol];

	for (int s = 0; s < inst->n_sol; s++)
	{
		// Initialize adjacency matrix and list of neighbors
		inst->X_ic[s] = new int* [inst->n_individuals];

		for (int i = 0; i < inst->n_individuals; i++)
		{
			inst->X_ic[s][i] = new int[inst->n_groups];
		}

		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int c = 0; c < inst->n_groups; c++)
			{
				inst->X_ic[s][i][c] = 0;
			}
		}
	}

	// Populate the matrix
	int ii, cc;
	int sol_number = 0;

	while (getline(graphfile, line)) {

		if (sol_number >= inst->n_sol)
		{
			break;
		}
		//cout << "sol = " << sol_number << endl;

		stringstream ss(line);

		ss >> ii >> cc;

		//cout << "sol = " << sol_number << " i = " << ii << " c = " << cc << " X = " << inst->X_ic[sol_number][ii - 1][cc - 1] << endl;

		inst->X_ic[sol_number][ii - 1][cc - 1] = 1;

		//inst->X_ic[ii-1][cc-1] = 1;

		//cout << "(" << ii << ", " << jj << ") = "<< inst->alpha_data[ii - 1][jj - 1] << endl;

		if (ii == inst->n_individuals)
		{
			sol_number = sol_number + 1;
		}

	}

	graphfile.close();

	inst->sol_val_EXTERNAL = new int[inst->n_sol];

	if (inst->problem <= 0) {

		for (int s = 0; s < inst->n_sol; s++)
		{
			inst->sol_val_EXTERNAL[s] = 0;

			for (int i = 0; i < inst->n_individuals; i++)
			{
				for (int j = 0; j < inst->n_individuals; j++)
				{
					for (int c = 0; c < inst->n_groups; c++)
					{
						inst->sol_val_EXTERNAL[s] = inst->sol_val_EXTERNAL[s] + inst->alpha_data[i][j] * inst->X_ic[s][i][c] * inst->X_ic[s][j][c];
					}
				}
			}
		}

	}

	if (inst->problem >= 1) {

		for (int s = 0; s < inst->n_sol; s++)
		{
			inst->sol_val_EXTERNAL[s] = 0;

			int sol_val0 = 0;
			int sol_val1 = 0;

			for (int i = 0; i < inst->n_individuals; i++)
			{
				for (int j = 0; j < inst->n_individuals; j++)
				{
					sol_val0 = sol_val0 + inst->alpha_data[i][j] * inst->X_ic[s][i][0] * inst->X_ic[s][j][0];
				}
			}

			for (int c = 1; c < inst->n_groups; c++)
			{
				sol_val1 = 0;

				for (int i = 0; i < inst->n_individuals; i++)
				{
					for (int j = 0; j < inst->n_individuals; j++)
					{
						sol_val1 = sol_val1 + inst->alpha_data[i][j] * inst->X_ic[s][i][c] * inst->X_ic[s][j][c];
					}
				}

				if (sol_val1 < sol_val0) {
					sol_val0 = sol_val1;
				}
			}

			inst->sol_val_EXTERNAL[s] = sol_val0;
		}

	}

}

void read_data_beta(instance* inst)
{
	// define the dimensions of the v_data matrix

	inst->beta_data = new int* [inst->n_individuals];

	for (int i = 0; i < inst->n_individuals; i++) {

		inst->beta_data[i] = new int[inst->n_features];

		for (int f = 0; f < inst->n_features; f++) {

			inst->beta_data[i][f] = 0;

		}
	}

	char filename[256];
	//strcpy_s(filename, sizeof(filename), "beta_test.txt");
	strcpy_s(filename, sizeof(filename), inst->beta_file_data);

	//char* filename = inst->beta_file_data;

	FILE* file;

	if (fopen_s(&file, filename, "r") != 0) {
		perror("Error opening beta file");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < inst->n_individuals; i++)
	{
		for (int f = 0; f < inst->n_features; f++)
		{
			if (fscanf_s(file, "%d", &inst->beta_data[i][f]) == EOF) {
				perror("Error reading beta file");
				exit(EXIT_FAILURE);
			}
		}
	}

	fclose(file);
	printf("Beta file closed\n");

	//for (int i = 0; i < inst->n_individuals; i++) {

	//	for (int f = 0; f < inst->n_features; f++) {
	//		cout << "Reading beta ss: " << "(" << i + 1 << ", " << f + 1 << "):  " << inst->beta_data[i][f] << endl;
	//	}

	//}// for i


}

void read_data_theta(instance* inst) {

	// define the dimensions of the v_data matrix

	inst->theta_data = new int[inst->n_features];

	char filename[256];
	//strcpy_s(filename, sizeof(filename), "theta_test.txt");
	strcpy_s(filename, sizeof(filename), inst->theta_file_data);

	ifstream thetafile(filename);

	if (!thetafile) {

		cerr << "can not open the file " << inst->theta_file_data << endl;
		exit(-1);
	}

	string line;

	for (int f = 0; f < inst->n_features; f++) {

		getline(thetafile, line);

		stringstream ss(line);

		ss >> inst->theta_data[f];

		ss.clear();
		//cout << "Reading theta" << "(" << f + 1 << "):  " << inst->theta_data[f] << endl;

	}

	getline(thetafile, line);

	stringstream ss1(line);

	ss1 >> inst->tau_min;

	ss1.clear();

	getline(thetafile, line);

	stringstream ss2(line);

	ss2 >> inst->tau_max;

	ss2.clear();

	thetafile.close();

	cout << "Tau_min: " << inst->tau_min << endl;
	cout << "Tau_max: " << inst->tau_max << endl;

	cout << "I have just closed the theta file " << endl;

}





//------------------------------------------------
// Variables positions
//------------------------------------------------

int position_eta(instance* inst)
{
	return ((int)(pow(inst->n_individuals, 2) * inst->n_groups) + inst->n_individuals * inst->n_groups);
}

int position_x(instance* inst, int i, int c)
{
	return  (c * (inst->n_individuals) + i);

}

int position_w(instance* inst, int i, int j, int c)
{

	int start_w = inst->n_groups * inst->n_individuals;

	return  ((int)(start_w + c * ((int)pow(inst->n_individuals, 2)) + i * inst->n_individuals + j));

}





//------------------------------------------------
// Constraints
//------------------------------------------------


void build_c1(instance* inst)
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
			inst->rmatind[c] = position_x(inst, i, c);
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

void build_constr_tau(instance* inst)
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

			inst->rmatind[i] = position_x(inst, i, c);
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

			inst->rmatind[i] = position_x(inst, i, c);
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

void build_constr_w(instance* inst)
{
	for (int c = 0; c < inst->n_groups; c++)
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int j = 0; j < inst->n_individuals; j++)
			{

				if (i != j) {
					//--------------------------------------------------
					// w_{i,j}^c <= x_{i,c}
					//--------------------------------------------------

					// * creating the knapsack time 0 constraint *
					inst->rcnt = 1;
					inst->nzcnt = 2;
					inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
					inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

					inst->rhs[0] = 0;
					inst->sense[0] = 'L';

					inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
					inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
					inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

					inst->rmatind[0] = position_w(inst, i, j, c);
					inst->rmatval[0] = 1;

					inst->rmatind[1] = position_x(inst, i, c);
					inst->rmatval[1] = -1;

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
					// w_{i,j}^c <= x_{j,c}
					//--------------------------------------------------

					// * creating the knapsack time 0 constraint *
					inst->rcnt = 1;
					inst->nzcnt = 2;
					inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
					inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

					inst->rhs[0] = 0;
					inst->sense[0] = 'L';

					inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
					inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
					inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

					inst->rmatind[0] = position_w(inst, i, j, c);
					inst->rmatval[0] = 1;

					inst->rmatind[1] = position_x(inst, j, c);
					inst->rmatval[1] = -1;

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
					// w_{i,j}^c >= x_{i,c} + x_{j,c} - 1
					//--------------------------------------------------

					// * creating the knapsack time 0 constraint *
					inst->rcnt = 1;
					inst->nzcnt = 3;
					inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
					inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

					inst->rhs[0] = -1;
					inst->sense[0] = 'G';

					inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
					inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
					inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

					inst->rmatind[0] = position_w(inst, i, j, c);
					inst->rmatval[0] = 1;

					inst->rmatind[1] = position_x(inst, i, c);
					inst->rmatval[1] = -1;

					inst->rmatind[2] = position_x(inst, j, c);
					inst->rmatval[2] = -1;

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
				else {

					//--------------------------------------------------
					// w_{i,i}^c = x_{i,c}
					//--------------------------------------------------

					// * creating the knapsack time 0 constraint *
					inst->rcnt = 1;
					inst->nzcnt = 2;
					inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
					inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

					inst->rhs[0] = 0;
					inst->sense[0] = 'E';

					inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
					inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
					inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

					inst->rmatind[0] = position_w(inst, i, i, c);
					inst->rmatval[0] = 1;

					inst->rmatind[1] = position_x(inst, i, c);
					inst->rmatval[1] = -1;

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
	}
};

void build_constr_homogeneiry(instance* inst)
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

					inst->rmatind[count_var] = position_x(inst, i, c0);
					inst->rmatval[count_var] = (double)inst->beta_data[i][f];
					//inst->rmatval[count_var] = 5.0;
					//cout << "positive" << count_var << endl;

					count_var++;

				}
				for (int i = 0; i < inst->n_individuals; i++) {

					inst->rmatind[count_var] = position_x(inst, i, c1);
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

void build_constr_opt_bound(instance* inst)
{
	int count_var;

	inst->rcnt = 1;
	inst->nzcnt = (inst->n_groups) * (inst->n_individuals) * (inst->n_individuals);

	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

	inst->rhs[0] = inst->opt_bound;
	inst->sense[0] = 'G';

	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	count_var = 0;

	for (int c = 0; c < inst->n_groups; c++) {
		for (int i = 0; i < inst->n_individuals; i++) {
			for (int j = 0; j < inst->n_individuals; j++) {

				inst->rmatind[count_var] = position_w(inst, i, j, c);
				inst->rmatval[count_var] = (double)inst->alpha_data[i][j];

				count_var++;

			}
		}
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

};

void build_constr_symmetry_breaking(instance* inst)
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

			inst->rmatind[count_var] = position_x(inst, i, c);
			inst->rmatval[count_var] = (double)i;

			count_var++;

		}
		for (int i = 0; i < inst->n_individuals; i++) {

			inst->rmatind[count_var] = position_x(inst, i, c + 1);
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

void build_constr_min_sum_first_cut(instance* inst)
{

	inst->rcnt = 1;
	inst->nzcnt = 1 + (inst->n_groups * inst->n_individuals);

	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

	inst->sense[0] = 'L';

	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	double dual0;
	double dual1;
	double dual2;

	int counter = 0;
	dual2 = 0.0;

	get_zeta_first_cut(inst, (double)0.5);

	inst->rmatind[counter] = position_eta(inst);
	inst->rmatval[counter] = (double)1.0;

	counter++;

	for (int c = 0; c < inst->n_groups; c++)
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			dual0 = (double)0.0;
			dual1 = (double)0.0;

			for (int j = 0; j < inst->n_individuals; j++)
			{
				if (j == i) {
					continue;
				}

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

			inst->rmatind[counter] = position_x(inst, i, c);
			inst->rmatval[counter] = (dual0 + dual1);

			counter++;

		}

	}//for c

	inst->rhs[0] = dual2;

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

};

void build_constr_min_max_first_cut(instance* inst)
{
	inst->rcnt = 1;
	inst->nzcnt = 1 + (inst->n_individuals);

	inst->rhs = (double*)calloc(inst->rcnt, sizeof(double));
	inst->sense = (char*)calloc(inst->rcnt, sizeof(double));

	inst->sense[0] = 'L';

	double dual0;
	double dual1;
	double dual2;

	int counter;

	get_zeta_first_cut(inst, 0.5);

	inst->rmatbeg = (int*)calloc(inst->rcnt, sizeof(int));
	inst->rmatind = (int*)calloc(inst->nzcnt, sizeof(int));
	inst->rmatval = (double*)calloc(inst->nzcnt, sizeof(double));

	for (int c = 0; c < inst->n_groups; c++)
	{
		dual0 = 0.0;
		dual1 = 0.0;
		dual2 = 0.0;

		counter = 0;

		inst->rmatind[counter] = position_eta(inst);
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

			inst->rmatind[counter] = position_x(inst, i, c);
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




//------------------------------------------------
// Lagrangian multipliers
//------------------------------------------------


void get_zeta_first_cut(instance* inst, double rr) {

	for (int c = 0; c < inst->n_groups; c++)
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int j = 0; j < inst->n_individuals; j++)
			{
				inst->zeta_i[0][i][j][c] = 0.0;
				inst->zeta_j[0][i][j][c] = 0.0;
				inst->zeta[0][i][j][c] = 0.0;

				if (j == i) {
					continue;
				}
				else {

					if ((inst->X_ic[0][i][c] + inst->X_ic[0][j][c]) <= 1)
					{
						//if (1 >= 0)
						if (inst->alpha_data[i][j] >= 0)
						{
							if ((inst->X_ic[0][i][c] + inst->X_ic[0][j][c]) >= 1) {
								inst->zeta_i[0][i][j][c] = 0.0;
								inst->zeta_j[0][i][j][c] = 0.0;
								inst->zeta[0][i][j][c] = inst->alpha_data[i][j];
							}
							else {
								inst->zeta_i[0][i][j][c] = inst->alpha_data[i][j];
								inst->zeta_j[0][i][j][c] = inst->alpha_data[i][j];
								inst->zeta[0][i][j][c] = 0.0;
							}
						}
						else
						{
							if ((inst->X_ic[0][i][c] + inst->X_ic[0][j][c]) >= 1) {
								inst->zeta_i[0][i][j][c] = 0.0;
								inst->zeta_j[0][i][j][c] = 0.0;
								inst->zeta[0][i][j][c] = -inst->alpha_data[i][j];
							}
							else {
								inst->zeta_i[0][i][j][c] = -inst->alpha_data[i][j];
								inst->zeta_j[0][i][j][c] = -inst->alpha_data[i][j];
								inst->zeta[0][i][j][c] = 0.0;
							}
						}

					}
					else {

						//rr = rand() % 1000;
						if (inst->alpha_data[i][j] >= 0)
						{
							inst->zeta_i[0][i][j][c] = (double)(rr)*inst->alpha_data[i][j];
							inst->zeta_j[0][i][j][c] = (double)(1 - rr) * inst->alpha_data[i][j];
							inst->zeta[0][i][j][c] = (double)0.0;
						}
						else {

							inst->zeta_i[0][i][j][c] = (double)0.0;
							inst->zeta_j[0][i][j][c] = (double)0.0;
							inst->zeta[0][i][j][c] = -inst->alpha_data[i][j];

						}
					}
				}
			}
		}
	}


}

void get_zeta_first_cut_multiple(instance* inst, double rr)
{
	for (int ell = 0; ell < inst->n_sol; ell++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			for (int i = 0; i < inst->n_individuals; i++)
			{
				for (int j = 0; j < inst->n_individuals; j++)
				{
					inst->zeta_i[ell][i][j][c] = 0.0;
					inst->zeta_j[ell][i][j][c] = 0.0;
					inst->zeta[ell][i][j][c] = 0.0;

					if (j == i) {
						continue;
					}
					else {

						if ((inst->X_ic[ell][i][c] + inst->X_ic[ell][j][c]) <= 1)
						{
							//if (1 >= 0)
							if (inst->alpha_data[i][j] >= 0)
							{
								if ((inst->X_ic[ell][i][c] + inst->X_ic[ell][j][c]) >= 1) {
									inst->zeta_i[ell][i][j][c] = 0.0;
									inst->zeta_j[ell][i][j][c] = 0.0;
									inst->zeta[ell][i][j][c] = inst->alpha_data[i][j];
								}
								else {
									inst->zeta_i[ell][i][j][c] = inst->alpha_data[i][j];
									inst->zeta_j[ell][i][j][c] = inst->alpha_data[i][j];
									inst->zeta[ell][i][j][c] = 0.0;
								}
							}
							else
							{
								if ((inst->X_ic[ell][i][c] + inst->X_ic[ell][j][c]) >= 1) {
									inst->zeta_i[ell][i][j][c] = 0.0;
									inst->zeta_j[ell][i][j][c] = 0.0;
									inst->zeta[ell][i][j][c] = -inst->alpha_data[i][j];
								}
								else {
									inst->zeta_i[ell][i][j][c] = -inst->alpha_data[i][j];
									inst->zeta_j[ell][i][j][c] = -inst->alpha_data[i][j];
									inst->zeta[ell][i][j][c] = 0.0;
								}
							}

						}
						else {

							//rr = rand() % 1000;
							if (inst->alpha_data[i][j] >= 0)
							{
								inst->zeta_i[ell][i][j][c] = (double)(rr)*inst->alpha_data[i][j];
								inst->zeta_j[ell][i][j][c] = (double)(1 - rr) * inst->alpha_data[i][j];
								inst->zeta[ell][i][j][c] = (double)0.0;
							}
							else {

								inst->zeta_i[ell][i][j][c] = (double)0.0;
								inst->zeta_j[ell][i][j][c] = (double)0.0;
								inst->zeta[ell][i][j][c] = -inst->alpha_data[i][j];

							}
						}
					}
				}
			}
		}
	}


}

void get_zeta(instance* inst, double rr)
{

	for (int c = 0; c < inst->n_groups; c++)
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int j = 0; j < inst->n_individuals; j++)
			{
				inst->zeta_i[0][i][j][c] = 0.0;
				inst->zeta_j[0][i][j][c] = 0.0;
				inst->zeta[0][i][j][c] = 0.0;

				if (j == i) {
					continue;
				}

				if ((inst->x_ic[i][c] + inst->x_ic[j][c]) <= 1.01)
				{
					if (inst->alpha_data[i][j] >= 0)
					{
						if ((inst->x_ic[i][c] + inst->x_ic[j][c]) >= 0.9) {
							inst->zeta_i[0][i][j][c] = 0.0;
							inst->zeta_j[0][i][j][c] = 0.0;
							inst->zeta[0][i][j][c] = inst->alpha_data[i][j];
						}
						else {
							inst->zeta_i[0][i][j][c] = inst->alpha_data[i][j];
							inst->zeta_j[0][i][j][c] = inst->alpha_data[i][j];
							inst->zeta[0][i][j][c] = 0.0;
						}
					}
					else
					{
						if ((inst->x_ic[i][c] + inst->x_ic[j][c]) >= 1) {
							inst->zeta_i[0][i][j][c] = 0.0;
							inst->zeta_j[0][i][j][c] = 0.0;
							inst->zeta[0][i][j][c] = -inst->alpha_data[i][j];
						}
						else {
							inst->zeta_i[0][i][j][c] = -inst->alpha_data[i][j];
							inst->zeta_j[0][i][j][c] = -inst->alpha_data[i][j];
							inst->zeta[0][i][j][c] = 0.0;
						}
					}

				}
				else {

					//rr = rand() % 1000;
					if (inst->alpha_data[i][j] >= 0.0)
					{
						inst->zeta_i[0][i][j][c] = (double)(rr)*inst->alpha_data[i][j];
						inst->zeta_j[0][i][j][c] = (double)(1 - rr) * inst->alpha_data[i][j];
						inst->zeta[0][i][j][c] = (double)0.0;
					}
					else {

						inst->zeta_i[0][i][j][c] = (double)0.0;
						inst->zeta_j[0][i][j][c] = (double)0.0;
						inst->zeta[0][i][j][c] = -inst->alpha_data[i][j];

					}
				}
			}
		}
	}


}

void get_zeta_MW(instance* inst) {

	for (int c = 0; c < inst->n_groups; c++)
	{
		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int j = 0; j < inst->n_individuals; j++)
			{
				inst->zeta_i[0][i][j][c] = (double)0.0;
				inst->zeta_j[0][i][j][c] = (double)0.0;
				inst->zeta[0][i][j][c] = (double)0.0;

				if (j == i) {
					continue;
				}

				if (((inst->x_ic[i][c] + inst->x_ic[j][c]) <= (double)1.1))
				{
					if (inst->alpha_data[i][j] >= (double)0.0)
					{
						// CASE 1: either x_ic = 1 or x_jc = 1 (but not both)
						if (inst->x_ic[i][c] >= (double)0.99) {

							inst->zeta_i[0][i][j][c] = (double)0.0;
							inst->zeta_j[0][i][j][c] = (double)inst->alpha_data[i][j];
							inst->zeta[0][i][j][c] = (double)0.0;
						}
						if (inst->x_ic[j][c] >= (double)0.99) {

							inst->zeta_i[0][i][j][c] = (double)inst->alpha_data[i][j];
							inst->zeta_j[0][i][j][c] = (double)0.0;
							inst->zeta[0][i][j][c] = (double)0.0;
						}

						// CASE 2: both x_ic = 0 and x_jc = 0
						if (((inst->x_ic[i][c] + inst->x_ic[j][c]) <= (double)0.01))
						{
							if (inst->X_ic[0][i][c] <= inst->X_ic[0][j][c]) {

								inst->zeta_i[0][i][j][c] = (double)inst->alpha_data[i][j];
								inst->zeta_j[0][i][j][c] = (double)0.0;
								inst->zeta[0][i][j][c] = (double)0.0;
							}
							else
							{
								inst->zeta_i[0][i][j][c] = (double)0.0;
								inst->zeta_j[0][i][j][c] = (double)inst->alpha_data[i][j];
								inst->zeta[0][i][j][c] = (double)0.0;
							}
						}
					}
					else
					{
						inst->zeta_i[0][i][j][c] = (double)0.0;
						inst->zeta_j[0][i][j][c] = (double)0.0;
						inst->zeta[0][i][j][c] = (double)0.0;
					}
				}

				// CASE 2: both x_ic = 1 and x_jc = 1

				if (((inst->x_ic[i][c] + inst->x_ic[j][c]) >= (double)1.99))
				{

					//rr = rand() % 1000;
					//if (1 >= 0)
					if (inst->alpha_data[i][j] >= (double)0.0)
					{
						if (inst->X_ic[0][i][c] <= inst->X_ic[0][j][c])
						{
							inst->zeta_i[0][i][j][c] = (double)inst->alpha_data[i][j];
							inst->zeta_j[0][i][j][c] = (double)0.0;
							inst->zeta[0][i][j][c] = (double)0.0;
						}
						else {
							inst->zeta_j[0][i][j][c] = (double)inst->alpha_data[i][j];
							inst->zeta_i[0][i][j][c] = (double)0.0;
							inst->zeta[0][i][j][c] = (double)0.0;
						}
					}
					else {

						inst->zeta_i[0][i][j][c] = (double)0.0;
						inst->zeta_j[0][i][j][c] = (double)0.0;
						inst->zeta[0][i][j][c] = -inst->alpha_data[i][j];

					}
				}
			}
		}
	}
}

void get_zeta_MW_multiple(instance* inst) 
{
	for (int ell = 0; ell < inst->n_sol; ell++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			for (int i = 0; i < inst->n_individuals; i++)
			{
				for (int j = 0; j < inst->n_individuals; j++)
				{
					inst->zeta_i[ell][i][j][c] = (double)0.0;
					inst->zeta_j[ell][i][j][c] = (double)0.0;
					inst->zeta[ell][i][j][c] = (double)0.0;

					if (j == i) {
						continue;
					}

					if (((inst->x_ic[i][c] + inst->x_ic[j][c]) <= (double)1.1))
					{
						if (inst->alpha_data[i][j] >= (double)0.0)
						{
							// CASE 1: either x_ic = 1 or x_jc = 1 (but not both)
							if (inst->x_ic[i][c] >= (double)0.99) {

								inst->zeta_i[ell][i][j][c] = (double)0.0;
								inst->zeta_j[ell][i][j][c] = (double)inst->alpha_data[i][j];
								inst->zeta[ell][i][j][c] = (double)0.0;
							}
							if (inst->x_ic[j][c] >= (double)0.99) {

								inst->zeta_i[ell][i][j][c] = (double)inst->alpha_data[i][j];
								inst->zeta_j[ell][i][j][c] = (double)0.0;
								inst->zeta[ell][i][j][c] = (double)0.0;
							}

							// CASE 2: both x_ic = 0 and x_jc = 0
							if (((inst->x_ic[i][c] + inst->x_ic[j][c]) <= (double)0.01))
							{
								if (inst->X_ic[ell][i][c] <= inst->X_ic[ell][j][c]) {

									inst->zeta_i[ell][i][j][c] = (double)inst->alpha_data[i][j];
									inst->zeta_j[ell][i][j][c] = (double)0.0;
									inst->zeta[ell][i][j][c] = (double)0.0;
								}
								else
								{
									inst->zeta_i[ell][i][j][c] = (double)0.0;
									inst->zeta_j[ell][i][j][c] = (double)inst->alpha_data[i][j];
									inst->zeta[ell][i][j][c] = (double)0.0;
								}
							}
						}
						else
						{
							inst->zeta_i[ell][i][j][c] = (double)0.0;
							inst->zeta_j[ell][i][j][c] = (double)0.0;
							inst->zeta[ell][i][j][c] = (double)0.0;
						}
					}

					// CASE 2: both x_ic = 1 and x_jc = 1

					if (((inst->x_ic[i][c] + inst->x_ic[j][c]) >= (double)1.99))
					{

						//rr = rand() % 1000;
						//if (1 >= 0)
						if (inst->alpha_data[i][j] >= (double)0.0)
						{
							if (inst->X_ic[ell][i][c] <= inst->X_ic[ell][j][c])
							{
								inst->zeta_i[ell][i][j][c] = (double)inst->alpha_data[i][j];
								inst->zeta_j[ell][i][j][c] = (double)0.0;
								inst->zeta[ell][i][j][c] = (double)0.0;
							}
							else {
								inst->zeta_j[ell][i][j][c] = (double)inst->alpha_data[i][j];
								inst->zeta_i[ell][i][j][c] = (double)0.0;
								inst->zeta[ell][i][j][c] = (double)0.0;
							}
						}
						else {

							inst->zeta_i[ell][i][j][c] = (double)0.0;
							inst->zeta_j[ell][i][j][c] = (double)0.0;
							inst->zeta[ell][i][j][c] = -inst->alpha_data[i][j];

						}
					}
				}
			}
		}
	}
}





//------------------------------------------------
// Auxiliary functions for the heuristics
//------------------------------------------------


// Function to find the index of the most similar point not already in a group
int* findMostSimilarPoint_not_assigned(instance* inst, int sol,  int group, int group_size, int group_from)
{
	//(inst, sol, group_size_whichmin, group_size[c], c)

	int* out = new int[2];
	int count_i;

	int maxSimilarity = -(int)pow(inst->n_individuals, 2);
	int mostSimilarPoint = 1;

	int average_sim;
	bool first_cond;
	first_cond = true;
	// Find the closest element to the group
	for (int i = 0; i < inst->n_individuals; ++i)
	{
		if (group_from < inst->n_groups)
		{
			first_cond = (inst->Heur_X_ic[sol][i][group_from] >= 1);
		}		

		if (first_cond & (inst->Heur_X_ic[sol][i][group] <=0))
		{
			//cout << "ZOZZO sol = " << sol <<", c = " << group_from  << ", i = "<< i << endl;

			average_sim = 0;

			for (int j = 0; j < inst->n_individuals; ++j)
			{
				average_sim = average_sim + inst->Heur_X_ic[sol][j][group] * (inst->alpha_data[j][i] + inst->alpha_data[i][j]);
			}

			average_sim = average_sim / group_size;

			if (average_sim > maxSimilarity)
			{
				count_i = 0;
				for (int ccc = 0; ccc < inst->n_groups; ccc++)
				{
					count_i = count_i + inst->Heur_X_ic[sol][i][ccc];
				}
				if (count_i <= 0) 
				{
					maxSimilarity = average_sim;
					mostSimilarPoint = i;
				}

			}
		}
	}


	out[0] = mostSimilarPoint;
	out[1] = maxSimilarity;

	return out;

}

// Function to find the index of the most similar point not already in a group
int* findMostSimilarPoint(instance* inst, int sol, int group, int group_size, int group_from)
{
	//(inst, sol, group_size_whichmin, group_size[c], c)

	int* out = new int[2];
	int count_i;

	int maxSimilarity = -(int)pow(inst->n_individuals, 2);
	int mostSimilarPoint = 1;

	int average_sim;
	bool first_cond;
	first_cond = true;
	// Find the closest element to the group
	for (int i = 0; i < inst->n_individuals; ++i)
	{
		if (group_from < inst->n_groups)
		{
			first_cond = (inst->Heur_X_ic[sol][i][group_from] >= 1);
		}

		if (first_cond & (inst->Heur_X_ic[sol][i][group] <= 0))
		{
			//cout << "ZOZZO sol = " << sol <<", c = " << group_from  << ", i = "<< i << endl;

			average_sim = 0;

			for (int j = 0; j < inst->n_individuals; ++j)
			{
				average_sim = average_sim + inst->Heur_X_ic[sol][j][group] * (inst->alpha_data[j][i] + inst->alpha_data[i][j]);
			}

			average_sim = average_sim / group_size;

			if (average_sim > maxSimilarity)
			{
				maxSimilarity = average_sim;
				mostSimilarPoint = i;
			}
		}
	}


	out[0] = mostSimilarPoint;
	out[1] = maxSimilarity;

	return out;

}

// Function to calculate the correlation between two vectors
double calculateCorrelation(const std::vector<int>& v1, const std::vector<int>& v2) {
	
	int n = v1.size();

	// Calculate the mean of each vector
	double mean1 = 0.0, mean2 = 0.0;
	for (int i = 0; i < n; ++i) {
		mean1 += v1[i];
		mean2 += v2[i];
	}
	mean1 /= n;
	mean2 /= n;

	// Calculate the correlation
	double numerator = 0.0, denominator1 = 0.0, denominator2 = 0.0;
	for (int i = 0; i < n; ++i) {
		numerator += (v1[i] - mean1) * (v2[i] - mean2);
		denominator1 += std::pow(v1[i] - mean1, 2);
		denominator2 += std::pow(v2[i] - mean2, 2);
	}

	double correlation = numerator / std::sqrt(denominator1 * denominator2);
	return correlation;
}

// Function to extract a specific row from the 2D vector
std::vector<int> getRow(const int* const* matrix, int row, int size) {
	std::vector<int> result;
	for (int i = 0; i < size; ++i) {
		result.push_back(matrix[row][i]);
	}
	return result;
}

// Function to extract a specific row from the 2D vector
std::vector<int> getCol(const int* const* matrix, int col, int size) {
	std::vector<int> result;
	for (int i = 0; i < size; ++i) {
		result.push_back(matrix[i][col]);
	}
	return result;
}

// Compute the inner products between rows of inst->alpha_data[][]
int computeInnerProduct_row(instance* inst, int i, int j) {
	// Ensure the indices are within bounds
	if (i < 0 || i >= inst->n_individuals || j < 0 || j >= inst->n_individuals) {
		throw std::out_of_range("Row indices out of bounds");
	}

	int innerProduct = 0;

	// Compute the inner product of row i and row j
	for (int k = 0; k < inst->n_individuals; ++k) {
		innerProduct += inst->alpha_data[i][k] * inst->alpha_data[j][k];
	}

	return innerProduct;
}

// Compute the inner products between columns of inst->alpha_data[][]
int computeInnerProduct_col(instance* inst, int i, int j) {
	// Ensure the indices are within bounds
	if (i < 0 || i >= inst->n_individuals || j < 0 || j >= inst->n_individuals) {
		throw std::out_of_range("Row indices out of bounds");
	}

	int innerProduct = 0;

	// Compute the inner product of row i and row j
	for (int k = 0; k < inst->n_individuals; ++k) {
		innerProduct += inst->alpha_data[k][i] * inst->alpha_data[k][j];
	}

	return innerProduct;
}

void init_group_feature(instance* inst, int sol) {

	for (int f = 0; f < inst->n_features; f++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			inst->group_features[f][c] = 0;

			for (int i = 0; i < inst->n_individuals; i++)
			{
				inst->group_features[f][c] = inst->group_features[f][c] + inst->Heur_X_ic[sol][i][c] * inst->beta_data[i][f];
			}
		}
	}
}

void group_feature(instance* inst, int c, int sol) {

	for (int f = 0; f < inst->n_features; f++)
	{
		inst->group_features[f][c] = 0;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			if (inst->Heur_X_ic[sol][i][c] == 1) 
			{
				inst->group_features[f][c] = inst->group_features[f][c] + inst->beta_data[i][f];
			}
		}

	}
}

void objective(instance* inst, int method, int sol)
{
	if (method == 0) {//Heuristic

		inst->sol_val_HEUR_DOWN[sol] = 0;

		if (inst->problem <= 0) {

			for (int i = 0; i < (inst->n_individuals-1); i++)
			{
				for (int j =  (i + 1); j < inst->n_individuals; j++)
				{
					for (int c = 0; c < inst->n_groups; c++)
					{
						if ((inst->Heur_X_ic[sol][i][c] == 1) & (inst->Heur_X_ic[sol][j][c] == 1))
						{
							inst->sol_val_HEUR_DOWN[sol] = inst->sol_val_HEUR_DOWN[sol] + inst->alpha_data[i][j] + inst->alpha_data[j][i];
						}
					}
				}
			}

		}

		if (inst->problem >= 1)
		{

			int sol_val0 = 0;
			int sol_val1 = 0;

			for (int i = 0; i < inst->n_individuals; i++)
			{
				for (int j = 0; j < inst->n_individuals; j++)
				{
					sol_val0 = sol_val0 + inst->alpha_data[i][j] * inst->Heur_X_ic[sol][i][0] * inst->Heur_X_ic[sol][j][0];
				}
			}

			for (int c = 1; c < inst->n_groups; c++)
			{
				sol_val1 = 0;

				for (int i = 0; i < inst->n_individuals; i++)
				{
					for (int j = 0; j < inst->n_individuals; j++)
					{
						sol_val1 = sol_val1 + inst->alpha_data[i][j] * inst->Heur_X_ic[sol][i][c] * inst->Heur_X_ic[sol][j][c];
					}
				}

				if (sol_val1 < sol_val0) {
					sol_val0 = sol_val1;
				}
			}

			inst->sol_val_HEUR_DOWN[sol] = sol_val0;

		}
	}
	if (method == 1) {//Exact

		inst->sol_val = 0;

		if (inst->problem <= 0) {

			for (int i = 0; i < inst->n_individuals; i++)
			{
				for (int j = 0; j < inst->n_individuals; j++)
				{
					for (int c = 0; c < inst->n_groups; c++)
					{
						inst->sol_val = inst->sol_val + inst->alpha_data[i][j] * inst->X_ic[0][i][c] * inst->X_ic[0][j][c];
					}
				}
			}

		}

		if (inst->problem >= 1)
		{

			int sol_val0 = 0;
			int sol_val1 = 0;

			for (int i = 0; i < inst->n_individuals; i++)
			{
				for (int j = 0; j < inst->n_individuals; j++)
				{
					sol_val0 = sol_val0 + inst->alpha_data[i][j] * inst->X_ic[sol][i][0] * inst->X_ic[sol][j][0];
				}
			}

			for (int c = 1; c < inst->n_groups; c++)
			{
				sol_val1 = 0;

				for (int i = 0; i < inst->n_individuals; i++)
				{
					for (int j = 0; j < inst->n_individuals; j++)
					{
						sol_val1 = sol_val1 + inst->alpha_data[i][j] * inst->X_ic[sol][i][c] * inst->X_ic[sol][j][c];
					}
				}

				if (sol_val1 < sol_val0) {
					sol_val0 = sol_val1;
				}
			}

			inst->sol_val = sol_val0;

		}
	}
}

int Check_all_theta_constr(instance* inst)
{
	int check = 0;

	for (int f = 0; f < inst->n_features; f++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			for (int cc = 0; cc < inst->n_groups; cc++)
			{
				if ((abs(inst->group_features[f][c] - inst->group_features[f][cc]) > inst->theta_data[f])) {
					check = check + 1;
				}
			}
		}
	}

	return check;

}

int Check_unique_assignment_constr(instance* inst, int sol)
{
	int check = 0;

	for (int i = 0; i < inst->n_individuals; i++)
	{
		int count_i = 0;
		for (int ccc = 0; ccc < inst->n_groups; ccc++)
		{
			count_i = count_i + inst->Heur_X_ic[sol][i][ccc];
		}
		if (count_i != 1)
		{
			check = 1;
		}
	}

	return check;

}

int Check_theta_constr(instance* inst, int c_exept, int cc_exept, int f_exept)
{
	int check = 0;

	for (int f = 0; f < inst->n_features; f++)
	{
		for (int c = 0; c < inst->n_groups; c++)
		{
			for (int cc = 0; cc < inst->n_groups; cc++)
			{
				if ((c != c_exept) & (cc != cc_exept) & (f != f_exept))
				{
					check = check + (abs(inst->group_features[f][c] - inst->group_features[f][cc]) - inst->theta_data[f]);

				}
			}
		}
	}

	return check;

}

int find_max_2(int** ARRAY_BETA, int** CONTROL_X, int SIZE1, int SIZE2, int DIM, int fixDIM, int group, int exept)
{

	int size;

	if (DIM == 1)
	{
		size = SIZE1;
	}
	else
	{
		size = SIZE2;
	}

	int index = 0;
	int test;
	int beta_test = -1000;

	for (int i = 1; i < size; i++)
	{
		if ((CONTROL_X[i][group] > 0) & (i != exept))
		{
			test = ARRAY_BETA[i][fixDIM];

			if (test > beta_test)
			{
				beta_test = test;
				index = i;
			}
		}
	}

	return index;
}

int find_min_2(int** ARRAY_BETA, int** CONTROL_X, int SIZE1, int SIZE2, int DIM, int fixDIM, int group, int exept)
{

	int size;

	if (DIM == 1)
	{
		size = SIZE1;

	}
	else {
		size = SIZE2;
	}

	int index = 0;
	int test;
	int beta_test = size;

	for (int i = 1; i < size; i++)
	{
		if ((CONTROL_X[i][group] > 0) & (i != exept))
		{
			test = ARRAY_BETA[i][fixDIM];

			if (test < beta_test)
			{
				beta_test = test;
				index = i;
			}
		}
	}

	return index;

}

// Function to generate a random index
int generateRandomIndex(instance* inst, int sol, int c) {

	// Vector to store valid indices
	std::vector<int> validIndices;

	// Populate the valid indices vector
	for (int i = 0; i < inst->n_individuals; ++i)
	{
		if (inst->Heur_X_ic == nullptr || inst->Heur_X_ic[sol] == nullptr || inst->Heur_X_ic[sol][i] == nullptr)
		{
			throw std::runtime_error("Error in generateRandomIndex(): Heur_X_ic is not properly initialized.");
		}

		if (c < inst->n_groups) 
		{
			if (inst->Heur_X_ic[sol][i][c] > 0) 
			{
				//cout << "index added: " << i << ", validIndices size: " << validIndices.size() << endl;
				validIndices.push_back(i);
				//cout << "Index added: " << i << ", validIndices size: " << validIndices.size() << endl;
			}
		}
		else
		{
			validIndices.push_back(i);
		}
	}

	//std::random_device rd;
	//std::mt19937 gen(rd());
	//std::mt19937 gen(21051983); // Fixed seed
	std::uniform_int_distribution<> dis(0, validIndices.size() - 1);

	int randomIdx = dis(gen);

	// Return the selected index
	return validIndices[randomIdx];
}

// Function to compare two int** arrays
bool areSolutionsEqual(instance* inst, int* solution1, int* solution2)
{	
	for (int j = 0; j < inst->n_individuals; ++j)
	{
		if (solution1[j] != solution2[j]) {
			return false;
		}
	}
	
	return true;
}

// Function to create a deep copy of a 2D array
int* deepCopyArray(instance* inst, int* source) {
	int* copy = new int [inst->n_individuals];
	for (int i = 0; i < inst->n_individuals; ++i) {
		copy[i]= source[i];
	}
	return copy;
}

int makeCopy(int x) { return x; }

void printExploredSolution(instance* inst, const std::vector<int*>& solutions)
{
	std::cout << "ExploredSolution contents:\n";

	for (size_t s = 0; s < solutions.size(); ++s) 
	{ // Iterate over each int** in the vector
		std::cout << "Solution " << s + 1 << ":\n";
		if (solutions[s] != nullptr) 
		{ // Check if the solution exists
			for (int i = 0; i < inst->n_individuals; ++i) { // Iterate over rows
					std::cout << solutions[s][i] << " ";
			}
			std::cout << "\n"; // End of row
		}
		else {
			std::cout << "(null)\n";
		}
		std::cout << "----\n"; // Separator between solutions
	}
	if (solutions.empty()) {
		std::cout << "(The vector is empty)\n";
	}
}

// Function to check if a new solution exists in the vector ExploredSolution
bool isSolutionInVector(instance* inst, const std::vector<int*>& ExploredSolution, int* newSolution) {
	
	for (size_t s = 0; s < ExploredSolution.size(); ++s)
	{
		bool check = true; // Assume a match until proven otherwise
		if (ExploredSolution[s] != nullptr)
		{
			for (int j = 0; j < inst->n_individuals; ++j)
			{
				//std::cout << ExploredSolution[s][j] << " ";
				//std::cout << "Check " << j << " component " << ExploredSolution[s][j] << " " << newSolution[j] << std::endl;
				if (ExploredSolution[s][j] != newSolution[j])
				{
					check = false;
					break; // No need to check further rows for this solution
				}
			}
		}
		else {
			cout << "Nullptr in function isSolutionInVector()" << endl;
			exit(-1);
		}
		//std::cout << "Check = " << check << std::endl;
		if (check) {
			return true; // Return immediately if a match is found
		}
	}

	return false; // No match found in the vector
}

void clearExploredSolution(instance* inst, std::vector<int**>& solutions) {
	
	cout << "Minkione: "<< solutions.size() << endl;

	// Iterate over each int** in the vector
	for (int ss = 0; ss < solutions.size(); ss ++)
	{

		cout << "Stronzone: " << ss << endl;
		if (solutions[ss] != nullptr) { // Check if the outer pointer is valid
			// Deallocate the memory for the int** (row by row)
			for (int i = 0; i < inst->n_individuals; ++i) {
				if (solutions[ss][i] != nullptr) { // Check if the row pointer is valid
					delete[] solutions[ss][i];    // Free the row
					//solutions[ss][i] = nullptr;  // Set to null for safety
				}
			}
			delete[] solutions[ss]; // Free the outer pointer
			//solutions[ss] = nullptr; // Set the pointer to null for safety
		}
	}

	cout << "Minkione again " << endl;

	solutions.clear(); // Clear the vector

	cout << "Minkione again " << endl;

}






//------------------------------------------------
// Heuristics
//------------------------------------------------


void R_MN_SAP_ConstructuveHeur(instance* inst)
{
	// Create a vector of integers
	std::vector<int> ExploredElements;

	int most_sim_first;

	int gg_c;
	int gg_cc;
	int gg_diff;
	int beta_c;
	int beta_cc;
	int ind_i;
	int ind_ii;
	int ind_c;
	int ind_cc;
	int infeasible;
	int infeasible_new;
	bool breakk;
	int increment;
	int check0;
	int check1;
	int obj0;
	int obj1;
	bool test;

	int ii;
	int jj;
	int val_ij;
	double cor_ij;
	double cor_ji;
	int group_size_c;
	int Group_beta_big;
	int Group_beta_small;
	int* Result_c = new int[2];

	int* group_size = new int[inst->n_groups];
	int group_size_min;
	int group_size_whichmin;
	int group_size_max;
	int group_size_whichmax;
	int MyTest;
	int Theta_iter_limit = 2*inst->n_groups * inst->n_individuals;

	//std::random_device rd; // Seed for randomness
	//std::mt19937 gen(rd()); // Mersenne Twister engine

	cout << endl << "CONSTRUCTING " << inst->Heu_sol  << " HEURISTIC SOLUTIONS " << endl << endl;

	#pragma omp parallel for schedule(guided)
	for (int sol = 0; sol < inst->Heu_sol; sol++)
	{
		std::uniform_real_distribution<> dis(0.0, 1.0); // Range [0, 1)
		double random_number;

		cout <<"............................................" << endl;
		//------------------------------------------------
		// Initialization
		//------------------------------------------------

		int first_element = 0;
		for (int c = 0; c < inst->n_groups; c++)
		{
			if (ExploredElements.empty()) // Attach the first two individuals to group 0
			{
				most_sim_first = -(int)pow(inst->n_individuals, 2);
				ii = first_element;
				jj = first_element + 1;

				for (int i = 0; i < inst->n_individuals - 1; ++i)
				{
					for (int j = i + 1; j < inst->n_individuals; ++j)
					{
						//cor_ij = calculateCorrelation(getRow(inst->alpha_data, i, inst->n_individuals), getRow(inst->alpha_data, j, inst->n_individuals));
						//cor_ji = calculateCorrelation(getCol(inst->alpha_data, i, inst->n_individuals), getCol(inst->alpha_data, j, inst->n_individuals));

						cor_ij = computeInnerProduct_row(inst, i, j);
						cor_ji = computeInnerProduct_col(inst, i, j);

						//val_ij = (inst->alpha_data[i][j] + inst->alpha_data[j][i]);
						val_ij = (inst->alpha_data[i][j] + inst->alpha_data[j][i] + cor_ij + cor_ji);

						if (val_ij > most_sim_first) {

							// Generate a random number
							//random_number = dis(gen);

							//if (random_number > 0) 
							//{
								ii = i;
								jj = j;

								most_sim_first = val_ij;
							//}
						}
					}
				}

				first_element = first_element + 2;

				inst->Heur_X_ic[sol][ii][0] = 1;
				inst->Heur_X_ic[sol][jj][0] = 1;

				ExploredElements.push_back(ii);
				ExploredElements.push_back(jj);

			}
			else {	// Initialize all the other groups c > 0

				most_sim_first = -(int)pow(inst->n_individuals, 2);

				ii = first_element;
				jj = first_element + 1;

				int isInList_i;
				int isInList_j;

				for (int i = 0; i < inst->n_individuals - 1; ++i)
				{
					isInList_i = 0;
					for (int iii : ExploredElements)
					{
						if (iii == i) {
							isInList_i = 1;
							break;
						}
					}
					if (isInList_i <= 0)
					{
						for (int j = i + 1; j < inst->n_individuals; ++j)
						{
							val_ij = (inst->alpha_data[i][j] + inst->alpha_data[j][i]);

							if (val_ij > most_sim_first)
							{
								isInList_j = 0;

								for (int jjj : ExploredElements) {
									if (jjj == j) {
										isInList_j = 1;
										break;
									}
								}

								if (isInList_j <= 0)
								{
									ii = i;
									jj = j;

									most_sim_first = val_ij;
								}
							}
						}
					}
				}

				first_element = first_element + 2;

				inst->Heur_X_ic[sol][ii][c] = 1;
				inst->Heur_X_ic[sol][jj][c] = 1;

				ExploredElements.push_back(ii);
				ExploredElements.push_back(jj);

			}
		}

		//-----------------------------------------------------------------------
		// Main loop
		//-----------------------------------------------------------------------
		// All groups have been initialized (they contain two individuals each)
		// The following loop associate the remaining n - 2*c individuals to groups
		//-----------------------------------------------------------------------

		int element = 0;
		int group_select = 0;
		int value = 0;

		// While there is still non-assigned individuals
		while (ExploredElements.size() < inst->n_individuals)
		{
			//for (int iiii : ExploredElements) {
			//	std::cout << iiii << " ";
			//}
			//std::cout << " " << endl;

			value = -(int)pow(inst->n_individuals, 4);

			for (int c = 0; c < inst->n_groups; c++)
			{
				Result_c = findMostSimilarPoint_not_assigned(inst, sol, c, ExploredElements.size(), inst->n_groups + 1);

				//cout << "Soy el mejor para el grupo " << c << ". Me llamo " << Result_c[0] << endl;
				//cout << "Estoy mal initializado :-( " << value << endl;
				//cout << "Culpa suya: " << Result_c[1] << endl;

				if (value < Result_c[1]) {

					element = Result_c[0];
					value = Result_c[1];
					group_select = c;

				}
			}

			//cout << "Explored size =" << ExploredElements.size() << ", element =" << element << ", group =" << group_select << endl;

			inst->Heur_X_ic[sol][element][group_select] = 1;
			ExploredElements.push_back(element);

		}

		//-----------------------------------------------------------------------
		// Compute and print objective function
		//-----------------------------------------------------------------------

		objective(inst, 0, sol);

		cout << "OBJECTIVE HEURISTIC BEFORE THE TAU SWAP = " << inst->sol_val_HEUR_DOWN[sol] << endl;

		//-----------------------------------------------------------------------
		// Tau-feasibility
		//-----------------------------------------------------------------------
		// The following loop move individuals from the largest groups to  
		// the smallest ones to guarantee that the group sizes are between 
		// tau_min and tau_max
		//-----------------------------------------------------------------------

		MyTest = -1;
		while (MyTest < 0)
		{
			group_size_min = inst->n_individuals;
			group_size_whichmin = 0;

			group_size_max = -1;
			group_size_whichmax = 0;
			
			for (int c = 0; c < inst->n_groups; c++)
			{
				group_size_c = 0;

				for (int i = 0; i < inst->n_individuals; ++i)
				{
					if (inst->Heur_X_ic[sol][i][c] >= 1)
					{
						group_size_c = group_size_c + 1;
					}
				}
				group_size[c] = group_size_c;

				if (group_size_c > group_size_max) {
					group_size_whichmax = c;
					group_size_max = group_size_c;
				}
				if (group_size_c < group_size_min) {
					group_size_whichmin = c;
					group_size_min = group_size_c;
				}
				//cout << " " << group_size_c ;
			}
			//cout << endl;

			for (int c = 0; c < inst->n_groups; c++)
			{
				//cout << "Size of groups: " << endl;
				//for (int ccc = 0; ccc < inst->n_groups; ccc++)
				//{
				//	cout << " " << group_size[ccc];
				//}
				//cout << endl;

				if (group_size[c] < inst->tau_min)
				{
					ExploredElements.clear();

					// Just compute the size of group c:
					for (int i = 0; i < inst->n_individuals; ++i)
					{
						if (inst->Heur_X_ic[sol][i][c] >= 1)
						{
							ExploredElements.push_back(i);
						}
						if (inst->Heur_X_ic[sol][i][group_size_whichmax] >= 1)
						{
							//cout << i << " ";
						}
					}
					//cout << endl;

					//-----------------------
					// Since the size of group c is smaller than tau_min, find the member in the 
					// largest group which is close (based on alpha) to group c
					//-----------------------

					Result_c = findMostSimilarPoint(inst, sol, c, group_size[c], group_size_whichmax);

					element = Result_c[0];
					//value = Result_c[1];
					
					//cout << "Selected element:  " << element << ", belonging to group =" << group_size_whichmax << endl;
					//cout << "Size of the biggest group = " << group_size[group_size_whichmax] << endl;

					if (group_size[group_size_whichmax] >= inst->tau_min)
					{
						//cout << "Move individual " << element << ", from group " << group_size_whichmax << " to group " << c << endl;

						inst->Heur_X_ic[sol][element][group_size_whichmax] = 0;
						inst->Heur_X_ic[sol][element][c] = 1;

						group_size[group_size_whichmax] = group_size[group_size_whichmax] - 1;
						group_size[c] = group_size[c] + 1;

						group_size_max = 0;
						group_size_min = inst->n_individuals;
						for (int ggg = 0; ggg < inst->n_groups; ggg++)
						{
							if (group_size[ggg] > group_size_max) {
								group_size_whichmax = ggg;
								group_size_max = group_size[ggg];
							}
							if (group_size[ggg] < group_size_min) {
								group_size_whichmin = ggg;
								group_size_min = group_size[ggg];
							}
						}
					}
					else
					{
						cout << "Warning: No feasible solution exist" << endl;
						exit(-1);
					}

				}

				if (group_size[c] > inst->tau_max)
				{

					ExploredElements.clear();

					// Just compute the size of group group_size_whichmax:
					for (int i = 0; i < inst->n_individuals; ++i)
					{
						if (inst->Heur_X_ic[sol][i][group_size_whichmin] >= 1)
						{
							ExploredElements.push_back(i);
						}
						//if (inst->Heur_X_ic[sol][i][c] >= 1)
						//{
						//	cout << i <<" " ;
						//}
						
					}
					//cout << endl;

					//-----------------------
					// Since the size of group c is larger than tau_max, find the member c which 
					// is close (based on alpha) to the smallest group (group_size_whichmin)
					//-----------------------

					Result_c = findMostSimilarPoint(inst, sol, group_size_whichmin, group_size[c], c);

					element = Result_c[0];
					//value = Result_c[1];

					//cout << "Selected element:  " << element << ", belonging to group =" << c << endl;
					//cout << "Size of the smallest group = " << group_size[group_size_whichmin] << endl;

					if (group_size[group_size_whichmin] <= inst->tau_max)
					{
						//cout << "Move individual " << element << ", from group " << c << " to group " << group_size[group_size_whichmin] << endl;

						inst->Heur_X_ic[sol][element][group_size_whichmin] = 1;
						inst->Heur_X_ic[sol][element][c] = 0;

						group_size[group_size_whichmin] = group_size[group_size_whichmin] + 1;
						group_size[c] = group_size[c] - 1;

						group_size_max = 0;
						group_size_min = inst->n_individuals;
						for (int ggg = 0; ggg < inst->n_groups; ggg++)
						{
							if (group_size[ggg] > group_size_max) {
								group_size_whichmax = ggg;
								group_size_max = group_size[ggg];
							}
							if (group_size[ggg] < group_size_min) {
								group_size_whichmin = ggg;
								group_size_min = group_size[ggg];
							}
						}
					}
					else
					{
						cout << "Warning: No feasible solution exist" << endl;
						exit(-1);
					}

				}
			}

			MyTest = 1;
			for (int ggg = 0; ggg < inst->n_groups; ggg++)
			{
				if (group_size[ggg] > inst->tau_max) {
					MyTest = -1;
					break;
				}
				if (group_size[ggg] < inst->tau_min) {
					MyTest = -1;
					break;
				}
			}
			if (MyTest > 0)
			{
				break;
			}
		}

		//-----------------------------------------------------------------------
		// Compute and print objective function
		//-----------------------------------------------------------------------

		objective(inst, 0, sol);

		cout << "OBJECTIVE HEURISTIC BEFORE THE THETA SWAP = " << inst->sol_val_HEUR_DOWN[sol] << endl;

		//-----------------------------------------------------------------------
		// Theta-feasibility
		//-----------------------------------------------------------------------
		// The following loop perform bilaterar exchanges between groups 
		//-----------------------------------------------------------------------

		init_group_feature(inst, sol);

		infeasible = Check_all_theta_constr(inst);
		increment = 0;

		while (infeasible > 0)
		{
			//cout << endl << "-----------------infeasible : " << infeasible << endl;

			for (int f = 0; f < inst->n_features; f++)
			{
				for (int c = 0; c < inst->n_groups-1; c++)
				{
					for (int cc = c+1; cc < inst->n_groups; cc++)
					{

						if (infeasible <= 0)
						{
							break;
						}

						if (increment > Theta_iter_limit)
						{
							infeasible = 1;
							break;
						}

						group_feature(inst, c, sol);
						group_feature(inst, cc, sol);

						gg_c = inst->group_features[f][c];
						gg_cc = inst->group_features[f][cc];

						gg_diff = abs(gg_c - gg_cc);
						
						//cout << "f = "<< f<< ", c = " << c << ", cc = " << cc << "  ----->  " << gg_diff - inst->theta_data[f] << endl;
						//cout << "size[" << c << "]= " << group_size[c] << endl;
						//cout << "size[" << cc << "]= " << group_size[cc] << endl;	
						//cout << endl << endl;

						if (inst->theta_data[f] < gg_diff)
						{
							if (gg_c > gg_cc) {
								Group_beta_big = c;
								Group_beta_small = cc;
							}
							else {
								Group_beta_big = cc;
								Group_beta_small = c;
							}

							//cout << "Theta[" << f << "] = " << inst->theta_data[f] << endl << endl;
							//cout << "Group " << Group_beta_big << endl << endl;
							//for (int y = 0; y < inst->n_individuals; y++)
							//{
							//	if (inst->Heur_X_ic[sol][y][Group_beta_big] == 1) {
							//		cout << "Beta[" << y << "] = " << inst->beta_data[y][f] << endl;
							//	}
							//}

							ind_i = generateRandomIndex(inst, sol, Group_beta_big);

							//ind_i = find_max_2(inst->beta_data, inst->Heur_X_ic[sol], inst->n_individuals, inst->n_features, 1, f, Group_beta_big, inst->n_individuals);

							//cout << "The total level of feature " << f << " in group " << Group_beta_big << " is " << inst->group_features[f][Group_beta_big] << endl;
							//cout << "Drop the following member :" << ind_i << endl;

							//cout << endl << "Group " << Group_beta_small << endl << endl;
							//for (int y = 0; y < inst->n_individuals; y++)
							//{
							//	if (inst->Heur_X_ic[sol][y][Group_beta_small] == 1) {
							//		cout << "Beta[" << y << "] = " << inst->beta_data[y][f] << endl;
							//	}
							//}

							ind_ii = generateRandomIndex(inst, sol, Group_beta_small);
							//ind_ii = find_min_2(inst->beta_data, inst->Heur_X_ic[sol], inst->n_individuals, inst->n_features, 1, f, Group_beta_small, ind_i);

							//cout << "The total level of feature " << f << " in group " << Group_beta_small << " is " << inst->group_features[f][Group_beta_small] << endl;
							//cout << "Drop the following member :" << ind_ii << endl;

							inst->Heur_X_ic[sol][ind_i][Group_beta_big] = 0;
							inst->Heur_X_ic[sol][ind_i][Group_beta_small] = 1;
							inst->Heur_X_ic[sol][ind_ii][Group_beta_big] = 1;
							inst->Heur_X_ic[sol][ind_ii][Group_beta_small] = 0;

							group_feature(inst, c, sol);
							group_feature(inst, cc, sol);

							infeasible_new = Check_all_theta_constr(inst);

							if (infeasible_new <= infeasible) 
							{
								infeasible = infeasible_new;
							}
							else {

								inst->Heur_X_ic[sol][ind_i][Group_beta_big] = 1;
								inst->Heur_X_ic[sol][ind_i][Group_beta_small] = 0;
								inst->Heur_X_ic[sol][ind_ii][Group_beta_big] = 0;
								inst->Heur_X_ic[sol][ind_ii][Group_beta_small] = 1;

								group_feature(inst, c, sol);
								group_feature(inst, cc, sol);

							}
						}
						else {
						
							infeasible = Check_all_theta_constr(inst);

						}

						//cout << endl << "-----------------infeasible : " << infeasible << endl;

						increment++;

					}
				}
			}


			if (infeasible <= 0)
			{
				break;
			}

			if (increment > Theta_iter_limit)
			{
				break;
			}

		}

		objective(inst, 0, sol);
		//cout << "OBJECTIVE HEURISTIC AFTER THE THETA SWAP = " << inst->sol_val_HEUR_DOWN[sol] << endl;

		for (int c = 0; c < inst->n_groups; c++)
		{
			group_feature(inst, c, sol);
		}

		for (int f = 0; f < inst->n_features; f++)
		{
			for (int c = 0; c < inst->n_groups; c++)
			{
				for (int cc = 0; cc < inst->n_groups; cc++)
				{
					gg_c = inst->group_features[f][c];
					gg_cc = inst->group_features[f][cc];

					if (gg_cc - gg_c > inst->theta_data[f]) {

						cout << "Warning: No feasible solution exist" << endl;
						exit(-1);

					}
				}
			}
		}

		objective(inst, 0, sol);

		cout << "SOLUTION " << sol << " CREATED. OBJECTIVE FUNCTION =" << inst->sol_val_HEUR_DOWN[sol] << endl << endl;

		// Clear the vector
		ExploredElements.clear(); 
	}

}


void R_MN_SAP_LocalSearch(instance* inst)
{
	// Initialize adjacency matrix and list of neighbors
	inst->Tabu_X_ic = new int* [inst->n_individuals];

	for (int i = 0; i < inst->n_individuals; i++)
	{
		inst->Tabu_X_ic[i] = new int[inst->n_groups];

		for (int c = 0; c < inst->n_groups; ++c)
		{
			inst->Tabu_X_ic[i][c] = 0;
		}
	}

	int most_sim_first;
	int new_objective;
	int best_objective;
	int Best_objective;
	int best_ii = 0;
	int best_jj = 0;

	int gg_c;
	int gg_cc;
	int gg_diff;
	int infeasible;
	int infeasible_new;

	int improving;

	int group_size_c;
	int Group_beta_big;
	int Group_beta_small;
	int* Result_c = new int[2];

	int* group_size = new int[inst->n_groups];
	int* group_of_i = new int[inst->n_individuals];
	bool check;

	int Check_theta_feasibility;
	int Theta_iter_limit = 2 * inst->n_groups * inst->n_individuals;

	//std::random_device rd; // Seed for randomness
	//std::mt19937 gen(rd()); // Mersenne Twister engine
	//std::uniform_real_distribution<> dis(0.0, 1.0); // Range [0, 1)
	//double random_number;

	cout << endl << "NUMBER OF SOLUTIONS: " << inst->Heu_sol << endl << endl;

	#pragma omp parallel for schedule(guided)
	for (int sol = 0; sol < inst->Heu_sol; sol++)
	{
		objective(inst, 0, sol);

		cout << "SOLUTION " << sol <<", STARTING VALUE " << inst->sol_val_HEUR_DOWN[sol]  << "............................................" << endl;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int c = 0; c < inst->n_groups; c++)
			{
				inst->Tabu_X_ic[i][c] = inst->Heur_X_ic[sol][i][c];
				//cout << inst->Heur_X_ic[sol][i][c] << " "; 

				if (inst->Heur_X_ic[sol][i][c] >= 1)
				{
					group_of_i[i] = c;
				}
				//cout << "The group of " << i << " is " << group_of_i[i] << endl;
			}
			//cout << endl;
		}

		init_group_feature(inst, sol);
		best_objective = inst->sol_val_HEUR_DOWN[sol];
		Best_objective = inst->sol_val_HEUR_DOWN[sol];

		improving = 1;

		while (improving > 0)
		{
			//infeasible = Check_all_theta_constr(inst);
			//cout << endl << "-----------------infeasible : " << infeasible << endl;

			check = false;

			int* copy_group_of_i = deepCopyArray(inst, group_of_i);

			for (int i = 0; i < (inst->n_individuals - 1); i++)
			{
				if (check) {
					break;
				}

				for (int j = i + 1; j < inst->n_individuals; j++)
				{
					if (copy_group_of_i[i] == copy_group_of_i[j]) { continue; }
					// Switch group

					inst->Heur_X_ic[sol][i][copy_group_of_i[i]] = 0;
					inst->Heur_X_ic[sol][i][copy_group_of_i[j]] = 1;
					inst->Heur_X_ic[sol][j][copy_group_of_i[i]] = 1;
					inst->Heur_X_ic[sol][j][copy_group_of_i[j]] = 0;

					objective(inst, 0, sol);

					new_objective = inst->sol_val_HEUR_DOWN[sol];

					inst->Heur_X_ic[sol][i][copy_group_of_i[i]] = 1;
					inst->Heur_X_ic[sol][i][copy_group_of_i[j]] = 0;
					inst->Heur_X_ic[sol][j][copy_group_of_i[i]] = 0;
					inst->Heur_X_ic[sol][j][copy_group_of_i[j]] = 1;

					//objective(inst, 0, sol);
					//cout << "RI-STAMPALA QUI: " << inst->sol_val_HEUR_DOWN[sol] << endl;
					//printExploredSolution(inst, ExploredSolution);
					//cout << "STRONZONE " << group_of_i[i] << " " << group_of_i[j] << endl;

					//group_of_i[i] = copy_group_of_i[j];
					//group_of_i[j] = copy_group_of_i[i];

					//cout << "individual " << i << "(" << group_of_i[i] <<") , and " << j << "(" << group_of_i[j] << ") " << "Best objective = " << best_objective << ", New objective = " << new_objective << endl;

					if (best_objective < new_objective)
					{
						group_of_i[i] = copy_group_of_i[j];
						group_of_i[j] = copy_group_of_i[i];

						group_feature(inst, copy_group_of_i[i], sol);
						group_feature(inst, copy_group_of_i[j], sol);

						//out << "individual " << i << "(" << group_of_i[i] << ") , and " << j << "(" << group_of_i[j] << ") " << "Best objective = " << best_objective << ", New objective = " << new_objective << endl;

						if (Check_all_theta_constr(inst) <= 0)
						{
							best_objective = new_objective;

							best_ii = i;
							best_jj = j;

							inst->Heur_X_ic[sol][i][copy_group_of_i[best_ii]] = 0;
							inst->Heur_X_ic[sol][i][copy_group_of_i[best_jj]] = 1;
							inst->Heur_X_ic[sol][j][copy_group_of_i[best_ii]] = 1;
							inst->Heur_X_ic[sol][j][copy_group_of_i[best_jj]] = 0;

							objective(inst, 0, sol);

							//cout << "[WHAT ABOUT THIS] " << inst->sol_val_HEUR_DOWN[sol] << endl;

							check = true;
							break;
						}
					}
				}
			}

			if (Best_objective < best_objective) 
			{
				Best_objective = best_objective;

				inst->Heur_X_ic[sol][best_ii][copy_group_of_i[best_ii]] = 0;
				inst->Heur_X_ic[sol][best_ii][copy_group_of_i[best_jj]] = 1;
				inst->Heur_X_ic[sol][best_jj][copy_group_of_i[best_jj]] = 0;
				inst->Heur_X_ic[sol][best_jj][copy_group_of_i[best_ii]] = 1;

				inst->Tabu_X_ic[best_ii][copy_group_of_i[best_ii]] = 0;
				inst->Tabu_X_ic[best_ii][copy_group_of_i[best_jj]] = 1;
				inst->Tabu_X_ic[best_jj][copy_group_of_i[best_jj]] = 0;
				inst->Tabu_X_ic[best_jj][copy_group_of_i[best_ii]] = 1;

				for (int ii = 0; ii < inst->n_individuals; ii++)
				{
					for (int cc = 0; cc < inst->n_groups; cc++)
					{
						if (inst->Heur_X_ic[sol][ii][cc] >= 1)
						{
							group_of_i[ii] = cc;
						}
						//cout << "The group of " << i << " is " << group_of_i[i] << endl;
					}
					//cout << endl;
				}

				objective(inst, 0, sol);

			}
			else
			{
				for (int i = 0; i < inst->n_individuals; i++)
				{
					for (int c = 0; c < inst->n_groups; c++)
					{
						inst->Heur_X_ic[sol][i][c] = inst->Tabu_X_ic[i][c];
						//cout << inst->Heur_X_ic[sol][i][c] << " "; 

						if (inst->Heur_X_ic[sol][i][c] >= 1)
						{
							group_of_i[i] = c;
						}
						//cout << "The group of " << i << " is " << group_of_i[i] << endl;
					}
					//cout << endl;
				}

				objective(inst, 0, sol);

				improving = -1;
			}			

			cout << "INCUMBENT SOLUTION IN THE LOCAL SEARCH = " << inst->sol_val_HEUR_DOWN[sol] << endl;

		}

		cout << "FINAL OBJECTIVE HEURISTIC IN THE LOCAL SEARCH = " << inst->sol_val_HEUR_DOWN[sol] << endl;

		//for (int c = 0; c < inst->n_groups; c++)
		//{
		//	group_feature(inst, c, sol);
		//}

		//for (int i = 0; i < inst->n_individuals - 1; i++)
		//{
		//	for (int c = 0; c < inst->n_groups; c++)
		//	{
		//		cout << inst->Heur_X_ic[sol][i][c] << " ";
		//	}

		//	cout << endl;
		//}
		//cout << "---" << endl;
	}

}


void R_MN_SAP_TabuSearch(instance* inst)
{
	// Initialize adjacency matrix and list of neighbors
	inst->Tabu_X_ic = new int* [inst->n_individuals];

	for (int i = 0; i < inst->n_individuals; i++)
	{
		inst->Tabu_X_ic[i] = new int[inst->n_groups];

		for (int c = 0; c < inst->n_groups; ++c)
		{
			inst->Tabu_X_ic[i][c] = 0;
		}
	}
	
	int new_objective;
	int best_objective;
	int Best_objective;

	int gg_c;
	int gg_cc;
	int gg_diff;
	int infeasible;
	int infeasible_new;

	int not_improving;

	int group_size_c;

	int* group_size = new int[inst->n_groups];
	int* group_of_i = new int[inst->n_individuals];

	int Check_theta_feasibility;
	int Theta_iter_limit = 2 * inst->n_groups * inst->n_individuals;

	int best_ii = 0;
	int best_jj = 0;

	cout << endl << "IMPROVING " << inst->Heu_sol << " SOLUTIONS BY TS" << endl << endl;

	#pragma omp parallel for schedule(guided)
	for (int sol = 0; sol < inst->Heu_sol; sol++)
	{
		objective(inst, 0, sol);

		cout << "SOLUTION " << sol << ", STARTING VALUE " << inst->sol_val_HEUR_DOWN[sol] << "............................................" << endl;

		// Create a vector of arrays
		std::vector<int*> ExploredSolution;
		
		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int c = 0; c < inst->n_groups; c++)
			{
				//cout << inst->Heur_X_ic[sol][i][c] << " "; 
				inst->Tabu_X_ic[i][c] = inst->Heur_X_ic[sol][i][c];

				if (inst->Heur_X_ic[sol][i][c] >= 1)
				{
					group_of_i[i] = c;
				}
				//cout << "The group of " << i << " is " << group_of_i[i] << endl;
			}
			//cout << endl;
		}

		// Add a copy of Tabu_X_ic to ExploredSolution
		if (ExploredSolution.empty()) {
			ExploredSolution.push_back(deepCopyArray(inst,group_of_i));
		}
		else {
			ExploredSolution.push_back(deepCopyArray(inst, group_of_i));
		}

		init_group_feature(inst, sol);
		best_objective = inst->sol_val_HEUR_DOWN[sol];
		Best_objective = inst->sol_val_HEUR_DOWN[sol];

		not_improving = 0;

		while(not_improving < inst->Tabu_not_improving_max)
		{
			infeasible = Check_all_theta_constr(inst);
			int infea = Check_unique_assignment_constr(inst, sol);
			cout << endl << "-----------------infeasible : " << infeasible << " " << infea << endl;

			int* copy_group_of_i = deepCopyArray(inst, group_of_i);

			for (int i = 0; i < (inst->n_individuals - 1); i++)
			{
				for (int j = i + 1; j < inst->n_individuals; j++)
				{
					if (group_of_i[i] == group_of_i[j]) { continue; }

					// Switch group

					inst->Heur_X_ic[sol][i][copy_group_of_i[i]] = 0;
					inst->Heur_X_ic[sol][i][copy_group_of_i[j]] = 1;
					inst->Heur_X_ic[sol][j][copy_group_of_i[i]] = 1;
					inst->Heur_X_ic[sol][j][copy_group_of_i[j]] = 0;

					objective(inst, 0, sol);

					new_objective = inst->sol_val_HEUR_DOWN[sol];

					//printExploredSolution(inst, ExploredSolution);
					//cout << "STRONZONE " << group_of_i[i] << " " << group_of_i[j] << endl;

					group_of_i[i] = copy_group_of_i[j];
					group_of_i[j] = copy_group_of_i[i];
					
					//printExploredSolution(inst, ExploredSolution);
					//cout << "STRONZONE " << group_of_i[i] << " " << group_of_i[j] << endl;

					//cout << "individual " << i << "(" << group_of_i[i] <<") , and " << j << "(" << group_of_i[j] << ") " << "Best objective = " << best_objective << ", New objective = " << new_objective << endl;

					//cout << "Is this solution present? " << isSolutionInVector(inst, ExploredSolution, group_of_i) << endl;

					// Check if the newSolution is in the vector
					if (isSolutionInVector(inst, ExploredSolution, group_of_i))
					{
						std::cout << "Solution is already in the tabu list.\n";

						inst->Heur_X_ic[sol][i][copy_group_of_i[i]] = 1;
						inst->Heur_X_ic[sol][i][copy_group_of_i[j]] = 0;
						inst->Heur_X_ic[sol][j][copy_group_of_i[i]] = 0;
						inst->Heur_X_ic[sol][j][copy_group_of_i[j]] = 1;

						//group_feature(inst, group_of_i[i], sol);
						//group_feature(inst, group_of_i[j], sol);

						//objective(inst, 0, sol);
						continue;
					}

					if (best_objective < new_objective)
					{
						group_feature(inst, group_of_i[i], sol);
						group_feature(inst, group_of_i[j], sol);

						//out << "individual " << i << "(" << group_of_i[i] << ") , and " << j << "(" << group_of_i[j] << ") " << "Best objective = " << best_objective << ", New objective = " << new_objective << endl;

						if (Check_all_theta_constr(inst) <= 0)
						{
							best_objective = new_objective;

							best_ii = i;
							best_jj = j;
						}
						//Check_theta_feasibility = true;

						//for (int f = 0; f < inst->n_features; f++)
						//{
						//	gg_c = inst->group_features[f][group_of_i[i]];
						//	gg_cc = inst->group_features[f][group_of_i[j]];

						//	gg_diff = abs(gg_c - gg_cc);							

						//	if (inst->theta_data[f] < gg_diff)
						//	{
						//		Check_theta_feasibility = false;

						//		break;
						//	}
						//}

						//if (Check_theta_feasibility)
						//{
						//	best_objective = new_objective;

						//	best_ii = i;
						//	best_jj = j;
						//}
					}

					inst->Heur_X_ic[sol][i][copy_group_of_i[i]] = 1;
					inst->Heur_X_ic[sol][i][copy_group_of_i[j]] = 0;
					inst->Heur_X_ic[sol][j][copy_group_of_i[i]] = 0;
					inst->Heur_X_ic[sol][j][copy_group_of_i[j]] = 1;

					group_of_i[i] = copy_group_of_i[i];
					group_of_i[j] = copy_group_of_i[j];

					//for (int ii = 0; ii < inst->n_individuals; ii++)
					//{
					//	for (int cc = 0; cc < inst->n_groups; cc++)
					//	{
					//		if (inst->Heur_X_ic[sol][ii][cc] >= 1)
					//		{
					//			group_of_i[ii] = cc;
					//		}
					//		//cout << "The group of " << i << " is " << group_of_i[i] << endl;
					//	}
					//	//cout << endl;
					//}

					group_feature(inst, copy_group_of_i[i], sol);
					group_feature(inst, copy_group_of_i[j], sol);

					objective(inst, 0, sol);

				}
			}

			cout << "UPDATE. Best_objective =  " << Best_objective << ", best_objective = "<< best_objective  << ", i = " << best_ii << " j = " << best_jj << endl;

			inst->Heur_X_ic[sol][best_ii][group_of_i[best_ii]] = 0;
			inst->Heur_X_ic[sol][best_ii][group_of_i[best_jj]] = 1;
			inst->Heur_X_ic[sol][best_jj][group_of_i[best_jj]] = 0;
			inst->Heur_X_ic[sol][best_jj][group_of_i[best_ii]] = 1;

			for (int ii = 0; ii < inst->n_individuals; ii++)
			{
				for (int cc = 0; cc < inst->n_groups; cc++)
				{
					if (inst->Heur_X_ic[sol][ii][cc] >= 1)
					{
						group_of_i[ii] = cc;
					}
					//cout << "The group of " << i << " is " << group_of_i[i] << endl;
				}
				//cout << endl;
			}

			ExploredSolution.push_back(deepCopyArray(inst, group_of_i));

			if (Best_objective < best_objective) 
			{
				Best_objective = best_objective;

				inst->Tabu_X_ic[best_ii][copy_group_of_i[best_ii]] = 0;
				inst->Tabu_X_ic[best_ii][copy_group_of_i[best_jj]] = 1;
				inst->Tabu_X_ic[best_jj][copy_group_of_i[best_jj]] = 0;
				inst->Tabu_X_ic[best_jj][copy_group_of_i[best_ii]] = 1;

				//cout << "[IMROVING] INCUMBENT SOLUTION IN THE TABU SEARCH = " << best_objective  << endl;

			}
			else
			{
				not_improving++;

				//cout << "[NOT IMROVING] INCUMBENT SOLUTION IN THE TABU SEARCH = " << best_objective << endl;

				if (not_improving < inst->Tabu_not_improving_max)
				{					
					break;				
				}

			}		
		}

		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int c = 0; c < inst->n_groups; c++)
			{
				//cout << inst->Heur_X_ic[sol][i][c] << " "; 
				inst->Heur_X_ic[sol][i][c] = inst->Tabu_X_ic[i][c];
			}
			//cout << endl;
		}

		objective(inst, 0, sol);

		// Cleanup (manually clear the vector and deallocate memory)
		for (int* solution : ExploredSolution)
		{
			delete[] solution;
		}
		ExploredSolution.clear();

		cout << endl << "FINAL OBJECTIVE HEURISTIC IN THE TABU SEARCH = " << inst->sol_val_HEUR_DOWN[sol] << endl;

	}

}


void R_MN_SAP_SimulatedAnnealing(instance* inst)
{
	int new_objective;
	int best_objective;
	int Best_objective;

	int gg_c;
	int gg_cc;
	int gg_diff;
	int infeasible;
	int infeasible_new;

	int not_improving;
	int group_size_c;

	int* group_size = new int[inst->n_groups];
	int* group_of_i = new int[inst->n_individuals];

	int Check_theta_feasibility;
	int Theta_iter_limit = 2 * inst->n_groups * inst->n_individuals;

	int best_ii = 0;
	int best_jj = 0;

	cout << endl << "IMPROVING " << inst->Heu_sol << " SOLUTIONS BY SA" << endl << endl;

	#pragma omp parallel for schedule(guided)
	for (int sol = 0; sol < inst->Heu_sol; sol++)
	{
		//std::random_device rd; // Seed for randomness
		//std::mt19937 gen(rd()); // Mersenne Twister engine
		//std::mt19937 gen(21051983); // Fixed seed
		std::uniform_real_distribution<> dis(0.0, 1.0); // Range [0, 1)
		double random_number;
		int Iter_counter = 0;
		double Temperature = 1.0;

		// Create a vector of arrays
		std::vector<int*> ExploredSolution;

		cout << "............................................" << endl;

		for (int i = 0; i < inst->n_individuals; i++)
		{
			for (int c = 0; c < inst->n_groups; c++)
			{
				if (inst->Heur_X_ic[sol][i][c] >= 1)
				{
					group_of_i[i] = c;
				}
				//cout << "The group of " << i << " is " << group_of_i[i] << endl;
			}
			//cout << endl;
		}

		objective(inst, 0, sol);

		init_group_feature(inst, sol);

		best_objective = inst->sol_val_HEUR_DOWN[sol];
		Best_objective = inst->sol_val_HEUR_DOWN[sol];

		int* copy_group_of_i;

		while (Iter_counter < inst->SA_Iter_max)
		{
			Temperature = (1.0 - (double)(Iter_counter + 1.0)/inst->SA_Iter_max);

			//infeasible = Check_all_theta_constr(inst);
			//int infea = Check_unique_assignment_constr(inst, sol);
			//cout << endl << "-----------------infeasible : " << infeasible << " " << infea << endl;

			copy_group_of_i = deepCopyArray(inst, group_of_i);

			int i = 0;
			int j = 0;
			bool same_group_ij = (group_of_i[i] == group_of_i[j]);

			while (same_group_ij) 
			{
				i = generateRandomIndex(inst, sol, inst->n_groups + 1);
				j = generateRandomIndex(inst, sol, inst->n_groups + 1);

				same_group_ij = (group_of_i[i] == group_of_i[j]);

				if (!same_group_ij)
				{
					break;
				}
			}

			// Switch group

			inst->Heur_X_ic[sol][i][copy_group_of_i[i]] = 0;
			inst->Heur_X_ic[sol][i][copy_group_of_i[j]] = 1;
			inst->Heur_X_ic[sol][j][copy_group_of_i[i]] = 1;
			inst->Heur_X_ic[sol][j][copy_group_of_i[j]] = 0;

			group_of_i[i] = copy_group_of_i[j];
			group_of_i[j] = copy_group_of_i[i];

			group_feature(inst, group_of_i[i], sol);
			group_feature(inst, group_of_i[j], sol);

			//out << "individual " << i << "(" << group_of_i[i] << ") , and " << j << "(" << group_of_i[j] << ") " << "Best objective = " << best_objective << ", New objective = " << new_objective << endl;

			if (Check_all_theta_constr(inst) <= 0)
			{
				objective(inst, 0, sol);

				new_objective = inst->sol_val_HEUR_DOWN[sol];

				double TransitionProb = (double)exp((((double)new_objective - (double)best_objective)/Temperature));

				random_number = dis(gen);

				if (TransitionProb > random_number)
				{
					//cout << "ITER = " << Iter_counter << endl;
					// Generate a random number

					//cout << "UPDATE. Best_obj =  " << best_objective << ", new_obj = " << new_objective << ", Trans prob = " << TransitionProb << ", Rand = " << random_number << ", i = " << i << " j = " <<j << endl;

					best_objective = new_objective;

					best_ii = i;
					best_jj = j;

					for (int ii = 0; ii < inst->n_individuals; ii++)
					{
						for (int cc = 0; cc < inst->n_groups; cc++)
						{
							if (inst->Heur_X_ic[sol][ii][cc] >= 1)
							{
								group_of_i[ii] = cc;
							}
							//cout << "The group of " << i << " is " << group_of_i[i] << endl;
						}
						//cout << endl;
					}
				}
				else 
				{
					inst->Heur_X_ic[sol][i][copy_group_of_i[i]] = 1;
					inst->Heur_X_ic[sol][i][copy_group_of_i[j]] = 0;
					inst->Heur_X_ic[sol][j][copy_group_of_i[i]] = 0;
					inst->Heur_X_ic[sol][j][copy_group_of_i[j]] = 1;

					group_of_i[i] = copy_group_of_i[i];
					group_of_i[j] = copy_group_of_i[j];

					group_feature(inst, copy_group_of_i[i], sol);
					group_feature(inst, copy_group_of_i[j], sol);

					objective(inst, 0, sol);

				}

			}
			else
			{
				inst->Heur_X_ic[sol][i][copy_group_of_i[i]] = 1;
				inst->Heur_X_ic[sol][i][copy_group_of_i[j]] = 0;
				inst->Heur_X_ic[sol][j][copy_group_of_i[i]] = 0;
				inst->Heur_X_ic[sol][j][copy_group_of_i[j]] = 1;

				group_of_i[i] = copy_group_of_i[i];
				group_of_i[j] = copy_group_of_i[j];

				group_feature(inst, copy_group_of_i[i], sol);
				group_feature(inst, copy_group_of_i[j], sol);

				objective(inst, 0, sol);

			}

			Iter_counter++;

		}//while


		if (Best_objective > best_objective)
		{
			for (int i = 0; i < inst->n_individuals; i++)
			{
				for (int c = 0; c < inst->n_groups; ++c)
				{
					inst->Heur_X_ic[sol][i][c] = inst->X_ic[sol][i][c];
				}
			}

			init_group_feature(inst, sol);
		}

		objective(inst, 0, sol);

		cout << endl << "FINAL OBJECTIVE HEURISTIC IN THE SA = " << inst->sol_val_HEUR_DOWN[sol] << endl;

	}

}




void free_data(instance* inst)
/*****************************************************************/
{
	// Delete n_neighbors_up and n_neighbors_down
	delete[] inst->n_neighbors_up;
	delete[] inst->n_neighbors_down;

	// Delete neighbors_up and neighbors_down
	for (int i = 0; i < inst->n_individuals; i++) {
		delete[] inst->neighbors_up[i];
		delete[] inst->neighbors_down[i];
	}
	delete[] inst->neighbors_up;
	delete[] inst->neighbors_down;

	// Delete alpha_data
	for (int i = 0; i < inst->n_individuals; i++) {
		delete[] inst->alpha_data[i];
	}
	delete[] inst->alpha_data;

	// Delete beta_data
	for (int i = 0; i < inst->n_individuals; i++) {
		delete[] inst->beta_data[i];
	}
	delete[] inst->beta_data;

}


