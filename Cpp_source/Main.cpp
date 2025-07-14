

#define _CRT_SECURE_NO_WARNINGS

#include <chrono>
#include "global_functions.h"
#include "global_variables.h"
#include "MILP_min_sum.h"
#include "MILP_min_max.h"
#include "BEN_min_max.h"
#include "BEN_min_sum.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>  // OpenMP header for parallel programming

/*****************************************************************/
int main(int argc, char** argv)
/*****************************************************************/
{
	instance inst;

	bool debugging = false;
	inst.n_cuts = 1;
	bool reoptimize = true;

	if (debugging) 
	{
		inst.n_individuals = 21;
		inst.n_features = 2;
		inst.n_groups = 3;

		inst.opt_bound = 60;
		inst.output_file_data = "Random_output_table.txt";
		inst.output_file_x = "Random_output_solution.txt";
		inst.iteration_BEN = 0;

		generate_random_data(&inst);

		/*Param8*/inst.method = 1;			// 0 if B&C, 1 if Benders
		/*Param8*/inst.problem = 1;			// 0 if MinSum, 1 if MinMax
		/*Param8*/inst.Heuristic_only = 0;			// 0 if no Taylor cut, 1 if Taylor cut 
		/*Param8*/inst.LB_heur = 0;			// 0 if no LB cut, 1 if LB cut
		/*Param8*/inst.UB_heur = 0;			// 0 if no UB cut, 1 if UB cut
		/*Param8*/inst.MagnantiW = 0;		// 0 if the first cut is included
		/*Param9*/inst.external_sol = 1;	// if this parameter is equal to zero, no 10-13 cut can be applied
		/*Param10*/inst.number_of_CPU = 8;
		/*Param11*/inst.timelimit = 60;
		/*Param12*/inst.TOLL_OPTIMALITY = 0.0;

		inst.option = 1; //CPLEX B&C

	}else{

		inst.opt_bound = 60;
		inst.output_file_data = "output_table.txt";
		inst.output_file_x = "output_solution.txt";
		inst.iteration_BEN = 0;
		inst.solution_time_HEUR = 0.0;
		inst.Heu_sol = 1;
		inst.Tabu_not_improving_max = 5;
		
		strcpy_s(inst.alpha_file_data, sizeof(inst.alpha_file_data), argv[1]);
		strcpy_s(inst.beta_file_data, sizeof(inst.beta_file_data), argv[2]);
		strcpy_s(inst.theta_file_data, sizeof(inst.theta_file_data), argv[3]);
		strcpy_s(inst.x_file_data, sizeof(inst.x_file_data), argv[4]);

		/*Param5*/inst.n_individuals = (int)atof(argv[5]);
		/*Param6*/inst.n_features = (int)atof(argv[6]);
		/*Param7*/inst.n_groups = (int)atof(argv[7]);
		/*Param8*/inst.method = (int)atof(argv[8]);					// 0 if B&C, 1 if Benders
		/*Param9*/inst.problem = (int)atof(argv[9]);				// 0 if MinSum, 1 if MinMax
		/*Param10*/inst.Heuristic_only = (int)atof(argv[10]);		// 1 if only the heuristics has to be applied 
		/*Param11*/inst.LB_heur = (int)atof(argv[11]);				// 0 if no LB cut, 1 if LB cut
		/*Param12*/inst.UB_heur = (int)atof(argv[12]);				// 0 if no UB cut, 1 if UB cut
		/*Param13*/inst.MagnantiW = (int)atof(argv[13]);			// 0 if the first cut is included
		/*Param14*/inst.external_sol = (int)atof(argv[14]);			// if this parameter is equal to zero, no 10-13 cut can be applied
		/*Param15*/inst.number_of_CPU = (int)atof(argv[15]);
		/*Param16*/inst.timelimit = (int)atof(argv[16]);
		/*Param17*/inst.TOLL_OPTIMALITY = (double)atof(argv[17]);
		/*Param18*/inst.warm_start = (bool)atof(argv[18]);
		/*Param14*/inst.get_feasible = (bool)atof(argv[19]);

	
		inst.SA_Iter_max = inst.n_individuals * inst.n_groups;

		//inst.warm_start = false;
		//inst.get_feasible = false;

		read_data_alpha(&inst);
		read_data_beta(&inst);
		read_data_theta(&inst);

		inst.option = 1; //CPLEX B&C
		
		if (inst.external_sol <= 0) // no initial solution is used
		{
			inst.Heuristic_only = 0;
			inst.LB_heur = 0;
			inst.UB_heur = 0;
			inst.MagnantiW = 0;

			inst.n_sol = inst.Heu_sol;

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = -1;
				inst.sol_val_UP[s] = (int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
				inst.sol_val_DOWN[s] = -(int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
			}

		}
		
		// -------------------------
		// inst.external_sol == 1
		// -------------------------
		//
		// Importing inst.n_sol external solutions from local files and
		// using them directly to create Benders cuts and warm starts
		//
		// -------------------------

		if (inst.external_sol == 1) //an external solution is included
		{
			cout << "INITIAL SOLUTION CONFIGURATION: " << inst.external_sol << endl << endl;

			read_data_x(&inst);

			//--------------------------------------------------
			// Initialize arrays
			//--------------------------------------------------

			inst.Heu_sol = inst.n_sol;

			// Initialize adjacency matrix and list of neighbors
			inst.Heur_X_ic = new int** [inst.Heu_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				// Initialize adjacency matrix and list of neighbors
				inst.Heur_X_ic[s] = new int* [inst.n_individuals];

				for (int i = 0; i < inst.n_individuals; i++)
				{
					inst.Heur_X_ic[s][i] = new int[inst.n_groups];

					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.Heur_X_ic[s][i][c] = 0;
					}
				}
			}	

			// Order the solutions, so that the first is always the best known solution

			int Best_sol_so_fat = -pow(inst.n_individuals, 3);
			int BestIndex = -1;

			for (int s = 0; s < inst.n_sol; s++)
			{
				if (Best_sol_so_fat <= inst.sol_val_EXTERNAL[s])
				{
					Best_sol_so_fat = inst.sol_val_EXTERNAL[s];
					BestIndex = s;
				}
			}

			inst.sol_val_EXTERNAL[BestIndex] = inst.sol_val_EXTERNAL[0];
			inst.sol_val_EXTERNAL[0] = Best_sol_so_fat;

			int in_the_meantime = 0;
			for (int i = 0; i < inst.n_individuals; i++)
			{
				for (int c = 0; c < inst.n_groups; ++c)
				{
					in_the_meantime = inst.X_ic[0][i][c];
					inst.X_ic[0][i][c] = inst.X_ic[BestIndex][i][c];
					inst.X_ic[BestIndex][i][c] = in_the_meantime;
				}
			}

			inst.sol_val_UP = new int [inst.n_sol];
			inst.sol_val_DOWN = new int [inst.n_sol];
			inst.sol_val_HEUR_DOWN = new int[inst.n_sol];

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_HEUR_DOWN[s] = inst.sol_val_EXTERNAL[s];
				inst.sol_val_DOWN[s] = inst.sol_val_EXTERNAL[s];
				inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_DOWN[s]) + 2;
				//inst.sol_val_UP = (int)round((1 + 0.1) * inst.sol_val_DOWN) + 2;
				inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_DOWN[s]) - 2;
				//inst.sol_val_DOWN = (int)round((1 - 0.7) * inst.sol_val_DOWN) - 2;
				cout << "SOLUTION " << s << ", STARTING VALUE " << inst.sol_val_EXTERNAL[s] << "............................................" << endl;
			}
		}

		// -------------------------
		// inst.external_sol == 2
		// -------------------------
		//
		// Importing inst.n_sol external solutions from local files 
		// and improving them using a local search previous to use   
		// these solutions to create Benders cuts and warm starts
		//
		// -------------------------

		if (inst.external_sol == 2) //an external solution is included but improved by the internal LS
		{
			cout << "INITIAL SOLUTION CONFIGURATION: " << inst.external_sol << endl << endl;

			read_data_x(&inst);

			inst.Heu_sol = inst.n_sol;

			//--------------------------------------------------
			// Initialize arrays
			//--------------------------------------------------

			// Initialize adjacency matrix and list of neighbors
			inst.Heur_X_ic = new int** [inst.Heu_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				// Initialize adjacency matrix and list of neighbors
				inst.Heur_X_ic[s] = new int* [inst.n_individuals];

				for (int i = 0; i < inst.n_individuals; i++)
				{
					inst.Heur_X_ic[s][i] = new int[inst.n_groups];

					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.Heur_X_ic[s][i][c] = inst.X_ic[s][i][c];
					}
				}
			}

			inst.group_features = new int* [inst.n_features];
			for (int f = 0; f < inst.n_features; f++)
			{
				inst.group_features[f] = new int[inst.n_groups];
			}

			inst.sol_val_UP = new int[inst.n_sol];
			inst.sol_val_DOWN = new int[inst.n_sol];
			inst.sol_val_HEUR_DOWN = new int[inst.n_sol];

			auto start_time = std::chrono::high_resolution_clock::now();

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				inst.sol_val_HEUR_DOWN[s] = inst.sol_val_EXTERNAL[s];
			}

			R_MN_SAP_LocalSearch(&inst);

			auto end_time = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time) / 1000000.0;

			inst.solution_time_HEUR = static_cast<double>(duration.count());

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				for (int i = 0; i < inst.n_individuals; i++)
				{
					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.X_ic[s][i][c] = inst.Heur_X_ic[s][i][c];
					}
				}
			}

			// Order the solutions, so that the first is always the best known solution

			int Best_sol_so_fat = -pow(inst.n_individuals,3);
			int BestIndex = -1;

			for (int s = 0; s < inst.n_sol; s++)
			{
				if (Best_sol_so_fat <= inst.sol_val_HEUR_DOWN[s])
				{
					Best_sol_so_fat = inst.sol_val_HEUR_DOWN[s];
					BestIndex = s;
				}
			}

			int in_the_meantime = 0;
			for (int i = 0; i < inst.n_individuals; i++)
			{
				for (int c = 0; c < inst.n_groups; ++c)
				{
					in_the_meantime = inst.X_ic[0][i][c];
					inst.X_ic[0][i][c] = inst.X_ic[BestIndex][i][c];
					inst.X_ic[BestIndex][i][c] = in_the_meantime;
				}
			}

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_DOWN[s] = inst.sol_val_HEUR_DOWN[s];
				inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_DOWN[s]) + 2;
				//inst.sol_val_UP = (int)round((1 + 0.1) * inst.sol_val_DOWN) + 2;
				inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_DOWN[s]) - 2;
				//inst.sol_val_DOWN = (int)round((1 - 0.7) * inst.sol_val_DOWN) - 2;
			}
		}

		// -------------------------
		// inst.external_sol == 3
		// -------------------------
		//
		// Importing inst.n_sol external solutions from local files 
		// and improving them based on a tabu search previous to use  
		// these solutions to create Benders cuts and warm starts
		//
		// -------------------------

		if (inst.external_sol == 3) 
		{
			cout << "INITIAL SOLUTION CONFIGURATION: " << inst.external_sol << endl << endl;

			read_data_x(&inst);

			inst.Heu_sol = inst.n_sol;

			//--------------------------------------------------
			// Initialize arrays
			//--------------------------------------------------

			// Initialize adjacency matrix and list of neighbors
			inst.Heur_X_ic = new int** [inst.Heu_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				// Initialize adjacency matrix and list of neighbors
				inst.Heur_X_ic[s] = new int* [inst.n_individuals];

				for (int i = 0; i < inst.n_individuals; i++)
				{
					inst.Heur_X_ic[s][i] = new int[inst.n_groups];

					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.Heur_X_ic[s][i][c] = inst.X_ic[s][i][c];
					}
				}
			}

			inst.group_features = new int* [inst.n_features];
			for (int f = 0; f < inst.n_features; f++)
			{
				inst.group_features[f] = new int[inst.n_groups];
			}

			inst.sol_val_UP = new int[inst.n_sol];
			inst.sol_val_DOWN = new int[inst.n_sol];
			inst.sol_val_HEUR_DOWN = new int[inst.n_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				inst.sol_val_HEUR_DOWN[s] = inst.sol_val_EXTERNAL[s];
			}

			auto start_time = std::chrono::high_resolution_clock::now();

			R_MN_SAP_TabuSearch(&inst);

			auto end_time = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time) / 1000000.0;

			inst.solution_time_HEUR = static_cast<double>(duration.count());

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				for (int i = 0; i < inst.n_individuals; i++)
				{
					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.X_ic[s][i][c] = inst.Heur_X_ic[s][i][c];
					}
				}
			}

			int Best_sol_so_fat = -pow(inst.n_individuals, 3);
			int BestIndex = -1;

			for (int s = 0; s < inst.n_sol; s++)
			{
				if (Best_sol_so_fat <= inst.sol_val_HEUR_DOWN[s])
				{
					Best_sol_so_fat = inst.sol_val_HEUR_DOWN[s];
					BestIndex = s;
				}
			}

			int in_the_meantime = 0;
			for (int i = 0; i < inst.n_individuals; i++)
			{
				for (int c = 0; c < inst.n_groups; ++c)
				{
					in_the_meantime = inst.X_ic[0][i][c];
					inst.X_ic[0][i][c] = inst.X_ic[BestIndex][i][c];
					inst.X_ic[BestIndex][i][c] = in_the_meantime;
				}
			}

			for (int s = 0; s < inst.n_sol; s++)
			{
				//objective(&inst, 0, s);

				inst.sol_val_DOWN[s] = inst.sol_val_HEUR_DOWN[s];
				inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_DOWN[s]) + 2;
				//inst.sol_val_UP = (int)round((1 + 0.1) * inst.sol_val_DOWN) + 2;
				inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_DOWN[s]) - 2;
				//inst.sol_val_DOWN = (int)round((1 - 0.7) * inst.sol_val_DOWN) - 2;
			}
		}
		
		// -------------------------
		// inst.external_sol == 4
		// -------------------------
		//
		// Importing inst.n_sol external solutions from local files 
		// and improving them based on a tabu search previous to use  
		// these solutions to create Benders cuts and warm starts
		//
		// -------------------------

		if (inst.external_sol == 4)
		{
			cout << "INITIAL SOLUTION CONFIGURATION: " << inst.external_sol << endl << endl;

			read_data_x(&inst);

			inst.Heu_sol = inst.n_sol;

			//--------------------------------------------------
			// Initialize arrays
			//--------------------------------------------------

			// Initialize adjacency matrix and list of neighbors
			inst.Heur_X_ic = new int** [inst.Heu_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				// Initialize adjacency matrix and list of neighbors
				inst.Heur_X_ic[s] = new int* [inst.n_individuals];

				for (int i = 0; i < inst.n_individuals; i++)
				{
					inst.Heur_X_ic[s][i] = new int[inst.n_groups];

					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.Heur_X_ic[s][i][c] = inst.X_ic[s][i][c];
					}
				}
			}

			inst.group_features = new int* [inst.n_features];
			for (int f = 0; f < inst.n_features; f++)
			{
				inst.group_features[f] = new int[inst.n_groups];
			}

			inst.sol_val_UP = new int[inst.n_sol];
			inst.sol_val_DOWN = new int[inst.n_sol];
			inst.sol_val_HEUR_DOWN = new int[inst.n_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				inst.sol_val_HEUR_DOWN[s] = inst.sol_val_EXTERNAL[s];
			}

			auto start_time = std::chrono::high_resolution_clock::now();

			R_MN_SAP_SimulatedAnnealing(&inst);

			auto end_time = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time) / 1000000.0;

			inst.solution_time_HEUR = static_cast<double>(duration.count());

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				for (int i = 0; i < inst.n_individuals; i++)
				{
					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.X_ic[s][i][c] = inst.Heur_X_ic[s][i][c];
					}
				}
			}

			int Best_sol_so_fat = -pow(inst.n_individuals, 3);
			int BestIndex = -1;

			for (int s = 0; s < inst.n_sol; s++)
			{
				if (Best_sol_so_fat <= inst.sol_val_HEUR_DOWN[s])
				{
					Best_sol_so_fat = inst.sol_val_HEUR_DOWN[s];
					BestIndex = s;
				}
			}

			int in_the_meantime = 0;
			for (int i = 0; i < inst.n_individuals; i++)
			{
				for (int c = 0; c < inst.n_groups; ++c)
				{
					in_the_meantime = inst.X_ic[0][i][c];
					inst.X_ic[0][i][c] = inst.X_ic[BestIndex][i][c];
					inst.X_ic[BestIndex][i][c] = in_the_meantime;
				}
			}

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_DOWN[s] = inst.sol_val_HEUR_DOWN[s];
				inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_DOWN[s]) + 2;
				//inst.sol_val_UP = (int)round((1 + 0.1) * inst.sol_val_DOWN) + 2;
				inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_DOWN[s]) - 2;
				//inst.sol_val_DOWN = (int)round((1 - 0.7) * inst.sol_val_DOWN) - 2;
			}
		}

		// -------------------------
		// inst.external_sol == 5
		// -------------------------
		//
		// Generating inst.n_sol solutions using a greedy method and
		// using them directly to create Benders cuts and warm starts
		//
		// -------------------------

		if (inst.external_sol == 5) //hierarchical heuristic is used as initial solution
		{
			cout << "INITIAL SOLUTION CONFIGURATION: " << inst.external_sol << endl << endl;

			inst.n_sol = 1; // Generate n_sol initial feasible solutions using HC
			inst.Heu_sol = inst.n_sol; // Generate n_sol initial feasible solutions using HC

			//--------------------------------------------------
			// Initialize arrays
			//--------------------------------------------------

			// Initialize adjacency matrix and list of neighbors
			inst.Heur_X_ic = new int** [inst.Heu_sol];
			inst.X_ic = new int** [inst.Heu_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				// Initialize adjacency matrix and list of neighbors
				inst.Heur_X_ic[s] = new int* [inst.n_individuals];
				inst.X_ic[s] = new int* [inst.n_individuals];

				for (int i = 0; i < inst.n_individuals; i++)
				{
					inst.X_ic[s][i] = new int[inst.n_groups];
					inst.Heur_X_ic[s][i] = new int[inst.n_groups];

					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.X_ic[s][i][c] = 0;
						inst.Heur_X_ic[s][i][c] = 0;
					}
				}
			}

			inst.group_features = new int* [inst.n_features];
			for (int f = 0; f < inst.n_features; f++)
			{
				inst.group_features[f] = new int[inst.n_groups];
			}

			inst.sol_val_UP = new int[inst.n_sol];
			inst.sol_val_DOWN = new int[inst.n_sol];
			inst.sol_val_EXTERNAL = new int[inst.n_sol];
			inst.sol_val_HEUR_DOWN = new int[inst.n_sol];

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = 0;
				inst.sol_val_HEUR_DOWN[s] = 0;
				inst.sol_val_UP[s] = (int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
				inst.sol_val_DOWN[s] = -(int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
			}

			auto start_time = std::chrono::high_resolution_clock::now();

			R_MN_SAP_ConstructuveHeur(&inst);

			auto end_time = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time) / 1000000.0;

			inst.solution_time_HEUR = static_cast<double>(duration.count());

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = inst.sol_val_HEUR_DOWN[s];

				inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_DOWN[s]) + 2;
				//inst.sol_val_UP = (int)round((1 + 0.1) * inst.sol_val_DOWN) + 2;
				inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_DOWN[s]) - 2;
				//inst.sol_val_DOWN = (int)round((1 - 0.7) * inst.sol_val_DOWN) - 2;		
			}

			for (int s = 0; s < inst.n_sol; s++)
			{
				for (int i = 0; i < inst.n_individuals; i++)
				{
					for (int c = 0; c < inst.n_groups; c++)
					{
						inst.X_ic[s][i][c] = inst.Heur_X_ic[s][i][c];
					}
				}

				if (inst.sol_val_HEUR_DOWN[s] > 0)
				{
					inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_HEUR_DOWN[s]) + 2;
					inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_HEUR_DOWN[s]) - 2;
				}
				if (inst.sol_val_HEUR_DOWN[s] < 0)
				{
					inst.sol_val_UP[s] = (int)round((1 - 0.2) * inst.sol_val_HEUR_DOWN[s]) + 2;
					inst.sol_val_DOWN[s] = (int)round((1 + 0.2) * inst.sol_val_HEUR_DOWN[s]) - 2;
				}
				if (inst.sol_val_HEUR_DOWN[s] == 0)
				{
					inst.sol_val_UP[s] = 2;
					inst.sol_val_DOWN[s] = -2;
				}
			}

			int Best_sol_so_fat = -pow(inst.n_individuals, 3);
			int BestIndex = -1;

			for (int s = 0; s < inst.n_sol; s++)
			{
				if (Best_sol_so_fat <= inst.sol_val_HEUR_DOWN[s])
				{
					Best_sol_so_fat = inst.sol_val_HEUR_DOWN[s];
					BestIndex = s;
				}
			}

			int in_the_meantime = 0;
			for (int i = 0; i < inst.n_individuals; i++)
			{
				for (int c = 0; c < inst.n_groups; ++c)
				{
					in_the_meantime = inst.X_ic[0][i][c];
					inst.X_ic[0][i][c] = inst.X_ic[BestIndex][i][c];
					inst.X_ic[BestIndex][i][c] = in_the_meantime;
				}
			}

		}

		// -------------------------
		// inst.external_sol == 6
		// -------------------------
		//
		// Generating inst.n_sol solutions using a greedy method and
		// and improving them based on a local search previous to use  
		// these solutions to create Benders cuts and warm starts
		//
		// -------------------------

		if (inst.external_sol == 6) //hierarchical heuristic is used as initial solution
		{
			cout << "INITIAL SOLUTION CONFIGURATION: " << inst.external_sol << endl << endl;

			inst.n_sol = 2;
			inst.Heu_sol = 2;

			//--------------------------------------------------
			// Initialize arrays
			//--------------------------------------------------

			// Initialize adjacency matrix and list of neighbors
			inst.Heur_X_ic = new int** [inst.Heu_sol];
			inst.X_ic = new int** [inst.Heu_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				// Initialize adjacency matrix and list of neighbors
				inst.Heur_X_ic[s] = new int* [inst.n_individuals];
				inst.X_ic[s] = new int* [inst.n_individuals];

				for (int i = 0; i < inst.n_individuals; i++)
				{
					inst.X_ic[s][i] = new int[inst.n_groups];
					inst.Heur_X_ic[s][i] = new int[inst.n_groups];

					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.X_ic[s][i][c] = 0;
						inst.Heur_X_ic[s][i][c] = 0;
					}
				}
			}

			inst.group_features = new int* [inst.n_features];
			for (int f = 0; f < inst.n_features; f++)
			{
				inst.group_features[f] = new int[inst.n_groups];
			}

			inst.sol_val_UP = new int[inst.n_sol];
			inst.sol_val_DOWN = new int[inst.n_sol];
			inst.sol_val_EXTERNAL = new int[inst.n_sol];
			inst.sol_val_HEUR_DOWN = new int[inst.n_sol];

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = 0;
				inst.sol_val_HEUR_DOWN[s] = 0;
				inst.sol_val_UP[s] = (int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
				inst.sol_val_DOWN[s] = -(int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
			}

			R_MN_SAP_ConstructuveHeur(&inst);

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = inst.sol_val_HEUR_DOWN[s];
			}

			auto start_time = std::chrono::high_resolution_clock::now();

			R_MN_SAP_LocalSearch(&inst);

			auto end_time = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time) / 1000000.0;

			inst.solution_time_HEUR = static_cast<double>(duration.count());

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_DOWN[s]) + 2;
				//inst.sol_val_UP = (int)round((1 + 0.1) * inst.sol_val_DOWN) + 2;
				inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_DOWN[s]) - 2;
				//inst.sol_val_DOWN = (int)round((1 - 0.7) * inst.sol_val_DOWN) - 2;

				//print_data(&inst);
				//cout << inst.sol_val_EXTERNAL[s] << endl;				
			}

			for (int s = 0; s < inst.n_sol; s++)
			{
				for (int i = 0; i < inst.n_individuals; i++)
				{
					for (int c = 0; c < inst.n_groups; c++)
					{
						inst.X_ic[s][i][c] = inst.Heur_X_ic[s][i][c];
					}
				}

				if (inst.sol_val_HEUR_DOWN[s] > 0)
				{
					inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_HEUR_DOWN[s]) + 2;
					inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_HEUR_DOWN[s]) - 2;
				}
				if (inst.sol_val_HEUR_DOWN[s] < 0)
				{
					inst.sol_val_UP[s] = (int)round((1 - 0.2) * inst.sol_val_HEUR_DOWN[s]) + 2;
					inst.sol_val_DOWN[s] = (int)round((1 + 0.2) * inst.sol_val_HEUR_DOWN[s]) - 2;
				}
				if (inst.sol_val_HEUR_DOWN[s] == 0)
				{
					inst.sol_val_UP[s] = 2;
					inst.sol_val_DOWN[s] = -2;
				}
			}

			int Best_sol_so_fat = -pow(inst.n_individuals, 3);
			int BestIndex = -1;

			for (int s = 0; s < inst.n_sol; s++)
			{
				if (Best_sol_so_fat <= inst.sol_val_HEUR_DOWN[s])
				{
					Best_sol_so_fat = inst.sol_val_HEUR_DOWN[s];
					BestIndex = s;
				}
			}

			int in_the_meantime = 0;
			for (int i = 0; i < inst.n_individuals; i++)
			{
				for (int c = 0; c < inst.n_groups; ++c)
				{
					in_the_meantime = inst.X_ic[0][i][c];
					inst.X_ic[0][i][c] = inst.X_ic[BestIndex][i][c];
					inst.X_ic[BestIndex][i][c] = in_the_meantime;
				}
			}

		}

		// -------------------------
		// inst.external_sol == 7
		// -------------------------
		//
		// Generating inst.n_sol solutions using a greedy method and
		// and improving them based on a tabu search previous to use  
		// these solutions to create Benders cuts and warm starts
		//
		// -------------------------

		if (inst.external_sol == 7) //hierarchical heuristic is used as initial solution
		{
			cout << "INITIAL SOLUTION CONFIGURATION: " << inst.external_sol << endl << endl;

			inst.n_sol = 2;
			inst.Heu_sol = 2;

			//--------------------------------------------------
			// Initialize arrays
			//--------------------------------------------------

			// Initialize adjacency matrix and list of neighbors
			inst.Heur_X_ic = new int** [inst.Heu_sol];
			inst.X_ic = new int** [inst.Heu_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				// Initialize adjacency matrix and list of neighbors
				inst.Heur_X_ic[s] = new int* [inst.n_individuals];
				inst.X_ic[s] = new int* [inst.n_individuals];

				for (int i = 0; i < inst.n_individuals; i++)
				{
					inst.X_ic[s][i] = new int[inst.n_groups];
					inst.Heur_X_ic[s][i] = new int[inst.n_groups];

					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.X_ic[s][i][c] = 0;
						inst.Heur_X_ic[s][i][c] = 0;
					}
				}
			}

			inst.group_features = new int* [inst.n_features];
			for (int f = 0; f < inst.n_features; f++)
			{
				inst.group_features[f] = new int[inst.n_groups];
			}

			inst.sol_val_UP = new int[inst.n_sol];
			inst.sol_val_DOWN = new int[inst.n_sol];
			inst.sol_val_EXTERNAL = new int[inst.n_sol];
			inst.sol_val_HEUR_DOWN = new int[inst.n_sol];

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = 0;
				inst.sol_val_HEUR_DOWN[s] = 0;
				inst.sol_val_UP[s] = (int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
				inst.sol_val_DOWN[s] = -(int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
			}
			auto start_time = std::chrono::high_resolution_clock::now();

			R_MN_SAP_ConstructuveHeur(&inst);

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = inst.sol_val_HEUR_DOWN[s];
			}

			R_MN_SAP_TabuSearch(&inst);

			auto end_time = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time) / 1000000.0;

			inst.solution_time_HEUR = static_cast<double>(duration.count());

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_DOWN[s]) + 2;
				//inst.sol_val_UP = (int)round((1 + 0.1) * inst.sol_val_DOWN) + 2;
				inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_DOWN[s]) - 2;
				//inst.sol_val_DOWN = (int)round((1 - 0.7) * inst.sol_val_DOWN) - 2;
		
			}

			for (int s = 0; s < inst.n_sol; s++)
			{
				for (int i = 0; i < inst.n_individuals; i++)
				{
					for (int c = 0; c < inst.n_groups; c++)
					{
						inst.X_ic[s][i][c] = inst.Heur_X_ic[s][i][c];
					}
				}

				if (inst.sol_val_HEUR_DOWN[s] > 0)
				{
					inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_HEUR_DOWN[s]) + 2;
					inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_HEUR_DOWN[s]) - 2;
				}
				if (inst.sol_val_HEUR_DOWN[s] < 0)
				{
					inst.sol_val_UP[s] = (int)round((1 - 0.2) * inst.sol_val_HEUR_DOWN[s]) + 2;
					inst.sol_val_DOWN[s] = (int)round((1 + 0.2) * inst.sol_val_HEUR_DOWN[s]) - 2;
				}
				if (inst.sol_val_HEUR_DOWN[s] == 0)
				{
					inst.sol_val_UP[s] = 2;
					inst.sol_val_DOWN[s] = -2;
				}
			}


			int Best_sol_so_fat = -pow(inst.n_individuals, 3);
			int BestIndex = -1;

			for (int s = 0; s < inst.n_sol; s++)
			{
				if (Best_sol_so_fat <= inst.sol_val_HEUR_DOWN[s])
				{
					Best_sol_so_fat = inst.sol_val_HEUR_DOWN[s];
					BestIndex = s;
				}
			}

			int in_the_meantime = 0;
			for (int i = 0; i < inst.n_individuals; i++)
			{
				for (int c = 0; c < inst.n_groups; ++c)
				{
					in_the_meantime = inst.X_ic[0][i][c];
					inst.X_ic[0][i][c] = inst.X_ic[BestIndex][i][c];
					inst.X_ic[BestIndex][i][c] = in_the_meantime;
				}
			}

		}

		// -------------------------
		// inst.external_sol == 8
		// -------------------------
		//
		// Generating inst.n_sol solutions using a greedy method and
		// and improving them based on a simulating annealing previous to use  
		// these solutions to create Benders cuts and warm starts
		//
		// -------------------------

		if (inst.external_sol == 8) //hierarchical heuristic is used as initial solution
		{
			cout << "INITIAL SOLUTION CONFIGURATION: " << inst.external_sol << endl << endl;

			inst.n_sol = 2;
			inst.Heu_sol = 2;

			//--------------------------------------------------
			// Initialize arrays
			//--------------------------------------------------

			// Initialize adjacency matrix and list of neighbors
			inst.Heur_X_ic = new int** [inst.Heu_sol];
			inst.X_ic = new int** [inst.Heu_sol];

			for (int s = 0; s < inst.Heu_sol; s++)
			{
				// Initialize adjacency matrix and list of neighbors
				inst.Heur_X_ic[s] = new int* [inst.n_individuals];
				inst.X_ic[s] = new int* [inst.n_individuals];

				for (int i = 0; i < inst.n_individuals; i++)
				{
					inst.X_ic[s][i] = new int[inst.n_groups];
					inst.Heur_X_ic[s][i] = new int[inst.n_groups];

					for (int c = 0; c < inst.n_groups; ++c)
					{
						inst.X_ic[s][i][c] = 0;
						inst.Heur_X_ic[s][i][c] = 0;
					}
				}
			}

			inst.group_features = new int* [inst.n_features];
			for (int f = 0; f < inst.n_features; f++)
			{
				inst.group_features[f] = new int[inst.n_groups];
			}

			inst.sol_val_UP = new int[inst.n_sol];
			inst.sol_val_DOWN = new int[inst.n_sol];
			inst.sol_val_EXTERNAL = new int[inst.n_sol];
			inst.sol_val_HEUR_DOWN = new int[inst.n_sol];

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = 0;
				inst.sol_val_HEUR_DOWN[s] = 0;
				inst.sol_val_UP[s] = (int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
				inst.sol_val_DOWN[s] = -(int)(1000 * pow(inst.n_individuals, 2) * inst.n_groups);
			}
			auto start_time = std::chrono::high_resolution_clock::now();

			R_MN_SAP_ConstructuveHeur(&inst);

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_EXTERNAL[s] = inst.sol_val_HEUR_DOWN[s];
			}

			R_MN_SAP_SimulatedAnnealing(&inst);

			auto end_time = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time) / 1000000.0;

			inst.solution_time_HEUR = static_cast<double>(duration.count());

			for (int s = 0; s < inst.n_sol; s++)
			{
				inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_DOWN[s]) + 2;
				//inst.sol_val_UP = (int)round((1 + 0.1) * inst.sol_val_DOWN) + 2;
				inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_DOWN[s]) - 2;
				//inst.sol_val_DOWN = (int)round((1 - 0.7) * inst.sol_val_DOWN) - 2;

				//print_data(&inst);
				//cout << inst.sol_val_EXTERNAL[s] << endl;				
			}

			for (int s = 0; s < inst.n_sol; s++)
			{
				for (int i = 0; i < inst.n_individuals; i++)
				{
					for (int c = 0; c < inst.n_groups; c++)
					{
						inst.X_ic[s][i][c] = inst.Heur_X_ic[s][i][c];
					}
				}

				if (inst.sol_val_HEUR_DOWN[s] > 0)
				{
					inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.sol_val_HEUR_DOWN[s]) + 2;
					inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.sol_val_HEUR_DOWN[s]) - 2;
				}
				if (inst.sol_val_HEUR_DOWN[s] < 0)
				{
					inst.sol_val_UP[s] = (int)round((1 - 0.2) * inst.sol_val_HEUR_DOWN[s]) + 2;
					inst.sol_val_DOWN[s] = (int)round((1 + 0.2) * inst.sol_val_HEUR_DOWN[s]) - 2;
				}
				if (inst.sol_val_HEUR_DOWN[s] == 0)
				{
					inst.sol_val_UP[s] = 2;
					inst.sol_val_DOWN[s] = -2;
				}
			}

			int Best_sol_so_fat = -pow(inst.n_individuals, 3);
			int BestIndex = -1;

			for (int s = 0; s < inst.n_sol; s++)
			{
				if (Best_sol_so_fat <= inst.sol_val_HEUR_DOWN[s])
				{
					Best_sol_so_fat = inst.sol_val_HEUR_DOWN[s];
					BestIndex = s;
				}
			}

			int in_the_meantime = 0;
			for (int i = 0; i < inst.n_individuals; i++)
			{
				for (int c = 0; c < inst.n_groups; ++c)
				{
					in_the_meantime = inst.X_ic[0][i][c];
					inst.X_ic[0][i][c] = inst.X_ic[BestIndex][i][c];
					inst.X_ic[BestIndex][i][c] = in_the_meantime;
				}
			}

		}

	}

	cout << "Number of groups:" << inst.n_groups << endl << endl;

	////////////////////////////////////////////////////////////////////////////////////////

	if (inst.Heuristic_only <= 0) {

		if (inst.method <= 0) // B&C
		{
			if (inst.problem <= 0) //MaxSum
			{
				build_MILP_min_sum(&inst);
				solve_MILP_min_sum(&inst);

				if (reoptimize)
				{
					clean_MILP_min_sum(&inst);

					double mytime = inst.solution_time;

					inst.n_sol = 1;
					if (inst.MagnantiW == 2) {
						inst.MagnantiW = 1;
					}

					for (int s = 0; s < inst.n_sol; s++)
					{
						inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.objval) + 2;
						inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.objval) - 2;
			
					}

					build_MILP_min_sum(&inst);
					solve_MILP_min_sum(&inst);

					inst.solution_time = mytime + inst.solution_time;
				}

			}

			if (inst.problem >= 1) //MaxMin
			{
				build_MILP_min_max(&inst);
				solve_MILP_min_max(&inst);

				if (reoptimize)
				{
					clean_MILP_min_max(&inst);

					double mytime = inst.solution_time;

					inst.n_sol = 1;
					if (inst.MagnantiW == 2) {
						inst.MagnantiW = 1;
					}

					for (int s = 0; s < inst.n_sol; s++)
					{
						inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.objval) + 2;
						inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.objval) - 2;

					}

					build_MILP_min_max(&inst);
					solve_MILP_min_max(&inst);

					inst.solution_time = mytime + inst.solution_time;
				}
			}

		}

		if (inst.method >= 1) // Benders
		{

			if (inst.problem <= 0) //MaxSum
			{
				build_BEN_min_sum(&inst);
				solve_BEN_min_sum(&inst);

				if (reoptimize)
				{
					clean_BEN_min_sum(&inst);

					inst.n_sol = 1;
					if (inst.MagnantiW == 2) {
						inst.MagnantiW = 1;
					}

					double mytime = inst.solution_time;

					for (int s = 0; s < inst.n_sol; s++)
					{
						inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.objval) + 1;
						inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.objval) - 1;

					}

					build_BEN_min_sum(&inst);
					solve_BEN_min_sum(&inst);

					inst.solution_time = mytime + inst.solution_time;
				}
			}

			if (inst.problem >= 1) //MaxMin
			{
				build_BEN_min_max(&inst);
				solve_BEN_min_max(&inst);

				if (reoptimize)
				{
					clean_BEN_min_max(&inst);
					double mytime = inst.solution_time;

					inst.n_sol = 1;
					if (inst.MagnantiW == 2) {
						inst.MagnantiW = 1;
					}

					for (int s = 0; s < inst.n_sol; s++)
					{
						inst.sol_val_UP[s] = (int)round((1 + 0.2) * inst.objval) + 1;
						inst.sol_val_DOWN[s] = (int)round((1 - 0.2) * inst.objval) - 1;
					}

					build_BEN_min_max(&inst);
					solve_BEN_min_max(&inst);

					inst.solution_time = mytime + inst.solution_time;
				}

			}

		}

	}

	///////////////////////////////////////////////////////////////////////////////////

	ofstream compact_file_aggregate;
	compact_file_aggregate.open(inst.output_file_data, ios::app);

	compact_file_aggregate << fixed
		<< inst.x_file_data << "\t"
		<< inst.n_individuals << "\t"
		<< inst.n_features << "\t"
		<< inst.n_groups << "\t"
		<< inst.method << "\t"
		<< inst.problem << "\t"
		<< inst.Heuristic_only << "\t"
		<< inst.LB_heur << "\t"
		<< inst.UB_heur << "\t"
		<< inst.MagnantiW << "\t"
		<< inst.external_sol << "\t"
		<< inst.warm_start << "\t"
		<< inst.get_feasible << "\t"
		<< inst.number_of_CPU << "\t"
		<< inst.Heu_sol << "\t";

		// Loop through the arrays and print each element
		for (size_t s = 0; s < inst.Heu_sol; ++s) {
			compact_file_aggregate << inst.sol_val_EXTERNAL[s] << "\t";
		}
		for (size_t s = 0; s < inst.Heu_sol; ++s) {
			compact_file_aggregate << inst.sol_val_HEUR_DOWN[s] << "\t";
		}
		// Print the remaining values
		compact_file_aggregate
		<< inst.objval << "\t"
		<< inst.bestobjval << "\t"
		<< inst.solution_time_HEUR << "\t"
		<< inst.solution_time << "\t"
		<< inst.nodecount << "\t"
		<< inst.TOLL_OPTIMALITY << "\t"
		<< inst.iteration_BEN << "\t"
		<< endl;
	compact_file_aggregate.close();

	inst.status = CPXfreeprob(inst.env_MILP, &(inst.lp_MILP));
	if (inst.status != 0) { printf("error in CPXfreeprob\n"); exit(-1); }

	inst.status = CPXcloseCPLEX(&(inst.env_MILP));
	if (inst.status != 0) { printf("error in CPXcloseCPLEX\n"); exit(-1); }

	free(inst.x);

	////////////////////////////////////////////////////////////////////////////////////////
		
	cout << endl;

	printf("\nDONE!\n\n");

	return 1;

}

