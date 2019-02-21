#include "GRBInterface.h"
#include "KGraph.h"
#include "Clique.h"
#include "wkplex.h"
//#include "wclique.h"
#include <sstream>
#include <string>
//#include <ilcplex/cplex.h>
#include <ctime>
#include <iomanip> 
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include "ConnectorEnumeration.h"
#include <unordered_map>



template <class T>
inline std::string to_string(const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

using namespace std;
long* g_degree;
int seed;
int u = 0;
KGraph G;


int main(int argc, char *argv[])
{

	// change the precision for outputting the running time. I want it to be rounded to 2 decimal places.
	cout.setf(ios::fixed);
	cout.precision(2);

	time_t start_time = clock();
	if (argc<2)
		cerr << "ERROR: Not enough arguments.";


	else if (strcmp(argv[1], "Preprocessing") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		time_t Preprocessing_start = clock();
		vector<long> solution = Preprocessing(g);
		cout << "Preprocessing time is " << (double)(clock() - Preprocessing_start) / CLOCKS_PER_SEC << " " << endl;
		cerr << "Size of I is " << solution.size() << endl;
		for (long i = 0; i < solution.size(); i++) cerr << solution[i] << " ";
	}


	else if (strcmp(argv[1], "Heuristic") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);

		cerr << "Number of nodes are " << Heuristic_sol.size() << endl;
		cerr << "node are ";
		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			cerr << Heuristic_sol[i] << " ";
		}
		cerr << "\n";

		vector<long> close_vertices;
		long num_close_vertices = 0;

		vector<bool> new_nodes(g.n, true);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			new_nodes[Heuristic_sol[i]] = false;
		}

		for (long v = 0; v < g.n; v++)
		{
			vector <long> dist_from_v = g.ShortestPathsUnweighted(v, new_nodes);
			for (long w = v + 1; w < g.n; w++)
			{
				if (dist_from_v[w] <= k)
				{
					num_close_vertices++;
				}
			}
		}
		cout << "heuristic solution is " << 2 * num_close_vertices << endl;
		cout << "Time that spent in Heuristic in sec is: " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
	}


	else if (strcmp(argv[1], "Thin") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();


		time_t Heuristic_start = clock();
		//vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);
		vector<long> Heuristic_sol;

		cerr << "Time that spent in Heuristic in sec is: " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;


		bool subOpt;
		vector<long> Deleted = solveDCNP_thin_formulation(g, k, b, Heuristic_sol, subOpt);

		/*cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << Deleted.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < Deleted.size(); i++)
		cout << Deleted[i] << " ";*/

		//Use the following for batch files:
		cout << g.name << " " << k << " " << b << " ";
		cout << Heuristic_time << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		for (long i = 0; i < Deleted.size(); i++)
		{
			cout << Deleted[i] << " ";
		}
		cout << "\n";
	}

	else if (strcmp(argv[1], "Thin_weighted") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		//vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);
		vector<long> Heuristic_sol;

		cerr << "Time that spent in Heuristic in sec is: " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;


		bool subOpt;
		vector<long> Deleted = solveDCNP_thin_formulation_weighted(g, k, b, Heuristic_sol, subOpt);

		cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << Deleted.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < Deleted.size(); i++)
		cout << Deleted[i] << " ";

		//Use the following for batch files:
		/*cout << g.name << " " << k << " " << b << " ";
		cout << Heuristic_time << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		for (long i = 0; i < Deleted.size(); i++)
		{
			cout << Deleted[i] << " ";
		}
		cout << "\n";*/
	}

	else if (strcmp(argv[1], "Thin_fractional") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);
		cerr << "Time that spent in Heuristic in sec is: " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;

		//vector<long> Heuristic_sol;

		bool subOpt;
		vector<long> Deleted = solveDCNP_thin_formulation_fractional(g, k, b, Heuristic_sol, subOpt);

		cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << Deleted.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < Deleted.size(); i++)
		cout << Deleted[i] << " ";

		//Use the following for batch files:
		//cout << g.name << " " << k << " " << b << " ";
		//cout << Heuristic_time << " ";
		//cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		///*for (long i = 0; i < Deleted.size(); i++)
		//{
		//	cout << Deleted[i] << " ";
		//}*/
		//cout << "\n";
	}


	else if (strcmp(argv[1], "Path_like") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);
		cerr << "Time that spent in Heuristic in sec is: " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		//vector<long> Heuristic_sol;

		bool subOpt;
		vector<long> Deleted = solveDCNP_path_like(g, k, b, Heuristic_sol, subOpt);

		/*cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << Deleted.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < Deleted.size(); i++)
		cout << Deleted[i] << " ";*/

		//Use the following for batch files:
		cout << g.name << " " << k << " " << b << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		cout << Heuristic_time << " ";
		for (long i = 0; i < Deleted.size(); i++)
		{
			cout << Deleted[i] << " ";
		}
		cout << "\n";
	}


	else if (strcmp(argv[1], "DCNP_Veremyev") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);
		cerr << "Time spent in Heuristic in sec is: " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		//vector<long> Heuristic_sol;

		bool subOpt;
		vector<long> Deleted = solveDCNP_Veremyev(g, k, b, Heuristic_sol, subOpt);

		cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes Deleted equal " << Deleted.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < Deleted.size(); i++)
		{
			cout << Deleted[i] << " ";
		}

		//Use the following for batch files:
		/*cout << g.name << " " << k << " " << b << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		cout << Heuristic_time << " ";
		for (long i = 0; i < Deleted.size(); i++)
		{
		cout << Deleted[i] << " ";
		}
		cout << "\n";*/
	}

	else if (strcmp(argv[1], "CNP") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long b = atol(argv[4]); //budget 
		time_t start = clock();

		long k = g.n - 1;


		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);


		cerr << "Time that spent in Heuristic in sec is: " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;


		bool subOpt;
		vector<long> Deleted = solveDCNP_thin_formulation_fractional(g, k, b, Heuristic_sol, subOpt);

		/*cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << Deleted.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < Deleted.size(); i++)
		cout << Deleted[i] << " ";*/

		//Use the following for batch files:
		cout << g.name << " " << k << " " << b << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		cout << Heuristic_time << " ";
		for (long i = 0; i < Deleted.size(); i++)
		{
			cout << Deleted[i] << " ";
		}
		cout << "\n";
	}

	//to have distance matrix in weighted graph G
	else if (strcmp(argv[1], "distance_matrix_G") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		vector < vector<long>> all_dist;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i_to = g.ShortestPathsWeighted(i);
			all_dist.push_back(dist_from_i_to);
		}
		for (long i = 0; i < g.n; i++)
		{
			for (long j = i + 1; j < g.n; j++)
			{
				cerr << "distance from vertex " << i << " to vertex " << j << " is " << all_dist[i][j] << endl;
			}
		}
	}

	//to have distance matrix in weighted graph G[S]
	else if (strcmp(argv[1], "distance_matrix_GS") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		vector < vector<long>> all_dist;
		vector<bool> S(g.n, true);
		/*S[1] = false;
		S[5] = false;*/


		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i_to = g.ShortestPathsWeighted(i, S);
			all_dist.push_back(dist_from_i_to);
		}
	
		for (long i = 0; i < g.n; i++)
		{
			for (long j = i + 1; j < g.n; j++)
			{
				cerr << "distance from vertex " << i << " to vertex " << j << " is " << all_dist[i][j] << endl;
			}
		}
	}

	
	else if (strcmp(argv[1], "BH_dijkstra") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		vector<bool>S(g.n, true);
		S[1] = false;
		S[2] = false;


		vector<long> path;
		g.BinaryHeapDijkstra(0, 5, S);
	}


	else
	{
		cout << "ERROR: Your command is not valid.";
	}
	return EXIT_SUCCESS;
}
