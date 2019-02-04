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


	else if (strcmp(argv[1], "Heuristic") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);

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
		cout << "heur is " << 2 * num_close_vertices << endl;
		cout << "Time that spent in Heuristic in sec is: " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
	}



	else if (strcmp(argv[1], "Thin") == 0)
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
		vector<long> Deleted = solveDCNP_thin_formulation(g, k, b, Heuristic_sol, subOpt);

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

	else
	{
		cout << "ERROR: Your command is not valid.";
	}
	return EXIT_SUCCESS;
}
