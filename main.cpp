#include "GRBInterface.h"
#include "KGraph.h"
#include <sstream>
#include <string>
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




int main(int argc, char *argv[])
{

	//change the precision for outputting the running time. I want it to be rounded to 2 decimal places.
	cout.setf(ios::fixed);
	cout.precision(2);

	time_t start_time = clock();
	if (argc<2)
		cerr << "ERROR: Not enough arguments.";


	else if (strcmp(argv[1], "FindNonCriticalNodes") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		time_t Preprocessing_start = clock();
		vector<long> set_I = Preprocessing(g);
		cout << "Preprocessing time is " << (double)(clock() - Preprocessing_start) / CLOCKS_PER_SEC << " " << endl;
		cerr << "Number of simplicial vertices fixed = " << set_I.size() << endl;
		//for (long i = 0; i < set_I.size(); i++) cerr << set_I[i] << " ";
	
		long leaves = 0;
		for (long i = 0; i<g.n; i++)
			if (g.degree[i] == 1)
			{
				long v = g.adj[i][0];
				if (g.degree[v] > 1 || i < v) leaves++;
			}
		cerr << "Number of leaves fixed = " << leaves << endl;
	}


	else if (strcmp(argv[1], "Heuristic") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		cout << "Heuristic time is " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		cerr << "nodes chosen by heuristic are ";
		for (long i = 0; i < Heuristic_sol.size(); i++) cerr << Heuristic_sol[i] << " ";
		cerr << endl;
		cerr << "heuristic solution is " << obj(g, Heuristic_sol, k);
	}

	/*To solve DCNP when distances are measured in terms of hops
	* Thin formulation is used.
	* Integer separation is used.*/
	else if (strcmp(argv[1], "Thin") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();


		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		//vector<long> Heuristic_sol;

		cerr << "Heuristic time is " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;


		vector<long> criticalNodes = solveDCNP_thin_formulation(g, k, b, Heuristic_sol);

		/*cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << criticalNodes.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < criticalNodes.size(); i++)
		cout << criticalNodes[i] << " ";*/

		//Use the following for batch files:
		cout << g.name << " " << k << " " << b << " ";
		cout << Heuristic_time << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		for (long i = 0; i < criticalNodes.size(); i++)
		{
			cout << criticalNodes[i] << " ";
		}
		cout << "\n";
	}


	/*To solve DCNP in edge-weighted graphs.
	* Thin formulation is used.
	* Integer separation is used*/
	else if (strcmp(argv[1], "Thin_weighted") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		//vector<long> Heuristic_sol = Greedy_Heuristic(g, k, b);
		vector<long> Heuristic_sol;

		cerr << "Heuristic time is " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;


		vector<long> criticalNodes = solveDCNP_thin_formulation_weighted(g, k, b, Heuristic_sol);

		cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << criticalNodes.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < criticalNodes.size(); i++)
		cout << criticalNodes[i] << " ";

		//Use the following for batch files:
		/*cout << g.name << " " << k << " " << b << " ";
		cout << Heuristic_time << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		for (long i = 0; i < criticalNodes.size(); i++)
		{
			cout << criticalNodes[i] << " ";
		}
		cout << "\n";*/
	}


	/*To solve DCNP with Thin formulation when distances 
	  are measured in terms of hops.
	  fractional separation is used. */
	else if (strcmp(argv[1], "Thin_fractional") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		cerr << "Heuristic time is " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;

		//vector<long> Heuristic_sol;

		vector<long> criticalNodes = solveDCNP_thin_formulation_fractional(g, k, b, Heuristic_sol);

		cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << criticalNodes.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < criticalNodes.size(); i++)
		cout << criticalNodes[i] << " ";

		//Use the following for batch files:
		/*cout << g.name << " " << k << " " << b << " ";
		cout << Heuristic_time << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		for (long i = 0; i < criticalNodes.size(); i++)
		{
			cout << criticalNodes[i] << " ";
		}
		cout << "\n";*/
	}

	/*To solve DCNP with Path-like formulation when distances 
	  are measured in terms of hops and k=3. */
	else if (strcmp(argv[1], "Path_like_k3") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long b = atol(argv[4]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, 3, b);
		cerr << "Heuristic time is " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		//vector<long> Heuristic_sol;

		vector<long> criticalNodes = solveDCNP_path_like_k3(g, b, Heuristic_sol);

		/*cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes deleted equal " << criticalNodes.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < criticalNodes.size(); i++)
		cout << criticalNodes[i] << " ";*/

		//Use the following for batch files:
		cout << g.name << " " << 3 << " " << b << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		cout << Heuristic_time << " ";
		for (long i = 0; i < criticalNodes.size(); i++)
		{
			cout << criticalNodes[i] << " ";
		}
		cout << "\n";
	}


	/*To solve DCNP when distances are measured in terms of hops
	* Recursive formulation is used. */
	else if (strcmp(argv[1], "DCNP_Veremyev") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		cerr << "Heuristic time is " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		//vector<long> Heuristic_sol;

		vector<long> criticalNodes = solveDCNP_Veremyev(g, k, b, Heuristic_sol);

		cout << "time in sec: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cout << "# nodes Deleted equal " << criticalNodes.size() << endl;
		cout << "** Label of deleted Node(s) is(are) ";
		for (long i = 0; i < criticalNodes.size(); i++)
		{
			cout << criticalNodes[i] << " ";
		}

		//Use the following for batch files:
		/*cout << g.name << " " << k << " " << b << " ";
		cout << (double)(clock() - start) / CLOCKS_PER_SEC << " ";
		cout << Heuristic_time << " ";
		for (long i = 0; i < criticalNodes.size(); i++)
		{
		cout << criticalNodes[i] << " ";
		}
		cout << "\n";*/
	}

	else
	{
		cout << "ERROR: Your command is not valid.";
	}
	return EXIT_SUCCESS;
}
