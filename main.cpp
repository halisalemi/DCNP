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
#include <unordered_map>
#include "GRBInterface.h"
#include "KGraph.h"
#include "ConnectorEnumeration.h"




int main(int argc, char *argv[])
{

	//Change the precision for outputting the running time. It is rounded to 2 decimal places.
	cout.setf(ios::fixed);
	cout.precision(2);
	cerr.setf(ios::fixed);
	cerr.precision(2);


	time_t start_time = clock();
	if (argc<2)
		cerr << "ERROR: Not enough arguments.";


	//Variable fixing to find a set of non-critical nodes in unweighted graphs.
	else if (strcmp(argv[1], "FindNonCriticalNodes") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		time_t Variable_fixing_start = clock();
		vector<long> set_I = FindNonCriticalNodes(g);
		cout << "Variable fixing time = " << (double)(clock() - Variable_fixing_start) / CLOCKS_PER_SEC << " " << endl;
	}


	//DCNP Heuristic for graphs with hop-based distances
	else if (strcmp(argv[1], "DCNPHeuristic") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		cout << "Heuristic time = " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		cerr << "Nodes chosen by heuristic = ";
		for (long i = 0; i < Heuristic_sol.size(); i++) cerr << Heuristic_sol[i] << " ";
		cerr << endl;
		cerr << "Heuristic solution = " << obj(g, Heuristic_sol, k) << endl;
	}


	//DCNP Heuristic for graphs with edge-weighted distances
	else if (strcmp(argv[1], "DCNPHeuristicWeighted") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic_Weighted(g, k, b);
		cout << "Heuristic time = " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		cerr << "Nodes chosen by heuristic = ";
		for (long i = 0; i < Heuristic_sol.size(); i++) cerr << Heuristic_sol[i] << " ";
		cerr << endl;
		cerr << "Heuristic solution = " << obj_weighted(g, Heuristic_sol, k) << endl;
	}


	/*To solve DCNP for graphs with hop-based distances.
	* Thin formulation with integer separation*/
	else if (strcmp(argv[1], "Thin_I") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		cout << g.name << " Thin_I " << k << " " << b << " ";
		time_t start = clock();

		cerr << "Finding heuristic solution" << endl;
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;

		bool subopt;
		vector<bool> Initial(5, false);
		Initial[1] = true; //minimal length-exactly-1 i,j-connectors are added initially.
		Initial[2] = true; //minimal length-exactly-2 i,j-connectors are added initially.
		Initial[3] = true; //minimal length-exactly-3 i,j-connectors are added initially (when k >= 3).
		Initial[4] = true; //minimal length-exactly-4 i,j-connectors are added initially (when k >= 4).


		vector<long> criticalNodes = Thin_I(g, k, b, Heuristic_sol, subopt, Initial);
		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << " ";

		if (subopt) cout << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "# critical nodes = " << criticalNodes.size() << endl;
			cerr << "** Label of critical node(s) = ";
			PrintVectorLong(criticalNodes);
			cout << "# close vertex pairs in G-D = " << obj(g, criticalNodes, k) << endl;
		}
	}

	/*To solve DCNP for graphs with hop-based distances.
	* Thin formulation with fractional separation*/
	else if (strcmp(argv[1], "Thin_F") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		cout << g.name << " Thin_F " << k << " " << b << " ";
		time_t start = clock();
		
		cerr << "Finding heuristic solution" << endl;
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;

		bool subopt;
		vector<bool> Initial(5, false);
		Initial[1] = false; //minimal length-exactly-1 i,i-connectors are added initially.
		Initial[2] = true; //minimal length-exactly-2 i,i-connectors are added initially.
		Initial[3] = true; //minimal length-exactly-3 i,i-connectors are added initially (when k>=3).
		Initial[4] = false; //minimal length-exactly-4 i,i-connectors are added initially (when k>=4).
		vector<long> criticalNodes = Thin_F(g, k, b, Heuristic_sol, subopt, Initial);
		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << " ";

		if(subopt) cout << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "# critical nodes = " << criticalNodes.size() << endl;
			cerr << "** Label of critical node(s) =  ";
			PrintVectorLong(criticalNodes);
			cout << "# close vertex pairs in G-D = " << obj(g, criticalNodes, k) << endl;
		}
	}


	/*To solve DCNP for graphs with hop-based distances and k=3 
	 * Path-like formulation*/
	else if (strcmp(argv[1], "Path_like_k3") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = 3;
		long b = atol(argv[4]); //budget 
		cout << g.name << " Path_like_k3 3 " << b << " ";
		time_t start = clock();

		cerr << "Finding heuristic solution" << endl;
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;
		
		bool subopt;
		vector<long> criticalNodes = Path_like_k3(g, b, Heuristic_sol, subopt);
		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << " ";

		if (subopt) cout << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "# critical nodes = " << criticalNodes.size() << endl;
			cerr << "** Label of critical node(s) =  ";
			PrintVectorLong(criticalNodes);
			cout << "# close vertex pairs in G-D = " << obj(g, criticalNodes, k) << endl;
		}
	}

	/*To solve DCNP for graphs with hop-based distances and k=4
	* Path-like formulation*/
	else if (strcmp(argv[1], "Path_like_k4") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = 4;
		long b = atol(argv[4]); //budget 
		cout << g.name << " Path_like_k 4 " << b << " ";
		time_t start = clock();

		cerr << "Finding heuristic solution" << endl;
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;

		bool subopt;
		vector<long> criticalNodes = Path_like_k4(g, b, Heuristic_sol, subopt);
		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << " ";

		if (subopt) cout << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "# critical nodes = " << criticalNodes.size() << endl;
			cerr << "** Label of critical node(s) =  ";
			PrintVectorLong(criticalNodes);
			cout << "# close vertex pairs in G-D = " << obj(g, criticalNodes, k) << endl;
		}
	}


	/*To solve DCNP for graphs with hop-based distances
	* Recursive formulation*/
	else if (strcmp(argv[1], "Recursive") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget
		cout << g.name << " Recursive " << k << " " << b << " ";
		time_t start = clock();

		cerr << "Finding heuristic solution" << endl;
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;
		
		bool subopt;
		vector<long> criticalNodes = Recursive(g, k, b, Heuristic_sol, subopt);
		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << " ";

		if (subopt) cout << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "# critical nodes = " << criticalNodes.size() << endl;
			cerr << "** Label of critical node(s) =  ";
			PrintVectorLong(criticalNodes);
			cout << "# close vertex pairs in G-D = " << obj(g, criticalNodes, k) << endl;
		}
	}


	/*To solve DCNP for graphs with edge-weighted distances.
	* Thin formulation with integer separation*/
	else if (strcmp(argv[1], "Thin_weighted") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		cout << g.name << " " << k << " " << b << " ";
		time_t start = clock();

		cerr << "Finding heuristic solution" << endl;
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic_Weighted(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cout << "Heuristic time = " << Heuristic_time << " ";

		bool subopt;
		vector<long> criticalNodes = Thin_Weighted(g, k, b, Heuristic_sol, subopt);
		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << " ";

		if (subopt) cout << "Model is infeasible or other problem." << endl;
		else
		{
			cerr << "# critical nodes = " << criticalNodes.size() << endl;
			cerr << "** Label of critical node(s) = ";
			PrintVectorLong(criticalNodes);
			cout << "# close vertex pairs in G-D = " << obj_weighted(g, criticalNodes, k) << endl;
		}
	}

	else
	{
		cout << "ERROR: Your command is not valid." << endl;
	}
	return EXIT_SUCCESS;
}
