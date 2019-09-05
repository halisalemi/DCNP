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

	//change the precision for outputting the running time. It is rounded to 2 decimal places.
	cout.setf(ios::fixed);
	cout.precision(2);
	cerr.setf(ios::fixed);
	cerr.precision(2);


	time_t start_time = clock();
	if (argc<2)
		cerr << "ERROR: Not enough arguments.";


	/*Preprocessing procedure to find the non-critical nodes
	/*It works for arbitrary k and b*/
	else if (strcmp(argv[1], "FindNonCriticalNodes") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		time_t Preprocessing_start = clock();
		vector<long> set_I = FindNonCriticalNodes(g);
		cout << "Preprocessing time = " << (double)(clock() - Preprocessing_start) / CLOCKS_PER_SEC << " " << endl;
	}


	//DCNP Heuristic to find distance-based critical nodes
	else if (strcmp(argv[1], "DCNPHeuristic") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		cout << "Heuristic time = " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		cerr << "Nodes chosen by heuristic = ";
		for (long i = 0; i < Heuristic_sol.size(); i++) cerr << Heuristic_sol[i] << " ";
		cerr << endl;
		cerr << "Heuristic solution = " << obj(g, Heuristic_sol, k) << endl;
	}

	else if (strcmp(argv[1], "DCNPHeuristicWeighted") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic_Weighted(g, k, b);
		cout << "Heuristic time = " << (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC << endl;
		cerr << "Nodes chosen by heuristic = ";
		for (long i = 0; i < Heuristic_sol.size(); i++) cerr << Heuristic_sol[i] << " ";
		cerr << endl;
		cerr << "Heuristic solution = " << obj_weighted(g, Heuristic_sol, k) << endl;
	}


	/*To solve DCNP for graphs with hop-based distances.
	* Thin formulation with integer separation is used. */
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
		Initial[1] = true; //minimal length-exactly-1 i,i-connectors are added initially
		Initial[2] = true; //minimal length-exactly-2 i,i-connectors are added initially
		Initial[3] = true; //minimal length-exactly-3 i,i-connectors are added initially (when k>=3)
		Initial[4] = false; //minimal length-exactly-4 i,i-connectors are added initially (when k>=2)

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
	* Thin formulation with fractional separation is used. */
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
		Initial[1] = true; //minimal length-exactly-1 i,i-connectors are added initially
		Initial[2] = true; //minimal length-exactly-2 i,i-connectors are added initially
		Initial[3] = true; //minimal length-exactly-3 i,i-connectors are added initially (when k>=3)
		Initial[4] = false; //minimal length-exactly-4 i,i-connectors are added initially (when k>=2)
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


	/*To solve DCNP when distances are measured in terms of hops and k=3. 
	 * Path-like formulation is used */
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

	/*To solve DCNP when distances are measured in terms of hops and k=4.
	* Path-like formulation is used */
	else if (strcmp(argv[1], "Path_like_k4") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long k = 4;
		long b = atol(argv[4]); //budget 
		cout << g.name << " Path_like_k3 3 " << b << " ";
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


	/*To solve DCNP when distances are measured in terms of hops.
	* Recursive formulation is used. */
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
	* Thin formulation with integer separation is used.*/
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


	else if (strcmp(argv[1], "FindPercentiles") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		double alpha = atol(argv[4]); //We'll consider alpha^th percentile
		vector<vector<long>> DistanceMatrix;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsWeighted(i);
			DistanceMatrix.push_back(dist_from_i);
		}
		vector<long> AllDist;
		for (long i = 0; i < g.n; i++)
		{
			for (long j = i + 1; j < g.n; j++)
			{
				AllDist.push_back(DistanceMatrix[i][j]);
			}
		}
		sort(AllDist.begin(), AllDist.end());
		int index = int(double(alpha / 100)*AllDist.size());
		cerr << "k = " << AllDist[index];
	}

	
	else
	{
		cout << "ERROR: Your command is not valid." << endl;
	}
	return EXIT_SUCCESS;
}
