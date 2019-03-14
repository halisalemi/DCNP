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



	time_t start_time = clock();
	if (argc<2)
		cerr << "ERROR: Not enough arguments.";


	/*Preprocessing procedure to find the non-critical nodes
	/*It works for arbitrary k and b and also for edge-weighted graphs*/
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


	/*To solve DCNP when distances are measured in terms of hops
	* Thin formulation with integer separation is used. */
	else if (strcmp(argv[1], "Thin") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;
		
		vector<long> criticalNodes = solveDCNP_thin_formulation(g, k, b, Heuristic_sol);

		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cerr << "# critical nodes = " << criticalNodes.size() << endl;
		cerr << "** Label of critical node(s) = ";
		PrintVectorLong(criticalNodes);
		cerr << "# close vertex pairs in G-D = " << obj(g, criticalNodes, k) << endl;
	}


	/*To solve DCNP in edge-weighted graphs.
	* Thin formulation with integer separation is used.*/
	else if (strcmp(argv[1], "Thin_weighted") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();
		cerr << "Finding heuristic solution" << endl;
		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;
		
		vector<long> criticalNodes = solveDCNP_thin_formulation_weighted(g, k, b, Heuristic_sol);

		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cerr << "# critical nodes = " << criticalNodes.size() << endl;
		cerr << "** Label of critical node(s) = ";
		PrintVectorLong(criticalNodes);
	}


	/*To solve DCNP when distances are measured in terms of hops
	* Thin formulation with fractional separation is used. */
	else if (strcmp(argv[1], "Thin_fractional") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long k = atol(argv[4]); //short distance
		long b = atol(argv[5]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, k, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;
		

		vector<long> criticalNodes = solveDCNP_thin_formulation_fractional(g, k, b, Heuristic_sol);

		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cerr << "# critical nodes = " << criticalNodes.size() << endl;
		cerr << "** Label of critical node(s) =  ";
		PrintVectorLong(criticalNodes);
		cerr << "# close vertex pairs in G-D = " << obj(g, criticalNodes, k) << endl;
	}


	/*To solve DCNP when distances are measured in terms of hops and k=3 
	 * Path-like formulation is used */
	else if (strcmp(argv[1], "Path_like_k3") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);

		long b = atol(argv[4]); //budget 
		time_t start = clock();

		time_t Heuristic_start = clock();
		vector<long> Heuristic_sol = DCNP_Heuristic(g, 3, b);
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cout << "Heuristic time = " << Heuristic_time << endl;
		

		vector<long> criticalNodes = solveDCNP_path_like_k3(g, b, Heuristic_sol);

		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cerr << "# critical nodes = " << criticalNodes.size() << endl;
		cerr << "** Label of critical node(s) =  ";
		PrintVectorLong(criticalNodes);
		cerr << "# close vertex pairs in G-D = " << obj(g, criticalNodes, 3) << endl;
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
		double Heuristic_time = (double)(clock() - Heuristic_start) / CLOCKS_PER_SEC;
		cerr << "Heuristic time = " << Heuristic_time << endl;
		
		vector<long> criticalNodes = solveDCNP_Veremyev(g, k, b, Heuristic_sol);

		cout << "Total time = " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
		cerr << "# critical nodes = " << criticalNodes.size() << endl;
		cerr << "** Label of critical node(s) =  ";
		PrintVectorLong(criticalNodes);
		cerr << "# close vertex pairs in G-D = " << obj(g, criticalNodes, k) << endl;
	}

	else if (strcmp(argv[1], "Diameter") == 0)
	{
		KGraph g(argv[3], argv[3], argv[2]);
		long diameter = g.DiameterWeighted();
		cerr << "diameter = " << diameter << endl;
		vector< vector< long> > clusters;
		vector<long> degreeZero;
		g.FindConnectedComponents(clusters, degreeZero);
		cerr << "# isolated vertices " << degreeZero.size() << endl;
		for (long i = 0; i < degreeZero.size(); i++)
			cerr << degreeZero[i] << endl;
		cerr << "# components " << clusters.size() << endl;
	}

	else
	{
		cout << "ERROR: Your command is not valid." << endl;
	}
	return EXIT_SUCCESS;
}
