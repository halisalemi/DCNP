#include "GRBInterface.h"
#include <sstream>
#include <ctime>
#include <iostream>
#include <vector>
#include <queue>
#include <list>
#include <algorithm>
#include <limits>
#include <stack>
#include "ConnectorEnumeration.h"
#include "Binary_Heap.h"
#include <unordered_map>
#include <random>


using namespace std;

vector<bool> boolify(vector<long> &S, long n)
{
	vector<bool> Sbool(n, false);
	for (long i = 0; i < S.size(); i++) Sbool[S[i]] = true;
	return Sbool;
}

vector<bool> boolify(double *y, long n)
{
	vector<bool> Sbool(n, false);
	for (long i = 0; i < n; i++) if (y[i]>0.5) Sbool[i] = true;
	return Sbool;
}

void PrintVectorLong(vector<long>&S)
{
	for (long i = 0; i < S.size(); i++) cerr << S[i] << " ";
	cerr << "\n";
}


long IntegerSeparation::numCallbacks = 0;
double IntegerSeparation::TotalCallbackTime = 0;
long IntegerSeparation::numLazyCutsInteger = 0;


long IntegerSeparationWeighted::numCallbacks = 0;
double IntegerSeparationWeighted::TotalCallbackTime = 0;
long IntegerSeparationWeighted::numLazyCutsInteger = 0;


long FractionalSeparation::numCallbacks = 0;
double FractionalSeparation::TotalCallbackTimeInteger = 0;
double FractionalSeparation::TotalCallbackTimeFractional = 0;
long FractionalSeparation::numLazyCutsInteger = 0;
long FractionalSeparation::numLazyCutsFractional = 0;

double ViolationThreshold = 0.05;

void IntegerSeparation::callback()
{
	try
	{
		time_t start = clock();
		if (where == GRB_CB_MIPSOL)
		{
			numCallbacks++;

			//get the solution vector (for x and y variables) from Gurobi 
			populate_x_MIPSOL();
			populate_y_MIPSOL();

			//now, make it boolean
			vector<bool> NonDeletedVertices(n, false);
			NonDeletedVertices = boolify(y, n);
			NonDeletedVertices.flip();

			
			//find a violated length-k i,j-connector inequality (if any exist)
			vector<long> predecessor;
			vector<long> distance_from_i;
			for (long i = 0; i < n; i++)
			{
				distance_from_i = g->ShortestPathsUnweighted(i, NonDeletedVertices, predecessor);
				for (long counter = 0; counter < gs->degree[i]; counter++)
				{
					long j = gs->adj[i][counter];
					if (j > i && distance_from_i[j] <= k && x[hashing[i*n + j]] < 0.5)
					{
						GRBLinExpr expr = xvars[hashing[i*n + j]] + yvars[i] + yvars[j];

						long q = predecessor[j];
						while (q != i)
						{
							expr += yvars[q];
							q = predecessor[q];
						}
						addLazy(expr >= 1);
						numLazyCutsInteger++;
					}
					predecessor.clear();
				}
			}
		}
		TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}

void FractionalSeparation::callback()
{
	try
	{
		//if we find a new MIP incumbent
		if (where == GRB_CB_MIPSOL)
		{
			time_t start = clock();
			numCallbacks++;

			//get the solution vector for x and y variables from Gurobi
			populate_x_MIPSOL();
			populate_y_MIPSOL();

			vector<bool> NonDeletedVertices(n, false);
			NonDeletedVertices = boolify(y, n);
			NonDeletedVertices.flip();

			vector<long> predecessor;
			vector<long> distance_from_i;
			for (long i = 0; i < n; i++)
			{
				distance_from_i = g->ShortestPathsUnweighted(i, NonDeletedVertices, predecessor);
				for (long counter = 0; counter < gs->degree[i]; counter++)
				{
					long j = gs->adj[i][counter];
					if (j > i && distance_from_i[j] <= k && x[hashing[i*n + j]] < 0.5)
					{
						GRBLinExpr expr = xvars[hashing[i*n + j]] + yvars[i] + yvars[j];
						long q = predecessor[j];
						while (q != i)
						{
							expr += yvars[q];
							q = predecessor[q];
						}
						addLazy(expr >= 1);
						numLazyCutsInteger++;
					}
					predecessor.clear();
				}
			}
			TotalCallbackTimeInteger += (double)(clock() - start) / CLOCKS_PER_SEC;
		}

		//if we explore a MIP node
		else if (where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
		{
			numCallbacks++;
			time_t start = clock();

			//get the solution vector for x and y variables from Gurobi
			populate_x_MIPNODE();
			populate_y_MIPNODE();

			for (long i = 0; i < n; i++)
			{
				d_and_p output = d_and_p_function(*g, *gs, i, k, y);
				for (long counter3 = 0; counter3 < gs->degree[i]; counter3++)
				{
					long j = gs->adj[i][counter3];
					if (j > i) //due to hashing
					{
						if (1 - x[hashing[i*n + j]] - output.d[j][k] > ViolationThreshold)
						{
							GRBLinExpr expr = xvars[hashing[i*n + j]] + yvars[j] + yvars[i];
							long q = output.p[j][k];
							long ss = k - 1;
							for (long counter4 = 2; counter4 <= k; counter4++)
							{
								if (q == i) break;
								expr += yvars[q];
								q = output.p[q][ss];
								ss--;
							}
							addCut(expr >= 1);
							numLazyCutsFractional++;
						}
					}
				}
			}
			TotalCallbackTimeFractional += (double)(clock() - start) / CLOCKS_PER_SEC;
		}
	}

	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}


void IntegerSeparationWeighted::callback()
{
	try
	{
		if (where == GRB_CB_MIPSOL)
		{
			numCallbacks++;
			time_t start = clock();

			populate_x_MIPSOL();
			populate_y_MIPSOL();


			//now, make it boolean
			vector<bool> NonDeletedVertices(n, false);
			NonDeletedVertices = boolify(y, n);
			NonDeletedVertices.flip();


			//find a violated length-k i,j-connector inequality (if any exist)
			vector<long> predecessor;
			vector<long> distance_from_i;
			for (long i = 0; i < n; i++)
			{
				distance_from_i = g->ShortestPathsWeighted(i, NonDeletedVertices, predecessor);
				for (long counter = 0; counter < gs->degree[i]; counter++)
				{
					long j = gs->adj[i][counter];
					if (j > i && distance_from_i[j] <= k && x[hashing[i*n + j]] < 0.5)
					{
						GRBLinExpr expr = xvars[hashing[i*n + j]] + yvars[i] + yvars[j];

						long q = predecessor[j];
						while (q != i)
						{
							expr += yvars[q];
							q = predecessor[q];
						}
						addLazy(expr >= 1);
						numLazyCutsInteger++;
					}
					predecessor.clear();
				}
			}
			TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}



/*The following function uses O(|E^k|) variables.
Here, we use hashing for power graph edges.
Integer separation is performed.*/
vector<long> Thin_I(KGraph &g, long k, long b, vector<long> &Heuristic_sol, bool &subopt, vector<bool> &Initial)
{
	vector<long> Deleted;
	subopt = true;
	cerr << "creating power graph " << endl;
	KGraph gs = g.CreatePowerGraph(k);
	cout << "|E^k| = " << gs.m << " ";
	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3); //use concurrent method to solve LP relaxation.
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		GRBModel model = GRBModel(env);
		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar *X = model.addVars(gs.m, GRB_CONTINUOUS);
		model.update();

		cerr << "Adding objective function." << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < gs.m; i++)
			Objective += X[i];
		model.setObjective(Objective, GRB_MINIMIZE);
		model.update();


		cerr << "Hashing edges" << endl;
		unordered_map<long, long> hash_edges;
		long cur = 0;
		for (long u = 0; u < gs.n; u++)
		{
			for (long i = 0; i < gs.adj[u].size(); i++)
			{
				long v = gs.adj[u][i];
				if (v > u) hash_edges.insert(make_pair(u*gs.n + v, cur++));
			}
		}

		cerr << "Fixing variables" << endl;
		vector<long> NonCriticalNodes = FindNonCriticalNodes(g);
		for (long i = 0; i < NonCriticalNodes.size(); i++)
			Y[NonCriticalNodes[i]].set(GRB_DoubleAttr_UB, 0);

		long len1 = 0;
		if (Initial[1])
		{
			cerr << "Adding minimal length-exactly-1 i,j-connectors" << endl;
			for (long i = 0; i < g.n; i++)
			{
				for (long it = 0; it < g.degree[i]; it++)
				{
					long j = g.adj[i][it];
					if (j > i)
					{
						model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[j] >= 1); //checking i>j is due to hashing
						len1++;
					}
				}
			}
		}

		long len2 = 0;
		if (Initial[2])
		{
			cerr << "Adding minimal length-exactly-2 i,j-connectors" << endl;
			for (long i = 0; i < g.n; i++)
			{
				vector<bool> i_neighbors = boolify(g.adj[i], g.n);
				for (long it1 = 0; it1 < g.degree[i]; it1++)
				{
					long u = g.adj[i][it1];
					for (long it2 = 0; it2 < g.degree[u]; it2++)
					{
						long j = g.adj[u][it2];
						if (i >= j) continue;
						if (i_neighbors[j]) continue;
						model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[j] + Y[u] >= 1);
						len2++;
					}
				}
			}
		}
	
		long len3 = 0;
		if (Initial[3] && k >= 3)
		{
			cerr << "Adding minimal length-exactly-3 i,j-connectors" << endl;
			for (long i = 0; i < g.n; i++)
			{
				vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
				for (long it1 = 0; it1 < g.degree[i]; it1++)
				{
					long u = g.adj[i][it1];
					vector<long> dist_from_u = g.ShortestPathsUnweighted(u);
					for (long it2 = 0; it2 < g.degree[u]; it2++)
					{
						long v = g.adj[u][it2];
						if (dist_from_i[v] == 1) continue;
						for (long it2 = 0; it2 < g.degree[v]; it2++)
						{
							long j = g.adj[v][it2];
							if (dist_from_i[j] == 1) continue;
							if (dist_from_u[j] == 1) continue;
							if (i >= j) continue;
							model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[u] + Y[v] + Y[j] >= 1);
							len3++;
						}
					}
				}
			}
		}

		cerr << "len1 = " << len1 << endl;
		cerr << "len2 = " << len2 << endl;
		cerr << "len3 = " << len3 << endl;


		if (Initial[4] && k >= 4)
		{
			cerr << "Adding minimal length-exactly-4 i,j-connectors" << endl;
			for (long i = 0; i < g.n; i++)
			{
				vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
				for (long j = i + 1; j < g.n; j++)
				{
					if (dist_from_i[j] == 1 || dist_from_i[j] > 4) continue;
					vector<long> dist_from_j = g.ShortestPathsUnweighted(j);
					for (long it1 = 0; it1 < g.degree[i]; it1++)
					{
						long p = g.adj[i][it1];
						vector<bool> p_neighbors = boolify(g.adj[p], g.n);
						if (dist_from_j[p] != 2) continue;
						for (long it2 = 0; it2 < g.degree[p]; it2++)
						{
							long v = g.adj[p][it2];
							if (dist_from_i[v] != 2 || dist_from_j[v] != 2) continue;
							for (long it3 = 0; it3 < g.degree[v]; it3++)
							{
								long q = g.adj[v][it3];
								if (dist_from_i[q] == 2 && dist_from_j[q] == 1 && p_neighbors[q]==false)
									model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[p] + Y[v] + Y[q] + Y[j] >= 1); //type1 V_12, V_22, V_21
								if (dist_from_i[q] == 3 && dist_from_j[q] == 1)
									model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[p] + Y[v] + Y[q] + Y[j] >= 1); //type2 V_12, V_22, V_31
							}
						}
					}
				}
			}
			for (long i = 0; i < g.n; i++) //type3 V_13, V_22, V_21
			{
				vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
				for (long j = i + 1; j < g.n; j++)
				{
					if (dist_from_i[j] == 1 || dist_from_i[j] > 4) continue;
					vector<long> dist_from_j = g.ShortestPathsUnweighted(j);
					for (long it1 = 0; it1 < g.degree[i]; it1++)
					{
						long u = g.adj[i][it1];
						if (dist_from_j[u] != 3) continue;
						for (long it2 = 0; it2 < g.degree[u]; it2++)
						{
							long v = g.adj[u][it2];
							if (dist_from_i[v] != 2 || dist_from_j[v] != 2) continue;
							for (long it3 = 0; it3 < g.degree[v]; it3++)
							{
								long q = g.adj[v][it3];
								if (dist_from_i[q] == 2 && dist_from_j[q] == 1)
									model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[u] + Y[v] + Y[q] + Y[j] >= 1); //type3 V_13, V_22, V_21
								if (dist_from_i[q] == 3 && dist_from_j[q] == 1)
									model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[u] + Y[v] + Y[q] + Y[j] >= 1); //type4 V_13, V_22, V_31
							}
						}
					}
				}
			}
		}


		cerr << "Adding budget constraints" << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
			expr1 += Y[i];
		model.addConstr(expr1 <= b);
		model.update();


		cerr << "Adding lazy constraints" << endl;
		IntegerSeparation cb = IntegerSeparation(X, Y, &g, &gs, hash_edges,k,b);
		model.setCallback(&cb);

		cerr << "Providing Initial Solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}
		vector<bool> NonDeletedNodes = boolify(Heuristic_sol, g.n);
		NonDeletedNodes.flip();

		for (long i = 0; i < g.n; i++)
		{
			vector <long> dist_from_i = g.ShortestPathsUnweighted(i, NonDeletedNodes);
			for (long j = i + 1; j < g.n; j++)
				if (dist_from_i[j] <= k)
					X[hash_edges[i*g.n + j]].set(GRB_DoubleAttr_Start, 1);
		}

		cerr << "Optimizing" << endl;
		model.optimize();

		double bestUB = model.get(GRB_DoubleAttr_ObjVal);
		double bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << "bestLB = " << bestLB << " bestUB = " << bestUB << " ";

		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B nodes = " << NumOfBandBNodes << endl;
		cerr << "# callbacks = " << IntegerSeparation::numCallbacks << endl;
		cerr << "# lazy cuts = " << IntegerSeparation::numLazyCutsInteger << endl;
		cerr << "Time spent in callback = " << IntegerSeparation::TotalCallbackTime << endl;


		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subopt = false;
			for (long i = 0; i < g.n; i++)
				if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
					Deleted.push_back(i);
		}

		delete[] Y;
		delete[] X;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return Deleted;
}


/*To solve DCNP when distances are measured in terms of hops
* Thin formulation with fractional separation is used. */
vector<long> Thin_F(KGraph &g, long k, long b, vector<long> &Heuristic_sol, bool &subopt, vector<bool> &Initial)
{
	vector<long> Deleted;
	subopt = true;
	cerr << "Creating power graph" << endl;
	KGraph gs = g.CreatePowerGraph(k);
	cerr << "|E^k| = " << gs.m << endl;
	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3); //use concurrent method to solve LP relaxation.
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar *X = model.addVars(gs.m, GRB_CONTINUOUS);
		model.update();

		cerr << "Adding objective function" << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < gs.m; i++)
			Objective += X[i];
		model.setObjective(Objective, GRB_MINIMIZE);
		model.update();

		cerr << "Hashing edges" << endl;
		unordered_map<long, long> hash_edges;
		long cur = 0;
		for (long u = 0; u < gs.n; u++)
		{
			for (long i = 0; i < gs.degree[u]; i++)
			{
				long v = gs.adj[u][i];
				if (v > u) hash_edges.insert(make_pair(u*gs.n + v, cur++));
			}
		}

		cerr << "Fixing variables" << endl;
		vector<long> NonCriticalNodes = FindNonCriticalNodes(g);
		for (long i = 0; i < NonCriticalNodes.size(); i++)
			Y[NonCriticalNodes[i]].set(GRB_DoubleAttr_UB, 0);


		if (Initial[1])
		{
			cerr << "Adding minimal length-exactly-1 i,j-connectors" << endl;
			for (long i = 0; i < g.n; i++)
			{
				for (long it = 0; it < g.degree[i]; it++)
				{
					long j = g.adj[i][it];
					if (j > i) model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[j] >= 1); //checking i>j is due to hashing
				}
			}
		}


		if (Initial[2])
		{
			cerr << "Adding minimal length-exactly-2 i,j-connectors" << endl;
			for (long i = 0; i < g.n; i++)
			{
				vector<bool> i_neighbors = boolify(g.adj[i], g.n);
				for (long it1 = 0; it1 < g.degree[i]; it1++)
				{
					long u = g.adj[i][it1];
					for (long it2 = 0; it2 < g.degree[u]; it2++)
					{
						long j = g.adj[u][it2];
						if (i >= j) continue;
						if (i_neighbors[j]) continue;
						model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[j] + Y[u] >= 1);
					}
				}
			}
		}


		if (Initial[3] && k >= 3)
		{
			cerr << "Adding minimal length-exactly-3 i,j-connectors" << endl;
			for (long i = 0; i < g.n; i++)
			{
				vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
				for (long it1 = 0; it1 < g.degree[i]; it1++)
				{
					long u = g.adj[i][it1];
					vector<long> dist_from_u = g.ShortestPathsUnweighted(u);
					for (long it2 = 0; it2 < g.degree[u]; it2++)
					{
						long v = g.adj[u][it2];
						if (dist_from_i[v] == 1) continue;
						for (long it2 = 0; it2 < g.degree[v]; it2++)
						{
							long j = g.adj[v][it2];
							if (dist_from_i[j] == 1) continue;
							if (dist_from_u[j] == 1) continue;
							if (i >= j) continue;
							model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[u] + Y[v] + Y[j] >= 1);
						}
					}
				}
			}
		}


		if (Initial[4] && k >= 4)
		{
			cerr << "Adding minimal length-exactly-4 i,j-connectors" << endl;
			for (long i = 0; i < g.n; i++)
			{
				vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
				for (long j = i + 1; j < g.n; j++)
				{
					if (dist_from_i[j] == 1 || dist_from_i[j] > 4) continue;
					vector<long> dist_from_j = g.ShortestPathsUnweighted(j);
					for (long it1 = 0; it1 < g.degree[i]; it1++)
					{
						long p = g.adj[i][it1];
						vector<bool> p_neighbors = boolify(g.adj[p], g.n);
						if (dist_from_j[p] != 2) continue;
						for (long it2 = 0; it2 < g.degree[p]; it2++)
						{
							long v = g.adj[p][it2];
							if (dist_from_i[v] != 2 || dist_from_j[v] != 2) continue;
							for (long it3 = 0; it3 < g.degree[v]; it3++)
							{
								long q = g.adj[v][it3];
								if (dist_from_i[q] == 3 && dist_from_j[q] == 1 && p_neighbors[q] == false)
									model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[p] + Y[v] + Y[q] + Y[j] >= 1); //type1 V_12, V_22, V_21
								if (dist_from_i[q] == 3 && dist_from_j[q] == 1)
									model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[p] + Y[v] + Y[q] + Y[j] >= 1); //type2 V_12, V_22, V_31
							}
						}
					}
				}
			}
			for (long i = 0; i < g.n; i++) 
			{
				vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
				for (long j = i + 1; j < g.n; j++)
				{
					if (dist_from_i[j] == 1 || dist_from_i[j] > 4) continue;
					vector<long> dist_from_j = g.ShortestPathsUnweighted(j);
					for (long it1 = 0; it1 < g.degree[i]; it1++)
					{
						long u = g.adj[i][it1];
						if (dist_from_j[u] != 3) continue;
						for (long it2 = 0; it2 < g.degree[u]; it2++)
						{
							long v = g.adj[u][it2];
							if (dist_from_i[v] != 2 || dist_from_j[v] != 2) continue;
							for (long it3 = 0; it3 < g.degree[v]; it3++)
							{
								long q = g.adj[v][it3];
								if (dist_from_i[q] == 2 && dist_from_j[q] == 1)
									model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[u] + Y[v] + Y[q] + Y[j] >= 1); //type3 V_13, V_22, V_21
								if (dist_from_i[q] == 3 && dist_from_j[q] == 1)
									model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[u] + Y[v] + Y[q] + Y[j] >= 1); //type4 V_13, V_22, V_31
							}
						}
					}
				}
			}
		}
		

		cerr << "Adding budget constraints" << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
			expr1 += Y[i];
		model.addConstr(expr1 <= b);
		model.update();


		cerr << "Adding lazy constraints" << endl;
		FractionalSeparation cb = FractionalSeparation(X, Y, &g, &gs, hash_edges, k, b);
		model.setCallback(&cb);


		cerr << "Providing initial solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}
		vector<bool> NonDeletedNodes = boolify(Heuristic_sol, g.n);
		NonDeletedNodes.flip();

		for (long i = 0; i < g.n; i++)
		{
			vector <long> dist_from_i = g.ShortestPathsUnweighted(i, NonDeletedNodes);
			for (long j = i + 1; j < g.n; j++)
				if (dist_from_i[j] <= k)
					X[hash_edges[i*g.n + j]].set(GRB_DoubleAttr_Start, 1);
		}


		cerr << "Optimizing" << endl;
		model.optimize();


		double bestUB = model.get(GRB_DoubleAttr_ObjVal);
		double bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << "bestLB = " << bestLB << " bestUB = " << bestUB << " ";

		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B nodes = " << NumOfBandBNodes << endl;
		cerr << "# callbacks = " << FractionalSeparation::numCallbacks << endl;
		cerr << "# lazy cuts in integer separation part = " << FractionalSeparation::numLazyCutsInteger << endl;
		cerr << "# lazy cuts in fractional separation part = " << FractionalSeparation::numLazyCutsFractional << endl;
		cerr << "Time spent in integer part of callback = " << FractionalSeparation::TotalCallbackTimeInteger << endl;
		cerr << "Time spent in fractional part of callback = " << FractionalSeparation::TotalCallbackTimeFractional << endl;

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subopt = false;
			for (long i = 0; i < g.n; i++)
				if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
					Deleted.push_back(i);
		}
		
		delete[] Y;
		delete[] X;

	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return Deleted;
}


vector<long> Path_like_k3(KGraph &g, long b, vector<long> &Heuristic_sol, bool &subopt)
{
	vector<long> Deleted;
	long k = 3;
	subopt = true;
	cerr << "Creating power graph " << endl;
	KGraph gk = g.CreatePowerGraph(k);
	cerr << "|E^k| = " << gk.m << endl;
	time_t start_EnumerateConnector = clock();
	map<vector<long>, long, classcomp>map = EnumerateLength3Connector(g);
	cerr << " Number of minimal length-3 connectors = " << map.size() << endl;
	cerr << "Connector Enumeration time = " << (double)(clock() - start_EnumerateConnector) / CLOCKS_PER_SEC << endl;

	try
	{
		GRBEnv env = GRBEnv();
		//env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3); //use concurrent method to solve LP relaxation.
		long time = 3600 - (clock() - start_EnumerateConnector) / CLOCKS_PER_SEC;
		cerr << "time = " << time << endl;
		env.set(GRB_DoubleParam_TimeLimit, time);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(gk.m, GRB_CONTINUOUS);
		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar *Z = model.addVars(map.size(), GRB_BINARY);
		model.update();


		cerr << "Adding objective function" << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < gk.m; i++)
			Objective += X[i];
		model.setObjective(Objective, GRB_MINIMIZE);
		model.update();


		cerr << "Hashing edges" << endl;
		unordered_map<long, long> hash_edges;
		long cur = 0;
		for (long u = 0; u < gk.n; u++)
		{
			for (long i = 0; i < gk.degree[u]; i++)
			{
				long v = gk.adj[u][i];
				if (v > u) hash_edges.insert(make_pair(u*gk.n + v, cur++));
			}
		}


		cerr << "Fixing variables" << endl;
		vector<long> NonCriticalNodes = FindNonCriticalNodes(g);
		for (long i = 0; i < NonCriticalNodes.size(); i++)
			Y[NonCriticalNodes[i]].set(GRB_DoubleAttr_UB, 0);


		cerr << "Adding constraints (2b) & (2e)" << endl;
		for (long u = 0; u < g.n; u++)
		{
			vector<long> dist_from_u = g.ShortestPathsUnweighted(u);
			for (long v = u + 1; v < g.n; v++)
			{
				if (dist_from_u[v] > 3) continue;
				vector<long> dist_from_v = g.ShortestPathsUnweighted(v);
				if (dist_from_u[v] == 1)
				{
					std::map<std::vector<long>, long>::iterator it = map.find({ u,v });
					model.addConstr(Z[it->second] + Y[u] + Y[v] >= 1); //2b
					model.addConstr(X[hash_edges[u*g.n + v]] >= Z[it->second]); //2e
					continue; //No other minimal length-3 i,j-connectors exist
				}
				for (long v_neighbors_iterator = 0; v_neighbors_iterator < g.degree[v]; v_neighbors_iterator++)
				{
					long ii = g.adj[v][v_neighbors_iterator];
					if (dist_from_u[ii] == 1)
					{
						vector<long> sortedsubset = sortnodes(u, v, ii);
						std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
						model.addConstr(Z[it->second] + Y[u] + Y[v] + Y[ii] >= 1); //2b
						model.addConstr(X[hash_edges[u*g.n + v]] >= Z[it->second]); //2e
					 }
					if (dist_from_u[ii] == 2)
					{
						for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.degree[ii]; ii_neighbors_iterator++)
						{
							long jj = g.adj[ii][ii_neighbors_iterator];
							if (dist_from_u[jj] == 1 && dist_from_v[jj] == 2)
							{
								vector<long> sortedsubset = sort4nodes(u, v, ii, jj);
								std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
								model.addConstr(Z[it->second] + Y[u] + Y[v] + Y[ii] + Y[jj] >= 1); //2b
								model.addConstr(X[hash_edges[u*g.n + v]] >= Z[it->second]); //2e
							}
						}
					}
				}
			}
		}
	

		cerr << "Adding budget constraints" << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
			expr1 += Y[i];
		model.addConstr(expr1 <= b);
		model.update();


		cerr << "Providing initial solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}
		vector<bool> NonDeletedNodes = boolify(Heuristic_sol, g.n);
		NonDeletedNodes.flip();

		for (long i = 0; i < g.n; i++)
		{
			vector <long> dist_from_i = g.ShortestPathsUnweighted(i, NonDeletedNodes);
			for (long j = i + 1; j < g.n; j++)
				if (dist_from_i[j] <= k)
					X[hash_edges[i*g.n + j]].set(GRB_DoubleAttr_Start, 1);
		}


		cerr << "Optimizing" << endl;
		model.optimize();

		double bestUB = model.get(GRB_DoubleAttr_ObjVal);
		double bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << "bestLB = " << bestLB << " bestUB = " << bestUB << " ";


		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B Nodes : " << NumOfBandBNodes << endl;


		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subopt = false;
			for (long i = 0; i < g.n; i++)
				if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
					Deleted.push_back(i);
		}

		delete[] Y;
		delete[] X;
		delete[] Z;

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return Deleted;
}


vector<long> Path_like_k4(KGraph &g, long b, vector<long> &Heuristic_sol, bool &subopt)
{
	vector<long> Deleted;
	long k = 4;
	subopt = true;
	cerr << "Creating power graph " << endl;
	KGraph gk = g.CreatePowerGraph(k);
	cerr << "|E^k| = " << gk.m << endl;
	time_t start_EnumerateConnector = clock();

	map<vector<long>, long, classcomp>map1 = EnumerateLength3Connector(g);
	map<vector<long>, long, classcomp>map2 = EnumerateLength4Connector(g);
	cerr << "map1 size = " << map1.size() << endl;
	cerr << "map2 size = " << map2.size() << endl;


	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3); //use concurrent method to solve LP relaxation.
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(gk.m, GRB_CONTINUOUS);
		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar *Z = model.addVars(map1.size(), GRB_BINARY);
		GRBVar *W = model.addVars(map2.size(), GRB_BINARY);
		model.update();


		cerr << "Adding objective function" << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < gk.m; i++)
			Objective += X[i];
		model.setObjective(Objective, GRB_MINIMIZE);
		model.update();


		cerr << "Hashing edges" << endl;
		unordered_map<long, long> hash_edges;
		long cur = 0;
		for (long u = 0; u < gk.n; u++)
		{
			for (long i = 0; i < gk.degree[u]; i++)
			{
				long v = gk.adj[u][i];
				if (v > u) hash_edges.insert(make_pair(u*gk.n + v, cur++));
			}
		}


		cerr << "Fixing variables" << endl;
		vector<long> NonCriticalNodes = FindNonCriticalNodes(g);
		for (long i = 0; i < NonCriticalNodes.size(); i++)
			Y[NonCriticalNodes[i]].set(GRB_DoubleAttr_UB, 0);



		cerr << "Adding constraints (2b) & (2e)" << endl;
		for (long u = 0; u < g.n; u++)
		{
			vector<long> dist_from_u = g.ShortestPathsUnweighted(u);
			for (long v = u + 1; v < g.n; v++)
			{
				if (dist_from_u[v] > 3) continue;
				vector<long> dist_from_v = g.ShortestPathsUnweighted(v);
				if (dist_from_u[v] == 1)
				{
					std::map<std::vector<long>, long>::iterator it = map1.find({ u,v });
					model.addConstr(Z[it->second] + Y[u] + Y[v] >= 1); //2b
					model.addConstr(X[hash_edges[u*g.n + v]] >= Z[it->second]); //2e
					continue; //No other minimal length-3 i,j-connectors exist
				}
				for (long v_neighbors_iterator = 0; v_neighbors_iterator < g.degree[v]; v_neighbors_iterator++)
				{
					long ii = g.adj[v][v_neighbors_iterator];
					if (dist_from_u[ii] == 1)
					{
						vector<long> sortedsubset = sortnodes(u, v, ii);
						std::map<std::vector<long>, long>::iterator it = map1.find(sortedsubset);
						model.addConstr(Z[it->second] + Y[u] + Y[v] + Y[ii] >= 1); //2b
						model.addConstr(X[hash_edges[u*g.n + v]] >= Z[it->second]); //2e
					}
					if (dist_from_u[ii] == 2)
					{
						for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.degree[ii]; ii_neighbors_iterator++)
						{
							long jj = g.adj[ii][ii_neighbors_iterator];
							if (dist_from_u[jj] == 1 && dist_from_v[jj] == 2)
							{
								vector<long> sortedsubset = sort4nodes(u, v, ii, jj);
								std::map<std::vector<long>, long>::iterator it = map1.find(sortedsubset);
								model.addConstr(Z[it->second] + Y[u] + Y[v] + Y[ii] + Y[jj] >= 1); //2b
								model.addConstr(X[hash_edges[u*g.n + v]] >= Z[it->second]); //2e
							}
						}
					}
				}
			}
		}

		cerr << "New" << endl;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist_from_i[j] == 1 || dist_from_i[j] > 4) continue;
				vector<long> dist_from_j = g.ShortestPathsUnweighted(j);
				for (long it1 = 0; it1 < g.degree[i]; it1++)
				{
					long p = g.adj[i][it1];
					vector<bool> p_neighbors = boolify(g.adj[p], g.n);
					if (dist_from_j[p] != 2) continue;
					for (long it2 = 0; it2 < g.degree[p]; it2++)
					{
						long v = g.adj[p][it2];
						if (dist_from_i[v] != 2 || dist_from_j[v] != 2) continue;
						for (long it3 = 0; it3 < g.degree[v]; it3++)
						{
							long q = g.adj[v][it3];
							if (dist_from_i[q] == 2 && dist_from_j[q] == 1 && p_neighbors[q] == false)
							{
								vector<long> sortedsubset = sort5nodes(i, p, v, q, j);
								std::map<std::vector<long>, long>::iterator it = map2.find(sortedsubset);
								model.addConstr(W[it->second] + Y[i] + Y[p] + Y[v] + Y[q] + Y[j] >= 1); //2b
								model.addConstr(X[hash_edges[i*g.n + j]] >= W[it->second]); //2e
							}

							if (dist_from_i[q] == 3 && dist_from_j[q] == 1)
							{
								vector<long> sortedsubset = sort5nodes(i, p, v, q, j);
								std::map<std::vector<long>, long>::iterator it = map2.find(sortedsubset);
								model.addConstr(W[it->second] + Y[i] + Y[p] + Y[v] + Y[q] + Y[j] >= 1); //2b
								model.addConstr(X[hash_edges[i*g.n + j]] >= W[it->second]); //2e
							}	

						}
					}
				}
			}
		}

		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist_from_i[j] == 1 || dist_from_i[j] > 4) continue;
				vector<long> dist_from_j = g.ShortestPathsUnweighted(j);
				for (long it1 = 0; it1 < g.degree[i]; it1++)
				{
					long u = g.adj[i][it1];
					if (dist_from_j[u] != 3) continue;
					for (long it2 = 0; it2 < g.degree[u]; it2++)
					{
						long v = g.adj[u][it2];
						if (dist_from_i[v] != 2 || dist_from_j[v] != 2) continue;
						for (long it3 = 0; it3 < g.degree[v]; it3++)
						{
							long q = g.adj[v][it3];
							if (dist_from_i[q] == 2 && dist_from_j[q] == 1)
							{
								vector<long> sortedsubset = sort5nodes(i, u, v, q, j);
								std::map<std::vector<long>, long>::iterator it = map2.find(sortedsubset);
								model.addConstr(W[it->second] + Y[i] + Y[u] + Y[v] + Y[q] + Y[j] >= 1); //2b
								model.addConstr(X[hash_edges[i*g.n + j]] >= W[it->second]); //2e
							}
								
							if (dist_from_i[q] == 3 && dist_from_j[q] == 1)
							{
								vector<long> sortedsubset = sort5nodes(i, u, v, q, j);
								std::map<std::vector<long>, long>::iterator it = map2.find(sortedsubset);
								model.addConstr(W[it->second] + Y[i] + Y[u] + Y[v] + Y[q] + Y[j] >= 1); //2b
								model.addConstr(X[hash_edges[i*g.n + j]] >= W[it->second]); //2e
							}			
						}
					}
				}
			}
		}

		cerr << "Adding budget constraints" << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
			expr1 += Y[i];
		model.addConstr(expr1 <= b);
		model.update();


		cerr << "Providing initial solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}
		vector<bool> NonDeletedNodes = boolify(Heuristic_sol, g.n);
		NonDeletedNodes.flip();

		for (long i = 0; i < g.n; i++)
		{
			vector <long> dist_from_i = g.ShortestPathsUnweighted(i, NonDeletedNodes);
			for (long j = i + 1; j < g.n; j++)
				if (dist_from_i[j] <= k)
					X[hash_edges[i*g.n + j]].set(GRB_DoubleAttr_Start, 1);
		}


		cerr << "Optimizing" << endl;
		model.optimize();

		double bestUB = model.get(GRB_DoubleAttr_ObjVal);
		double bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << "bestLB = " << bestLB << " bestUB = " << bestUB << " ";


		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B Nodes : " << NumOfBandBNodes << endl;


		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subopt = false;
			for (long i = 0; i < g.n; i++)
				if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
					Deleted.push_back(i);
		}

		delete[] Y;
		delete[] X;
		delete[] Z;
		delete[] W;

	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return Deleted;
}


vector<long> Recursive(KGraph &g, long k, long b, vector<long> &Heuristic_sol, bool &subopt)
{
	vector<long> Deleted;
	subopt = true;
	cerr << "Creating power graph " << endl;
	KGraph gs = g.CreatePowerGraph(k);

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3); //use concurrent method to solve LP relaxation.
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);

		time_t Building_Model_Start = clock();
		GRBVar *X = model.addVars(gs.m, GRB_CONTINUOUS);
		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar ***U = new GRBVar**[g.n]; //U[i][j][s] denotes whether distance between i and j is at most s in the interdicted graph
		for (long i = 0; i<g.n; i++)	
		{
			GRBVar **U_temp = new GRBVar*[g.n];
			for (long j = i+1; j < g.n; j++)
				U_temp[j] = model.addVars(k, GRB_BINARY);

			U[i] = U_temp;
		}
		model.update();

		cerr << "Hashing edges" << endl;
		unordered_map<long, long> hash_edges;
		long cur = 0;
		for (long u = 0; u < gs.n; u++)
		{
			for (long i = 0; i < gs.adj[u].size(); i++)
			{
				long v = gs.adj[u][i];
				if (v > u) hash_edges.insert(make_pair(u*gs.n + v, cur++));
			}
		}

		cerr << "Adding objective function" << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < gs.m; i++)
			Objective += X[i];
		model.setObjective(Objective, GRB_MINIMIZE);
		model.update();

		cerr << "Fixing variables" << endl;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist_from_i[j] > k)
					for (long s = 0; s < k; s++)
						U[i][j][s].set(GRB_DoubleAttr_UB, 0);

				else
					for (long s = 0; s < dist_from_i[j] - 1; s++)
						U[i][j][s].set(GRB_DoubleAttr_UB, 0);
			}
		}

		vector<long> NonCriticalNodes = FindNonCriticalNodes(g);
		for (long i = 0; i < NonCriticalNodes.size(); i++)
			Y[NonCriticalNodes[i]].set(GRB_DoubleAttr_UB, 0);

		cerr << "Adding constraints (1b)" << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
			expr1 += Y[i];
		model.addConstr(expr1 <= b);


		cerr << "Adding constraints (1c)" << endl;
		for (long i = 0; i < g.n; i++)
		{
			for (long counter = 0; counter < g.degree[i]; counter++)
			{
				long j = g.adj[i][counter];
				if (j > i) model.addConstr(U[i][j][0] + Y[i] + Y[j] >= 1);
			}
		}

		cerr << "Adding constraint (1d)" << endl;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist_from_i[j] > k) continue;

				for (long s = dist_from_i[j]; s < k; s++)
					model.addConstr(U[i][j][s] + Y[i] <= 1);
			}
		}

		cerr << "Adding constraint (1e)" << endl;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
			for (long counter = 0; counter < g.degree[i]; counter++)
			{
				long j = g.adj[i][counter];
				if (dist_from_i[j] > k) continue;
				if (j > i)
					for (long s = dist_from_i[j]; s < k; s++)
						model.addConstr(U[i][j][s] == U[i][j][0]);
			}
		}

		cerr << "Adding constraints (1f)" << endl;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist_from_i[j] > k) continue;
				if (dist_from_i[j] == 1) continue;
				for (long s = dist_from_i[j]; s < k; s++)
				{
					GRBLinExpr expr = 0;
					for (long counter = 0; counter < g.degree[i]; counter++)
					{
						long t = g.adj[i][counter];
						if(j > t) expr += U[t][j][s - 1];
						if (t > j) expr += U[j][t][s - 1];
					}
					model.addConstr(U[i][j][s] <= expr);
				}
			}
		}

		cerr << "Adding constraints (1g)" << endl;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist_from_i[j] == 1) continue;
				for (long s = dist_from_i[j] - 1; s < k; s++)
				{
					for (long counter = 0; counter < g.degree[i]; counter++)
					{
						long t = g.adj[i][counter];
						if (j > t) model.addConstr(U[t][j][s - 1] <= U[i][j][s] + Y[i]);
						if (t > j) model.addConstr(U[j][t][s - 1] <= U[i][j][s] + Y[i]);
					}
				}
			}
		}

		cerr << "Adding constraints (1i)" << endl;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist_from_i[j] > k) continue;
				model.addConstr(U[i][j][k - 1] <= X[hash_edges[i*g.n + j]]);
			}
		}
		model.update();
		cerr << "Time spent for building the model = " << (double)(clock() - Building_Model_Start) / CLOCKS_PER_SEC << endl;


		cerr << "Providing initial solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}

		vector<bool> NonDeletedNodes = boolify(Heuristic_sol, g.n);
		NonDeletedNodes.flip();


		for (long v = 0; v < g.n; v++)
		{
			vector <long> dist_from_v = g.ShortestPathsUnweighted(v, NonDeletedNodes);
			for (long w = v + 1; w < g.n; w++)
			{
				if (dist_from_v[w] <= k)
				{
					X[hash_edges[v*g.n + w]].set(GRB_DoubleAttr_Start, 1);
					for (long s = dist_from_v[w]-1; s < k; s++)
						U[v][w][s].set(GRB_DoubleAttr_Start, 1);
				}
					
			}
		}


		cerr << "Optimizing" << endl;
		model.optimize();

		double bestUB = model.get(GRB_DoubleAttr_ObjVal);
		double bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << "bestLB = " << bestLB << " bestUB = " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subopt = false;
			for (long i = 0; i < g.n; i++)
				if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
					Deleted.push_back(i);
		}

		delete[] X;
		delete[] Y;
		delete[] U;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return Deleted;
}


vector<long> Thin_Weighted(KGraph &g, long k, long b, vector<long> &Heuristic_sol, bool &subopt)
{
	vector<long> Deleted;
	subopt = true;
	cerr << "Creating power graph " << endl;
	KGraph gs = g.CreatePowerGraphWeighted(k);
	cerr << "|E^k| = " << gs.m << endl;

	try
	{
		GRBEnv env = GRBEnv();
		//env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_Method, 3); //use concurrent method to solve LP relaxation.
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		GRBModel model = GRBModel(env);
		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar *X = model.addVars(gs.m, GRB_CONTINUOUS);
		model.update();

		cerr << "Adding objective function" << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < gs.m; i++)
		{
			Objective += X[i];
		}
		model.setObjective(Objective, GRB_MINIMIZE);
		model.update();


		cerr << "Hashing edges" << endl;
		unordered_map<long, long> hash_edges;
		long cur = 0;
		for (long u = 0; u < gs.n; u++)
		{
			for (long i = 0; i < gs.adj[u].size(); i++)
			{
				long v = gs.adj[u][i];
				if (v > u) hash_edges.insert(make_pair(u*gs.n + v, cur++));
			}
		}


		cerr << "Adding constraints when hop distance = 1" << endl;
		for (long i = 0; i < g.n; i++)
		{
			for (long counter = 0; counter < g.degree[i]; counter++)
			{
				long j = g.adj[i][counter];
				if (g.weight[i][counter] <= k)
					model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[j] >= 1);
			}
		}


		cerr << "Adding constraints for length-k i,j-connectors when hop distance = 2" << endl;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> hop_dist_from_i = g.ShortestPathsUnweighted(i);
			for (long counter1 = 0; counter1 < g.degree[i]; counter1++)
			{
				long u = g.adj[i][counter1];
				if (g.weight[i][counter1] > k) continue;
				for (long counter2 = 0; counter2 < g.degree[u]; counter2++)
				{
					long j = g.adj[u][counter2];
					if (i >= j) continue;
					if (hop_dist_from_i[j] == 1) continue;
					if (g.weight[i][counter1] + g.weight[u][counter2] > k) continue;
					model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[j] + Y[u] >= 1);
				}
			}
		}


		cerr << "Adding budget constraints" << endl;
		GRBLinExpr expr2 = 0;
		for (long i = 0; i < g.n; i++)
			expr2 += Y[i];
		model.addConstr(expr2 <= b);
		model.update();


		cerr << "Adding lazy constraints" << endl;
		IntegerSeparationWeighted cb = IntegerSeparationWeighted(X, Y, &g, &gs, hash_edges, k, b);
		model.setCallback(&cb);


		cerr << "Providing initial solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}
		vector<bool> NonDeletedNodes = boolify(Heuristic_sol, g.n);
		NonDeletedNodes.flip();

		for (long i = 0; i < g.n; i++)
		{
			vector <long> dist_from_i = g.ShortestPathsWeighted(i, NonDeletedNodes);
			for (long j = i + 1; j < g.n; j++)
			{
				if (dist_from_i[j] <= k)
					X[hash_edges[i*g.n + j]].set(GRB_DoubleAttr_Start, 1);
			}
		}

		cerr << "Optimizing" << endl;
		model.optimize();

		double bestUB = model.get(GRB_DoubleAttr_ObjVal);
		double bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << "bestLB = " << bestLB << ", bestUB = " << bestUB << " ";

		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B nodes = " << NumOfBandBNodes << endl;
		cerr << "# callbacks = " << IntegerSeparationWeighted::numCallbacks << endl;
		cerr << "# Lazy Cuts = " << IntegerSeparationWeighted::numLazyCutsInteger << endl;
		cerr << "Time spent in callback = " << IntegerSeparationWeighted::TotalCallbackTime << endl;

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subopt = false;
			for (long i = 0; i < g.n; i++)
				if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
					Deleted.push_back(i);
		}

		delete[] Y;
		delete[] X;

	}


	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}

	return Deleted;
}


//Check if a vertex i is simplicial
bool IsSimplicial(KGraph &g, long i)
{
	bool simplicial = true;
	for (long counter1 = 0; counter1 < g.degree[i]; counter1++)
	{
		long j = g.adj[i][counter1];
		if (g.CommonNeighbors(j, g.adj[i]) != g.degree[i] - 1)
			simplicial = false;
	}
	return simplicial;
}

/*DCNP preprocessing function
*It works for arbitrary k and b*/
vector<long> FindNonCriticalNodes(KGraph &g)
{
	vector<long> simplicial;

	for (long i = 0; i < g.n; i++)
		if (IsSimplicial(g, i))
			simplicial.push_back(i);

	vector<long> Rmap1;
	KGraph gnew = g.CreateInducedGraph(simplicial, Rmap1);

	vector<long> set_I;
	vector< vector< long> > components;
	vector<long> degreeZero;
	gnew.FindConnectedComponents(components, degreeZero);

	for (long i = 0; i < degreeZero.size(); i++)
		set_I.push_back(simplicial[Rmap1[simplicial[degreeZero[i]]]]);

	for (long i = 0; i < components.size(); i++)
		set_I.push_back(simplicial[Rmap1[simplicial[components[i][0]]]]);

	cerr << "# simplicial vertices fixed = " << set_I.size() << endl;

	long leaves = 0;
	for (long i = 0; i<g.n; i++)
		if (g.degree[i] == 1)
		{
			long v = g.adj[i][0];
			if (g.degree[v] > 1 || i < v) leaves++;
		}
	cerr << "# leaves fixed = " << leaves << endl;

	return set_I;
}

//Calculate number of vertex pairs with distance at most k in graph G-D (obj(G-D))
long obj(KGraph &g, vector<long> &deleted, long k)
{
	vector<bool> nondeleted(g.n, true);
	for (long i = 0; i < deleted.size(); i++) nondeleted[deleted[i]] = false;
	return obj(g, nondeleted, k);
}

long obj(KGraph &g, vector<bool> &nondeleted, long k)
{
	long num_close_vertices = 0;
	for (long v = 0; v < g.n; v++)
	{
		vector <long> dist_from_v = g.ShortestPathsUnweighted(v, nondeleted);
		for (long w = v + 1; w < g.n; w++)
			if (nondeleted[w] && dist_from_v[w] <= k)
				num_close_vertices++;
	}
	return num_close_vertices;
}

long obj_weighted(KGraph &g, vector<long> &deleted, long k)
{
	vector<bool> nondeleted(g.n, true);
	for (long i = 0; i < deleted.size(); i++) nondeleted[deleted[i]] = false;
	return obj_weighted(g, nondeleted, k);
}

long obj_weighted(KGraph &g, vector<bool> &nondeleted, long k)
{
	long num_close_vertices = 0;
	for (long v = 0; v < g.n; v++)
	{
		vector <long> dist_from_v = g.ShortestPathsWeighted(v, nondeleted);
		for (long w = v + 1; w < g.n; w++)
			if (nondeleted[w] && dist_from_v[w] <= k)
				num_close_vertices++;			
	}
	return num_close_vertices;
}


/*Based on Algorithm 1 of paper "A faster algorithm for betweenness centerality",
by Ulrik Brandes. The running time of this algorithm is O(mn) for unweighted graphs.*/
vector<long> FindTopTBetweennessCentralityNodes(KGraph &g, long T)
{
	vector<double> C_B(g.n, 0);
	for (long s = 0; s < g.n; s++)
	{
		stack<int> Stack;
		vector< vector<long> > P(g.n);
		vector<double> sigma(g.n, 0); //sigma[t] = number of shortest paths from s (source) to t
		sigma[s] = 1;
		vector<long> d(g.n, -1); //d[t] = distance from s to t
		d[s] = 0;
		vector<long> Q; //A queue
		Q.push_back(s);

		while (!Q.empty())
		{
			long v = Q[0];
			Q.erase(Q.begin());
			Stack.push(v);

			for (long i = 0; i < g.degree[v]; i++)
			{
				long w = g.adj[v][i];

				//w found for the first time?
				if (d[w] < 0)
				{
					Q.push_back(w);
					d[w] = d[v] + 1;
				}

				//shortest path to w via v?
				if (d[w] == d[v] + 1)
				{
					sigma[w] = sigma[w] + sigma[v];
					P[w].push_back(v);
				}
			}
		}

		vector<double> delta(g.n, 0.0);
		//Stack returns vertices in order of non-increasing distance from s
		while (!Stack.empty())
		{
			long w = Stack.top();
			Stack.pop();

			for (long j = 0; j < P[w].size(); j++)
			{
				long v = P[w][j];
				delta[v] = double(delta[v] + double((sigma[v] / sigma[w]))*(1 + delta[w]));
			}

			if (w != s)
				C_B[w] = C_B[w] + delta[w];
		}
	}

	//Find top T betweenness centrality nodes
	vector<long> TopT_BC;
	priority_queue<pair<double, long>> q;
	for (long i = 0; i < C_B.size(); i++)
		q.push(pair<double, long>(C_B[i], i));

	for (long i = 0; i < T; i++)
	{
		long index = q.top().second;
		TopT_BC.push_back(index);
		q.pop();
	}

	return TopT_BC;
}

vector<long> FindTopTBetweennessCentralityNodesWeighted(KGraph &g, long T, long k)
{
	vector<double> C_B(g.n, 0);
	vector<long> d;
	for (long s = 0; s < g.n; s++)
	{
		vector<bool> visited(g.n, false);
		d = g.ShortestPathsWeighted(s);
		stack<int> Stack;
		vector< vector<long> > P(g.n);
		vector<double> sigma(g.n, 0); //sigma[t] = number of shortest paths from s (source) to t
		sigma[s] = 1;
		vector<long> Q; //A queue
		Q.push_back(s);

		while (!Q.empty())
		{
			long v = Q[0];
			Q.erase(Q.begin());
			if (!visited[v] && d[v] <= k)
			{
				Stack.push(v);
				visited[v] = true;
				for (long i = 0; i < g.degree[v]; i++)
				{
					long w = g.adj[v][i];
					if (d[w] == d[v] + g.weight[v][i])
					{
						sigma[w] = sigma[w] + sigma[v];
						P[w].push_back(v);
						Q.push_back(w);
					}
				}
			}	
		}
		vector<double> delta(g.n, 0.0);

		//Stack returns vertices in order of non-increasing distance from s
		while (!Stack.empty())
		{
			long w = Stack.top();
			Stack.pop();

			for (long j = 0; j < P[w].size(); j++)
			{
				long v = P[w][j];
				delta[v] = double(delta[v] + double((sigma[v] / sigma[w]))*(1 + delta[w]));
			}

			if (w != s)
				C_B[w] = C_B[w] + delta[w];
		}
	}

	//Find top T betweenness centrality nodes
	vector<long> TopT_BC;
	priority_queue<pair<double, long>> q;
	for (long i = 0; i < C_B.size(); i++)
		q.push(pair<double, long>(C_B[i], i));

	for (long i = 0; i < T; i++)
	{
		long index = q.top().second;
		TopT_BC.push_back(index);
		q.pop();
	}

	return TopT_BC;
}


//DCNP Heuristic to find distance-based critical nodes
vector<long> DCNP_Heuristic(KGraph &g, long k, long b)
{
	vector<long> D_Star;
	vector<long> D = FindTopTBetweennessCentralityNodes(g, 2 * b); //D contains top 2b betweenness centrality nodes
	vector<bool> NonDeleted(g.n,true);
	for (long i = 0; i < D.size(); i++)
		NonDeleted[D[i]] = false;

	for (long counter = 0; counter < b; counter++)
	{
		long min = LONG_MAX;
		long argmin = -1;
		for (long i = 0; i < g.n; i++)
		{
			if (NonDeleted[i]) continue;
			NonDeleted[i] = true;
			long TempObj = obj(g, NonDeleted, k);
			if (TempObj < min)
			{
				min = TempObj;
				argmin = i;
			}
			NonDeleted[i] = false;
		}
		NonDeleted[argmin] = true;
	}
	NonDeleted.flip();
	for (long i = 0; i < g.n; i++)
		if (NonDeleted[i] == true)
			D_Star.push_back(i);

	return D_Star;
}

vector<long> DCNP_Heuristic_Weighted(KGraph &g, long k, long b)
{
	vector<long> D_Star;
	vector<long> D = FindTopTBetweennessCentralityNodesWeighted(g, 2 * b, k); //D contains top 2b betweenness centrality nodes
	vector<bool> NonDeleted(g.n, true);
	for (long i = 0; i < D.size(); i++)
		NonDeleted[D[i]] = false;

	for (long counter = 0; counter < b; counter++)
	{
		long min = LONG_MAX;
		long argmin = -1;
		for (long i = 0; i < g.n; i++)
		{
			if (NonDeleted[i]) continue;
			else NonDeleted[i] = true;
			long TempObj = obj_weighted(g, NonDeleted, k);
			if (TempObj < min)
			{
				min = TempObj;
				argmin = i;
			}
			NonDeleted[i] = false;
		}
		NonDeleted[argmin] = true;
	}
	NonDeleted.flip();
	for (long i = 0; i < g.n; i++)
		if (NonDeleted[i] == true)
			D_Star.push_back(i);

	return D_Star;
}


//Function to calculate variables d(v,s) and p(v,s) for fractional separation
d_and_p d_and_p_function(KGraph &goriginal, KGraph &gpower, long i, long k, double* y)
{
	long infty = numeric_limits<long>::max();
	vector<double> temp1(k + 1, infty);
	vector<long> temp2(k + 1, infty);
	vector < vector<double> > d(goriginal.n, temp1);
	vector< vector<long> > p(goriginal.n, temp2);

	for (long s = 0; s <= k; s++)
	{
		d[i][s] = y[i];
		p[i][s] = i;
	}

	for (long s = 1; s <= k; s++)
	{
		for (long counter1 = 0; counter1 < gpower.degree[i]; counter1++)
		{
			long v = gpower.adj[i][counter1];
			for (long counter2 = 0; counter2 < goriginal.degree[v]; counter2++)
			{
				long u = goriginal.adj[v][counter2];
				if (d[v][s] > d[u][s - 1] + y[v])
				{
					d[v][s] = d[u][s - 1] + y[v];
					p[v][s] = u;
				}
			}
		}
	}
	d_and_p result = { d,p };
	return result;
}



