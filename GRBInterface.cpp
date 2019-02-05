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
#include <unordered_map>



using namespace std;


string itos(int i) { stringstream s; s << i; return s.str(); }


long integer_separation::numCallbacks = 0;
double integer_separation::TotalCallbackTime = 0;
long integer_separation::numLazyCuts = 0;


long fractional_separation::numCallbacks = 0;
double fractional_separation::TotalCallbackTime = 0;
long fractional_separation::numLazyCuts = 0;


vector<vector<long>> k_hop_SP(KGraph &g, long source, long s, vector<long>dist_from_source)
{
	//initializing edge weights
	//g.Edgeweight(g);
	double infty = 1000000;

	vector<long> y(g.n);
	y[0] = 0; y[1] = 0.4; y[2] = 0.1;

	//vector<long> dist_from_source = g.ShortestPathsUnweighted(source);


	vector<double> temp1(s + 1, infty);
	vector<long> temp2(s + 1, 1000);
	vector < vector<double> > delta(g.n, temp1);
	vector< vector<long> > P(g.n, temp2);


	for (long k = 0; k <= s; k++)
	{
		delta[source][k] = y[source];
		P[source][k] = source;
	}

	/*for (long i = 0; i < g.degree[source]; i++)
	{
	long u = g.adj[source][i];
	delta[u][1] = g.weight[source][i];
	P[u][1] = source;
	}*/

	for (long k = 1; k <= s; k++)
	{
		for (long v = 0; v < g.n; v++)
		{
			if (v == source || dist_from_source[v] > k) continue;

			for (long i = 0; i < g.degree[v]; i++)
			{
				long u = g.adj[v][i];
				if (delta[u][k - 1] + y[v] < delta[v][k])
				{
					P[v][k] = u;
					delta[v][k] = delta[u][k - 1] + y[v];
				}
			}
		}
	}
	for (long k = 0; k <= s; k++)
	{
		for (long v = 0; v < g.n; v++)
		{
			cerr << "d[" << v << "][" << k << "] is " << delta[v][k] << endl;
		}
	}
	return P;
}

void integer_separation::callback()
{
	try
	{
		if (where == GRB_CB_MIPSOL)
		{
			/*cerr << "Found a new MIP incumbent" << endl;

			double objbest = getDoubleInfo(GRB_CB_MIPSOL_OBJBST);
			cerr << "Current best objective = " << objbest << endl;

			double nodecnt = getDoubleInfo(GRB_CB_MIPSOL_NODCNT);
			double objbnd = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
			cerr << "For node " << nodecnt << ", current best objective bound " << objbnd << endl;*/


			numCallbacks++;
			time_t start = clock();

			//get the solution vector (for y variables) from Gurobi 
			//g1 is duplicate of power graph
			double *y = new double[g1.n];
			y = getSolution(vars, g1.n);

			//now, make it and its complement boolean
			vector<bool> DeletedVertices(g1.n, false);
			vector<bool> NonDeletedVertices(g1.n, false);
			for (long i = 0; i < g1.n; i++)
			{
				if (y[i] > 0.5) DeletedVertices[i] = true;
				else NonDeletedVertices[i] = true;
			}

			//get the solution vector (for x variables) from Gurobi
			//g1 is duplicate of power graph
			double *x = new double[g1.m];
			x = getSolution(vars1, g1.m);

			//find a violated length-k i,j-connector inequality (if any exist)
			vector<long> predecessor;
			vector<long> distance;
			for (long u = 0; u < g2.n; u++)
			{
				vector<long>original_distance = g2.ShortestPathsUnweighted(u);
				distance = g2.ShortestPathsUnweighted(u, NonDeletedVertices, predecessor); //g2 is duplicate of original graph.
				for (long v = u + 1; v < g2.n; v++)
				{
					long counter = 0;
					if (original_distance[v] == 1) continue; //all connector inequalities are satisfied
					if (distance[v] <= k1 && x[hashing[u*g2.n + v]] == 0)
					{
						GRBLinExpr expr = vars1[hashing[u*g2.n + v]];
						long i = v;
						while (i != u)
						{
							expr += vars[i];
							i = predecessor[i];
						}
						expr += vars[i];
						addLazy(expr >= 1);
						counter++;
						numLazyCuts++;
						if (counter == 10 * b1) break; //with this we can change number of cuts. 
					}
					predecessor.clear();
				}
			}
		}
		/*	else if (where == GRB_CB_MIPNODE)
		{
		cerr << "Currently exploring a MIP node" << endl;

		double objbest = getDoubleInfo(GRB_CB_MIPNODE_OBJBST);
		cerr << "Current best objective = " << objbest << endl;

		double nodecnt = getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
		double objbnd = getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
		cerr << "For node " << nodecnt  << ", current best objective bound " << objbnd << endl;

		}

		else if (where == GRB_CB_MIP)
		{
		cerr << "Currently in MIP" << endl;

		double objbest = getDoubleInfo(GRB_CB_MIP_OBJBST);
		cerr << "Current best objective = " << objbest << endl;

		double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
		double objbnd = getDoubleInfo(GRB_CB_MIP_OBJBND);
		cerr << "For node " << nodecnt << ", current best objective bound " << objbnd << endl;
		}


		else if (where == GRB_CB_SIMPLEX)
		{
		cerr << "Currently in SIMPLEX." << endl;
		double objvalue = getDoubleInfo(GRB_CB_SPX_OBJVAL);
		cerr << "Current simplex obj value is " << objvalue << endl;
		}

		else if (where == GRB_CB_BARRIER)
		{

		cerr << "Currently in BARRIER." << endl;
		double primalobjval = getDoubleInfo(GRB_CB_BARRIER_PRIMOBJ);
		cerr << "Primal objective value for current barrier iterate is " << primalobjval <<endl;
		}*/
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}


void fractional_separation::callback()
{
	try
	{
		if (where == GRB_CB_MIPSOL)
		{
			/*cerr << "Found a new MIP incumbent" << endl;

			double objbest = getDoubleInfo(GRB_CB_MIPSOL_OBJBST);
			cerr << "Current best objective = " << objbest << endl;

			double nodecnt = getDoubleInfo(GRB_CB_MIPSOL_NODCNT);
			double objbnd = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
			cerr << "For node " << nodecnt << ", current best objective bound " << objbnd << endl;*/

			numCallbacks++;
			time_t start = clock();

			//get the solution vector (for y variables) from Gurobi 
			//g1 is duplicate of power graph
			double *y = new double[g1.n];
			y = getSolution(vars, g1.n);

			//now, make it and its complement boolean
			vector<bool> DeletedVertices(g1.n, false);
			vector<bool> NonDeletedVertices(g1.n, false);
			for (long i = 0; i < g1.n; i++)
			{
				if (y[i] > 0.5) DeletedVertices[i] = true;
				else NonDeletedVertices[i] = true;
			}

			//get the solution vector (for x variables) from Gurobi
			//g1 is duplicate of power graph
			double *x = new double[g1.m];
			x = getSolution(vars1, g1.m);

			//find a violated length-k i,j-connector inequality (if any exist)
			vector<long> predecessor;
			vector<long> distance;
			for (long u = 0; u < g2.n; u++)
			{
				vector<long>original_distance = g2.ShortestPathsUnweighted(u);
				distance = g2.ShortestPathsUnweighted(u, NonDeletedVertices, predecessor); //g2 is duplicate of original graph.
				for (long v = u + 1; v < g2.n; v++)
				{
					long counter = 0;
					if (original_distance[v] == 1) continue; //all connector inequalities are satisfied
					if (distance[v] <= k1 && x[hashing[u*g2.n + v]] == 0)
					{
						GRBLinExpr expr = vars1[hashing[u*g2.n + v]];
						long i = v;
						while (i != u)
						{
							expr += vars[i];
							i = predecessor[i];
						}
						expr += vars[i];
						addLazy(expr >= 1);
						counter++;
						numLazyCuts++;
						if (counter == 10 * b1) break; //with this we can change number of cuts. 
					}
					predecessor.clear();
				}
			}
		}

		else if (where == GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
		{
			//cerr << "In MIPNODE" << endl;
			numCallbacks++;
			time_t start = clock();


			//get the solution vector (for y variables) from Gurobi 
			//g1 is duplicate of power graph
			//g2 is duplicate of original graph 
			double *y = new double[g1.n];
			y = getNodeRel(vars, g1.n);

			//get the solution vector (for x variables) from Gurobi
			//double *x = getNodeRel(vars1, g1.m);
			double *x = new double[g1.m];
			x = getNodeRel(vars1, g1.m);

			for (long i = 0; i < g2.n; i++)
			{
				double infty = numeric_limits<double>::max();
				vector<double> temp1(k1 + 1, infty);
				vector<double> temp2(k1 + 1, infty);
				vector < vector<double> > d(g2.n, temp1);
				vector< vector<double> > p(g2.n, temp2);

				for (long s = 0; s <= k1; s++)
				{
					d[i][s] = y[i];
					p[i][s] = i;
				}

				for (long s = 1; s <= k1; s++)
				{
					for (long counter1 = 0; counter1 < g1.degree[i]; counter1++)
					{
						long v = g1.adj[i][counter1];
						for (long counter2 = 0; counter2 < g2.degree[v]; counter2++)
						{
							long u = g2.adj[v][counter2];
							if (d[v][s] > d[u][s - 1] + y[v])
							{
								d[v][s] = d[u][s - 1] + y[v];
								p[v][s] = u;
							}
						}
					}
				}

				double min = 1;
				long t = -1;
				for (long counter3 = 0; counter3 < g1.degree[i]; counter3++)
				{
					long j = g1.adj[i][counter3];
					if (j > i) //because of hashing 
					{
						if (x[hashing[i*g2.n + j]] + d[j][k1] < min)
						{
							min = x[hashing[i*g2.n + j]] + d[j][k1];
							t = j;
						}
					}
				}
				if (min < 0.9)
				{
					cerr << "min is " << min << endl;
					GRBLinExpr expr = vars1[hashing[i*g2.n + t]] + vars[t] + vars[i];
					long q = p[t][k1];
					long ss = k1 - 1;
					for (long counter = 2; counter <= k1; counter++)
					{
						expr += vars[q];
						q = p[q][ss];
						ss--;
					}
					addCut(expr >= 1);
					//cerr << "add cut" << endl;
				}
			}
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
vector<long> solveDCNP_thin_formulation(KGraph &g, long k, long b, vector<long> Heuristic_sol, bool &subOpt)
{
	vector<long> Deleted;
	subOpt = true;

	cerr << "creating power graph " << endl;
	KGraph gs = g.CreatePowerGraph(k);
	cerr << "|E^k| is " << gs.m << endl;

	try
	{
		GRBEnv env = GRBEnv();
		//env.set(GRB_IntParam_OutputFlag, 0);
		//env.set(GRB_IntParam_Method, 3); //use barrier method to solve LP relaxation.
		env.set(GRB_DoubleParam_TimeLimit, 14400);
		GRBModel model = GRBModel(env);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0.0);
		model.getEnv().set(GRB_IntParam_LazyConstraints, 1);


		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar *X = model.addVars(gs.m, GRB_CONTINUOUS);
		model.update();

		cerr << "Adding objective function." << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < gs.m; i++)
		{
			Objective += X[i];
		}
		Objective *= 2;
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


		cerr << "Adding constraints for edges." << endl;
		for (long u = 0; u < g.n; u++)
		{
			for (long i = 0; i < g.adj[u].size(); i++)
			{
				long v = g.adj[u][i];
				if (v > u) model.addConstr(X[hash_edges[u*g.n + v]] + Y[u] + Y[v] >= 1);
			}
		}


		long length_2_connecotrs = 0;
		cerr << "Adding constraints for length-2 connectors" << endl;
		for (long u = 0; u < g.n; u++)
		{
			for (long counter1 = 0; counter1 < g.degree[u]; counter1++)
			{
				vector<bool> neighbors(g.n, 0);
				long i = g.adj[u][counter1];
				for (long counter2 = 0; counter2 < g.degree[i]; counter2++)
				{
					neighbors[g.adj[i][counter2]] = true;
				}
				for (long counter3 = 0; counter3 < g.degree[u]; counter3++)
				{
					long j = g.adj[u][counter3];
					if (neighbors[j] == false && j>i)
					{
						model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[j] + Y[u] >= 1);
						//length_2_connecotrs++;
					}
				}
			}
		}
		cerr << "Number of minimal length-2 connecter inequalities is " << length_2_connecotrs << endl;


		vector < vector<long>> all_dist;
		for (long i = 0; i < g.n; i++)
		{
			vector<long> dist_from_i_to = g.ShortestPathsUnweighted(i);
			all_dist.push_back(dist_from_i_to);
		}

		//cerr << "Adding constraints for length-3 connectors" << endl; //what discussed in meeting, all the connectors
		//long length_3_connecotrs = 0;
		//for (long i = 0; i < g.n; i++)
		//{
		//	for (long u = 0; u < g.n; u++) //u=i+1
		//	{
		//		if (all_dist[i][u] == 1)
		//		{
		//			for (long v = 0; v < g.n; v++) //v=u+1
		//			{
		//				if (all_dist[u][v] == 1)
		//				{
		//					for (long j = 0; j < g.n; j++) //j=v+1
		//					{
		//						if (all_dist[v][j] == 1)
		//						{
		//							if (all_dist[i][v] == 2 & all_dist[u][j] == 2 &&j>i)
		//							{
		//								model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[u] + Y[v] + Y[j] >= 1);
		//								length_3_connecotrs++;
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}
		//	}
		//}
		//cerr << "Number of minimal length-3 connecter inequalities is " << length_3_connecotrs << endl;


		cerr << "Adding some of the length-3 i,j-connectors" << endl;
		for (long u = 0; u < g.n; u++)
		{
			for (long v = u + 1; v < g.n; v++)
			{
				if (all_dist[u][v] == 1)
				{
					for (long w = v + 1; w < g.n; w++)
					{
						if (all_dist[u][w] == 2 & all_dist[v][w] == 1)
						{
							for (long x = w + 1; x < g.n; x++)
							{
								if (all_dist[u][x] == 3 & all_dist[v][x] == 2 & all_dist[w][x] == 1 & x > u)
								{
									model.addConstr(X[hash_edges[u*g.n + x]] + Y[u] + Y[v] + Y[w] + Y[x] >= 1);
								}
							}
						}
					}
				}
			}
		}

		cerr << "Adding budget constraints" << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
		{
			expr1 += Y[i];
		}
		model.addConstr(expr1 <= b);
		model.update();

		cerr << "Adding Lazy Constraints." << endl;
		integer_separation cb = integer_separation(X, hash_edges, Y, gs, g, k, b);
		model.setCallback(&cb);


		cerr << "Providing Initial Solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}
		vector<bool> S(g.n, true);
		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			S[Heuristic_sol[i]] = false;
		}

		for (long i = 0; i < g.n; i++) //i is a vertex in original graph.
		{
			vector <long> dist_from_i = g.ShortestPathsUnweighted(i, S);
			for (long j = i + 1; j < g.n; j++) //j is another vertex in original graph.
			{
				if (dist_from_i[j] <= k) //i and j are close (dist <= s) in original graph.
				{
					X[hash_edges[i*g.n + j]].set(GRB_DoubleAttr_Start, 1);
				}
			}
		}

		vector<long> Rmap1;
		vector<long> Rmap2;
		vector<long> Neighbors_of_i;
		vector<long> simplicial;
		for (long i = 0; i < g.n; i++)
		{
			Neighbors_of_i.push_back(i);
			for (long counter = 0; counter < g.degree[i]; counter++)
			{
				Neighbors_of_i.push_back(g.adj[i][counter]);
			}
			sort(Neighbors_of_i.begin(), Neighbors_of_i.end());
			KGraph gs = g.CreateInducedGraph(Neighbors_of_i, Rmap1);
			if (gs.m == (gs.n*(gs.n - 1) / 2))
			{
				simplicial.push_back(i);
			}
			Neighbors_of_i.clear();
		}
		sort(simplicial.begin(), simplicial.end());
		KGraph gnew = g.CreateInducedGraph(simplicial, Rmap2);

		vector<long> final;
		vector< vector< long> > components;
		vector<long> degreeZero;
		gnew.FindConnectedComponents(components, degreeZero);

		/*cerr << "Number of components is " << components.size()+degreeZero.size() << endl;
		cerr << "vertices with degree equals 0 are ";*/
		for (long i = 0; i < degreeZero.size(); i++)
		{
			//cerr << simplicial[Rmap2[simplicial[degreeZero[i]]]] << " ";
			final.push_back(simplicial[Rmap2[simplicial[degreeZero[i]]]]);
		}
		//cerr << "\n";

		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				//cerr << "Here, we have " << simplicial[Rmap2[simplicial[components[i][j]]]] << " ";
				final.push_back(simplicial[Rmap2[simplicial[components[i][j]]]]);
				break;
			}
		}

		//cerr << "final size is " << final.size() << endl;
		for (long i = 0; i < final.size(); i++)
		{
			model.addConstr(Y[final[i]] == 0);
		}



		cerr << "Optimizing." << endl;
		model.optimize();

		long bestUB = model.get(GRB_DoubleAttr_ObjVal);
		long bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";



		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
			subOpt = false;
		else return Deleted;


		for (long i = 0; i < g.n; i++)
			if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
				Deleted.push_back(i);

		delete[] Y;
		//for (long i = 0; i < gs.m; i++)
		delete[] X;


		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B Nodes : " << NumOfBandBNodes << endl;
		cerr << "Number of callbacks is: " << integer_separation::numCallbacks << endl;
		cerr << "Number of Lazy Cuts is: " << integer_separation::numLazyCuts << endl;
		cerr << "Time spent in callback is: " << integer_separation::TotalCallbackTime << endl;

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



/*The following function uses O(|E^k|) variables. We use hashing for power graph edges. Fractional separation is performed.*/
vector<long> solveDCNP_thin_formulation_fractional(KGraph &g, long k, long b, vector<long> Heuristic_sol, bool &subOpt)
{
	vector<long> Deleted;
	subOpt = true;

	cerr << "creating power graph " << endl;
	KGraph gs = g.CreatePowerGraph(k);
	//cerr << "|E^k| is " << gs.m << endl;

	try
	{
		GRBEnv env = GRBEnv();
		//env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0.0);
		model.getEnv().set(GRB_IntParam_LazyConstraints, 1);


		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar *X = model.addVars(gs.m, GRB_CONTINUOUS);
		model.update();

		cerr << "Adding objective function." << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < gs.m; i++)
		{
			Objective += X[i];
		}
		Objective *= 2;
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

		cerr << "Adding constraints for edges." << endl;
		for (long u = 0; u < g.n; u++)
		{
			for (long i = 0; i < g.adj[u].size(); i++)
			{
				long v = g.adj[u][i];
				if (v > u) model.addConstr(X[hash_edges[u*g.n + v]] + Y[u] + Y[v] >= 1);
			}
		}


		//long length_2_connecotrs = 0;
		//cerr << "Adding constraints for length-2 connectors" << endl;
		//for (long u = 0; u < g.n; u++)
		//{
		//	for (long counter1 = 0; counter1 < g.degree[u]; counter1++)
		//	{
		//		vector<bool> neighbors(g.n, 0);
		//		long i = g.adj[u][counter1];
		//		for (long counter2 = 0; counter2 < g.degree[i]; counter2++)
		//		{
		//			neighbors[g.adj[i][counter2]] = true;
		//		}
		//		for (long counter3 = 0; counter3 < g.degree[u]; counter3++)
		//		{
		//			long j = g.adj[u][counter3];
		//			if (neighbors[j] == false && j > i)
		//			{
		//				model.addConstr(X[hash_edges[i*g.n + j]] + Y[i] + Y[j] + Y[u] >= 1);
		//				length_2_connecotrs++;
		//			}
		//		}
		//	}
		//}
		//cerr << "Number of minimal length-2 connecter inequalities is " << length_2_connecotrs << endl;

		//cerr << "Adding some of the length-3 i,j-connectors" << endl;
		//   vector < vector<long>> all_dist;
		//for (long i = 0; i < g.n; i++)
		//{
		//	vector<long> dist_from_i_to = g.ShortestPathsUnweighted(i);
		//	all_dist.push_back(dist_from_i_to);
		//}
		//for (long u = 0; u < g.n; u++)
		//{
		//	for (long v = u+1; v < g.n; v++) //v=0
		//	{
		//		if (all_dist[u][v] == 1)
		//		{
		//			for (long w = v+1; w < g.n; w++) //w=0
		//			{
		//				if (all_dist[u][w] == 2 & all_dist[v][w] == 1)
		//				{
		//					for (long x = w+1; x < g.n; x++) //x=0 
		//					{
		//						if (all_dist[u][x] == 3 & all_dist[v][x] == 2 & all_dist[w][x] == 1 & x > u)
		//						{
		//							model.addConstr(X[hash_edges[u*g.n + x]] + Y[u] + Y[v] + Y[w] + Y[x] >= 1);
		//						}
		//					}
		//				}
		//			}
		//		}
		//	}
		//}


		cerr << "Adding budget constraints" << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
		{
			expr1 += Y[i];
		}
		model.addConstr(expr1 <= b);
		model.update();

		cerr << "Adding Lazy Constraints." << endl;
		fractional_separation cb = fractional_separation(X, hash_edges, Y, gs, g, k, b);
		model.setCallback(&cb);



		cerr << "Providing Initial Solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}
		vector<bool> S(g.n, true);
		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			S[Heuristic_sol[i]] = false;
		}

		for (long i = 0; i < g.n; i++) //i is a vertex in original graph.
		{
			vector <long> dist_from_i = g.ShortestPathsUnweighted(i, S);
			for (long j = i + 1; j < g.n; j++) //j is another vertex in original graph.
			{
				if (dist_from_i[j] <= k) //i and j are close (dist <= s) in original graph.
				{
					X[hash_edges[i*g.n + j]].set(GRB_DoubleAttr_Start, 1);
				}
			}
		}



		cerr << "Optimizing." << endl;
		model.optimize();

		long bestUB = model.get(GRB_DoubleAttr_ObjVal);
		long bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
			subOpt = false;
		else return Deleted;


		for (long i = 0; i < g.n; i++)
			if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
				Deleted.push_back(i);

		delete[] Y;
		delete[] X;


		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B Nodes : " << NumOfBandBNodes << endl;
		cerr << "Number of callbacks is: " << fractional_separation::numCallbacks << endl;
		cerr << "Number of Lazy Cuts is: " << fractional_separation::numLazyCuts << endl;
		cerr << "Time spent in callback is: " << fractional_separation::TotalCallbackTime << endl;

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


vector<long> solveDCNP_path_like(KGraph &g, long k, long b, vector<long> Heuristic_sol, bool &subOpt)
{
	vector<long> Deleted;
	subOpt = true;

	cerr << "creating power graph " << endl;
	KGraph gk = g.CreatePowerGraph(k);

	map<vector<long>, long, classcomp>map = EnumerateLength3Connector(g);

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0.0);
		//model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		GRBVar **X = new GRBVar*[g.n];
		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar *Z = model.addVars(map.size(), GRB_BINARY);

		for (long i = 0; i < g.n; i++)
		{
			X[i] = model.addVars(g.n, GRB_BINARY);
		}
		model.update();


		cerr << "Adding objective function." << endl;
		for (long i = 0; i < g.n; i++)
		{
			for (long j = 0; j < g.n; j++)
			{
				X[i][j].set(GRB_DoubleAttr_Obj, 1);
			}
		}
		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.update();

		cerr << "Adding constraint z_C + \sum_{v \in C}{y_v} >= 1" << endl;
		for (long v = 0; v < g.n; v++)
		{
			vector<long> dist_from_v_to = g.ShortestPathsUnweighted(v);
			for (long u = v + 1; u < g.n; u++)
			{
				if (dist_from_v_to[u] == 1)
				{
					long minuv = min(u, v);
					long maxuv = max(u, v);
					std::map<std::vector<long>, long>::iterator it = map.find({ minuv, maxuv });
					model.addConstr(Z[it->second] + Y[u] + Y[v] >= 1);
				}

				if (dist_from_v_to[u] == 2 || dist_from_v_to[u] == 3)
				{
					vector<long> dist_from_u_to = g.ShortestPathsUnweighted(u);
					for (long v_neighbors_iterator = 0; v_neighbors_iterator < g.adj[v].size(); v_neighbors_iterator++)
					{
						long ii = g.adj[v][v_neighbors_iterator];
						if (dist_from_u_to[ii] == 1)
						{
							vector<long> sortedsubset = sortnodes(u, v, ii);
							std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
							model.addConstr(Z[it->second] + Y[u] + Y[v] + Y[ii] >= 1);
						}

						if (dist_from_u_to[ii] == 2)
						{
							for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.adj[ii].size(); ii_neighbors_iterator++)
							{
								long jj = g.adj[ii][ii_neighbors_iterator];
								if (dist_from_u_to[jj] == 1 && dist_from_v_to[jj] == 2)
								{
									vector<long> sortedsubset = sort4nodes(u, v, ii, jj);
									std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
									model.addConstr(Z[it->second] + Y[u] + Y[v] + Y[ii] + Y[jj] >= 1);
								}
							}
						}
					}
				}
			}
		}

		cerr << "Adding x_[i][j] >= Z_c constraints." << endl;
		for (long v = 0; v < g.n; v++)
		{
			vector<long> dist_from_v_to = g.ShortestPathsUnweighted(v);
			for (long u = v + 1; u < g.n; u++)
			{
				if (dist_from_v_to[u] > 3) continue;
				if (dist_from_v_to[u] == 1)
				{
					long minuv = min(u, v);
					long maxuv = max(u, v);
					std::map<std::vector<long>, long>::iterator it = map.find({ minuv, maxuv });
					model.addConstr(X[u][v] >= Z[it->second]);
				}
				if (dist_from_v_to[u] == 2 || dist_from_v_to[u] == 3)
				{
					vector<long> dist_from_u_to = g.ShortestPathsUnweighted(u);
					for (long v_neighbors_iterator = 0; v_neighbors_iterator < g.adj[v].size(); v_neighbors_iterator++)
					{
						long ii = g.adj[v][v_neighbors_iterator];
						if (dist_from_u_to[ii] == 1)
						{
							vector<long> sortedsubset = sortnodes(u, v, ii);
							std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
							model.addConstr(X[u][v] >= Z[it->second]);
						}
						if (dist_from_u_to[ii] == 2)
						{
							for (long ii_neighbors_iterator = 0; ii_neighbors_iterator < g.adj[ii].size(); ii_neighbors_iterator++)
							{
								long jj = g.adj[ii][ii_neighbors_iterator];
								if (dist_from_u_to[jj] == 1 && dist_from_v_to[jj] == 2)
								{
									vector<long> sortedsubset = sort4nodes(u, v, ii, jj);
									std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
									model.addConstr(X[u][v] >= Z[it->second]);
								}
							}
						}
					}
				}
			}
		}

		cerr << "Adding X[i][j] == X[j][i] constraints." << endl;
		for (long i = 0; i < g.n; i++)
		{
			for (long j = 0; j < g.n; j++)
			{
				model.addConstr(X[i][j] == X[j][i]);
			}
		}

		cerr << "Adding budget constraints" << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
		{
			expr1 += Y[i];
		}
		model.addConstr(expr1 <= b);
		model.update();



		cerr << "Providing Initial Solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 1);
		}
		vector<bool> S(g.n, true);
		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			S[Heuristic_sol[i]] = false;
		}
		for (long i = 0; i < g.n; i++) //i is a vertex in original graph.
		{
			vector <long> dist_from_i = g.ShortestPathsUnweighted(i, S);
			for (long j = i + 1; j < g.n; j++) //j is another vertex in original graph.
			{
				if (dist_from_i[j] <= k) //i and j are close (dist <= s) in original graph.
				{
					for (long q = 0; q < gk.degree[i]; q++) //q is a counter for neighbors of i in power graph. 
					{
						//If qth neighbor of i is j in power graph, then X[i][q] = 1.
						if (gk.adj[i][q] == j) X[i][q].set(GRB_DoubleAttr_Start, 1);
					}
				}
			}
		}


		cerr << "Optimizing." << endl;
		model.optimize();

		long bestUB = model.get(GRB_DoubleAttr_ObjVal);
		long bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
			subOpt = false;
		else return Deleted;

		for (long i = 0; i < g.n; i++)
			if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
				Deleted.push_back(i);

		delete[] Y;
		for (long i = 0; i < gk.n; i++)
		{
			delete[] X[i];
		}
		delete[] X;


		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B Nodes : " << NumOfBandBNodes << endl;

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

vector<long> solveDCNP_Veremyev(KGraph &g, long k, long b, vector<long> Heuristic_sol, bool &subOpt)
{
	vector<long> Deleted;
	subOpt = true;

	try
	{
		GRBEnv env = GRBEnv();
		//env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3);  //use barrier method to solve LP relaxation.
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0.0);

		time_t Building_Model_Start = clock();
		GRBVar *Y = model.addVars(g.n, GRB_BINARY);
		GRBVar ***U = new GRBVar**[g.n]; //U[i][j][s] denotes whether distance between i and j is at most s in
		for (long i = 0; i<g.n; i++)	 // the interdicted graph
		{
			GRBVar **U_temp = new GRBVar*[g.n];
			for (long j = 0; j < g.n; j++)
			{
				U_temp[j] = model.addVars(k, GRB_BINARY);
			}
			U[i] = U_temp;
		}
		model.update();

		vector<long> dist_from_i_to;
		vector < vector<long>> all_dist;
		for (long i = 0; i < g.n; i++)
		{
			dist_from_i_to = g.ShortestPathsUnweighted(i);
			all_dist.push_back(dist_from_i_to);
		}


		cerr << "Adding objective function." << endl;
		GRBLinExpr Objective = 0;
		for (long i = 0; i < g.n; i++)
		{
			for (long j = i + 1; j < g.n; j++)
			{
				Objective += U[i][j][k - 1];
			}
		}
		Objective *= 2;
		model.setObjective(Objective, GRB_MINIMIZE);
		model.update();


		//Fixing variables
		cerr << "Fixing variables." << endl;
		for (long i = 0; i < g.n; i++)
		{
			for (long j = 0; j != i && j < g.n; j++)
			{
				if (all_dist[i][j] > k)
				{
					for (long s = 0; s < k; s++)
					{
						model.addConstr(U[i][j][s] == 0);
					}
				}
				else
				{
					for (long s = 0; s < all_dist[i][j] - 1; s++)
					{
						model.addConstr(U[i][j][s] == 0);
					}
				}
			}
		}


		//Adding budget constraint.
		cerr << "Adding constraint (1b)." << endl;
		GRBLinExpr expr1 = 0;
		for (long i = 0; i < g.n; i++)
		{
			expr1 += Y[i];
		}
		model.addConstr(expr1 <= b);


		//U[i][j][q] equals 1 iff i is adjacent to j and both of them are not deleted.
		cerr << "Adding constraint (1c)." << endl;
		for (long i = 0; i < g.n; i++)
		{
			for (long u = 0; u < g.degree[i]; u++) //u is a counter of neighbors of i.
			{
				long j = g.adj[i][u];
				model.addConstr(U[i][j][0] + Y[i] + Y[j] >= 1);
			}
		}

		//Enforce U[i][j][q] to be zero if i or j is deleted.
		cerr << "Adding constraint (1d)." << endl;
		for (long i = 0; i < g.n; i++)
		{
			//vector<long> dist_from_i_to = g.ShortestPathsUnweighted(i);
			for (long j = 0; j != i && j < g.n; j++)
			{
				//long q = dist_from_i_to[j];
				for (long q = all_dist[i][j]; q < k; q++) //q is like l in the original formulation.
				{
					model.addConstr(U[i][j][q] + Y[i] <= 1);
					model.addConstr(U[i][j][q] + Y[j] <= 1);
				}
			}
		}


		//U[i][j][q] = U[i][j][0]
		cerr << "Adding constraint (1e)." << endl;
		for (long i = 0; i < g.n; i++)
		{
			for (long u = 0; u < g.degree[i]; u++) //u is a counter of neighbors of i.
			{
				long j = g.adj[i][u];
				for (long q = all_dist[i][j]; q < k; q++)
				{
					model.addConstr(U[i][j][q] == U[i][j][0]);
				}
			}
		}

		cerr << "Adding constraints (1f) & (1g)." << endl;
		for (long i = 0; i < g.n; i++)
		{
			for (long q = 1; q < k; q++)
			{
				GRBLinExpr expr2 = 0;
				for (long j = 0; j < g.n; j++)
				{
					if (k >= all_dist[i][j] && all_dist[i][j] > 1)
					{
						for (long p = 0; p < g.degree[i]; p++)
						{
							long t = g.adj[i][p];
							model.addConstr(U[t][j][q - 1] <= U[i][j][q] + Y[i]); //(1g)
							expr2 += U[t][j][q - 1];
						}
						model.addConstr(U[i][j][q] <= expr2); //(1f)
					}
				}
			}
		}

		cerr << "Adding constraint (1h)." << endl;
		for (long q = 0; q < k; q++)
		{
			for (long i = 0; i < g.n; i++)
			{
				for (long j = 0; j != i && j < g.n; j++)
				{
					model.addConstr(U[i][j][q] == U[j][i][q]);
				}
			}
		}
		model.update();
		cerr << "Time spent for building the model is: " << (double)(clock() - Building_Model_Start) / CLOCKS_PER_SEC << endl;


		cerr << "Providing Initial Solution" << endl;
		for (long i = 0; i<g.n; i++)
			Y[i].set(GRB_DoubleAttr_Start, 1);

		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			long v = Heuristic_sol[i];
			Y[v].set(GRB_DoubleAttr_Start, 0);
		}
		vector<bool> S(g.n, true);
		for (long i = 0; i < Heuristic_sol.size(); i++)
		{
			S[Heuristic_sol[i]] = false;
		}
		for (long v = 0; v < g.n; v++)
		{
			vector <long> dist_from_v = g.ShortestPathsUnweighted(v, S);
			for (long w = v + 1; w < g.n; w++)
			{
				for (long q = 0; q < k; q++)
				{
					U[v][w][q].set(GRB_DoubleAttr_Start, 0);
				}
				for (long q = dist_from_v[w] - 1; q < k; q++)
				{
					U[v][w][q].set(GRB_DoubleAttr_Start, 1);
				}
			}
		}



		cerr << "Optimizing." << endl;
		model.optimize();

		long bestUB = model.get(GRB_DoubleAttr_ObjVal);
		long bestLB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
			subOpt = false;
		else return Deleted;

		for (long i = 0; i < g.n; i++)
			if (Y[i].get(GRB_DoubleAttr_X) > 0.5)
				Deleted.push_back(i);

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



//Based on algorithm in paper "A faster algorithm for betweenness centerality".
vector<long> BC(KGraph &g, long b)
{

	vector<float> C_B(g.n, 0);
	for (long s = 0; s < g.n; s++)
	{
		vector<long> Stack;
		vector< vector<long> > P(g.n);
		vector<float> sigma(g.n, 0); //sigma[t] = number of shortest paths from s (source) to t
		sigma[s] = 1;
		vector<long> d(g.n, -1); //d[t] = distance from s to t
		d[s] = 0;
		vector<long> Q; //A queue
		Q.push_back(s);

		while (!Q.empty())
		{
			long v = Q[0];
			Q.erase(Q.begin());
			Stack.push_back(v);

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
				if (d[w] == d[v] + 1 /*&& d[w] <=3*/)
				{
					sigma[w] = sigma[w] + sigma[v];
					P[w].push_back(v);
				}
			}
		}

		vector<float> delta(g.n, 0.0);
		//Stack returns vertices in order of non-increasing distance from s
		while (!Stack.empty())
		{

			long w = Stack.back();
			Stack.pop_back();

			for (long j = 0; j < P[w].size(); j++)
			{
				long v = P[w][j];
				delta[v] = float(delta[v] + float((sigma[v] / sigma[w]))*(1 + delta[w]));
			}

			if (w != s)
			{
				C_B[w] = C_B[w] + delta[w];
			}
		}
	}


	vector<long> top2b_BC; //It saves top 2b vertices in terms of the most BC values.

	vector<float> C_B_copy = C_B;
	sort(C_B_copy.rbegin(), C_B_copy.rend()); //C_B_copy is sorted

	vector<float> top2B;
	for (long i = 0; i < 2 * b; i++)
	{
		top2B.push_back(C_B_copy[i]);
	}

	long counter = 0;
	for (long j = 0; j < C_B.size() && counter <= 2 * b; j++)
	{
		for (long i = 0; i < top2B.size() && counter <= 2 * b; i++)
		{
			if (C_B[j] == top2B[i])
			{
				top2b_BC.push_back(j);
				counter++;
				break;
			}
		}
	}
	return top2b_BC;
}



vector<long> Greedy_Heuristic(KGraph &g, long k, long b)
{

	long N = 1; //number of total iterations of algorithm
	vector<long> D = BC(g, b);  
	vector<long> D_Star;
	float q_star = INFINITY;  
	vector<long> D_copy = D; //D_copy is just a copy of D
	long delta_b = b / 2;
	long obj = 0;

	do
	{
		while (D.size() > b - delta_b)
		{
			if (D.size() == b & obj <= q_star)
			{
				D_Star = D;
				q_star = obj;
				N++;
			}

			vector<long> close_vertices;
			vector<bool> current_nondeleted_vertices(g.n, true);

			//Deleting all vertices of set D from graph
			for (long i = 0; i < D.size(); i++)
			{
				current_nondeleted_vertices[D[i]] = false;
			}

			for (long j = 0; j < D.size(); j++)
			{
				//Adding back the vertices of set D to graph G one by one
				current_nondeleted_vertices[D[j]] = true;

				for (long v = 0; v < g.n; v++)
				{
					vector <long> dist_from_v = g.ShortestPathsUnweighted(v, current_nondeleted_vertices);
					for (long w = v + 1; w < g.n; w++)
					{
						if (dist_from_v[w] <= k)
						{
							obj++;
						}
					}
				}
				close_vertices.push_back(obj);
				obj = 0;
				current_nondeleted_vertices[D[j]] = false;
			}

			long min_size = *min_element(close_vertices.begin(), close_vertices.end());
			for (long i = 0; i < D.size(); i++)
			{
				if (close_vertices[i] == min_size)
				{
					D.erase(D.begin() + i);
				}
			}
		}


		while (D.size() < b + delta_b)
		{
			if (D.size() == b & obj <= q_star)
			{
				D_Star = D;
				q_star = obj;
				N++;
			}

			vector<long>far_apart;
			long num_far_apart_vertices = 0;

			for (long j = 0; j < D_copy.size(); j++)
			{
				vector<bool> new_nodes(g.n, true);
				for (long i = 0; i < D.size(); i++)
				{
					new_nodes[D[i]] = false;
				}
				new_nodes[D_copy[j]] = false;

				for (long v = 0; v < g.n; v++)
				{
					vector <long> dist_from_v = g.ShortestPathsUnweighted(v, new_nodes);
					for (long w = v + 1; w < g.n; w++)
					{
						if (dist_from_v[w] > k)
						{
							num_far_apart_vertices++;
						}
					}
				}
				far_apart.push_back(num_far_apart_vertices);
				num_far_apart_vertices = 0;
			}

			long max_size = *max_element(far_apart.begin(), far_apart.end());
			for (long i = 0; i < D_copy.size(); i++)
			{
				if (far_apart[i] == max_size)
				{
					D.push_back(D_copy[i]);
				}
			}

			vector<bool> remained_nodes(g.n, true);
			for (long i = 0; i < D.size(); i++)
			{
				remained_nodes[D[i]] = false;
			}

			obj = 0;

			for (long v = 0; v < g.n; v++)
			{
				vector <long> dist_from_v = g.ShortestPathsUnweighted(v, remained_nodes);
				for (long w = v + 1; w < g.n; w++)
				{
					if (dist_from_v[w] <= k)
					{
						obj++;
					}
				}
			}
		}

	} while (N <= 1);

	return D_Star;
}
