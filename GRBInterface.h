#ifndef GRBINTERFACE_H
#define GRBINTERFACE_H
#include "gurobi_c++.h"
#include "KGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include<climits>
#include <unordered_map>

//Boolify function
vector<bool> boolify(vector<long> &S, long n);
vector<bool> boolify(double *y, long n);

//Printing function
void PrintVectorLong(vector<long>&S);

/*DCNP preprocessing function
*It works for arbitrary k and b*/
vector<long> FindNonCriticalNodes(KGraph &g);

//Check if a vertex is simplicial
bool IsSimplicial(KGraph &g, long i);

//Calculate number of vertex pairs with distance at most k in graph G-D (obj(G-D)) where G is unweighted
long obj(KGraph &g, vector<long> &deleted, long k);
long obj(KGraph &g, vector<bool> &nondeleted, long k);

//Calculate number of vertex pairs with distance at most k in graph G-D (obj(G-D)) where G is weighted
long obj_weighted(KGraph &g, vector<long> deleted, long k);

//Data structure to store variables d(v,s) and p(v,s) for fractional separation
struct d_and_p
{
	vector < vector<double> > d;
	vector< vector<double> > p;
};
//Function to calculate variables d(v,s) and p(v,s) for fractional separation
d_and_p d_and_p_function(KGraph &goriginal, KGraph &gpower, long i, long k, double *y);


//Betweenneess Centrality function
vector<long> FindTopTBetweennessCentralityNodes(KGraph &g, long T);

//DCNP Heuristic to find distance-based critical nodes
vector<long> DCNP_Heuristic(KGraph &g, long s, long B);

/*To solve DCNP when distances are measured in terms of hops
* Thin formulation with integer separation is used. */
vector<long> solveDCNP_thin_formulation(KGraph &g, long k, long b, vector<long> &Heuristic_sol);

/*To solve DCNP in edge-weighted graphs.
* Thin formulation with integer separation is used.*/
vector<long> solveDCNP_thin_formulation_weighted(KGraph &g, long k, long b, vector<long> &Heuristic_sol);

/*To solve DCNP when distances are measured in terms of hops
* Thin formulation with fractional separation is used. */
vector<long> solveDCNP_thin_formulation_fractional(KGraph &g, long k, long b, vector<long> &Heuristic_sol);

/*To solve DCNP when distances are measured in terms of hops and k=3
* Path-like formulation is used */
vector<long> solveDCNP_path_like_k3(KGraph &g, long b, vector<long> &Heuristic_sol);

/*To solve DCNP when distances are measured in terms of hops
* Recursive formulation is used. */
vector<long> solveDCNP_Veremyev(KGraph &g, long k, long b, vector<long> &Heuristic_sol);


/*** callback functions for DCNP***/

//Integer separation for unweighted graphs
class integer_separation : public GRBCallback
{
public:
	GRBVar *vars;
	GRBVar *vars1;
	KGraph g1;
	KGraph g2;
	long k1;
	long b1;
	unordered_map<long, long> hashing;
	double *y;
	double *x;


	integer_separation(GRBVar *xvars, unordered_map<long, long> hash_edges, GRBVar *yvars, KGraph &gs, KGraph &g, long k, long b)
	{
		vars = yvars;
		vars1 = xvars;
		g1.Duplicate(gs);
		g2.Duplicate(g);
		k1 = k;
		b1 = b;
		hashing = hash_edges;
	}
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCutsInteger;

	void populate_y()
	{
		y = new double[g1.n];
		y = getSolution(vars, g1.n);
	}

	void populate_x()
	{
		x = new double[g1.m];
		x = getSolution(vars1, g1.m);
	}
};

//Integer separation for edge-weighted graphs
class integer_separation_weighted : public GRBCallback
{
public:
	GRBVar *vars;
	GRBVar *vars1;
	KGraph g1;
	KGraph g2;
	long k1;
	long b1;
	unordered_map<long, long> hashing;

	integer_separation_weighted(GRBVar *xvars, unordered_map<long, long> hash_edges, GRBVar *yvars, KGraph &gs, KGraph &g, long k, long b)
	{
		vars = yvars;
		vars1 = xvars;
		g1.DuplicateForWeighted(gs);
		g2.DuplicateForWeighted(g);
		k1 = k;
		b1 = b;
		hashing = hash_edges;
	}
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCutsInteger;
};


//Fractional separation for unweighted graphs 
class fractional_separation : public GRBCallback
{
public:
	GRBVar *vars;
	GRBVar *vars1;
	KGraph g1;
	KGraph g2;
	long k1;
	long b1;
	unordered_map<long, long> hashing;
	double *y;
	double *x;


	fractional_separation(GRBVar *xvars, unordered_map<long, long> hash_edges, GRBVar *yvars, KGraph &gs, KGraph &g, long k, long b)
	{
		vars = yvars;
		vars1 = xvars;
		//g1 = gs;
		//g2 = g;
		g1.Duplicate(gs);
		g2.Duplicate(g);
		k1 = k;
		b1 = b;
		hashing = hash_edges;
	}
	void callback();
	static long numCallbacks;
	static double TotalCallbackTimeInteger;
	static double TotalCallbackTimeFractional;
	static long numLazyCutsInteger;
	static long numLazyCutsFractional;

	void populate_y()
	{
		y = new double[g1.n];
		y = getSolution(vars, g1.n);
	}

	void populate_x()
	{
		x = new double[g1.m];
		x = getSolution(vars1, g1.m);
	}

	void populate1_y()
	{
		y = new double[g1.n];
		y = getNodeRel(vars, g1.n);
	}
	void populate1_x()
	{
		x = new double[g1.m];
		x = getNodeRel(vars1, g1.m);
	}

};



#endif
