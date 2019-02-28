#ifndef GRBINTERFACE_H
#define GRBINTERFACE_H
#include "gurobi_c++.h"
#include "KGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include<climits>
#include <unordered_map>

string itos(int i);

//Check if a vertex is simplicial
bool IsSimplicial(KGraph &g, long i);

//Calculate number of vertex pairs with distance at most k in graph G-D (obj(G-D)) where G is unweighted
long obj(KGraph &g, vector<long> deleted, long k);
long obj(KGraph &g, vector<bool> nondeleted, long k);

//Calculate number of vertex pairs with distance at most k in graph G-D (obj(G-D)) where G is weighted
long obj_weighted(KGraph &g, vector<long> deleted, long k);

//Data structure to store variables d(v,s) and p(v,s) for fractional separation
struct d_and_p
{
	vector < vector<double> > d;
	vector< vector<double> > p;
};
//Function to calculate variables d(v,s) and p(v,s) for fractional separation
d_and_p d_and_p_function (KGraph &goriginal, KGraph &gpower, long i, long k, double *y);


//Preprocessing function
vector<long> Preprocessing(KGraph &g);

//Betweenneess Centrality function
vector<long> FindTopTBetweenessCentralityNodes(KGraph &g, long T);

//DCNP heuristic
vector<long> DCNP_Heuristic(KGraph &g, long s, long B);


//Thin formulation using power graph and hashing
vector<long> solveDCNP_thin_formulation(KGraph &g, long k, long b, vector<long> &Heuristic_sol);

//Thin formulation using power graph and hashing for weighted instances
vector<long> solveDCNP_thin_formulation_weighted(KGraph &g, long k, long b, vector<long> &Heuristic_sol);

//Thin formulation using power graph and hashing with fractional seoaration
vector<long> solveDCNP_thin_formulation_fractional(KGraph &g, long k, long b, vector<long> &Heuristic_sol);

//To solve DCNP with path-like formulation, using power graph 
vector<long> solveDCNP_path_like_k3(KGraph &g, long b, vector<long> &Heuristic_sol);


//Veremyev DCNP function
vector<long> solveDCNP_Veremyev(KGraph &g, long k, long b, vector<long> &Heuristic_sol);



// callback functions for DCNP


//Integer separation when we have O(|E^k|) variables and we use hashing
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
	vector < vector<long>> all_distances;

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
};


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



//Fractional separation when we have O(E|G^k|) variables and we use hashing
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

	fractional_separation(GRBVar *xvars, unordered_map<long, long> hash_edges, GRBVar *yvars, KGraph &gs, KGraph &g, long k, long b)
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
	static double TotalCallbackTimeInteger;
	static double TotalCallbackTimeFractional;
	static long numLazyCutsInteger;
	static long numLazyCutsFractional;
};




#endif
