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

//Preprocessing function
vector<long> Preprocessing(KGraph &g);

//Betweenneess Centrality (BC) function
vector<long> BC(KGraph &g, long B);

//Greedy4 algorithm
vector<long> Greedy_Heuristic(KGraph &g, long s, long B);


//Thin formulation using power graph and hashing
vector<long> solveDCNP_thin_formulation(KGraph &g, long k, long b, vector<long> Heuristic_sol, bool &subOpt);

//Thin formulation using power graph and hashing for weighted instances
vector<long> solveDCNP_thin_formulation_weighted(KGraph &g, long k, long b, vector<long> Heuristic_sol, bool &subOpt);

//Thin formulation using power graph and hashing with fractional seoaration
vector<long> solveDCNP_thin_formulation_fractional(KGraph &g, long k, long b, vector<long> Heuristic_sol, bool &subOpt);

//To solve DCNP with path-like formulation, using power graph 
vector<long> solveDCNP_path_like(KGraph &g, long k, long b, vector<long> Heuristic_sol, bool &subOpt);


//Veremyev DCNP function
vector<long> solveDCNP_Veremyev(KGraph &g, long s, long B, vector<long> Heuristic_sol, bool &subOpt);

// callback function for DCNP


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
	static double TotalCallbackTime;
	static long numLazyCutsInteger;
	static long numLazyCutsFractional;
};


#endif
