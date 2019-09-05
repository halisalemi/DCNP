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
long obj_weighted(KGraph &g, vector<long> &deleted, long k);
long obj_weighted(KGraph &g, vector<bool> &nondeleted, long k);

//Data structure to store variables d(v,s) and p(v,s) for fractional separation
struct d_and_p
{
	vector < vector<double> > d;
	vector< vector<long> > p;
};
//Function to calculate variables d(v,s) and p(v,s) for fractional separation
d_and_p d_and_p_function(KGraph &goriginal, KGraph &gpower, long i, long k, double *y);


//Betweenneess Centrality function
vector<long> FindTopTBetweennessCentralityNodes(KGraph &g, long T);

//Betweenneess Centrality function for weighted graphs
vector<long> FindTopTBetweennessCentralityNodesWeighted(KGraph &g, long T, long k);

//DCNP Heuristic to find distance-based critical nodes
vector<long> DCNP_Heuristic(KGraph &g, long s, long B);

//DCNP Heuristic to find distance-based critical nodes in weighted graphs
vector<long> DCNP_Heuristic_Weighted(KGraph &g, long s, long B);

/*To solve DCNP when distances are measured in terms of hops
* Thin formulation with integer separation is used. */
vector<long> Thin_I(KGraph &g, long k, long b, vector<long> &Heuristic_sol, bool &subopt, vector<bool> &Initial);

/*To solve DCNP when distances are measured in terms of hops
* Thin formulation with fractional separation is used. */
vector<long> Thin_F(KGraph &g, long k, long b, vector<long> &Heuristic_sol, bool &subopt, vector<bool> &Initial);

/*To solve DCNP when distances are measured in terms of hops and k=3
* Path-like formulation is used */
vector<long> Path_like_k3(KGraph &g, long b, vector<long> &Heuristic_sol, bool &subopt);

/*To solve DCNP when distances are measured in terms of hops and k=4
* Path-like formulation is used */
vector<long> Path_like_k4(KGraph &g, long b, vector<long> &Heuristic_sol, bool &subopt);

/*To solve DCNP when distances are measured in terms of hops
* Recursive formulation is used. */
vector<long> Recursive(KGraph &g, long k, long b, vector<long> &Heuristic_sol, bool &subopt);

/*To solve DCNP in edge-weighted graphs.
* Thin formulation with integer separation is used.*/
vector<long> Thin_Weighted(KGraph &g, long k, long b, vector<long> &Heuristic_sol, bool &subopt);


/*** callback functions for DCNP***/

class GeneralCallbackClass : public GRBCallback
{
protected:
	GRBVar* xvars; // x variables
	double* x; //x values
	GRBVar* yvars; //y variables
	double* y; //y values
	KGraph* g; //original graph
	KGraph* gs; //power graph pointer
	long n; //number of nodes of power graph = number of nodes of original graph
	long gs_m; // number of edges of power graph
	
public:
	GeneralCallbackClass(GRBVar* grb_x_, GRBVar* grb_y_, KGraph* g_,  KGraph* gs_) : xvars(grb_x_), yvars(grb_y_), g(g_), gs(gs_) 
	{
		n = g->n;
		gs_m = gs->m;
		x = new double[gs_m];
		y = new double[n];
	}
	virtual ~GeneralCallbackClass()
	{
		delete[] x;
		delete[] y;
	}
protected:
	void populate_x_MIPSOL()
	{
		x = getSolution(xvars, gs_m);
	}
	void populate_y_MIPSOL()
	{
		y = getSolution(yvars, n);
	}
	void populate_x_MIPNODE()
	{
		x = getNodeRel(xvars, gs_m);
	}
	void populate_y_MIPNODE()
	{
		y = getNodeRel(yvars, n);
	}
};

class IntegerSeparation : public GeneralCallbackClass
{
private:
	unordered_map<long, long> hashing;
	long k;
	long b;

public:
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCutsInteger;
	IntegerSeparation(GRBVar* grb_x, GRBVar* grb_y, KGraph* g, KGraph* gs, unordered_map<long, long> hash_edges, long k1, long b1) : GeneralCallbackClass(grb_x, grb_y, g, gs)
	{
		hashing = hash_edges;
		k = k1;
		b = b1;
	}
protected:
	void callback();
};

class IntegerSeparationWeighted : public GeneralCallbackClass
{
private:
	unordered_map<long, long> hashing;
	long k;
	long b;

public:
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCutsInteger;
	IntegerSeparationWeighted(GRBVar* grb_x, GRBVar* grb_y, KGraph* g, KGraph* gs, unordered_map<long, long> hash_edges, long k1, long b1) : GeneralCallbackClass(grb_x, grb_y, g, gs)
	{
		hashing = hash_edges;
		k = k1;
		b = b1;
	}
protected:
	void callback();
};

class FractionalSeparation : public GeneralCallbackClass
{	
private:
	unordered_map<long, long> hashing;
	long k;
	long b;

public: 
	static long numCallbacks;
	static double TotalCallbackTimeInteger;
	static double TotalCallbackTimeFractional;
	static long numLazyCutsInteger;
	static long numLazyCutsFractional;
	FractionalSeparation(GRBVar* grb_x, GRBVar* grb_y, KGraph* g, KGraph* gs, unordered_map<long, long> hash_edges, long k1, long b1) : GeneralCallbackClass(grb_x, grb_y, g, gs)
	{
		hashing = hash_edges;
		k = k1;
		b = b1;
	}
protected:
void callback();

};




#endif
