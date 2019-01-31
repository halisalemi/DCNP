#ifndef KGRAPH_H
#define KGRAPH_H
//#include <ilcplex/cplex.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>


extern double copy_time;
using namespace std;

vector<long> UnionOfSortedLists(vector<long>&u, vector<long>&v);
// find B\S, where B and S are sorted (and B=bigger, S=smaller). Assumes S\subseteq B.
vector<long> SetMinusOfSortedLists(vector<long>&B, vector<long>&S);
vector<long> SetIntersectionOfSortedLists(vector<long>&u, vector<long>&v);

class KGraph
{
public:
	// start new functions by Austin
	bool IsConnected(vector<long> S);
	bool IsConnected(vector<bool> S);
	bool IsConnected();

	long DiameterUnweighted();
	long DiameterUnweighted(vector<long> S);
	long DiameterUnweighted(vector<bool> S1);
	long LongestShortestPathUnweighted(long origin);
	vector<long> ShortestPathsUnweighted(long origin);
	vector<long> ShortestPathsUnweighted(long origin, vector<bool> &S);
	//Added by Ali
	vector<long> ShortestPathsUnweighted(long origin, vector<bool> &S, vector<long> &Pre);

	vector< vector<bool> > CreateAdjacenyMatrix();
	bool DeleteNode(long i);
	void ComplementGraph(const KGraph &rhs);
	// end new functions by Austin

	vector<long> *adj;  // stores the adj lists. Each list is maintained as a sorted vector.


	vector<double> *weight; // stores the edge weight.



	long *degree;       // stores the degree seq
	double *wt;
	long n;             // num of nodes
	long m;             // num of edges
	long Delta;         // higest degree. As of now, is instanciated when graph is read or copied, not maintained afterwards
	string name;        // name of the graph. Could be anything.

						/* Constructors */
	KGraph();
	KGraph(string nm);
	KGraph(long n);
	KGraph(string nm, string file, string type);
	KGraph(const KGraph &rhs);

	/* File IO utility functions */
	void ReadDIMACSGraph(string file);
	void ReadDIMACSColorGraph(string file);
	void ReadDIMACSGraphParallel(string file);
	void ReadDATGraph(string file);
	void ReadSNAPGraph(string file);
	bool WriteGVizGraph(string outfile);
	void WriteSNAPGraph(string outfile);
	void WriteDIMACSGraph(string outfile);

	/* General purpose utility functions */
	bool CheckValid();
	void Duplicate(const KGraph &rhs);
	void DuplicateConnected(const KGraph &rhs, map<long, long> &node_map);
	bool AddEdge(long i, long j, bool reverseToo, bool safe);
	bool DeleteEdge(long i, long j, bool reverseToo, bool safe);
	long ConnectedVertices();

	/* k-core and k-community functions (and related) */
	bool KCore(long k);
	bool KCommunity(long k);
	bool KCommunityParallel(long k);
	long CommonNeighbors(long u, vector<long> &v);
	long CommonNeighbors(long u, long v);
	bool HasKCommonNeighbors(long u, long v, long k);
	vector<long> CommonNeighborsList(long u, long v);
	vector<long> CommonNeighborsList(long u, vector<long> &v);

	/* functions for induced subgraphs */
	KGraph CreateInducedGraph(vector<long> &S);
	KGraph CreateInducedGraph(vector<long> &S, vector<long> &ReverseMap);
	KGraph CreateInducedGraph(vector<bool> &S);
	KGraph CreateInducedGraph(vector<bool> &S, vector<long> &ReverseMap);
	void FindConnectedComponents(vector< vector< long> > &components, vector<long> &degreeZero);
	KGraph CreateFarPairsGraph(long s);

	/* Residual - a small sized graph obtained by finding KCore or kcomm using a suitable k. */
	long FindResidualBinary(string method, long u, long b, long ostergard_limit, const KGraph &backup, bool &residual, bool debug = false);
	long FindResidualLinear(string method, long u, long b, long ostergard_limit, const KGraph &backup, bool &residual, bool debug = false);

	/* Destructors */
	void clear();
	~KGraph();

	KGraph CreatePowerGraph(long s);
	vector<long> FindHeuristicClique(vector<long> &degeneracyorder, vector<long> &rightdegree);
	vector<long> FindDegeneracyOrdering(vector<long> &rightdegree);
	vector<long> FindVerticesOfKCore(vector<long> &degeneracyorder, vector<long> &rightdegree, long k);

	void FindInducedGraph(vector<bool> &S);
	void FindInducedGraph(vector<long> &S);
};

#endif
