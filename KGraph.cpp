#include "KGraph.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <queue>
#include <set>
#include <map>
#include <omp.h>
#include "Binary_Heap.h"
using namespace std;


vector<long> UnionOfSortedLists(vector<long>&u, vector<long>&v)
{
	vector<long> soln;
	long p1 = 0;	// position in vector u
	long p2 = 0;	// position in vector v

	while (p1<u.size() || p2<v.size())
	{
		if (p2 >= v.size())
		{
			soln.push_back(u[p1]);
			p1++;
		}
		else if (p1 >= u.size())
		{
			soln.push_back(v[p2]);
			p2++;
		}
		else if (u[p1] < v[p2])
		{
			soln.push_back(u[p1]);
			p1++;
		}
		else if (u[p1] > v[p2])
		{
			soln.push_back(v[p2]);
			p2++;
		}
		else // they are equal, u[p1]=v[p2]
		{
			soln.push_back(u[p1]);
			p1++;
			p2++;
		}
	}
	return soln;
}

// find B\S, where B and S are sorted (and B=bigger, S=smaller). Assumes S\subseteq B.
vector<long> SetMinusOfSortedLists(vector<long>&B, vector<long>&S)
{
	vector<long> soln;
	long p = 0; // position in S.
	for (long i = 0; i < B.size(); i++)	// i is position in B.
	{
		if (p >= S.size() || B[i] < S[p]) soln.push_back(B[i]);
		else p++;
	}
	return soln;
}

vector<long> SetIntersectionOfSortedLists(vector<long>&u, vector<long>&v)
{
	vector<long> soln;
	long p1 = 0;	// position in vector u
	long p2 = 0;	// position in vector v

	while (p1<u.size() || p2<v.size())
	{
		if (p2 >= v.size())
		{
			p1++;
		}
		else if (p1 >= u.size())
		{
			p2++;
		}
		else if (u[p1] < v[p2])
		{
			p1++;
		}
		else if (u[p1] > v[p2])
		{
			p2++;
		}
		else // they are equal, u[p1]=v[p2]
		{
			soln.push_back(u[p1]);
			p1++;
			p2++;
		}
	}
	return soln;
}

bool KGraph::IsConnected(vector<bool> S)
{
	/* Is the subgraph induced by S connected? */
	long u, v;
	vector<bool> reached(n, false);
	vector<long> children, parents;

	for (long i = 0; i<n; i++) //find a root node
		if (S[i]) {
			children.push_back(i);
			reached[i] = true;
			break;
		}

	for (; !children.empty(); ) { //do BFS
		parents = children;
		children.clear();
		for (long i = 0; i<parents.size(); i++) { //for each parent, examine the children
			u = parents[i];
			for (long j = 0; j<degree[u]; j++) {
				v = adj[u][j];
				if (S[v] && !reached[v]) { //can only use vertices in S
					reached[v] = true;
					children.push_back(v);
				}
			}
		}
	}
	for (long i = 0; i<n; i++) //if a vertex in S hasn't been reached, return false
		if (S[i] && !reached[i])
			return false;

	return true;
}
bool KGraph::IsConnected(vector<long> S1)
{
	vector<bool> S(n, false); //convert S1 to bool
	for (long i = 0; i<S1.size(); i++)
		S[S1[i]] = true;
	return IsConnected(S);
}

bool KGraph::IsConnected()
{
	//Is G connected?
	vector<bool> S(n, true);
	return IsConnected(S);
}

long KGraph::DiameterUnweighted()
{
	/* Returns the (unweighted) diameter of a connected graph.
	If the graph is not connected it returns n.
	Solves a series of shortest paths problems. */
	vector<long> ShortestPaths = ShortestPathsUnweighted(0); //check to see if the graph is connected
	for (long i = 0; i<n; i++)
		if (ShortestPaths[i] == n)
			return n;

	long diameter = 0, temp_longest;
	for (long i = 0; i<n; i++) { //solve shortest paths problem, originating from node i
		temp_longest = LongestShortestPathUnweighted(i);
		if (temp_longest > diameter)
			diameter = temp_longest;
	}
	return diameter;
}

long KGraph::DiameterUnweighted(vector<long> S)
{
	/* returns the diameter of the graph induced by S */
	KGraph g = CreateInducedGraph(S);
	return g.DiameterUnweighted();
}

long KGraph::DiameterUnweighted(vector<bool> S1)
{
	vector<long> S;
	for (long i = 0; i<n; i++)
		if (S1[i])
			S.push_back(i);
	return DiameterUnweighted(S);
}

long KGraph::DiameterWeighted()
{
	/* Returns the (unweighted) diameter of a connected graph.
	If the graph is not connected it returns n.
	Solves a series of shortest paths problems. */
	vector<long> ShortestPaths = ShortestPathsUnweighted(0); //check to see if the graph is connected
	for (long i = 0; i<n; i++)
		if (ShortestPaths[i] == n)
			return n;

	long diameter = 0, temp_longest;
	for (long i = 0; i<n; i++)  //solve shortest paths problem, originating from node i
	{
		temp_longest = LongestShortestPathWeighted(i);
		if (temp_longest > diameter)
			diameter = temp_longest;
	}
	return diameter;
}

long KGraph::DiameterWeighted(vector<long> S)
{
	KGraph g = CreateInducedGraph(S);
	return g.DiameterWeighted();
}

vector<long> KGraph::ShortestPathsUnweighted(long origin, vector<bool> &S, vector<long> &predecessor)
{
	/*Finds the shortest paths from node v to all other nodes in graph G[S].
	Assumes the graph is connected.
	Performs BFS.
	Finds predecessor*/
	long u, v;
	vector<long> YourPredecessor(n);
	vector<long> dist(n, n); //shortest distance from origin node to each other node. dist[i] = n means i not reachable
	if (!S[origin]) return dist;  // if origin not in S, return infinities.
	vector<bool> reached(n, false);
	vector<long> children, parents;

	children.push_back(origin);
	dist[origin] = 0; //the origin node is distance 0 from itself
	YourPredecessor[origin] = origin;
	reached[origin] = true;

	for (long d = 1; !children.empty(); d++) { //for each distance
		parents = children;
		children.clear();
		for (long i = 0; i<parents.size(); i++) { //for each parent, examine the children
			u = parents[i];
			for (long j = 0; j<degree[u]; j++) {
				v = adj[u][j];
				if (!reached[v] && S[v]) {
					reached[v] = true;
					YourPredecessor[v] = u;
					dist[v] = d;
					children.push_back(v);
				}
			}
		}
	}
	predecessor = YourPredecessor;
	return dist;
}

vector<long> KGraph::ShortestPathsWeighted(long origin)
{
	vector<long> dist(n);
	vector<long> pre(n);
	vector<bool> visited(n, false);

	//h is our heap
	struct BinaryHeap* h = createMinHeap(n);

	//Initialize heap h with all vertices.
	for (long v = 0; v < n; v++)
	{
		dist[v] = INT_MAX;
		h->array[v] = newMinHeapNode(v, dist[v]);
		h->pos[v] = v;
	}

	// Make dist value of origin vertex as 0 so that it is extracted first
	h->array[origin] = newMinHeapNode(origin, dist[origin]);
	h->pos[origin] = origin;
	dist[origin] = 0;
	decreaseKey(h, origin, dist[origin]);
	pre[origin] = origin;

	// Initially size of min heap is equal to n
	h->size = n;

	/*In the following loop, min heap contains all nodes
	/whose shortest distance is not yet finalized.*/
	for (long counter = 0; counter < n; counter++)
	{
		//Extract the vertex with minimum distance value
		struct BinaryHeapNode* minHeapNode = extractMin(h);
		long u = minHeapNode->v; // Store the extracted vertex number
		visited[u] = true;

		/*Traverse through all adjacent vertices of u (the extracted
		vertex) and update their distance values*/
		for (long i = 0; i < degree[u]; i++)
		{
			long v = adj[u][i];

			if (!visited[v] && dist[u] + weight[u][i] < dist[v])
			{
				dist[v] = dist[u] + weight[u][i];
				pre[v] = u;
				decreaseKey(h, v, dist[v]);
			}
		}
	}
	return dist;
}

vector<long> KGraph::ShortestPathsWeighted(long origin, vector<bool> &S)
{
	vector<long> dist(n,INT_MAX);
	vector<long> pre(n);
	vector<bool> visited(n, false);

	//if origin not in S, return infinities.
	if (!S[origin]) return dist;  

	//h is our heap
	struct BinaryHeap* h = createMinHeap(n);

	//Initialize heap h with all vertices.
	for (long v = 0; v < n; v++)
	{
		dist[v] = INT_MAX;
		h->array[v] = newMinHeapNode(v, dist[v]);
		h->pos[v] = v;
	}

	// Make dist value of origin vertex as 0 so that it is extracted first
	h->array[origin] = newMinHeapNode(origin, dist[origin]);
	h->pos[origin] = origin;
	dist[origin] = 0;
	decreaseKey(h, origin, dist[origin]);
	pre[origin] = origin;

	// Initially size of min heap is equal to n
	h->size = n;

	/*In the following loop, min heap contains all nodes
	/whose shortest distance is not yet finalized.*/
	for (long counter = 0; counter < n; counter++)
	{
		//Extract the vertex with minimum distance value
		struct BinaryHeapNode* minHeapNode = extractMin(h);

		long u = minHeapNode->v; // Store the extracted vertex number

		if (!S[u]) continue;

		visited[u] = true;

		/*Traverse through all adjacent vertices of u (the extracted
		vertex) and update their distance values*/
		for (long i = 0; i < degree[u]; i++)
		{
			long v = adj[u][i];

			if (S[v] && !visited[v] && dist[u] + weight[u][i] < dist[v])
			{
				dist[v] = dist[u] + weight[u][i];
				pre[v] = u;
				decreaseKey(h, v, dist[v]);
			}
		}
	}

	return dist;
}

vector<long> KGraph::ShortestPathsWeighted(long origin, vector<bool> &S, vector<long> &Predecessor)
{
	vector<long> dist(n, INT_MAX);
	vector<long> pre(n);
	vector<bool> visited(n, false);

	//if origin not in S, return infinities.
	if (!S[origin]) return dist;

	//h is our heap
	struct BinaryHeap* h = createMinHeap(n);

	//Initialize heap h with all vertices.
	for (long v = 0; v < n; v++)
	{
		dist[v] = INT_MAX;
		h->array[v] = newMinHeapNode(v, dist[v]);
		h->pos[v] = v;
	}

	// Make dist value of origin vertex as 0 so that it is extracted first
	h->array[origin] = newMinHeapNode(origin, dist[origin]);
	h->pos[origin] = origin;
	dist[origin] = 0;
	decreaseKey(h, origin, dist[origin]);
	pre[origin] = origin;

	// Initially size of min heap is equal to n
	h->size = n;

	/*In the following loop, min heap contains all nodes
	/whose shortest distance is not yet finalized.*/
	for (long counter = 0; counter < n; counter++)
	{
		//Extract the vertex with minimum distance value
		struct BinaryHeapNode* minHeapNode = extractMin(h);

		long u = minHeapNode->v; // Store the extracted vertex number

		if (!S[u]) continue;

		visited[u] = true;

		/*Traverse through all adjacent vertices of u (the extracted
		vertex) and update their distance values*/
		for (long i = 0; i < degree[u]; i++)
		{
			long v = adj[u][i];

			if (S[v] && !visited[v] && dist[u] + weight[u][i] < dist[v])
			{
				dist[v] = dist[u] + weight[u][i];
				pre[v] = u;
				decreaseKey(h, v, dist[v]);
			}
		}
	}
	Predecessor = pre;
	return dist;
}

vector<long> KGraph::BinaryHeapDijkstra(long origin, long sink, vector<bool> &S)
{
	vector<long> dist(n, INT_MAX);
	vector<long> pre(n);
	vector<bool> visited(n, false);

	//if origin not in S, return infinities.
	if (!S[origin]) return dist;

	//h is our heap
	struct BinaryHeap* h = createMinHeap(n);

	//Initialize heap h with all vertices.
	for (long v = 0; v < n; v++)
	{
		dist[v] = INT_MAX;
		h->array[v] = newMinHeapNode(v, dist[v]);
		h->pos[v] = v;
	}

	// Make dist value of origin vertex as 0 so that it is extracted first
	h->array[origin] = newMinHeapNode(origin, dist[origin]);
	h->pos[origin] = origin;
	dist[origin] = 0;
	decreaseKey(h, origin, dist[origin]);
	pre[origin] = origin;

	// Initially size of min heap is equal to n
	h->size = n;

	/*In the following loop, min heap contains all nodes
	/whose shortest distance is not yet finalized.*/
	for (long counter = 0; counter < n; counter++)
	{
		//Extract the vertex with minimum distance value
		struct BinaryHeapNode* minHeapNode = extractMin(h);

		long u = minHeapNode->v; // Store the extracted vertex number

		if (!S[u]) continue;

		//cerr << "u = " << u << endl;
		visited[u] = true;

		/*Traverse through all adjacent vertices of u (the extracted
		vertex) and update their distance values*/
		for (long i = 0; i < degree[u]; i++)
		{
			long v = adj[u][i];
			//cerr << "v = " << v << endl;

			if (S[v] && !visited[v] && dist[u] + weight[u][i] < dist[v])
			{
				dist[v] = dist[u] + weight[u][i];
				pre[v] = u;
				decreaseKey(h, v, dist[v]);
			}
		}
	}
	vector<long> path;
	vector<long> empty;

	/*for (long i = 0; i < pre.size(); i++)
	{
		cerr << pre[i] << endl;
	}*/

	
	/*cerr << "Shortest path from vertex " << origin << " to vertex " << sink << " includes vertices: " << endl;
	for (long i = 0; i < path.size(); i++)
	{
		cerr << path[i] << " ";
	}*/

	long target = sink;
	path.push_back(target);
	while (target != origin)
	{
		target = pre[target];
		path.push_back(target);
	}

	return path;


}

vector<long> KGraph::ShortestPathsUnweighted(long origin)
{
	vector<bool> S(n, true);
	return ShortestPathsUnweighted(origin, S);
}
vector<long> KGraph::ShortestPathsUnweighted(long origin, vector<bool> &S)
{
	/*Finds the shortest paths from node v to all other nodes in graph G[S].
	Assumes the graph is connected.
	Performs BFS.*/
	long u, v;
	vector<long> dist(n, n); //shortest distance from origin node to each other node. dist[i] = n means i not reachable
	if (!S[origin]) return dist;  // if origin not in S, return infinities.
	vector<bool> reached(n, false);
	vector<long> children, parents;

	children.push_back(origin);
	dist[origin] = 0; //the origin node is distance 0 from itself
	reached[origin] = true;

	for (long d = 1; !children.empty(); d++) { //for each distance
		parents = children;
		children.clear();
		for (long i = 0; i<parents.size(); i++) { //for each parent, examine the children
			u = parents[i];
			for (long j = 0; j<degree[u]; j++) {
				v = adj[u][j];
				if (!reached[v] && S[v]) {
					reached[v] = true;
					dist[v] = d;
					children.push_back(v);
				}
			}
		}
	}
	return dist;
}
long KGraph::LongestShortestPathUnweighted(long origin)
{
	vector<long> SP = ShortestPathsUnweighted(origin);
	return (long)*max_element(SP.begin(), SP.end());
}

long KGraph::LongestShortestPathWeighted(long origin)
{
	vector<long> SP = ShortestPathsWeighted(origin);
	return (long)*max_element(SP.begin(), SP.end());
}

vector< vector<bool> > KGraph::CreateAdjacenyMatrix()
{
	/* Creates the associated adjacency matrix. Not recommended for large graphs. */
	long v;
	vector<bool> tempAdjZeros(n, false);
	vector< vector<bool> > adjMatrix(n, tempAdjZeros);
	for (long i = 0; i<n; i++)
	{
		for (long j = 0; j<degree[i]; j++)
		{
			v = adj[i][j];
			adjMatrix[i][v] = true;
		}
	}
	return adjMatrix;
}

bool KGraph::DeleteNode(long i)
{
	/* Deletes all edges which are incident to node i.
	This does NOT actually remove the node per se.
	The value g.n remains the same. */
	for (long j = 0; j<degree[i]; j++)
	{
		long v = adj[i][j];
		DeleteEdge(v, i, false, false);
	}
	m -= degree[i];
	adj[i].clear();
	degree[i] = 0;
	return true;
}


void KGraph::ComplementGraph(const KGraph &rhs)
{
	/* Makes the graph the complement of rhs graph.*/
	if (n>0)
		clear();
	n = rhs.n;
	m = n*(n - 1) / 2 - rhs.m;
	name = rhs.name;
	degree = new long[n];
	adj = new vector<long>[n];
	vector<long>::iterator it; //position within adjacency list
	for (long i = 0; i<n; i++)
	{
		degree[i] = n - rhs.degree[i] - 1;	//update degrees
		if (degree[i] == 0) continue;
		it = rhs.adj[i].begin();
		for (long j = 0; j<n; j++) //add edges.
		{
			if (j == i) continue; //no self-loops
			if (it == rhs.adj[i].end())
				adj[i].push_back(j);
			else if (*it == j)
				it++;
			else adj[i].push_back(j);
		}
	}
}

/* Default constructor. Does nothing but initializing n and m. */
KGraph::KGraph()
{
	n = 0;
	m = 0;
	Delta = 0;
}

/* Default constructor, but names the graph in addition. */
KGraph::KGraph(string nm)
{
	name = nm;
	n = 0;
	m = 0;
	Delta = 0;
}

/* Useful constructor. Names the graph, and reads it from a files. */
KGraph::KGraph(string nm, string file, string type)
{
	n = 0;
	m = 0;
	Delta = 0;
	name = nm;
	if (type == "dimacs")
		ReadDIMACSGraph(file);
	else if (type == "snap_d")
		ReadSNAPGraph(file);
	else if (type == "weighted_graph")
		ReadWeightedGraph(file);
	else if (type == "dimacs_color")
		ReadDIMACSColorGraph(file);
	else if (type == "DAT")
		ReadDATGraph(file);
	else
		cerr << "Format " << type << " not found.\n";
}

/* Default constructor, but initializes the number of nodes and the vectors in addition */
KGraph::KGraph(long nodes)
{
	n = nodes;
	m = 0;
	Delta = 0;
	degree = new long[n];
	memset(degree, 0, sizeof(long)*n);
	adj = new vector<long>[n];
}

/* Copy constructor */
KGraph::KGraph(const KGraph &rhs)
{
	Duplicate(rhs);
}

/* Copying function. Makes the calling graph same as the passed graph. */
void KGraph::Duplicate(const KGraph &rhs)
{
	if (n>0)
		clear();
	n = rhs.n;
	m = rhs.m;
	name = rhs.name;
	Delta = rhs.Delta;
	degree = new long[n];
	adj = new vector<long>[n];
	long i = 0;
	//#pragma omp parallel for	
	for (i = 0; i<n; i++)
	{
		degree[i] = rhs.degree[i];
		adj[i] = rhs.adj[i];
	}
}

void KGraph::DuplicateForWeighted(const KGraph &rhs)
{
	if (n>0)
		clear();
	n = rhs.n;
	m = rhs.m;
	name = rhs.name;
	Delta = rhs.Delta;
	degree = new long[n];
	adj = new vector<long>[n];
	long i = 0;
	//#pragma omp parallel for	
	for (i = 0; i<n; i++)
	{
		degree[i] = rhs.degree[i];
		adj[i] = rhs.adj[i];
	}
	weight = rhs.weight;
}

/* Copying function. Makes the calling graph same as the passed graph, but removes isolated vertices
and renumbers the remaining vertices, making n smaller. */
void KGraph::DuplicateConnected(const KGraph &rhs, map<long, long> &node_map)
{
	if (n>0)
		clear();

	map<long, long> node;
	for (long i = 0; i<rhs.n; i++)
		if (rhs.degree[i] > 0)
		{
			node[i] = n;
			node_map[n] = i;
			n++;
		}
	m = rhs.m;
	name = rhs.name;
	degree = new long[n];
	adj = new vector<long>[n];


	for (long i = 0; i<rhs.n; i++)
	{
		long t = 0;
		if (rhs.degree[i] <= 0)
			continue;
		t = node.find(i)->second;
		degree[t] = rhs.degree[i];
		adj[t].resize(rhs.degree[i]);
		for (long s = 0; s<rhs.degree[i]; s++)
		{
			adj[t][s] = node.find(rhs.adj[i][s])->second;
		}
	}
}

long KGraph::ConnectedVertices()
{
	/* Returns the number of nodes which have positive degree. */
	long connectedVertices = 0;
	for (long i = 0; i<n; i++)
		if (degree[i] > 0)
			connectedVertices++;
	return connectedVertices;
}

/* reverseToo: if false, j will be added to i's adj, but not the other way round.
* safe: if true, a check will be performed to make sure the edge does not already exists. */
bool KGraph::AddEdge(long i, long j, bool reverseToo, bool safe)
{
	vector<long>::iterator it;
	if (degree[i] == 0 || j>adj[i][degree[i] - 1])
		adj[i].push_back(j);
	else
	{
		it = lower_bound(adj[i].begin(), adj[i].end(), j);
		if (!safe)
			adj[i].insert(it, j);
		else
		{
			if (*it != j)
				adj[i].insert(it, j);
			else
				return false;
		}
	}
	degree[i]++;

	if (reverseToo)
	{
		if (degree[j] == 0)
			adj[j].push_back(i);
		else
		{
			it = lower_bound(adj[j].begin(), adj[j].end(), i);
			if (!safe)
				adj[j].insert(it, i);
			else
			{
				if (*it != i)
					adj[j].insert(it, i);
				else return false;
			}
		}
		degree[j]++;
	}
	return true;
}

/* reverseToo: if false, j will be removed from i's adj, but not the other way round.
* safe: if true, a check will be performed to make sure the edge exists. */
bool KGraph::DeleteEdge(long i, long j, bool reverseToo, bool safe)
{
	vector<long>::iterator it = lower_bound(adj[i].begin(), adj[i].end(), j);
	if (!safe)
		adj[i].erase(it);
	else
	{
		if (it != adj[i].end() && *it == j)
			adj[i].erase(it);
		else return false;
	}
	degree[i]--;
	if (reverseToo)
	{
		it = lower_bound(adj[j].begin(), adj[j].end(), i);
		if (!safe)
			adj[j].erase(it);
		else
		{
			if (it != adj[j].end() && *it == i)
				adj[j].erase(it);
			else return false;
		}
		degree[j]--;
	}
	return true;
}

/* The calling graphs becomes the max k-core of the original graph.
* If you need to maintain a copy of the original graph too, then Duplicate the graph before calling this. */
bool KGraph::KCore(long k)
{
	long v;
	long totalDeleted = 0;

	bool updated = false;
	do {
		updated = false;
		long i;
		//#pragma omp parallel for
		for (i = 0; i<n; i++)
		{
			if (degree[i]<k && degree[i] != 0)
			{
				long j;
				//#pragma omp critical
				{
					for (j = 0; j<degree[i]; j++)
					{
						long v = adj[i][j];
						DeleteEdge(v, i, false, false);
					}
					totalDeleted += degree[i];
					adj[i].clear();
					degree[i] = 0;
				}
				updated = true;
			}
		}
	} while (updated);
	m -= totalDeleted;
	return false;
}

bool KGraph::CheckValid()
{
	long m1 = 0;
	for (long i = 0; i<n; i++)
		m1 += degree[i];
	if (m1 != 2 * m)
	{
		cerr << "ERROR: " << m1 << " != " << 2 * m << endl;
		return false;
	}
	for (long i = 0; i<n; i++)
	{
		for (long j = 0; j<degree[i] - 1; j++)
		{
			if (adj[i][j] >= adj[i][j + 1])
			{
				cerr << "ERROR: " << "adj[i][j]=" << adj[i][j] << " >= " << adj[i][j + 1] << "adj[i][j+1]" << endl;
				return false;
			}
		}
	}
	return true;
}

/* The calling graphs becomes the max k-community of the original graph.
* If you need to maintain a copy of the original graph too, then Duplicate the graph before calling this. */
bool KGraph::KCommunity(long k)
{
	long totalDeleted = 0;
	KCore(k + 1);
	if (m == 0) return true;
	bool *check = new bool[n];
	for (long i = 0; i<n; i++)
		check[i] = true;

	bool updated = false;
	long v, v1, deleted;
	do {
		updated = false;
		for (long i = n - 1; i>-1; i--)
		{
			if (!check[i])
				continue;
			check[i] = false;
			deleted = 0;
			for (long j = 0; j<degree[i]; j++)
			{
				v = adj[i][j];
				if (i>v && check[v])
					continue;
				if (!HasKCommonNeighbors(i, v, k))
				{
					updated = true;
					check[v] = true;
					check[i] = true;
					adj[i].erase(adj[i].begin() + j); // DeleteEdge not called because this is faster as location of j is known					
					DeleteEdge(v, i, false, false);
					degree[i]--;
					deleted++;
					j--;
					if (degree[i]<k + 1)
					{
						for (long t = 0; t<degree[i]; t++)
						{
							v1 = adj[i][t];
							check[v1] = true;
							DeleteEdge(v1, i, false, false);
						}
						deleted += degree[i];
						adj[i].clear();
						degree[i] = 0;
					}
					if (degree[v]<k + 1)
					{
						for (long t = 0; t<degree[v]; t++)
						{
							v1 = adj[v][t];
							check[v1] = true;
							DeleteEdge(v1, v, false, false);
						}
						deleted += degree[v];
						adj[v].clear();
						degree[v] = 0;
					}
				}
			}
			totalDeleted += deleted;
			if (m - totalDeleted < ((k + 1)*k) / 2)
			{
				for (long t = 0; t<n; t++)
				{
					adj[t].clear();
					degree[t] = 0;
					m = 0;
				}
				delete check;
				return true;
			}
		}
	} while (updated);
	delete check;
	m -= totalDeleted;
	return true;
}

/* The calling graphs becomes the max k-community of the original graph.
* If you need to maintain a copy of the original graph too, then Duplicate the graph before calling this. */
bool KGraph::KCommunityParallel(long k)
{
	long totalDeleted = 0;
	KCore(k + 1);
	if (m == 0) return true;
	bool *check = new bool[n];
	for (long i = 0; i<n; i++)
		check[i] = true;

	bool updated = false;
	long v, deleted;
	do {
		updated = false;
		long i;
#pragma omp parallel for// schedule(static,1)
		for (i = 0; i<n; i++)
		{
			//cout<<omp_get_num_threads()<<endl;
			if (!check[i])
				continue;
			check[i] = false;
			deleted = 0;
			for (long j = 0; j<degree[i]; j++)
			{
				bool hasK;
#pragma omp critical
				{
					v = adj[i][j];
					//if(i>v && check[v])
					//	continue;                                                 
					hasK = HasKCommonNeighbors(i, v, k);
				}
				if (!hasK)
				{
#pragma omp critical										
					{
						updated = true;
						check[v] = true;
						check[i] = true;
						adj[i].erase(adj[i].begin() + j); // DeleteEdge not called because this is faster as location of j is known
						DeleteEdge(v, i, false, true);
						degree[i]--;
						j--;
					}
				}
			}
		}
	} while (updated);
	delete check;
	m = 0;
	for (long i = 0; i<n; i++)
		m += degree[i];
	m = m / 2;
	cout << "m=" << m << endl;
	return true;
}

/* Returns the number of common neighbors node u has with v */
long KGraph::CommonNeighbors(long u, long v)
{
	long t = 0;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = degree[v];
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &adj[v].front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1<q1end && q<qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t++;
			q++;
			q1++;
		}
		else if (t1>t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Returns the number of common neighbors node u has with v */
bool KGraph::HasKCommonNeighbors(long u, long v, long k)
{
	long t = 0;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = degree[v];
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &adj[v].front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1<q1end && q<qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t++;
			if (t >= k)
				return true;
			q++;
			q1++;
		}
		else if (t1>t2)
			q1++;
		else
			q++;
	}
	if (t >= k)
		return true;
	return false;
}

/* Returns number of common neighbors node u has with v */
long KGraph::CommonNeighbors(long u, vector<long> &v)
{
	long t = 0;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = v.size();
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &v.front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1<q1end && q<qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t++;
			q++;
			q1++;
		}
		else if (t1>t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Returns the sorted list of common neighbors node u has with v */
vector<long> KGraph::CommonNeighborsList(long u, long v)
{
	vector<long> t;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = degree[v];
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &adj[v].front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1<q1end && q<qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t.push_back(t1);
			q++;
			q1++;
		}
		else if (t1>t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Returns the sorted list of common neighbors node u has with v */
vector<long> KGraph::CommonNeighborsList(long u, vector<long> &v)
{
	vector<long> t;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = v.size();
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &v.front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1<q1end && q<qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t.push_back(t1);
			q++;
			q1++;
		}
		else if (t1>t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Finds the subgraph induced by all the nodes in S. S is sorted */
void KGraph::FindInducedGraph(vector<bool> &S)
{
	for (long i = 0; i<n; i++)
		if (!S[i])
			DeleteNode(i);
}
KGraph KGraph::CreateFarPairsGraph(long s)
{
	// initialize some graph g that the function will return 
	KGraph g(n);

	// find the "far" nodes (i.e., those of distance greater than s)
	vector<long> dist;
	for (long i = 0; i < n; i++)
	{
		dist = ShortestPathsUnweighted(i);
		for (long j = 0; j < n; j++)
		{
			if (dist[j] > s)
			{
				g.adj[i].push_back(j);
			}
		}
		// update degree of i, max degree of graph, and edge count
		g.degree[i] = g.adj[i].size();
		if (g.degree[i] > g.Delta) g.Delta = g.degree[i];
		g.m += g.degree[i];
	}
	// we double counted the edges above, so we need to halve it.
	g.m /= 2;
	return g;
}

KGraph KGraph::CreateInducedGraph(vector<bool> &S)
{
	vector<long> S1;
	for (long i = 0; i < n; i++)
	{
		if (S[i]) S1.push_back(i);
	}
	return CreateInducedGraph(S1);
}

KGraph KGraph::CreateInducedGraph(vector<long> &S)
{
	vector<long> ReverseMap;
	return CreateInducedGraph(S, ReverseMap);
}
KGraph KGraph::CreateInducedGraph(vector<bool> &S, vector<long> &ReverseMap)
{
	vector<long> vecS;
	for (long i = 0; i < n; i++)
	{
		if (S[i])
		{
			vecS.push_back(i);
		}
	}
	return CreateInducedGraph(vecS, ReverseMap);
}
KGraph KGraph::CreateInducedGraph(vector<long> &S, vector<long> &ReverseMap)
{
	/* Finds the subgraph induced by all the nodes in S. S is sorted */
	unsigned long S_size = S.size();
	KGraph g(S_size);
	ReverseMap.resize(n, -1);
	for (long i = 0; i<S_size; i++)
		ReverseMap[S[i]] = i;
	for (unsigned long i = 0; i<S_size; i++)
	{
		g.adj[i] = CommonNeighborsList(S[i], S);
		g.degree[i] = g.adj[i].size();
		for (long j = 0; j<g.degree[i]; j++) //relabel the vertices for the new, smaller graph
			g.adj[i][j] = ReverseMap[g.adj[i][j]];
		g.m += g.degree[i];
	}
	g.m /= 2;
	return g;
}

/* This function finds all the connected components in a graph and places them in separate clusters.
* Nodes with degree 0 are placed in a vector called degreeZero.
* Arguments:
* clusters: the vector in which the components are stored. Could be non-empty, signifying
previously found components (For consistency, vertices in old clusters should have degree 0).
* degreeZero: array where singleton nodes are placed. Cleared before being filled. */
void KGraph::FindConnectedComponents(vector< vector< long> > &clusters, vector<long> &degreeZero)
{
	//cerr<<m<<" edges to start off in FindConnectedComponents."<<endl;
	long v;
	degreeZero.clear();
	bool* label = new bool[n];
	for (long i = 0; i<n; i++)
		label[i] = false;
	for (long i = 0; i<n; i++)
	{
		if (degree[i] != 0 && label[i] == false)
		{
			vector<long> cluster;
			cluster.push_back(i);
			label[i] = true;
			long c = 0;
			while (c != cluster.size())
			{
				long j = cluster[c];
				if (label[j] == false)
				{
					cluster.push_back(j);
					label[j] = true;
				}
				for (long t = 0; t<degree[j]; t++)
				{
					v = adj[j][t];
					if (label[v] == false)
					{
						cluster.push_back(v);
						label[v] = true;
					}
				}
				c++;
			}
			sort(cluster.begin(), cluster.end());
			clusters.push_back(cluster);
		}
	}
	for (long i = 0; i<n; i++)
		if (label[i] == false)
			degreeZero.push_back(i);
	delete[] label;
}

KGraph KGraph::CreatePowerGraph(long s)
{
	KGraph g(n);
	g.m = 0;

	for (long i = 0; i<n; i++)
	{
		vector<long> dist = ShortestPathsUnweighted(i);
		for (long j = 0; j<n; j++)
		{
			if (i != j && dist[j] <= s)
			{
				g.adj[i].push_back(j);
				g.degree[i]++;
				g.m++;
			}
		}
	}
	g.m /= 2;
	return g;
}

KGraph KGraph::CreatePowerGraphWeighted(long s)
{
	KGraph g(n);
	g.m = 0;

	for (long i = 0; i<n; i++)
	{
		vector<long> dist = ShortestPathsWeighted(i);
		for (long j = 0; j<n; j++)
		{
			if (i != j && dist[j] <= s)
			{
				g.adj[i].push_back(j);
				g.degree[i]++;
				g.m++;
			}
		}
	}
	g.m /= 2;
	return g;
}


vector<long> KGraph::FindHeuristicClique(vector<long> &degeneracyorder, vector<long> &rightdegree)
{
	vector<long> clique;
	for (long i = 0; i<n && clique.empty(); i++)
	{
		long v = degeneracyorder[i];
		// if v neighbors all vertices after it in the ordering, 
		//		i.e., rightdegree[v] = (n-1) - (i+1) + 1 = n-i-1, 
		//		then v and all vertices after it form a clique.
		if (rightdegree[v] == n - i - 1)
		{
			clique.resize(n - i);
			for (long j = i; j<n; j++)	clique[j - i] = degeneracyorder[j];
			sort(clique.begin(), clique.end());
		}
	}
	return clique;
}
vector<long> KGraph::FindDegeneracyOrdering(vector<long> &rightdegree)
{
	long degeneracy = 0;
	rightdegree.resize(n);

	// initialize deg. Also update max degree Delta just in case.
	Delta = 0;
	for (long i = 0; i<n; i++)
	{
		rightdegree[i] = degree[i];
		Delta = max(Delta, rightdegree[i]);
	}

	// prepare the bins
	vector<long> bin(Delta + 1, (long)0);
	for (long i = 0; i<n; i++)	bin[rightdegree[i]]++;
	long start = 0;
	for (long d = 0; d <= Delta; d++)
	{
		long num = bin[d];
		bin[d] = start;
		start += num;
	}

	// initialize the ordering & position vectors
	vector<long> pos(n);	// pos[v] is position of vertex v in vert
	vector<long> vert(n);	// vert[v] is the v-th vertex in the ordering
	for (long i = 0; i<n; i++)
	{
		pos[i] = bin[rightdegree[i]];
		vert[pos[i]] = i;
		bin[rightdegree[i]]++;
	}

	// reset the bin starting points
	for (long d = Delta; d >= 1; d--) bin[d] = bin[d - 1];
	bin[0] = 0;

	// start peeling away minimum degree nodes
	for (long i = 0; i<n; i++)
	{
		long minv = vert[i];	// this is a min-degree vertex in the remaining graph
		bin[rightdegree[minv]]++;
		degeneracy = max(degeneracy, rightdegree[minv]);

		for (long j = 0; j<degree[minv]; j++) // adjust the degrees of the neighbors of v
		{
			long u = adj[minv][j];
			if (pos[u]>pos[minv])	// this means vertex u is still "in the graph" so we need to update its degree and its bucket
			{
				if (rightdegree[u] == rightdegree[minv])
				{
					long pw = bin[rightdegree[minv]];	// the first position of the bin that contains vertex minv
					long w = vert[pw];					// the vertex in that position
					if (u != w)						// if u is not the first vertex in the bin, swap u and w
					{
						vert[pw] = u;
						vert[pos[u]] = w;
						pos[w] = pos[u];
						pos[u] = pw;
					}
					bin[rightdegree[minv] - 1] = pos[minv] + 1;
					bin[rightdegree[u]]++;
					rightdegree[u]--;
				}
				else
				{
					long pw = bin[rightdegree[u]];
					long w = vert[pw];

					if (u != w)
					{
						vert[pw] = u;
						vert[pos[u]] = w;
						pos[w] = pos[u];
						pos[u] = pw;
					}
					bin[rightdegree[u]]++;
					rightdegree[u]--;
				}
			}
		}
	}
	//cerr << "\n Degeneracy = " << degeneracy << endl;
	return vert;
}


void KGraph::FindInducedGraph(vector<long> &S)
{
	vector<bool> Sbool(n, false);
	for (long i = 0; i < S.size(); i++) Sbool[S[i]] = true;
	return FindInducedGraph(Sbool);
}

vector<long> KGraph::FindVerticesOfKCore(vector<long> &degeneracyorder, vector<long> &rightdegree, long k)
{
	vector<long> vertices;
	for (long i = 0; i<n && vertices.empty(); i++)
	{
		long v = degeneracyorder[i];
		if (rightdegree[v] >= k)	// k-core is this vertex and all to the right in degeneracy order.
		{
			vertices.resize(n - i);
			for (long j = i; j<n; j++)	vertices[j - i] = degeneracyorder[j];
			sort(vertices.begin(), vertices.end());
		}
	}
	return vertices;
}

/* Reads a graph from a file in the DIMACS format.
* The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::ReadDIMACSGraph(string file)
{
	if (n>0)
		clear();
	cerr << "ReadDIMACSGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line[0] != '%')
			lineread = true;
	}
	istringstream nLine = istringstream(line);
	nLine >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i<n; i++)
	{
		lineread = false;
		line.clear();
		while (!lineread)
		{
			getline(input, line);
			if (line[0] != '%')
				lineread = true;
		}
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		if (line == "") continue;
		istringstream iLine = istringstream(line);
		v = -1;
		while (!iLine.eof())
		{
			iLine >> u;
			if (u != v)
			{
				adj[i].push_back(u - 1);
				degree[i]++;
				m++;
				//if(i>(u-1) && !binary_search(adj[u-1].begin(), adj[u-1].end(), i))
				//	cerr<<i<<"-"<<u<<" found, but "<<u<<"-"<<i<<"wasn't\n";
			}
			v = u;
		}
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		if (degree[i] != adj[i].size()) { cerr << " Error on line number " << i << endl; exit(0); }
		Delta = max(degree[i], Delta);
	}
	cerr << endl;
	m = m / 2;
	if (m1 != m)
	{
		cerr << "WRONG DATA!!!!!!!!!!!!!!!!!!!!! " << m << " != " << m1 << endl;
		exit(0);
	}
}

/* Reads a graph from a file in the DIMACS-2 format (clique/coloring challenges).
* The format is described at  */
void KGraph::ReadDIMACSColorGraph(string file)
{
	if (n>0)
		clear();
	cerr << "ReadDIMACSColoringGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line != "" && line[0] != 'c')
			lineread = true;
	}
	istringstream nLine = istringstream(line);
	string p;
	nLine >> p >> temp >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m1 / 1000000 << " dots: ";
	for (long i = 0; i<m1; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		getline(input, line);
		istringstream nLine = istringstream(line);
		nLine >> temp;
		if (temp == "e")
		{
			nLine >> u >> v;
			if (u == v)
				continue;
			adj[u - 1].push_back(v - 1);
			adj[v - 1].push_back(u - 1);
		}
		else i--;
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i<n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	if (m1 != m)
	{
		cerr << "Possible error in ReadDirectedGraphFromFile: " << m1 << "!=" << m << "\n";
		//	exit(0);
	}
}

/* Reads a graph from a file in the DIMACS format.
* The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::ReadDIMACSGraphParallel(string file)
{
	if (n>0)
		clear();
	cerr << "ReadDIMACSGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line[0] != '%')
			lineread = true;
	}
	istringstream nLine = istringstream(line);
	nLine >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << n / 100000 << " dots: ";
	int nthreads = 0;
	int tid = 0;
	long lineNum = 0;

#pragma omp parallel shared(input,lineNum) private(line,tid)
	{
		long u, v;
		long share = 0;
		tid = nthreads++;
		cout << "Number of threads = " << tid << endl;
		string line;
		bool lineread;
		long i = 0;

		while (i<n && lineNum<n)
		{
			lineread = false;
			line.clear();
#pragma omp critical
			if (lineNum<n)
			{
#pragma omp flush (lineNum)
				while (!lineread && lineNum<n)
				{
					getline(input, line);
					if (line[0] != '%')
						lineread = true;
				}
				i = lineNum;
				lineNum++;
			}

			if ((i + 1) % 100000 == 0)
				cerr << ".";

			share++;
			if (i<n)
			{
				//cerr<<tid<<"------Using \t"<<i<<"\t"<<lineNum<<endl;
				if (line == "") continue;
				istringstream iLine = istringstream(line);
				v = -1;

				while (!iLine.eof())
				{
					iLine >> u;
					if (u != v)
					{
						adj[i].push_back(u - 1);
						degree[i]++;
#pragma omp atomic
						m++;
						//if(i>(u-1) && !binary_search(adj[u-1].begin(), adj[u-1].end(), i))
						//	cerr<<i<<"-"<<u<<" found, but "<<u<<"-"<<i<<"wasn't\n";
					}
					v = u;
				}
				sort(adj[i].begin(), adj[i].end());
				adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
				if (degree[i] != adj[i].size()) { cerr << " Error on line number " << i << " " << degree[i] << " " << adj[i].size() << endl; exit(0); }
#pragma omp critical(delta)
				Delta = max(degree[i], Delta);
			}
		}
		cerr << tid << " share = " << (double)share / (double)n << endl;
	}
	cerr << endl;
	m = m / 2;
	if (m1 != m)
	{
		cerr << "WRONG DATA!!!!!!!!!!!!!!!!!!!!! " << m << " != " << m1 << endl;
		exit(0);
	}
}
// reads the .dat graphs from Simonetti et al (2011) The Minimum Connected Dominating Set Problem: Formulation, Valid Inequalities and a Branch-and-Cut Algorithm
void KGraph::ReadDATGraph(string file)
{
	cerr << "ReadDATGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	input >> n >> m;
	cerr << n << " nodes, " << m << " edges suggested. ";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m / 1000000 << " dots: ";
	for (long i = 0; i<m; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		input >> u >> v;
		if (u == v)
			continue;
		v--; u--;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i<n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	cerr << m << " edges read\n";
}
/* Reads a graph from a file in the SNAP format.
* The format is described at http://snap.stanford.edu/data/index.html
* Note: if (i,j) is read from file, both (i,j) and (j,i) are added to the graph to make sure its undirected*/
void KGraph::ReadSNAPGraph(string file)
{
	cerr << "ReadSNAPGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	input >> n >> temp >> m >> temp;

	cerr << n << " nodes, " << m << " edges suggested. ";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m / 1000000 << " dots: ";
	for (long i = 0; i<m; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		input >> u >> v;
		if (u == v)
			continue;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i<n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	cerr << m << " edges read\n";
}

/* Reads a weighted graph from a file.
It assumes that there are no duplicates
and the vertices of an edge are sorted
like 1 2 rather than 2 1 */
void KGraph::ReadWeightedGraph(string file)
{
	cerr << "ReadWeightedGraph ";
	m = 0;
	n = 0;
	long w;
	Delta = 0;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	input >> n >> temp >> m >> temp;


	cerr << n << " nodes, " << m << " edges suggested. ";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	weight = new vector<double>[n];
	long mcopy = m;
	long ncopy = n;

	for (long i = 0; i<m; i++)
	{
		input >> u >> v >> w;
		if (u == v)
			continue;
		adj[u].push_back(v);
		adj[v].push_back(u);
		weight[u].push_back(w);
		weight[v].push_back(w);
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i<n; i++)
	{
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	cerr << m << " edges read\n";
}


/* Writes the graph to a file in the SNAP format.
* The format is described at http://snap.stanford.edu/data/index.html
* Note: if (i,j) is read from file, both (i,j) and (j,i) are added to the graph to make sure its undirected*/
void KGraph::WriteSNAPGraph(string file)
{
	cerr << "WriteSNAPGraph ";
	ofstream output;
	output.open(file.c_str(), ios::out);
	if (!output.is_open())
	{
		cout << "File " << file << " could not be opened!!!\n";
		return;
	}

	output << n << " nodes, " << m << " edges.\n";

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i<n; i++)
	{
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		for (long j = 0; j<degree[i]; j++)
			output << i << "\t" << adj[i][j] << endl;
	}
}

/* Writes the graph to a file in the DIMACS-10 format.
* The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::WriteDIMACSGraph(string file)
{
	cerr << "WriteDIMACSGraph: " << file << endl;
	ofstream output;
	output.open(file.c_str(), ios::out);
	if (!output.is_open())
	{
		cout << "File " << file << " could not be opened!!!\n";
		return;
	}

	output << n << " " << m << endl;

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i<n; i++)
	{
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		for (long j = 0; j<degree[i]; j++)
			output << adj[i][j] + 1 << " ";
		output << endl;
	}
}

/* Writes the graph in the format such that GraphViz can plot it.
* No position is specified, so graphviz determines whats best */
bool KGraph::WriteGVizGraph(string outfile)
{
	ofstream gviz;
	gviz.open(outfile.c_str(), ios::out);
	gviz << "graph test\n{\nnode [shape=point, pin=true, fontsize=1];\n";
	double scale = 10;
	for (long i = 0; i<n; i++)
		gviz << i + 1 << ";\n";

	long temp;
	for (long i = 0; i<n; i++)
	{
		for (long j = 0; j<degree[i]; j++)
		{
			temp = adj[i][j];
			if (i<temp)
				gviz << i + 1 << " -- " << temp + 1
				<< "[" << "color = gray55" << "]"
				<< ";\n";
		}
	}
	gviz << "}";
	gviz.close();
	return true;
}

/* The argument 'method' refers to k-core or k-community etc. The function tries to find
* the largest k such that the k-core or k-community is non-empty, and changes the graph
* to the corresponding residual structure. */
long KGraph::FindResidualBinary(string method, long u, long b, long ostergard_limit, const KGraph &backup, bool &residual, bool debug)
{
	KGraph temp;
	u = u - 2;
	b = b - 2;
	long k = u;
	long last_k = k;
	long l = max(b, (long)0);
	vector<double> connectedVertices(u + 1, -1);

	do {
		if (last_k >= k && temp.m != 0)
			Duplicate(temp);
		else if (last_k >= k)
			Duplicate(backup);

		if (method == "kcomm")
			KCommunity(k);
		else if (method == "kcore")
			KCore(k + 1);
		else if (method == "hybrid")
			KCore(k + 1);

		connectedVertices[k] = ConnectedVertices();
		if (debug && method == "kcomm") cerr << "\n---- k=" << k << ", nodes left = " << connectedVertices[k] << ".";
		else if (debug) cerr << "\n---- k=" << k + 1 << ", nodes left = " << connectedVertices[k] << ".";

		if (!residual && connectedVertices[k] < ostergard_limit && connectedVertices[k] > 0)
		{
			cerr << endl;
			return k + 2;
		}

		if (m != 0 && temp.m != m)
			temp.Duplicate(*this);

		last_k = k;
		if (m == 0)
			u = k;
		else
			l = k;
		k = l + (u - l) / 2;
	} while (u - l > 1);
	if (debug) cerr << endl;

	if (connectedVertices[k] == 0) // At the very end, we might have the wrong k and reduced graphs. So resolve with the right k
		k = k - 1;

	if (connectedVertices[k] >= ostergard_limit && method == "hybrid")
	{
		if (debug) cerr << "Trying kcomm binary now\n";
		return FindResidualBinary("kcomm", k + 2, b, 0, backup, residual, true);
	}

	if (temp.m != 0)
		Duplicate(temp);
	else
		Duplicate(backup);

	if (method == "kcomm")
		KCommunity(k);
	else if (method == "kcore")
		KCore(k + 1);
	else if (method == "hybrid")
		KCore(k + 1);

	connectedVertices[k] = ConnectedVertices();
	if (debug) cerr << "---- k=" << k << ", nodes left = " << connectedVertices[k] << endl;
	residual = true;
	return k + 2;
}

/* The argument 'method' refers to k-core or k-community etc. The function tries to find
* the largest k such that the k-core or k-community is non-empty, and changes the graph
* to the corresponding residual structure. */
long KGraph::FindResidualLinear(string method, long u, long b, long ostergard_limit, const KGraph &backup, bool &residual, bool debug)
{
	long k = b - 2;
	long connectedVertices;
	long nleft = n;
	time_t start_time = clock();
	while (m>0)
	{
		if (method == "kcomm")
			KCommunity(k);
		else if (method == "kcore")
			KCore(k + 1);
		else if (method == "hybrid")
			KCore(k + 1);

		connectedVertices = ConnectedVertices();
		if (debug && method == "kcomm") cerr << "\n---- k=" << k << ", nodes left = " << connectedVertices << ".";
		else if (debug) cerr << "\n---- k=" << k + 1 << ", nodes left = " << connectedVertices << ".";
		if (connectedVertices != 0) nleft = connectedVertices;
		if (!residual && connectedVertices < ostergard_limit && connectedVertices > 0)
		{
			residual = false;
			cerr << endl;
			return k + 2;
		}
		k++;
	}
	if (debug) cerr << endl;
	k = k - 2;
	Duplicate(backup);
	if (nleft >= ostergard_limit && method == "hybrid")
	{
		if (debug) cerr << "Trying kcomm linear now\n";
		return FindResidualLinear("kcommE", k + 2, b, 0, backup, residual, true);
	}
	/*else if(nleft >= ostergard_limit && method=="kcomm")
	{
	if(debug) cerr<<"Trying kcommE linear now\n";
	return FindResidualLinear("kcommE",k+2,b,0,backup,residual,true);
	}*/

	if (method == "kcomm")
		KCommunity(k);
	else if (method == "kcore")
		KCore(k + 1);
	else if (method == "hybrid")
		KCore(k + 1);
	connectedVertices = 0;

	residual = true;
	return k + 2;
}

void KGraph::clear()
{
	//cout<<"===Clearing DG "<<name<<"\n";
	if (n>0)
	{
		for (long i = 0; i<n; i++)
			adj[i].clear();
		delete[]adj;
		delete[]degree;
		n = 0;
		m = 0;
	}
}

KGraph::~KGraph()
{
	//cout<<"===Destructing DG "<<name<<"\n";
	if (n>0)
	{
		for (long i = 0; i<n; i++)
			adj[i].clear();
		delete[]adj;
		delete[]degree;
		n = 0;
		m = 0;
	}
}
