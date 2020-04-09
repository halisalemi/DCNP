#include "ConnectorEnumeration.h"
#include "GRBInterface.h"
#include <map>


using namespace std;

vector<long> sortnodes(long i, long j, long k)
{
	vector<long> sorted;
	sorted.push_back(i);
	sorted.push_back(j);
	sorted.push_back(k);

	sort(sorted.begin(), sorted.end());

	return sorted;
}

vector<long> sort4nodes(long i, long j, long k, long t)
{
	vector<long> sorted;
	sorted.push_back(i);
	sorted.push_back(j);
	sorted.push_back(k);
	sorted.push_back(t);

	sort(sorted.begin(), sorted.end());

	return sorted;
}

vector<long> sort5nodes(long i, long p, long v, long q, long j)
{
	vector<long> sorted;
	sorted.push_back(i);
	sorted.push_back(p);
	sorted.push_back(v);
	sorted.push_back(q);
	sorted.push_back(j);

	sort(sorted.begin(), sorted.end());

	return sorted;
}

map<vector<long>, long, classcomp> EnumerateLength3Connector(KGraph &g)
{
	//create map, mapsize, and path
	map<vector<long>, long, classcomp>map;
	long mapsize = 0;
	vector<long> path;

	for (long a = 0; a < g.n; a++)
	{
		//compute SP from node a 
		vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
		for (long b = a + 1; b < g.n; b++)
		{
			// If dist(a,b)=1 then nothing is needed to connect them.
			if (dist_from_a[b] == 1)
			{
				long minab = min(a, b);
				long maxab = max(a, b);
				std::map<std::vector<long>, long>::iterator it = map.find({minab, maxab});
				if (it == map.end())
				{
					map.insert(pair<vector<long>, long>({ minab, maxab }, mapsize));
					mapsize++;
				}
			}

			// If dist(a,b)>3 then they cannot be connected with short path.
			if (dist_from_a[b] == 2 || dist_from_a[b] == 3)
			{
				//compute SP from node b 
				vector<long> dist_from_b = g.ShortestPathsUnweighted(b);

				for (long a_neighbors_iterator = 0; a_neighbors_iterator < g.adj[a].size(); a_neighbors_iterator++)
				{
					long i = g.adj[a][a_neighbors_iterator];
					if (dist_from_b[i] == 1) //if i belongs to V_{11}
					{
						//Now, we have found the a,i,b path
						vector<long> sortedsubset = sortnodes(a, i, b);
						std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
						if (it == map.end())
						{
							map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
							mapsize++;
						}
					}
					if (dist_from_b[i] == 2)  // if i belongs to V_{12}
					{
						for (long i_neighbors_iterator = 0; i_neighbors_iterator < g.adj[i].size(); i_neighbors_iterator++)
						{
							long j = g.adj[i][i_neighbors_iterator];
							if (dist_from_a[j] == 2 && dist_from_b[j] == 1) // if j\in N(i) belongs to V_{21}
							{
								//Now, we have found the a,i,j,b path
								//checking whether a,i,j,b path does not belong to map and if not add it to the map	
								vector<long> sortedsubset = sort4nodes(a, b, i, j);
								std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
								if (it == map.end())
								{
									//cerr << minij << " " << maxij << endl;
									map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
									mapsize++;
								}
							}
						}
					}
				}
			}
		}
	}
	return map;
}

std::map<std::vector<long>, long, classcomp> EnumerateLength4Connector(KGraph &g)
{
	//create map, mapsize, and path
	map<vector<long>, long, classcomp>map;
	long mapsize = 0;
	vector<long> path;

	for (long i = 0; i < g.n; i++)
	{
		//compute SP from node i 
		vector<long> dist_from_i = g.ShortestPathsUnweighted(i);
		for (long j = i + 1; j < g.n; j++)
		{
			if (dist_from_i[j] == 1 || dist_from_i[j] > 4) continue;
			vector<long> dist_from_j = g.ShortestPathsUnweighted(j);
			for (long it1 = 0; it1 < g.degree[i]; it1++)
			{
				long p = g.adj[i][it1];
				vector<bool> p_neighbors = boolify(g.adj[p], g.n);
				if (dist_from_j[p] != 2) continue; //Now p belongs to V_12
				for (long it2 = 0; it2 < g.degree[p]; it2++)
				{
					long v = g.adj[p][it2];
					if (dist_from_i[v] != 2 || dist_from_j[v] != 2) continue; //Now v belongs to V_22
					for (long it3 = 0; it3 < g.degree[v]; it3++)
					{
						long q = g.adj[v][it3];
						if (dist_from_i[q] == 2 && dist_from_j[q] == 1 && p_neighbors[q] == false) //type1 V_12, V_22, V_21
						{
							vector<long> sortedsubset = sort5nodes(i, p, v, q, j);
							std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
							if (it == map.end())
							{
								map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
								mapsize++;
							}
							sortedsubset.clear();
						}
						if (dist_from_i[q] == 3 && dist_from_j[q] == 1) //type2 V_12, V_22, V_31
						{
							vector<long> sortedsubset = sort5nodes(i, p, v, q, j);
							std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
							if (it == map.end())
							{
								map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
								mapsize++;
							}
							sortedsubset.clear();
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
				if (dist_from_j[u] != 3) continue; //Now u belongs to V_13
				for (long it2 = 0; it2 < g.degree[u]; it2++)
				{
					long v = g.adj[u][it2];
					if (dist_from_i[v] != 2 || dist_from_j[v] != 2) continue; //Now v belongs to V_22
					for (long it3 = 0; it3 < g.degree[v]; it3++)
					{
						long q = g.adj[v][it3];
						if (dist_from_i[q] == 2 && dist_from_j[q] == 1) //type3 V_13, V_22, V_21
						{
							vector<long> sortedsubset = sort5nodes(i, u, v, q, j);
							std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
							if (it == map.end())
							{
								map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
								mapsize++;
							}
						}
						if (dist_from_i[q] == 3 && dist_from_j[q] == 1) //type4 V_13, V_22, V_31
						{
							vector<long> sortedsubset = sort5nodes(i, u, v, q, j);
							std::map<std::vector<long>, long>::iterator it = map.find(sortedsubset);
							if (it == map.end())
							{
								map.insert(pair<vector<long>, long>(sortedsubset, mapsize));
								mapsize++;
							}
						}
					}
				}
			}
		}
	}
	return map;
}
