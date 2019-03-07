#include "ConnectorEnumeration.h"
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
								vector<long> sortedsubset = sort4nodes(a, b, i,j);
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
