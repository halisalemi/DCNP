#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include "KGraph.h"

using namespace std;
// Use to impose a total ordering on paths
struct classcomp
{
	// if lhs comes before rhs return true
	bool operator()(const std::vector<long>& lhs, const std::vector<long>& rhs) const
	{
		//smaller path comes first. As an example {7,9} comes before {4,6,7} 
		if (lhs.size() < rhs.size()) return 1;
		if (lhs.size() > rhs.size()) return 0;

		// if paths have same length lexicographically smaller ones comes first
		for (int i = 0; i < lhs.size(); i++)
		{
			if (lhs[i] < rhs[i]) return 1;
			if (lhs[i] > rhs[i]) return 0;
		}
		return 0;
	}
};

map<vector<long>, long, classcomp> EnumerateLength3Connector(KGraph &g);
map<vector<long>, long, classcomp> EnumerateLength4Connector(KGraph &g);


vector<long> sortnodes(long i, long j, long k);
vector<long> sort4nodes(long i, long j, long k, long t);
vector<long> sort5nodes(long i, long p, long v, long q, long j);

