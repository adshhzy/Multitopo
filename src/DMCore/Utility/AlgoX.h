#ifndef _ALGOX_H_
#define _ALGOX_H_

#include "../Utility/BigInt.h"
#include <assert.h>
#include <iostream>
#include <queue>
#include <vector>


using namespace std;

class Sub{
public:
	BigInt curSets;
	BigInt nxtSets;
	BigInt coveredItems;
	float cost;

	Sub();
	Sub(int _nSet, int _nItem);
	Sub(const Sub &_sub);
	Sub & operator=(const Sub &_sub);

	void print();

	~Sub();
};

class mySubCompare{
public:
	bool operator() (const Sub& lhs, const Sub& rhs ) const {
		return lhs.cost > rhs.cost;
	}
};

typedef priority_queue<Sub, vector<Sub>, mySubCompare> pqueue;

class AlgoX{
public:
	static void AlgoXExactSetCover(int _nItem, int _nSet, 
		vector<vector<int> > & _sets, 
		vector<float> & weights, 
		vector<vector<int> > & _conflicts,
		vector<vector<int> > & grouping, 
		vector<float> & groupingCost);

};

#endif  //_ALGOX_H_
