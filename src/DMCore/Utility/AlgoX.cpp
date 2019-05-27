
#include "AlgoX.h"

using namespace std;

Sub::Sub(){
}
Sub::~Sub(){
}
Sub::Sub(int _nSet, int _nItem):
	curSets(BigInt(_nSet)), 
	nxtSets(BigInt(_nSet)),
	coveredItems(BigInt(_nItem)),
	cost(0.0){}

Sub::Sub(const Sub & _sub){
	curSets = _sub.curSets;
	nxtSets = _sub.nxtSets;
	coveredItems = _sub.coveredItems;
	cost    = _sub.cost;
}
Sub & Sub::operator=(const Sub &_sub){
	curSets = _sub.curSets;
	nxtSets = _sub.nxtSets;
	coveredItems = _sub.coveredItems;
	cost    = _sub.cost;
	return (*this);
}
void Sub::print(){
	cout << "curSets:" << endl;
	curSets.print();
	cout << "nxtSets:" << endl;
	nxtSets.print();
	cout << "coveredItems:" << endl;
	coveredItems.print();
}

// the interface specified for Grouping algorithm, packed all the functions in one.
void AlgoX::AlgoXExactSetCover(int nItem, int nSet, vector<vector<int> > & _sets, vector<float> & weights, vector<vector<int> > & _conflicts, vector<vector<int> > &grouping, vector<float> &groupingCost)
{
	vector<BigInt> sets;
	vector<BigInt> conflicts;
	vector<BigInt> item2sets;

	// convert vector<vector<int>> sets to BigInt sets
	// construct item2sets
	sets = vector<BigInt>(nSet, BigInt(nItem));
	item2sets = vector<BigInt>(nItem, BigInt(nSet));
	for (int i = 0; i<nSet; ++i){
		sets[i].set(_sets[i]);
		for (unsigned int j = 0; j<_sets[i].size(); ++j){
			item2sets[_sets[i][j]].set(i);
		}
	}
	// construct confliction graph
	conflicts = vector<BigInt>(nSet, BigInt(nSet));
	for (int i = 0; i<nSet; ++i){
		conflicts[i].set(i);
		for (int j = i + 1; j<nSet; ++j){
			if (sets[i].shareBit(sets[j])){
				conflicts[i].set(j);
				conflicts[j].set(i);
			}
		}
	}

	// ------------ solve ------------ //
	// construct the initial Sub: 
	// curSets = {}; nxtSets = All; cost =0.0;
	Sub initSub(nSet, nItem);
	initSub.nxtSets.flip();

	//BFS
	pqueue fringe;
	fringe.push(initSub);
	Sub curSub(nSet, nItem);
	while (!fringe.empty()) {
		curSub = fringe.top();
		fringe.pop();

		// ----- handleCurSub(curSub, fringe); ----- //
		// if _sub is a solution, claim one success
		if (curSub.coveredItems.isAllOne()){
			// success, all items are covered
			// push sol to solutions
			grouping.push_back(curSub.curSets.getOnesInd());
			groupingCost.push_back(curSub.cost);
			continue;
		}
		
		int nxtI = curSub.coveredItems.find1st(0);
		if (nxtI == -1){
			continue;
		}
		BigInt candidates(nSet);
		candidates = item2sets[nxtI] & curSub.nxtSets;
		while (!candidates.isZero()){
			// pop the pos of first 1, remove that 1 from candidate
			int seti = candidates.popOne();
			Sub nbSub(curSub);
			nbSub.curSets.set(seti);
			nbSub.nxtSets -= conflicts[seti];
			nbSub.coveredItems |= sets[seti];
			nbSub.cost += weights[seti];

			// if deadend, do not push to the pqueue 
			//!isDeadend(nbSub)
			if (nbSub.coveredItems.isAllOne() || !nbSub.nxtSets.isZero()){
				fringe.push(nbSub);
			}
		}
	}
}

