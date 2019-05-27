#include "Graph.h"
#include <stack>
#include <queue>
#include <vector>
#include <iostream>

using namespace std;

void Graph::printG(){
	for (int i=0; i<nV; i++){
		cout << "< " << i << " >:";
		for (unsigned int j=0; j<adj[i].size(); j++){
			cout << " " << adj[i][j];
		}
		cout << endl;
	}
}

void UnDirGraph::DFS(int srcV){
	stack<int> fringe;
	fringe.push(srcV);

	int curV;
	vector<bool> visited(nV, false);
	cout << "DFS: ";
	while (!fringe.empty()){
		curV = fringe.top();
		fringe.pop();
		if (visited[curV]) continue;
		visited[curV] = true;
		cout << curV << " ";
		for (auto nxtV : adj[curV]){
			if(!visited[nxtV]){
				fringe.push(nxtV);
			}
		}

	}
	cout << endl;

}
void UnDirGraph::labelConnectedComponents(vector<int> &labels, vector<int> &compSize){
	BigInt visited(nV);
	labels.resize(nV, -1);

	int compID = 0;
	while(!visited.isAllOne()){
		queue<int> fringe;
		fringe.push(visited.find1st(0));
		visited.set(fringe.front());
		int compN = 0;
		while(!fringe.empty()){
			int curV = fringe.front();
			fringe.pop();
			labels[curV] = compID;
			compN++;
			for(auto nbV : adj[curV]){
				if(!visited.get(nbV)){
					fringe.push(nbV);
					visited.set(nbV);
				}
			}
		}
		compSize.push_back(compN);
		compID++;
	}
}

void UnDirGraph::findConnectedComponents(vector<vector<int> > &comps){
	BigInt visited(nV);

	while(!visited.isAllOne()){
		queue<int> fringe;
		vector<int> oneComp;
		int curV = visited.find1st(0);
		fringe.push(curV);
		visited.set(curV);
		int compN = 0;
		while(!fringe.empty()){
			curV = fringe.front();
			oneComp.push_back(curV);
			fringe.pop();
			compN++;
			for(auto nbV : adj[curV]){
				if(!visited.get(nbV)){
					fringe.push(nbV);
					visited.set(nbV);
				}
			}
		}
		comps.push_back(oneComp);
	}
}