#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include "BigInt.h"

using namespace std;

class Graph{
protected:
	int nV;

public:
	vector<vector<int> > adj;
	Graph(){};
	Graph(int _nV) : nV(_nV){
		adj.resize(nV);
	}
	Graph(Graph &g){
		nV = g.nV;
		adj = g.adj;
	}
	Graph& operator=(const Graph &g) {
		if (this != &g) {
			nV = g.nV;
			adj = g.adj;
		}
		return *this;
	}

	virtual void inline addEdge(int v1, int v2) = 0;
	int inline size(){
		return nV;
	}
	void inline resize(int _nV){
		nV = _nV;
		adj.resize(nV, vector<int>());
	}
	void inline reset(){
		adj.resize(nV, vector<int>());
	}
	void inline removeDuplicates(){
		for (int i = 0; i < adj.size(); ++i){
			sort(adj[i].begin(), adj[i].end());
			adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		}
	}
	void printG();
};

class DirGraph : public Graph{
public:
	DirGraph() : Graph() {}
	DirGraph(int _nV) : Graph(_nV) {}
	void inline addEdge(int v1, int v2);
};
void inline DirGraph::addEdge(int v1, int v2){
	adj[v1].push_back(v2);
}

class UnDirGraph : public Graph{
public:
	UnDirGraph() : Graph() {}
	UnDirGraph(int _nV) : Graph(_nV) {}
	void inline addEdge(int v1, int v2);
	void DFS(int srcV);
	void labelConnectedComponents(vector<int> &labels, vector<int> &compSize);
	void findConnectedComponents(vector<vector<int> > &comps);
};
void inline UnDirGraph::addEdge(int v1, int v2){
	adj[v1].push_back(v2);
	adj[v2].push_back(v1);
}

#endif

