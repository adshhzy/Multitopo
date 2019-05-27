
#include "UnionFind_o.h"

UnionFind_o::UnionFind_o(int nEle){
	//init: rank[i] = 0; parent[i] = i;
	parent.resize(nEle);
	rank.resize(nEle);
	for (int i = 0; i < nEle; ++i){
		parent[i] = i;
	}
}

int UnionFind_o::Find(int x){
	if (parent[x] != x){
		// path compression
		parent[x] = Find(parent[x]);
	}
	return parent[x];
}

int UnionFind_o::Union(int x, int y){
	int xRoot = Find(x);
	int yRoot = Find(y);
	if (xRoot == yRoot){
		return xRoot;
	}
	if (rank[xRoot] < rank[yRoot]){
		parent[xRoot] = yRoot;
		return yRoot;
	}
	else if (rank[xRoot] > rank[yRoot]){
		parent[yRoot] = xRoot;
		return xRoot;
	}
	else{
		parent[yRoot] = xRoot;
		rank[xRoot]++;
		return xRoot;
	}
}

UnionFind_o::UnionFind_o(const UnionFind_o &_uf){
	parent = _uf.parent;
	rank = _uf.rank;
}
UnionFind_o & UnionFind_o::operator = (const UnionFind_o &_uf){
	parent = _uf.parent;
	rank = _uf.rank;
	return *this;
}

