#ifndef _UnionFind_o_H_
#define _UnionFind_o_H_

#include <vector>

using namespace std;

class UnionFind_o{
public:
    vector<int> parent;
    vector<int> rank;

    UnionFind_o(int nEle);
    int Find(int x);
    int Union(int x, int y);

    UnionFind_o(const UnionFind_o &_uf);
    UnionFind_o & operator= (const UnionFind_o &_uf);
};
/*
class UnionFind_o {

public:
    vector<int> _parent;
    vector<int> _elements;
    unordered_map<int, int> _ele2index;


    UnionFind_o(int nEle)
    {
        //init: rank[i] = 0; parent[i] = i;
        _parent.resize(nEle);
        _elements.resize(nEle);
        for (int i = 0; i < nEle; ++i){
            _parent[i] = i;
        }

    }
    UnionFind_o() {}
    ~UnionFind_o() {}
    // Add the initial elements for union-find. This function is supposed to be
    // called only once. If called twice, everything will be reset.
    void SetElements(const vector<int>& elements);

    // Add additional element for union-find. It can be useful if initial elements
    // are unkown in the beginning.
    void AddElement(int ele);
    int Find(int ele) { return Find_index(_ele2index[ele]); }
    int Find_index(int ind);

    // Union with elements.
    //void Union(int ele0, int ele1) { Union_roots(Find(ele0), Find(ele1)); }
    int Union(int ele0, int ele1) {
        int a = Find(ele0),b=Find(ele1);
        Union_roots(a, b);
        return b;
    }
    // Union with index instead of elements.
    void Union_indexes(int ind0, int ind1) {
        _parent[Find_index(ind0)] = Find_index(ind1);
    }
    // Union with root instead of elements or index.
    void Union_roots(int ind0, int ind1) { _parent[ind0] = ind1; }

    int getIndex(int ele) {
        // if (_ele2index.find(ele) == _ele2index.end()) return -1;
        return _ele2index[ele];
    }
    int getElement(int ind) {
        // if (ind >= _elements.size()) return -1;
        return _elements[ind];
    }
    void ExtractComponents(vector<vector<int>>& comps);
    int GetNumOfComponents();
};*/
#endif  //_UnionFind_o_H_
