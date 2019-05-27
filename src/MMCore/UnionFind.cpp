#include "UnionFind.h"
#include <unordered_map>
using namespace std;

void UnionFind::SetElements(const vector<int>& elements) {
  _elements = elements;
  _ele2index.clear();
  _parent.resize(_elements.size(), 0);
  for (int i = 0; i < _elements.size(); ++i) {
    _ele2index[_elements[i]] = i;
    _parent[i] = i;
  }
}
void UnionFind::AddElement(int ele) {
  _ele2index[ele] = _elements.size();
  _elements.push_back(ele);
  _parent.push_back(_parent.size());
}
int UnionFind::Find_index(int ind) {
  if (ind != _parent[ind]) {
    _parent[ind] = Find_index(_parent[ind]);
  }
  return _parent[ind];
}

void UnionFind::ExtractComponents(vector<vector<int>>& comps) {
  comps.clear();
  unordered_map<int, int> root2comp;
  for (int i = 0; i < _elements.size(); ++i) {
    int r = Find_index(i);
    // A new comp is found.
    if (root2comp.find(r) == root2comp.end()) {
      root2comp[r] = comps.size();
      comps.push_back({_elements[i]});
    } else {
      comps[root2comp[r]].push_back(_elements[i]);
    }
  }
}

int UnionFind::GetNumOfComponents() {
  int nComp = 0;
  for (int i=0; i<_parent.size(); ++i) {
    if (i == _parent[i]) {
      nComp++;
    }
  }
  return nComp;
}
