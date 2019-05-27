#ifndef _CELL_H_
#define _CELL_H_
#include <vector>
using namespace std;

class Interface {
 public:
  vector<vector<float>> _sufVs;
  vector<vector<int>> _sufFs;
  vector<vector<int>> _sufFMats;
  Interface() {}
  Interface(const vector<vector<float>>& sufVs,
            const vector<vector<int>>& sufFs,
            const vector<vector<int>>& sufFMats) {
    _sufVs = sufVs;
    _sufFs = sufFs;
    _sufFMats = sufFMats;
  }
  Interface(const Interface& itf) {
    _sufVs = itf._sufVs;
    _sufFs = itf._sufFs;
    _sufFMats = itf._sufFMats;
  }
  Interface& operator=(const Interface& itf) {
    _sufVs = itf._sufVs;
    _sufFs = itf._sufFs;
    _sufFMats = itf._sufFMats;
    return *this;
  }
  ~Interface() {}
};
class Cell {
 public:
  vector<vector<int>> _activeFs;
  vector<Interface> _sufs;
  Cell() {}
  Cell(const vector<vector<int>>& activeFs, const vector<Interface>& sufs) {
    _activeFs = activeFs;
    _sufs = sufs;
  }
  Cell(const Cell& cell) {
    _activeFs = cell._activeFs;
    _sufs = cell._sufs;
  }
  Cell& operator=(const Cell& cell) {
    _activeFs = cell._activeFs;
    _sufs = cell._sufs;
    return *this;
  }
  ~Cell() {}
};
#endif // _CELL_H_
