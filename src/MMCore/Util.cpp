#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "Util.h"
#include <stdlib.h>

using namespace std;

int MingUtility::GCD(int a, int b) {
  a = abs(a);
  b = abs(b);
  int tmp;
  while (a) {
    tmp = b % a;
    b = a;
    a = tmp;
  }
  return b;
}
int MingUtility::GCD(const vector<int>& nums) {
  vector<int> data = nums;
  return GCD_helper(data);
}
int MingUtility::GCD_helper(vector<int>& nums) {
  if (nums.empty()) {
    cout << "ERROR! Check GCD! Calling GCD for empty array." << endl;
    return -1;
  }
  if (nums.size() == 1) return nums[0];
  int a = nums.back();
  nums.pop_back();
  int b = nums.back();
  nums.pop_back();
  nums.push_back(GCD(a, b));
  return GCD_helper(nums);
}
void MingUtility::printRange(int start, int end, int step) {
  for (int i = start; i <= end; i += step) {
    cout << i << " ";
  }
  cout << endl;
}
void MingUtility::readVandT(const char* filename, vector<vector<float>>& Vs,
                        vector<vector<int>>& Ts, int dim) {
  ifstream ifs(filename, ios::in);
  if (!ifs) {
    cout << "Error! " << filename << " cannot be found!" << endl;
    return;
  }
  int nV = 0, nT = 0;
  if (ifs >> nV >> nT) {
    Vs.resize(nV, vector<float>(dim));
    Ts.resize(nT, vector<int>(3));
    int i = 0;
    while (i < nV && ifs >> Vs[i][0] >> Vs[i][1] >> Vs[i][2]) i++;
    i = 0;
    while (i < nT && ifs >> Ts[i][0] >> Ts[i][1] >> Ts[i][2]) i++;
  }
  ifs.close();
}
void MingUtility::readVandT(const char* filename, vector<vector<float>>& Vs,
                        vector<vector<int>>& Ts) {
  ifstream ifs(filename, ios::in);
  if (!ifs) {
    cout << "Error! " << filename << " cannot be found!" << endl;
    return;
  }
  string line;
  if (getline(ifs, line)) {
    // Read #V and #T.
    int nV = 0, nT = 0;
    istringstream in(line);
    in >> nV >> nT;

    // Read Vs.
    Vs.resize(nV);
    Ts.resize(nT, vector<int>(3));
    int i = 0;
    while (i < nV && getline(ifs, line)) {
      istringstream in(line);
      float val;
      while (in >> val) {
        Vs[i].push_back(val);
      }
      i++;
    }

    // Read Ts.
    i = 0;
    while (i < nT && ifs >> Ts[i][0] >> Ts[i][1] >> Ts[i][2]) i++;
  }
  ifs.close();
}
void MingUtility::processMesh(const vector<vector<float>>& Vs,
                          const vector<vector<int>>& Ts,
                          vector<vector<int>>& V2nbV) {
  typedef pair<int, int> edge;
  set<edge> Es;
  V2nbV.resize(Vs.size());
  for (auto oneT : Ts) {
    oneT.push_back(oneT[0]);
    for (int i = 0; i < 3; ++i) {
      int v1 = min(oneT[i], oneT[i + 1]);
      int v2 = max(oneT[i], oneT[i + 1]);
      edge e = make_pair(v1, v2);
      // New edge is found, add the neighbor V information.
      if (Es.find(e) == Es.end()) {
        Es.insert(e);
        V2nbV[v1].push_back(v2);
        V2nbV[v2].push_back(v1);
      }
    }
  }
}
int MingUtility::unionfind_find(int id, vector<int>& parent) {
  if (id != parent[id]) {
    parent[id] = unionfind_find(parent[id], parent);
  }
  return parent[id];
}

void MingUtility::combination(int n, int k, vector<vector<int>>& res) {
  res.clear();
  vector<int> buf;
  combinationHelper(0, n - 1, min(n, k), res, buf);
}
void MingUtility::combinationHelper(int s, int e, int k, vector<vector<int>>& res,
                                vector<int>& buf) {
  if (k <= 0) {
    res.push_back(buf);
    return;
  }
  for (int i = s; i <= e; ++i) {
    buf.push_back(i);
    combinationHelper(i + 1, e, k - 1, res, buf);
    buf.pop_back();
  }
}

bool MingUtility::ExtractLineFromSegments(const vector<vector<int>>& conns,
                                      vector<int>& res) {
  bool isClose = true;
  res.clear();
  if (conns.empty()) return isClose;
  unordered_map<int, vector<int>> v2nbv;
  unordered_set<int> vSet;
  // Find the neighboring Vs of each v.
  for (const auto& conn : conns) {
    v2nbv[conn[0]].push_back(conn[1]);
    v2nbv[conn[1]].push_back(conn[0]);
    vSet.insert(conn[0]);
    vSet.insert(conn[1]);
  }
  vector<int> vs(vSet.begin(), vSet.end());
  // Find the starting v: for close loop, start v can be another v; for open
  // curve, start v must be a v with only one neighbor.
  int src = vs[0];
  for (auto vi : vs) {
    if (v2nbv[vi].size() == 1) {
      isClose = false;
      src = vi;
      break;
    }
  }
  // Extract the loop starting from src.
  res.clear();
  res.push_back(src);
  int cur = src;
  int nxt = v2nbv[cur][0];
  int count = 0;  // If input is wrong, use this to quit anyway.
  int limit = 1000;
  while (nxt != src && count < limit) {
    count++;
    res.push_back(nxt);
    const auto& nbv = v2nbv[nxt];
    if (nbv.size() != 2) {
      break;
    }
    int tmp = nxt;
    nxt = nbv[0] + nbv[1] - cur;
    cur = tmp;
  }
  return isClose;
}
void MingUtility::WriteLevelSet(const vector<vector<float>>& sufVs,
                            const vector<vector<int>>& sufFs,
                            const vector<vector<int>>& sufFMats,
                            const vector<vector<int>>& segEs,
                            const char* filename) {
  ofstream ofs(filename, ofstream::out);
  ofs << sufVs.size() << " " << sufFs.size() << "\n";
  // SufVs.
  for (const auto& oneV : sufVs) {
    ofs << oneV[0] << " " << oneV[1] << " " << oneV[2] << "\n";
  }
  // SufFs and SufFMats.
  for (int fi = 0; fi < sufFs.size(); ++fi) {
    ofs << sufFs[fi][0] << " " << sufFs[fi][1] << " " << sufFs[fi][2] << " "
        << sufFMats[fi][0] << " " << sufFMats[fi][1] << "\n";
  }
  // Contour edges.
  ofs << segEs.size() << "\n";
  for (const auto& oneE : segEs) {
    ofs << oneE[0] << " " << oneE[1] << "\n";
  }
  ofs.close();
}
void MingUtility::LoadFileToVectorOfStrings(const char* filename,
                                        vector<string>& lines) {
  ifstream ifs;
  ifs.open(filename, ios::in);
  if (!ifs) {
    cerr << "Can't open input file " << filename << endl;
    return;
  }
  string line;
  while (ifs >> line) {
    lines.push_back(line);
  }
  ifs.close();
}
