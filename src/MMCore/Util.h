#ifndef _UTIL_H_
#define _UTIL_H_

#ifndef db1
#define db1 false
#endif

#ifndef db2
#define db2 false
#endif

#ifndef TIMING
#define TIMING false
#endif

#ifndef REUSESOLVER
#define REUSESOLVER true
#endif

#ifndef PERTURBATION
#define PERTURBATION 0.001
#endif

#ifndef ROUNDFACTOR
#define ROUNDFACTOR 100000.0f
#endif

#ifndef NMAXSUFOUTPUT
#define NMAXSUFOUTPUT 100
#endif

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <boost/functional/hash.hpp>
#include <math.h>
#include <algorithm>
#include <functional>

using namespace std;

// Hash function for any container, defined with boost::hash_range.
template <typename Container>
struct container_hash {
  std::size_t operator()(Container const& c) const {
    return boost::hash_range(c.begin(), c.end());
  }
};

class MingUtility {
 public:
  static void printRange(int start, int end, int step);

  template <class T>
  static void printVector(const vector<vector<T>>& in);

  template <class T>
  static void printVector(const vector<T>& in);

  template <class T>
  static void writeVector(const vector<vector<vector<T>>>& vec, ofstream& ofs);

  template <class T>
  static void writeVector(const vector<vector<T>>& vec, ofstream& ofs);

  template <class T>
  static void writeVector(const vector<T>& vec, ofstream& ofs);

  template <class T>
  static void writeVector_wB(const vector<T>& vec, ofstream& ofs);

  template <class T>
  static void writeVector_wB_tbc(const vector<T>& vec, ofstream& ofs);

  template <class T>
  static void writeVector(const T& vec, const char* filename);

  template <class T>
  static void writeVectorForm(const vector<vector<T>>& vec, ofstream& ofs);

  template <class T>
  static void writeVectorForm(const vector<T>& vec, ofstream& ofs);

  template <class T>
  static void writeVectorForm(const T& vec, const char* filename);

  static void readVandT(const char* filename, vector<vector<float>>& Vs,
                        vector<vector<int>>& Ts, int dim);
  static void readVandT(const char* filename, vector<vector<float>>& Vs,
                        vector<vector<int>>& Ts);

  static void processMesh(const vector<vector<float>>& Vs,
                          const vector<vector<int>>& Ts,
                          vector<vector<int>>& V2nbV);

  // The find() function for union find algorithm.
  static int unionfind_find(int id, vector<int>& parent);
  // The union() function for union find algorithm.
  static inline void unionfind_union(int r0, int r1, vector<int>& parent);
  // The find() + union() function for union find algorithm.
  static inline void unionfind_FindnUnion(int i0, int i1, vector<int>& parent);

  template <class T>
  static void fillVectorWithRange(T start, T step, vector<T>& vec);

  // template <class T>
  // static bool vectorSmaller(const vector<T>& vec0, const vector<T>& vec1);
  template <class T>
  static inline void sort2Dvector(vector<vector<T>>& vecs);

  static void combination(int n, int k, vector<vector<int>>& res);
  static void combinationHelper(int s, int e, int k, vector<vector<int>>& res,
                                vector<int>& buf);
  static int GCD(int a, int b);
  static int GCD(const vector<int>& nums);
  static int GCD_helper(vector<int>& nums);
  template <class T>
  static inline T L2Norm(const vector<T>& nums);
  template <class T>
  static void SumUpVectors(const vector<vector<T>>& vecs, vector<T>& res);
  template <class T>
  static void SumUpTwoVectors(const vector<T>& invec, vector<T>& outvec);
  template <class T>
  static void MultiplyVectorsByValue(T value, vector<T>& vec);

  // Return true if the line is a closed loop.
  static bool ExtractLineFromSegments(const vector<vector<int>>& conns,
                                      vector<int>& res);
  static void WriteLevelSet(const vector<vector<float>>& sufVs,
                            const vector<vector<int>>& sufFs,
                            const vector<vector<int>>& sufFMats,
                            const vector<vector<int>>& segEs,
                            const char* filename);
  static void LoadFileToVectorOfStrings(const char* filename,
                                        vector<string>& lines);
};
template <class T>
void MingUtility::fillVectorWithRange(T start, T step, vector<T>& vec) {
  int val = start;
  for (int i = 0; i < vec.size(); ++i) {
    vec[i] = val;
    val += step;
  }
}
template <class T>
void MingUtility::printVector(const vector<vector<T>>& in) {
  for (auto i = 0; i < in.size(); ++i) {
    cout << " <" << i << "> : ";
    printVector(in[i]);
  }
}

template <class T>
void MingUtility::printVector(const vector<T>& in) {
  for (auto ele : in) {
    cout << ele << " ";
  }
  cout << endl;
}

template <class T>
void MingUtility::writeVector(const vector<T>& vec, ofstream& ofs) {
  for (unsigned int j = 0; j < vec.size(); ++j) {
    if (j == 0)
      ofs << vec[j];
    else
      ofs << "," << vec[j];
  }
}

template <class T>
void MingUtility::writeVector_wB(const vector<T>& vec, ofstream& ofs) {
  ofs << "{";
  writeVector(vec, ofs);
  ofs << "}";
}

template <class T>
void MingUtility::writeVector_wB_tbc(const vector<T>& vec, ofstream& ofs) {
  writeVector_wB(vec, ofs);
  ofs << ",";
}

template <class T>
void MingUtility::writeVector(const vector<vector<T>>& vec, ofstream& ofs) {
  for (unsigned int i = 0; i < vec.size(); ++i) {
    if (i == 0)
      ofs << "{";
    else
      ofs << ",{";
    for (unsigned int j = 0; j < vec[i].size(); ++j) {
      if (j == 0)
        ofs << vec[i][j];
      else
        ofs << "," << vec[i][j];
    }
    ofs << "}";
  }
}

template <class T>
void MingUtility::writeVector(const vector<vector<vector<T>>>& vec, ofstream& ofs) {
  for (unsigned int i = 0; i < vec.size(); ++i) {
    if (i == 0)
      ofs << "{";
    else
      ofs << ",{";
    writeVector(vec[i], ofs);
    ofs << "}";
  }
}

template <class T>
void MingUtility::writeVector(const T& vec, const char* filename) {
  ofstream ofs(filename, ofstream::out);
  ofs << "{";
  writeVector(vec, ofs);
  ofs << "}";
  ofs.close();
}

template <class T>
void MingUtility::writeVectorForm(const vector<vector<T>>& vec, ofstream& ofs) {
  ofs << vec.size() << "\n";
  for (unsigned int i = 0; i < vec.size(); ++i) {
    writeVectorForm(vec[i], ofs);
  }
}

template <class T>
void MingUtility::writeVectorForm(const vector<T>& vec, ofstream& ofs) {
  // ofs << vec.size() << "\n";
  for (unsigned int j = 0; j < vec.size(); ++j) {
    ofs << vec[j] << " ";
  }
  ofs << "\n";
}

template <class T>
void MingUtility::writeVectorForm(const T& vec, const char* filename) {
  ofstream ofs(filename, ofstream::out);
  writeVectorForm(vec, ofs);
  ofs.close();
}
//
// template <class T>
// bool MingUtility::vectorSmaller(const vector<T>& vec0, const vector<T>& vec1) {
//   if (vec0.empty()) return true;
//   if (vec1.empty()) return false;
//   int i = 0, j = 0;
//   while (i < vec0.size() && j < vec1.size()) {
//     if (vec0[i] < vec1[j]) return true;
//     if (vec0[i] > vec1[j]) return false;
//     i++;
//     j++;
//   }
//   if (j == vec1.size()) return false;
//   return true;
// }

template <class T>
void MingUtility::sort2Dvector(vector<vector<T>>& vecs) {
  for (auto& vec : vecs) {
    sort(vec.begin(), vec.end());
  }
  sort(vecs.begin(), vecs.end());
}
void MingUtility::unionfind_union(int r0, int r1, vector<int>& parent) {
  parent[r0] = r1;
}

void MingUtility::unionfind_FindnUnion(int i0, int i1, vector<int>& parent) {
  unionfind_union(unionfind_find(i0, parent), unionfind_find(i1, parent),
                  parent);
}

template <class T>
T MingUtility::L2Norm(const vector<T>& nums) {
  T res = T(0.0);
  T exp = T(2.0);
  for (auto val : nums) {
    res += pow(val, exp);
  }
  return sqrt(res);
}

template <class T>
void MingUtility::SumUpVectors(const vector<vector<T>>& vecs, vector<T>& res) {
  if (vecs.empty()) return;
  int n = vecs[0].size();
  res.resize(n, 0);
  for (int i = 0; i < vecs.size(); ++i) {
    std::transform(res.begin(), res.end(), vecs[i].begin(), res.begin(),
                   std::plus<T>());
  }
}
template <class T>
void MingUtility::SumUpTwoVectors(const vector<T>& invec, vector<T>& outvec) {
  std::transform(outvec.begin(), outvec.end(), invec.begin(), outvec.begin(),
                 std::plus<T>());
}

template <class T>
void MingUtility::MultiplyVectorsByValue(T value, vector<T>& vec) {
  for (auto& ele : vec) {
    ele *= value;
  }
}

#endif  // _UTIL_H_
