#include "Mesh.h"
#include <unordered_map>
#include <unordered_set>

typedef unordered_map<vector<int>, int, container_hash<vector<int>>> efHash;

// For tri mesh, the vertex is NOT critical only if the vertices on its 1-ring
// neighborhood form exactly one single component for each of the two materials
// in interest.
bool TriMesh::isVCritical(int vi, const vector<int>& label, int mat0,
                          int mat1) const {
  // Vi's 1-ring neighbors.
  const vector<int>& voi = Vs[vi].nbV;
  const vector<int>& eoi = Vs[vi].rE;

  vector<int> mat2nV(2, 0);

  // Initialize UnionFind's parent data structure.
  vector<int> parent(voi.size());
  unordered_map<int, int> vi2localInd;
  for (int i = 0; i < voi.size(); ++i) {
    int v0 = voi[i];
    parent[i] = i;
    vi2localInd[v0] = i;
    if (label[v0] == mat0) {
      mat2nV[0]++;
    } else if (label[v0] == mat1) {
      mat2nV[1]++;
    }
  }

  // Special case: one of material has 0 component. It is the single-vertex
  // island situation, making vi critical.
  // TODO: In real world application, we don't like the single-vertex islands.
  // If you want to do something smart with those, here could be a good place.
  if (mat2nV[0] == 0 || mat2nV[1] == 0) return true;

  for (auto ei : eoi) {
    int v0 = Es[ei].v[0];
    int v1 = Es[ei].v[1];

    // Only union-find if the two vertices have the same label and the label is
    // one of the two labels we care about.
    if (label[v0] == label[v1] && (label[v0] == mat0 || label[v0] == mat1)) {
      // UnionFind : Find.
      int p1 = MingUtility::unionfind_find(vi2localInd[v0], parent);
      int p2 = MingUtility::unionfind_find(vi2localInd[v1], parent);
      // UnionFind : Union.
      parent[p1] = p2;
    }
  }
  vector<int> mats = {mat0, mat1};
  // Map a label to a representive vertex (root of union find) to check # of
  // components for each of the two materials.
  unordered_map<int, int> label2root;
  for (int i = 0; i < voi.size(); ++i) {
    int v0 = voi[i];
    if (label[v0] != mat0 && label[v0] != mat1) continue;
    int p1 = MingUtility::unionfind_find(i, parent);
    int l1 = label[v0];
    // If this label is unseen before, do nothing.
    if (label2root.find(l1) == label2root.end()) {
      label2root[l1] = p1;
    }
    // Otherwise, this label is seen before, check whether the root is the
    // same.
    else {
      if (label2root[l1] != p1) return true;
    }
  }
  return false;
}
void TriMesh::ExportMeshToMathematica(const char* filename) const {
  ofstream ofs(filename, ofstream::out);
  ofs << "{";
  // Vs
  ofs << "{";
  for (int i = 0; i < Vs.size(); ++i) {
    ofs << "{";
    Vs[i].write(ofs);
    if (i + 1 == Vs.size()) {
      ofs << "}";
    } else {
      ofs << "},";
    }
  }
  ofs << "},";
  // Es
  ofs << "{";
  for (int i = 0; i < Es.size(); ++i) {
    ofs << "{";
    Es[i].write(ofs);
    if (i + 1 == Es.size()) {
      ofs << "}";
    } else {
      ofs << "},";
    }
  }
  ofs << "},";
  // Fs
  ofs << "{";
  for (int i = 0; i < Fs.size(); ++i) {
    ofs << "{";
    Fs[i].write(ofs);
    if (i + 1 == Fs.size()) {
      ofs << "}";
    } else {
      ofs << "},";
    }
  }
  ofs << "}";

  ofs << "}";
  ofs.close();
}
void TriMesh::Load(const char* filename) {
  vector<vector<float>> inVs;
  vector<vector<int>> inFs;
  ifstream ifs(filename, ios::in);
  if (!ifs) {
    cout << "Error! " << filename << " cannot be found!" << endl;
    return;
  }
  int d = 0, nV = 0, nF = 0;
  if (ifs >> d >> nV >> nF) {
    inVs.resize(nV, vector<float>(2));
    inFs.resize(nF, vector<int>(3));
    int i = 0;
    while (i < nV && ifs >> inVs[i][0] >> inVs[i][1]) i++;
    i = 0;
    while (i < nF && ifs >> inFs[i][0] >> inFs[i][1] >> inFs[i][2]) i++;
  }
  ifs.close();
  Load(inVs, inFs);
}
void TriMesh::Load(const vector<vector<float>>& inVs,
                   const vector<vector<int>>& inFs) {
  int nV = inVs.size();
  int nF = inFs.size();

  // xid: enumerators, eg. 1,2,3
  // xi: index of x in the mesh
  // x: the vertex list of the element x
  // Initialize VCoords, Vs, Es, Fs.
  VCoords = inVs;
  Vs.resize(nV);
  Fs.resize(nF);
  for (int fi = 0; fi < nF; ++fi) {
    auto& oneF = Fs[fi];
    for (int vid = 0; vid < 3; ++vid) {
      oneF.v[vid] = inFs[fi][vid];
    }
  }

  // Define hashes for edges and faces.
  efHash ehash;

  // Combinations of taking pairs(edges) out from 3 vertices.
  vector<vector<int>> f2eid;
  MingUtility::combination(3, 2, f2eid);
  // Candidate Es for each Face.
  vector<vector<int>> candEs(3, vector<int>(2));

  // Read each tet to find edges and faces and fill the data structures.

  for (int fi = 0; fi < nF; ++fi) {
    // cout << "fid=" << fid << endl;
    const auto& f = inFs[fi];
    for (auto vi : f) {
      Vs[vi].nbF.push_back(fi);
    }

    // Gather candidate edges, sort each edge's vid.
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 2; ++j) {
        candEs[i][j] = f[f2eid[i][j]];
      }
      sort(candEs[i].begin(), candEs[i].end());
    }

    for (int eid = 0; eid < candEs.size(); ++eid) {
      // cout << "eid=" << eid << endl;
      auto& e = candEs[eid];
      int v0 = e[0];
      int v1 = e[1];
      int ei = -1;
      // cout << "enen" << endl;

      // cout << "v0=" << v0 << " v1=" << v1 << endl;
      // A new edge is found.
      if (ehash.find(e) == ehash.end()) {
        // cout << "enenIF" << endl;
        ei = Es.size();
        MeshEdge oneE;
        oneE.v[0] = v0;
        oneE.v[1] = v1;
        oneE.nbF.push_back(fi);
        Es.push_back(oneE);
        ehash[e] = ei;
        Vs[v0].nbV.push_back(v1);
        Vs[v1].nbV.push_back(v0);
        for (auto vi : e) {
          Vs[vi].nbE.push_back(ei);
        }
      }
      // And existing edge is met.
      else {
        // cout << "enenELSE" << endl;
        ei = ehash[e];
        Es[ei].nbF.push_back(fi);
      }
      // cout << "keke" << endl;
      // Ts[ti].edge[eid] = ei;
      Fs[fi].edge[eid] = ei;
    }
  }

  // Find the one ring neighborhood of each vertex.
  vector<int> oneE(2);
  for (int vi = 0; vi < nV; ++vi) {
    vector<int>& nbF = Vs[vi].nbF;
    vector<int>& rE = Vs[vi].rE;
    // For each its neighboring tet, find the single face that do not use vi.
    for (int fi : nbF) {
      int pos = 0;
      for (int vid = 0; vid < 3; ++vid) {
        int fvi = Fs[fi].v[vid];
        if (fvi == vi) continue;
        oneE[pos++] = fvi;
      }
      sort(oneE.begin(), oneE.end());
      rE.push_back(ehash[oneE]);
    }
    // TODO isBD is not set yet because right now it is not necessary. Add the
    // code later when isBD is needed.
  }
}
