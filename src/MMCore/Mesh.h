#ifndef _MESH_H_
#define _MESH_H_

#include "Util.h"
#include <vector>
#include <fstream>
#include <unordered_map>

using namespace std;
class MeshVertex {
 public:
  bool isBD;
  vector<int> nbV;
  vector<int> nbE;
  vector<int> nbF;
  // Edges in the one ring neighborhood.
  vector<int> rE;
  MeshVertex() : isBD(false) {}
  void writeBasic(ofstream& ofs) const {
    ofs << "{";
    MingUtility::writeVector(nbV, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(nbE, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(nbF, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(rE, ofs);
    ofs << "}";
  }
  void write(ofstream& ofs) const { writeBasic(ofs); }
};

class MeshEdge {
 public:
  int v[2];
  bool isBD;
  vector<int> nbF;//tet neighborhood faces
  // The edges of the nbF except current edge.
  // Used for levelset extraction.
  // TODO only done for tetmesh, haven't done for 2D trimesh.
  vector<int> nbE;

  MeshEdge() : isBD(false) {}
  void writeBasic(ofstream& ofs) const {
    ofs << "{" << v[0] << "," << v[1] << "},";
    ofs << "{";
    MingUtility::writeVector(nbF, ofs);
    ofs << "}";
  }
  void write(ofstream& ofs) const { writeBasic(ofs); }
  int getOtherV(int v1) const { return v[0] + v[1] - v1; }
};

class MeshFace {
 public:
  int v[3];
  int edge[3];
  bool isBD;

  MeshFace() : isBD(false) {}
  void writeBasic(ofstream& ofs) const {
    ofs << "{" << v[0] << "," << v[1] << "," << v[2] << "},";
    ofs << "{" << edge[0] << "," << edge[1] << "," << edge[2] << "}";
  }
  void write(ofstream& ofs) const { writeBasic(ofs); }
};

class MeshVertex3D : public MeshVertex {
 public:
  vector<int> nbT;
  // Faces in the one ring neighborhood.
  vector<int> rF;
  void write(ofstream& ofs) const {
    writeBasic(ofs);
    ofs << ",{";
    MingUtility::writeVector(nbT, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(rF, ofs);
    ofs << "}";
  }
};

class MeshEdge3D : public MeshEdge {
 public:
  vector<int> nbT;
  void write(ofstream& ofs) const {
    writeBasic(ofs);
    ofs << ",{";
    MingUtility::writeVector(nbT, ofs);
    ofs << "}";
  }
};

class MeshFace3D : public MeshFace {
 public:
  vector<int> nbT;
  void write(ofstream& ofs) const {
    writeBasic(ofs);
    ofs << ",{";
    MingUtility::writeVector(nbT, ofs);
    ofs << "}";
  }
};

class MeshTet3D {
 public:
  int v[4];
  int edge[6];
  int face[4];
  vector<int> nbT;

  MeshTet3D() {}
  void write(ofstream& ofs) const {
    ofs << "{" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << "},";
    ofs << "{" << edge[0] << "," << edge[1] << "," << edge[2] << "," << edge[3]
        << "," << edge[4] << "," << edge[5] << "},";
    ofs << "{" << face[0] << "," << face[1] << "," << face[2] << "," << face[3]
        << "},";
    ofs << "{";
    MingUtility::writeVector(nbT, ofs);
    ofs << "}";
  }
};

class _Mesh {
 public:
  vector<vector<float>> VCoords;
  vector<int> VMarkers;
  int _dim;
  _Mesh() : _dim(3) {}
  virtual ~_Mesh() {}
  inline void setDimension(int d) { _dim = d; }
  virtual void Load(const char* filename) = 0;
  virtual void Load(const vector<vector<float>>& inVs,
                    const vector<vector<int>>& inFs) = 0;
  virtual void Load(const vector<vector<float>>& inVs,
                    const vector<int>& inVMarkers,
                    const vector<vector<int>>& inFs) = 0;

  virtual void GenerateIndicatorFunction(int nMat, const vector<int>& labeling,
                                         vector<vector<float>>& v2f) = 0;
  // Given the labeling of its neighbor vertices, check whether vertex vi is
  // critical regarding to material 0 and 1.
  // Input :
  //  vi - index of the vertex of interest.
  //  label - the labeling of all the vertices.
  //  mat0/mat1 - the two labels we care about.
  virtual bool isVCritical(int vi, const vector<int>& label, int mat0,
                           int mat1) const = 0;

  // Function to export the complete mesh to a file, in Mathematica form.
  // The data can be loaded into Mathematica directly by calling
  // "mesh=ToExpression[Import["filename"]];" in Mathematica.
  virtual void ExportMeshToMathematica(const char* filename) const = 0;

  // Return the 1-ring neighbor Vs of vi.
  virtual inline const vector<int>& Get1RingNeighborVs(int vi) const = 0;

  virtual void PreProcessBD(const vector<int>& label, int nLabel) = 0;

  // My implementation of finding components starting from each connected
  // bdComp. Seems to be wrong. But I will keep it here for now.
  virtual bool ExtractLevelSetsBdComp(const vector<int>& label,
                                      vector<vector<int>>& comps,
                                      vector<vector<int>>& conns) const = 0;

  // Supposed to be my final implementation. Find all the connected surface
  // components and their connections. If successful, return 0; if the genus is
  // not 0, return -1; if islands exist, return -2.
  virtual int ExtractLevelSets(bool allowHighGenus, const vector<int>& label,
                               vector<vector<int>>& comps,
                               vector<vector<int>>& conns) const = 0;

  virtual int ExtractLevelSets(bool allowHighGenus, const vector<int>& label,
                       vector<vector<int>>& comps,
                       vector<vector<int>>& conns,
                       vector<int> &comp2eulerC,
                       int &junctionPoints) const = 0;
  // My another implementation of finding connected surface components, using
  // Union-Find. Not complete yet, and seems to be less efficient. Will probably
  // ditch it.
  virtual bool ExtractLevelSets(int nMat, const vector<int>& label,
                                vector<vector<int>>& comps,
                                vector<vector<int>>& conns) const = 0;


  virtual void ExtractBoundaryCurves(const vector<int>& label,
                                     vector<vector<float>>& crvVs,
                                     vector<vector<int>>& crvEs) const = 0;
  virtual void ExtractSmoothBoundaryCurves(
      const vector<int>& label, vector<vector<float>>& crvVs,
      vector<vector<int>>& crvEs) const = 0;
  // nbActiveFs is filled with the index of faces in activeFs, not the index of
  // faces in Fs.
  virtual void GetActiveBDEsNFs(const vector<int>& label) = 0;

  virtual void GenerateLevelSetFromLabeling(
      const vector<int>& label, vector<vector<float>>& sufVs,
      vector<vector<int>>& sufFs, vector<vector<int>>& sufFMats) const = 0;

  virtual void GetCurveInfoForArrangement(vector<vector<int>>& activeFs,
                                          vector<vector<float>>& crvVs,
                                          vector<vector<int>>& crvEs) const = 0;

  virtual void WriteLevelSet(const vector<vector<float>>& sufVs,
                             const vector<vector<int>>& sufFs,
                             const vector<vector<int>>& sufFMats,
                             const char* filename) const = 0;

  virtual void WriteLevelSet(const vector<vector<float>>& sufVs,
                             const vector<vector<int>>& sufFs,
                             const vector<vector<int>>& sufFMats,
                             const vector<int>& MappingMat,
                             const char* filename) const = 0;
  // Boundary vertices.
  vector<int> _bdVs;
  // Boundary edges.
  vector<int> _bdEs;
  // Vertex IDs for each boudary connected component.
  vector<vector<int>> _bdComps;
  // The label for each bdComp.
  vector<int> _bdComp2label;
  // A vector of bdComps for each label.
  vector<vector<int>> _label2bdComps;
  // // The frontier vertices for each bdComp.
  // // A frontier vertex is the vertex connecting a frontier edge in the
  // bdComp.
  // vector<vector<int>> _bdComp2ftVs;
  // The frontier edges for each bdComp.
  // A frontier edge is an edge connecting current bdComp to another bdComp.
  vector<vector<int>> _bdComp2ftEs;
  // The loops on the boundary, represented with the indexes of active edges.
  // Each edge is actually used twice, for its two different components of the
  // ends.
  vector<vector<int>> _bdLoops;
  // relapace edges with its two vertices
  vector<vector<int>> _bdLoopsV;
  // The bdComp that a bdLoop belongs to.
  vector<int> _bdLoop2Comp;
  // The bdLoops that a bdComp has.
  vector<vector<int>> _bdComp2Loops;
  // A vector of bdLoops for each label.
  vector<vector<int>> _label2bdLoops;
  vector<float> _vNbWeights;

  vector<vector<int>> _bdLoops2label;

};

// Class for triangular mesh.
class TriMesh : public _Mesh {
 public:
  vector<MeshVertex> Vs;
  vector<MeshEdge> Es;
  vector<MeshFace> Fs;
  ~TriMesh() {}

  void Load(const vector<vector<float>>& inVs, const vector<vector<int>>& inFs);
  void Load(const vector<vector<float>>& inVs,const vector<int>& inVMarkers, const vector<vector<int>>& inFs){}
  void Load(const char* filename);

  void GenerateIndicatorFunction(int nMat, const vector<int>& labeling,
                                 vector<vector<float>>& v2f) {}
  bool isVCritical(int vi, const vector<int>& label, int mat0, int mat1) const;
  void ExportMeshToMathematica(const char* filename) const;
  const vector<int>& Get1RingNeighborVs(int vi) const { return Vs[vi].nbV; }
  void PreProcessBD(const vector<int>& label, int nLabel) {
    // TODO
  }

  bool ExtractLevelSetsBdComp(const vector<int>& label,
                              vector<vector<int>>& comps,
                              vector<vector<int>>& conns) const {
    // TODO
    return false;
  }

  int ExtractLevelSets(bool allowHighGenus, const vector<int>& label,
                       vector<vector<int>>& comps,
                       vector<vector<int>>& conns) const {
    // TODO
    return false;
  }

  bool ExtractLevelSets(int nMat, const vector<int>& label,
                        vector<vector<int>>& comps,
                        vector<vector<int>>& conns) const {
    // TODO
    return false;
  }

  int ExtractLevelSets(bool allowHighGenus, const vector<int>& label,
                       vector<vector<int>>& comps,
                       vector<vector<int>>& conns,
                       vector<int> &comp2eulerC,
                       int &junctionPoints) const{
      // TODO
      return false;
  }
  void ExtractBoundaryCurves(const vector<int>& label,
                             vector<vector<float>>& crvVs,
                             vector<vector<int>>& crvEs) const {}
  void ExtractSmoothBoundaryCurves(const vector<int>& label,
                                   vector<vector<float>>& crvVs,
                                   vector<vector<int>>& crvEs) const {}
  // nbActiveFs is filled with the index of faces in activeFs, not the index of
  // faces in Fs.
  void GetActiveBDEsNFs(const vector<int>& label) {}
  void GenerateLevelSetFromLabeling(const vector<int>& label,
                                    vector<vector<float>>& sufVs,
                                    vector<vector<int>>& sufFs,
                                    vector<vector<int>>& sufFMats) const {}
  void GetCurveInfoForArrangement(vector<vector<int>>& activeFs,
                                  vector<vector<float>>& crvVs,
                                  vector<vector<int>>& crvEs) const {}
  void WriteLevelSet(const vector<vector<float>>& sufVs,
                     const vector<vector<int>>& sufFs,
                     const vector<vector<int>>& sufFMats,
                     const char* filename) const {}

  void WriteLevelSet(const vector<vector<float>>& sufVs,
                              const vector<vector<int>>& sufFs,
                              const vector<vector<int>>& sufFMats,
                              const vector<int>& MappingMat,
                              const char* filename) const{}
};

// Class for tetrahedron mesh.
class TetMesh : public _Mesh {
 public:
  vector<MeshVertex3D> Vs;
  vector<MeshEdge3D> Es;
  vector<MeshFace3D> Fs;
  vector<MeshTet3D> Ts;
  ~TetMesh() {}

  void Load(const vector<vector<float>>& inVs, const vector<vector<int>>& inTs);
  void Load(const vector<vector<float>>& inVs,
                      const vector<int>& inVMarkers,
                      const vector<vector<int>>& inTs);
  void Load(const char* filename);
  void GenerateIndicatorFunction(int nMat, const vector<int>& labeling,
                                 vector<vector<float>>& v2f);
  void GenerateIndicatorFunction(const int nMat, const vector<int>& labeling,vector<int>&mVs,
                                          vector<vector<float>>& v2f) ;
  bool isVCritical(int vi, const vector<int>& label, int mat0, int mat1) const;
  void ExportMeshToMathematica(const char* filename) const;
  const vector<int>& Get1RingNeighborVs(int vi) const { return Vs[vi].nbV; }
  void PreProcessBD(const vector<int>& label, int nLabel);
  bool ExtractLevelSetsBdComp(const vector<int>& label,
                              vector<vector<int>>& comps,
                              vector<vector<int>>& conns) const;
  int ExtractLevelSets(bool allowHighGenus, const vector<int>& label,
                       vector<vector<int>>& comps,
                       vector<vector<int>>& conns) const;
  bool ExtractLevelSets(int nMat, const vector<int>& label,
                        vector<vector<int>>& comps,
                        vector<vector<int>>& conns) const;
  int ExtractLevelSets(bool allowHighGenus, const vector<int>& label,
                       vector<vector<int>>& comps,
                       vector<vector<int>>& conns,
                       vector<int> &comp2VEF,
                       int &junctionPoints) const;

  // List of boundary faces.
  vector<int> _bdFs;
  // // The frontier faces for each bdComp.
  // // A frontier face is the face containing at least one (actually, exactly
  // two)
  // // frontier edge(s).
  // vector<vector<int>> _bdComp2ftFs;
  // Ming: just for quick try out of curve smoothing
  void ExtractBoundaryCurves(const vector<int>& label,
                             vector<vector<float>>& crvVs,
                             vector<vector<int>>& crvEs) const;
  void ExtractSmoothBoundaryCurves(const vector<int>& label,
                                   vector<vector<float>>& crvVs,
                                   vector<vector<int>>& crvEs) const;
  // nbActiveFs is filled with the index of faces in activeFs, not the index of
  // faces in Fs.
  void GetActiveBDEsNFs(const vector<int>& label);
  void SmoothBoundaryCurves();
  void GenerateLevelSetFromLabeling(const vector<int>& label,
                                    vector<vector<float>>& sufVs,
                                    vector<vector<int>>& sufFs,
                                    vector<vector<int>>& sufFMats) const;

  void GetCurveInfoForArrangement(vector<vector<int>>& activeFs,
                                  vector<vector<float>>& crvVs,
                                  vector<vector<int>>& crvEs) const;

  void WriteLevelSet(const vector<vector<float>>& sufVs,
                     const vector<vector<int>>& sufFs,
                     const vector<vector<int>>& sufFMats,
                     const char* filename) const;

  void WriteLevelSet(const vector<vector<float>>& sufVs,
                              const vector<vector<int>>& sufFs,
                              const vector<vector<int>>& sufFMats,
                              const vector<int>& MappingMat,
                              const char* filename) const;


  // The centroid of each T.
  vector<vector<float>> _tVCoords;
  vector<int> _bdActiveEs;
  vector<int> _bdActiveFs;
  vector<vector<int>> _bdCurveEs;
  vector<vector<float>> _bdCurveVs;
  unordered_map<int, int> _fi2bdActFInd;

  vector<int> _bdActiveF2Vs;  //for global mapping
  vector<int> _bdActiveE2Vs_fornonplcs; //for global mapping

private:
  vector<int> _bdCurveEsFlag;
  vector<int> _bdCurveEsInd;
  vector<vector<float>>_bdEdgeMidPoints;
  vector<int>_bdEMInd;

  void AddbdEdgeMidpoint2bdCurve();

public:
  void OutputCellBoundaryOriInd(vector<int>&Vori_Ind);
  void OutputCellBoundary(vector<double>&out_Vs, vector<uint>&out_Fs, vector<int>&Vori_Ind);
  void OutputCellBoundary(vector<double>&out_Vs, vector<uint>&out_Fs, vector<int>&Vori_Ind, vector<int>&outVmarkers);

  void CutTet(vector<double>&v2planeDist, vector<double> &outVpos, vector<uint> &outF2V);
};
#endif  //_MESH_H_
