
#ifndef _MM_LEVELSET_H_
#define _MM_LEVELSET_H_

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <cfloat>
// #include <boost/functional/hash.hpp>
#include "LP_solver.h"
#include "Util.h"
#include "Mesh.h"

using namespace std;

typedef vector<float> OffsetVector;
typedef vector<int> Range;

// // An Ordering is a sorted list of vertices by the ascending order of
// fi(v)-fj(v).
// // The 1st element in the pair is fi(v)-fj(v), the 2nd is the index of
// vertex.
// typedef vector<pair<float,int>> Ordering;
typedef vector<int> Ordering;
typedef vector<int> Labeling;
typedef vector<float> Ray;
typedef vector<int> RayInt;
//
// // Hash function for any container, defined with boost::hash_range.
// template <typename Container>
// struct container_hash {
//   std::size_t operator()(Container const& c) const {
//     return boost::hash_range(c.begin(), c.end());
//   }
// };
//
typedef unordered_map<Range, Labeling, container_hash<Range>> R2Lmap;
// For critical edge analysis, map a range to position in allRanges and its
// labeling.
typedef unordered_map<Range, pair<int, Labeling>, container_hash<Range>>
    R2ILmap;

class EventPoint {
 public:
  float _lamda;
  int _vi;
  int _mat0;
  int _mat1;
  float _score;
  float _span;
  OffsetVector offset;
  //Labeling label;
  EventPoint()
      : _lamda(0.0), _vi(0), _mat0(0), _mat1(0), _score(0.0), _span(1.0) {}
  EventPoint(float lamda, int vi, int mat0, int mat1)
      : _lamda(lamda),
        _vi(vi),
        _mat0(mat0),
        _mat1(mat1),
        _score(0.0),
        _span(1.0) {}
  EventPoint(const EventPoint& pt) {
    _lamda = pt._lamda;
    _vi = pt._vi;
    _mat0 = pt._mat0;
    _mat1 = pt._mat1;
    _score = pt._score;
    _span = pt._span;
    offset = pt.offset;
    //label = pt.label;
  }
  EventPoint& operator=(const EventPoint& pt) {
    _lamda = pt._lamda;
    _vi = pt._vi;
    _mat0 = pt._mat0;
    _mat1 = pt._mat1;
    _score = pt._score;
    _span = pt._span;
    offset = pt.offset;
    //label = pt.label;
    return *this;
  }
  void write(ofstream& ofs) const {
    ofs << "{" << _lamda << "," << _vi << "," << _mat0 << "," << _mat1 << ","
        << _score << "," << _span << ",";
    ofs << "{";
    MingUtility::writeVector(offset, ofs);
    ofs << "},{";
    //MingUtility::writeVector(label, ofs);
    ofs << "}";
    ofs << "}";
  }
};


class VEFG{
    int n_V,n_E,n_F;
    int eulerC;
};
class CellTopology {//use to detect duplication of red point (representation point) in one common topology region
    // store the contouring surface of a red point
    // == -> ray shooting cpp
 public:
  // Connected surface components.
  vector<vector<int>> _comps; //int->loop index
  vector<vector<int>> _conns;//surface touching int->index of _comps (pair)

  // To be removed in the future.
  OffsetVector _offset;
  Labeling _label;
  float _score;
  float _span;
  bool isContainIso = false;
  int _junctionpoints;
  // TODO ming
  // add vector<int> comp2genus;
  vector<int> comp2VEF;
  vector<int> mSeq;
  vector<int> isoComps2VEF;
  vector<int> isoComps2label;

  vector<vector<int>> comps2allloops;
  vector<vector<int>> isoComps2allloops;

  CellTopology() : _score(0.0), _span(0.0),_junctionpoints(0) {}
  CellTopology(const CellTopology& ct) {
    _comps = ct._comps;
    _conns = ct._conns;
    _offset = ct._offset;
    _label = ct._label;
    _score = ct._score;
    _span = ct._span;
    isContainIso = ct.isContainIso;
    comp2VEF = ct.comp2VEF;
    mSeq = ct.mSeq;
    isoComps2VEF = ct.isoComps2VEF;
    isoComps2label = ct.isoComps2label;
    _junctionpoints = ct._junctionpoints;

    comps2allloops = ct.comps2allloops;
    isoComps2allloops = ct.isoComps2allloops;
  }
  CellTopology& operator=(const CellTopology& ct) {
    _comps = ct._comps;
    _conns = ct._conns;
    _offset = ct._offset;
    _label = ct._label;
    _score = ct._score;
    _span = ct._span;
    isContainIso = ct.isContainIso;
    comp2VEF = ct.comp2VEF;
    mSeq = ct.mSeq;
    isoComps2VEF = ct.isoComps2VEF;
    isoComps2label = ct.isoComps2label;
    _junctionpoints = ct._junctionpoints;

    comps2allloops = ct.comps2allloops;
    isoComps2allloops = ct.isoComps2allloops;
    return *this;
  }
  void write(ofstream& ofs) const {
    ofs << "{" << _score << "," << _span << ",";
    ofs << "{";
    MingUtility::writeVector(_offset, ofs);
    ofs << "},{";
    MingUtility::writeVector(_label, ofs);
    ofs << "},{";
    MingUtility::writeVector(_comps, ofs);
    ofs << "},{";
    MingUtility::writeVector(_conns, ofs);
    ofs << "}";
    ofs << "}";
  }
  void writebinary(ofstream& ofs)  {



      auto writeVVector = [&ofs](vector<vector<int>>&vec){
          int num;
          num = vec.size();
          ofs.write((char*)(&num),sizeof(int));
          for(int i=0;i<vec.size();++i){
              num = vec[i].size();
              ofs.write((char*)(&num),sizeof(int));
          }
          for(int i=0;i<vec.size();++i){
              ofs.write((char*)(vec[i].data()),sizeof(int)*vec[i].size());
          }
      };
      auto writeVector = [&ofs](vector<int>&vec){
          int num;
          num = vec.size();
          ofs.write((char*)(&num),sizeof(int));

          ofs.write((char*)(vec.data()),sizeof(int)*vec.size());

      };
      auto writeVectorF = [&ofs](vector<float>&vec){
          int num;
          num = vec.size();
          ofs.write((char*)(&num),sizeof(int));

          ofs.write((char*)(vec.data()),sizeof(float)*vec.size());

      };
      ofs.write((char*)(&_score),sizeof(float));
      ofs.write((char*)(&_junctionpoints),sizeof(int));

      writeVectorF(_offset);
      writeVVector(_comps);
      writeVVector(_conns);
      writeVector(comp2VEF);

  }
  void readbinary(ifstream& ifs)  {



      auto readVVector = [&ifs](vector<vector<int>>&vec){
          int num;
          ifs.read((char*)(&num),sizeof(int));
          vec.resize(num);
          for(int i=0;i<vec.size();++i){
              ifs.read((char*)(&num),sizeof(int));
              vec[i].resize(num);
          }
          for(int i=0;i<vec.size();++i){
              ifs.read((char*)(vec[i].data()),sizeof(int)*vec[i].size());
          }
      };
      auto readVector = [&ifs](vector<int>&vec){
          int num;
          ifs.read((char*)(&num),sizeof(int));
          vec.resize(num);
          ifs.read((char*)(vec.data()),sizeof(int)*vec.size());

      };
      auto readVectorF = [&ifs](vector<float>&vec){
          int num;
          ifs.read((char*)(&num),sizeof(int));
          vec.resize(num);
          ifs.read((char*)(vec.data()),sizeof(float)*vec.size());

      };
      ifs.read((char*)(&_score),sizeof(float));
      ifs.read((char*)(&_junctionpoints),sizeof(int));
      readVectorF(_offset);
      readVVector(_comps);
      readVVector(_conns);
      readVector(comp2VEF);
  }

  bool writetxt(string filename)  {

      filename = filename + ".txt";
      ofstream fout(filename.data());
      if(fout.fail()){
          cout<<"Fail to create output file: "<<filename<<endl;
          return false;
      }
      fout<<_score<<endl;

      MingUtility::writeVector(_comps, fout);

      fout<<endl;



  }

  void Unique() {
    // TODO after comps is sort, conns is incorrect.
    MingUtility::sort2Dvector(_comps);
    // MingUtility::sort2Dvector(_conns);
  }
};

// MMLevelSet : Level Set on Multiple Material.
class MMLevelSet {
 public:
  MMLevelSet() : _N(0), _K(0), _upbound(1.0), _mesh(NULL) {}
  MMLevelSet(int N, int K) : _N(N), _K(K), _upbound(1.0), _mesh(NULL) {}
  MMLevelSet(const vector<vector<float>>& v2f)
      : _N(0), _K(0), _upbound(1.0), _mesh(NULL) {
    _V2F = v2f;
    if (_V2F.empty()) {
      cout << "STOP! Are you sure? 0 vertex is found." << endl;
      return;
    }
    _N = _V2F.size();
    _K = _V2F[0].size();
    // _score = 0.0;
  }
  ~MMLevelSet() {
    if (_mesh != NULL) delete _mesh;
  }

  // Take a offset vector (and its corresponding range) as the seed,
  // explore the neighboring ranges.
  // Find all the unique topologies withint those ranges.
  void ExploreTopologies(const OffsetVector& seedOffsetVector,
                         int maxRangeToExplore, int maxTopoToFind,
                         vector<Range>& topologies);

  // Some preprocessing before exploration.
  // 1. Find the boundary vertices. (_bdVs)
  // 2. Add perturbation to avoid overlapping hyperplane. (_V2F)
  // 3. Find the ordering of vertices for each pair of material. (Orders)
  // 4. Initialize the LP solver. (_solver)
  // 5. Pre-calculate log(fi) for each V, for score calculation. (_V2LogF)

  void PreProcess();

  // Given a range, calculate the score of the labeling represented by the
  // range.
  float Range2Score(const Range& range);

  // Given an offset, calculate the score of the labeling represented by the
  // offset.
  float OffsetVector2Score(const OffsetVector& offset);

  // Given a labeling, calculate the score of the labeling.
  float Labeling2Score(const Labeling& label);

  // Convert a range to a labeling.
  void Range2Labeling(const Range& range, Labeling& label);

  // Convert a range to a labeling.
  void Range2Offset(const Range& range, OffsetVector& offset);

  // Convert a range to a LP problem and solve.
  bool SolveRangeWithLP(const Range& range);

  // Convert a range to a LP problem and solve.
  bool SolveRangeWithLP(const Range& range, LPSolver& solver);

  // Given a range, find the neighbors of this range.
  void Range2NeighborRanges(const Range& range, vector<Range>& nbs);

  // check whether a range is valid (not empty).
  bool IsValidRange(const Range& range);

  // Convert an offset vector to a range.
  void OffsetVector2Range(const OffsetVector& offsetvector, Range& range) const;

  // Convert an offset vector to a labeling.
  void OffsetVector2Labeling(const OffsetVector& offsetvector,
                             Labeling& label) const;

  // Find the labeling for a subset of vertices with an offset vector specified.
  void OffsetVector2LabelingForVs(const OffsetVector& offsetvector,
                                  const vector<int>& Vs, Labeling& label) const;

  // Convert a labeling to a range.
  void Labeling2Range(const Labeling& label, Range& range) const;

  // Vetex to scalar function fi(v_p).
  // V2F[p][i] is the function value of vertex p on material i.
  // V2F.size() = N = # vertices
  // V2F[0].size() = K = # material
  vector<vector<float>> _V2F;
  vector<vector<float>> _V2LogF;

  vector<int>_initLabel;

  // The ordering of vertices for each pair of material.
  // Orders[i][j] is an ordering of vertices by ascending order of fi(v) -
  // fj(v).
  // Orders[j][i] is reverse of Orders[i][j]. Only Orders[i][j] is stored
  // (i<j).
  vector<vector<Ordering>> _Orders;

  // Number of vertices.
  int _N;

  // Number of material.
  // TODO: sometimes a cell does not use all the materials, so the _K we need to
  // set is actually smaller than nLabel. Handle it and make it consistent with
  // Mesh process and V2F reading.
  int _K;

  // The upper bound of si.
  float _upbound;

  // LP solver.
  LPSolver _solver;

  // For testing, fake score
  // float _score;

  // Flag used to keep track of which range is visited already.
  // R2Lmap visited;
  R2ILmap visited;

  // Vertices on the boundary that server as the boundary conditions.
  // The label on those vertices cannot be changed during range exploration.
  // TODO: instead of setting that additionaly to V2F, can we prove that the
  // vertices that cannot change the label are equevalent to the points with
  // scalar values to be (0,..,0,1,0,0..0) (has one and only one 1 value)?
  // Jan 4: instead of using unordered_set, use vector<bool> to indicate bdVs.
  // Improve efficiency by trading set.find() with more space.
  // unordered_set<int> _bdVs;
  vector<bool> _bdVs;

  //_Mesh* _mesh;
  TetMesh* _mesh;

  // Rays shot from the origin, for sampling.
  // Note each ray is NOT a unit-lenght vector. The other end touches the
  // boundary of the range space, which means min[ray]=0.0 and max[ray]=1.0.
  vector<Ray> _rays;

  // Critical points on each ray.
  vector<vector<EventPoint>> _criticalPtsPerRay;

  // Topology points on each ray.
  vector<vector<EventPoint>> _topologyPtsPerRay;

  // Unique cell topologies.
  vector<CellTopology> _cellTopologies;
  // For debug purpose.
  vector<vector<int>> _cellTopo2topoPtID;
  vector<vector<vector<int>>> _cellTopo2AllTopoPtID;
  vector<vector<vector<int>>> _badCellTopoPtID;

  vector<int>_cellNexplored2Topogroup;

  // Generate the sampling rays using integer based rejection-sampling.
  void GenerateRays(float density);
  void GenerateRays_helper(float density, int pos, RayInt& ray,
                           set<RayInt>& rays);

  // TODO: for now, I put the ray-sampling idea in this big function.
  // When the framework is further detailed, modify data-structure/functions
  // to
  // make the code better organized, and make the interface for the DP
  // algorithm.
  int ExploreTopologiesAlongRays(bool allowHighGenus, float density,
                                  float segLenLowBound);

  // Set the upper bound of si.
  inline void setVarUpbound(float upbound) { _upbound = upbound; }
  bool IsBoundaryValid(const Range& range, int mat1, int mat2);

  // For critical edge analysis in range space.
  // Find all the critical edges in the dual graph of the range pockets.
  vector<vector<int>> _V2nbVs;
  inline void loadV2nbVs(const char* filename) {
    vector<vector<float>> Vs;
    vector<vector<int>> Ts;
    MingUtility::readVandT(filename, Vs, Ts);
    MingUtility::processMesh(Vs, Ts, _V2nbVs);
  }

  void LoadMesh(const char* filename);
  void LoadCell(const vector<vector<float>>& cellVs,
                const vector<vector<int>>& cellTs,
                const vector<int>& cellInitLabel, int nMat);
  void LoadCell(const vector<vector<float>>& cellVs,
                const vector<int>&cellVMarkers,
                const vector<vector<int>>& cellTs,
                const vector<int>& cellInitLabel, int nMat);
  void PreProcessMeshBD();

  vector<Range> allRanges;
  vector<pair<int, int>> allEdges;
  map<pair<int, int>, int> edge2id;
  set<pair<int, int>> criticalEandV;


/************************************************************************/
//add for user control
private:
  vector<int>CellBoundary_VoriInd;
  vector<double>CellBoundary_Vs;
  vector<uint>CellBoundary_Fs;
  vector<int>CellBoundary_Vmarkers;
  vector<int>CellBoundary_OriLabel;

  void GetCellBoundaryInfo();
public:
  void clearupVectorSpace();
public:
  void OutputCellBoundary(vector<double>&Vs,vector<uint>&Fs,vector<int>&VMats);

  void CutTet(vector<double>&v2planeDist,vector<double>&outVpos,vector<uint>&outF2V);

  int AddConstraintPointInsideCell(vector<vector<float>> &aVs, vector<int> &aVMs, string outdir, int ci);

/************************************************************************/
public:

  void ReadTopologyFromFile(ifstream &ifs);
  void WriteTopologyToFile(ofstream &ofs);


/************************************************************************/
};

#endif  //_MM_LEVELSET_H_
