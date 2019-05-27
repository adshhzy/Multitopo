#ifndef _ARRANGEMENT_H_
#define _ARRANGEMENT_H_

#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../Utility/Utility.h"
#include "../Utility/BigInt.h"
#include "../Utility/Graph.h"
#include "../Utility/UnionFind_o.h"
#include "../Utility/VolReader.h"
#include "../Utility/AlgoX.h"
#include "tetgen.h"
#include "TagMapping.h"
#include "CellTet.h"
#include "CellSuf.h"
#include "../LuCtr2Suf/JuFair/meshFair.h"
//#include "../MingTriMultPoly/Algorithm/DMWT.h"
//#include <omp.h>
//#include "SmoothPatch.h"
#include "../Utility/MVC.h"

using namespace std;

class Arrangement{
public:
    /*add for connection, program would be reconstructed later*/
    tetgenio out;
    vector<int>_framef2Cs;
    vector<float>_frameVstack;
    vector<vector<float>> _frameVs;
    vector<vector<int>> _frameFs;
    vector<vector<int>> _frameCells;
    vector<vector<float>> _frameF2PlaneParas;
    vector<vector<float>>_Cs2PlaneParas;
    vector<float>_globalVpos;
    vector<vector<int>>_globalEdges;
    vector<vector<int>>_globalEdgesMat;
    int _nCell;
    int _nMat;
    int _nCrossSection;




    void connectParallelArrangement(vector<int>_framef2Cs,
                                    vector<float>_frameVstack,
                                    vector<vector<float>> _frameVs,
                                    vector<vector<int>> _frameFs,
                                    vector<vector<int>> _frameCells,
                                    vector<vector<float>> _frameF2PlaneParas,
                                    int _nCell,
                                    int _nCrossSection);


    /*Vs*/
    vector<vector<float> > Vs;	//{{x,y,z},...}

    /*Segs*/
    vector<vector<int> > seg2V;	//{{v1,...},...}
    vector<int> seg2F;
    vector<int> seg2lmat;
    vector<int> seg2cyc;
    vector<vector<int> > segGraph;
    vector<BigInt> segGraph_b;
    vector<int> V2seggV; // map a vertex of the curve to a vertex of the seg graph

    /*Fs*/
    vector<vector<int> > F2seg;	//{{seg1,...},...}
    vector<int> F2PlaneID;

    /*Cells*/
    vector<vector<int> > cell2F;	//{{f1,...},...}
    vector<vector<int> > cell2Fside;
    vector<vector<int> > cellGraph;
    vector<vector<int> > cell2cyc;
    vector<int> cell2ncyc;
    vector<vector<int> > cell2seg; // record segs in each cell
    vector<BigInt> cell2seg_b;
    vector<bool> cellEmpty;

    /*Cycs*/
    vector<vector<int> > cyc2seg;	//{{seg1,...},...}
    vector<BigInt> cyc2seg_b;	//same as cyc2seg, but represent segs as bits
    vector<vector<int> > cyc2V;

    /*Fronts*/
    vector<vector<int> > front2cyc;//{{cyc1,...},...}
    vector<BigInt> front2seg_b;

    bitIntHash cycsegb2cycid;
    vector<vector<float> > planeParas;
    vector<int> Order;
    vector<int> nonEmpty_Order;

    // the indexes in the 4 maps are local indexes
    vector<vector<int> > front_L_cellcyc2nf0cos;
    vector<vector<BigInt> > front_B_cellcyc2f1cyc;
    vector<vector<BigInt> > front_B_cellcyc2f0cyc;
    vector<vector<BigInt> > front_B_f0cyc2f1cyc;

    int nCell;
    int nF;
    int nSeg;
    int nSegV;
    int nCycs;
    int nFront;
    int maxNCyc;

    bool inputVol;
    bool inputBBox;

    OneCellSuf finalSuf;

    vector<vector<int> > allSets;
    vector<vector<vector<int> > > allGroupings; // for each cell, for each group, the set ids
    vector<vector<float> > allGroupingScore;
    vector<vector<float> > allGroupingIso;
    vector<vector<float> > cell2Loc_setIso;
    vector<vector<BigInt> > cell2Loc_setActEdge;
    vector<vector<int> > cell2Loc_Global_setMap;
    int nSets;

    AllCellTet allCellTets;
    vector<float> bbox_org; //{lo_x, lo_y, lo_z, hi_x, hi_y, hi_z}

    float cellFrameLowerCorner [3];
    float cellFrameScale [3];
    Volume volInfo;
    float t_randwb;

    float volLowerCorner [3];
    float volScale [3];

    vector<float> unifyTransVec;
    float unifyScaleRatio;

    vector<int> matMarks;

    AllCellSuf allCellSufs;

    Arrangement();
    ~Arrangement();

    void getCyc2V();
    void loadVol(const char * filename, const vector<float> & transVec, const float scaleRatio);
    void loadBBox(const char * filename);
    void constructArrangement(
            int ssvernum, int ssedgenum, int ssfacenum, int ssspacenum,
            float *ssver, int *ssedge, int **ssface, int **ssspace,
            int *ssfaceedgenum, int *ssspacefacenum, int **ssspace_planeside, int *ssface_planeindex,

            int planenum, float *pparam, float pbbox [],

            int *ctrfvernum, int *ctrfedgenum,
            float **ctrfverpos, int **ctrfvertype, int **ctrfverval,
            int **ctrfedge, int **ctrfedgetype, int **ctrfedgeval
            );

    void constructArrangementMM(
            int ssvernum, int ssedgenum, int ssfacenum, int ssspacenum,
            float *ssver, int *ssedge, int **ssface, int **ssspace,
            int *ssfaceedgenum, int *ssspacefacenum, int **ssspace_planeside, int *ssface_planeindex,

            int planenum, float *pparam, float pbbox [],

            int *ctrfvernum, int *ctrfedgenum,
            float **ctrfverpos, int **ctrfvertype, int **ctrfverval,
            int **ctrfedge, int **ctrfedgetype, int **ctrfedgeval
            );

    void moveIntersectionPtsBackToPlane();
    void moveSegPtsBackToPlane();
    void moveCellPtsBackToPlane();
    void writeArrangement(const char* filename);
    void sufCell_genSurface(vector<int> &resGrouping);


    void cleanTwoEndsOfOneSeg(vector<int> & oneSeg, vector<int> & toDelSegVs);
    void dealDegeneratedOneSegWithTwoV(vector<int> & oneSeg);
    void findCloseCyclesInFronts(vector<int>& _order); //given an search order, find the fronts, and xor'ed cycles
    void findCloseCyclesInCells();
    int countContinousOpenSeg(BigInt &cmnSegs_b);
    int findMatchVer(vector<float>& allVs, vector<float>& oneV);
    void traceOneSeg(vector<vector<int> > &nbVs,vector<vector<int> > &nbVsLeftMaterial, vector<int> &ind2ind, vector<bool> &unvisited, vector<int> &oneSeg, int end, int &leftMaterial);
    void traceCycsFromSegs(BigInt &segs_b, vector<int>& resCycs);
    inline bool isSegClosed(int segi);
    inline bool cycShareSeg(int cyc1, int cyc2);
    void preCalcFront_L_B();
    void explorTopologyInCell(int ci, bitIntHash& allSets_b);
    void tetCell_genCellTet(int ci, tetgenio &tetout);
    void tetCell_genCellTet2(int ci, tetgenio &tetout);

    void tetCell_tagActiveEdge(int ci, float iso, BigInt &actEdge);
    void tetCell_tagActiveEle(int ci, float iso, vector<float> &VPro, BigInt &actTet, BigInt &actFace, BigInt &actEdge);
    void tetCell_tetgen(int ci,
                        tetgenio &tetin,
                        tetgenio &tetout,
                        tetgenio &tetaddin,
                        unordered_set<int> &cellV,
                        unordered_set<int> &segV,
                        unordered_map<int, int> &segV2segi,
                        unordered_map<int, int> &cellV2tetV,
                        unordered_map<int, int> &segV2tetV);
    void tetCell_findComponent(int ci,
                               tetgenio &tetout,
                               unordered_map<int, int> &segV2tetV,
                               UnDirGraph &graph,
                               vector<vector<int> > &tetV2nbTet,
                               vector<int> &labelInVs,
                               vector<int> &labelOutVs,
                               vector<int> &labelXVs,
                               vector<queue<int> > &isoSufFringe,
                               vector<int> &tetV2cyci);
    void tetCell_Vaxman(const int ci, const vector<int> &labelInVs, const vector<int> &labelOutVs, const vector<int> &labelXVs, vector<float> &VPro);
    void tetCell_randomWalk(
            UnDirGraph &graph,
            vector<int> &labelInVs,
            vector<int> &labelOutVs,
            vector<int> &labelXVs,
            vector<float> &VPro,
            int tmpCi);
    void tetCell_volWeight(int ci, vector<float> &VProVol, vector<int> &labelXVs);
    void tetCell_criticalPoint(
            vector<vector<int> > &tetV2nbTet, vector<float> &VPro, tetgenio &tetout, vector<int> &labelXVs, vector<float> &criticalIsoVal, int tmpCi
            );
    int tetCell_criticalPoint_nComp(int ci, int vID, vector<int> &nbTets, vector<float> &VPro, tetgenio &tetout);
    void tetCell_findGroupings(int ci, vector<float> &VPro, vector<float> &VProVol, vector<int> &labelXVs, UnDirGraph &graph, vector<float> &criticalIsoVal, vector<queue<int> > &isoSufFringe, vector<int> &tetV2cyci, bitIntHash& allSets_b);
    float tetCell_findGroupings_GroupScore(vector<float> &VPro, vector<float> &VProVol, vector<int> &labelXVs, float iso);
    float tetCell_findGroupings_GroupScore_noVol(vector<float> &VPro, vector<float> &VProVol, vector<int> &labelXVs, float iso);
    float tetCell_findGroupings_mix_GroupScore(int ci, const vector<int> &curMixIsoGrouping, const vector<int> &cycFirstV, const vector<BigInt> &loc_set2actEdge, int nLocSet, const vector<float> &VPro, const vector<float> &VProVol, const vector<int> &labelXVs);

    int tetCell_getMatOfFaceNextToSeg(int ci, int fi, int ei, vector<int> &nxtSegVID);

    void sufCell_smoothSurface();
    void sufCell_CoarseSuf_marchingTet(int cellId, int groupId);
    void sufCell_marchingTet_calcInterpolateV(int cellId, int ei, float iso, vector<int> &tetE2sufV, vector<vector<float> > &sufV, unordered_map<int, int> &segV2sufV/*, int vi*/);
    void sufCell_marchingTet_calcInterpolateV(int cellId, int ei, float iso, unordered_map<int, int> &tetE2sufV, vector<vector<float> > &sufV, unordered_map<int, int> &segV2sufV);
    bool sufCell_isTetAct(int ci, int ti, float iso);
    void sufCell_findSegsInSufV(int cellId, unordered_map<int, int> &segV2sufV);
    void sufCell_mergeSegsAndSubdivideSufs(vector<vector<int> > &seg2cells, vector<vector<int> > &segInCellPos);
    void sufCell_mergeSubOneSegSuf(int cell0, int cell1, int pos0, int pos1);
    void sufCell_subOneSegSuf(vector<vector<float> > &sufV0,
                              vector<vector<float> > &sufV1,
                              vector<vector<int> > &sufF0,
                              vector<int> &sufSeg2F0,
                              vector<int> &sufSeg0,
                              vector<int> &sufSeg1,
                              vector<vector<int> > &merge0,
                              vector<int> &mergedSufSeg0
                              );
    inline void sufCell_hashSegFaces(vector<int> &oneTriangle, BigInt &isSufVASegV, unordered_map<vector<int>, int, myVecHash> &SegE2Fhash, vector<vector<int> > &sufF){
        int v0, v1;
        for(int i=0; i<3; ++i){
            if(i==2){
                v0 = oneTriangle[0]; v1 = oneTriangle[2];
            }else{
                v0 = oneTriangle[i]; v1 = oneTriangle[i+1];
            }
            if(isSufVASegV.get(v0)&&isSufVASegV.get(v1)){
                vector<int> oneHashEdge(2);
                oneHashEdge[0] = v0;
                oneHashEdge[1] = v1;
                SegE2Fhash[oneHashEdge] = sufF.size();
            }
        }
    }
    void sufCell_genSurface_Triangulation(vector<int> &resGrouping);
    void sufCell_TriangulateEachCell(int cellId, int gourpId);
    void sufCell_CoarseSuf_TMP(int cellId, int gourpId);
    int sufCell_MingTriMultPoly(
            // input
            vector<vector<vector<float> > > &inCurves,
            vector<vector<vector<float> > > &inNorms,
            bool useDelaunay,
            bool useMinSet,
            bool useNormal,
            float areaWeight,
            float edgeWeight,
            float dihedralWeight,
            float boundaryNormalWeight,
            // output
            int &_numofpoints,
            int &_numoftilingtris,
            float * _pPositions,
            float * _pNormals,
            int * _pFaceIndices,
            int cellId,
            int setId
            );

    void sufCell_saveFinalSufForSmooth(const int label);

    int sufCell_MingTriMultPoly(
            // input
            vector<vector<vector<float> > > &inCurves,
            vector<vector<vector<float> > > &inNorms,
            bool useDelaunay,
            bool useMinSet,
            bool useNormal,
            float areaWeight,
            float edgeWeight,
            float dihedralWeight,
            float boundaryNormalWeight,
            // output
            vector<vector<int> > &sufF,
            int cellId,
            int setId
            );
    void sufCell_JuFair( int t_nSmBefLoop, int t_nLoop, int t_nSmInLoop, bool m_doSwap);
    void sufCell_SmoothPatch( const int _nSub, const int _lapOrder, const int _lapIter );
    void sufCell_getRidOfNonManifold_simpleCase(int nSufV, vector<vector<int> > &sufF, vector<bool> &ignoreVs,vector<bool> &ignoreFs);
    void sufCell_getRidOfNonManifold(int nSufV, vector<vector<int> > &sufF, vector<bool> &ignoreVs,vector<bool> &ignoreFs);
    void sufCell_findAnchorV(int cellId, int curTet, float iso, int &anchorV, bool &isAnchorVPositive);
    bool sufCell_isTriOriented(int cellId, int curTet, float iso, vector<float> &v0, vector<float> &v1, vector<float> &v2);

    void tetCell_tetgen_scaleUp(tetgenio &tetin, float scale, float translate []);
    void tetCell_tetgen_scaleBack(tetgenio &tetout, float scale, float translate []);

    float getVolWeight(float xyzi [], float xyzj []);
    float getVolValAtPos(float xyz[]);

    inline float tetCell_calcGenus(int nC, int nB, int nV, int nE, int nF){
        return 0.5*(2*nC - nB - nV + nE - nF);
    }

    void writeLuCtrGraph(const char* filename);

    vector<vector<float> > cellVs;	//{{x,y,z},...}
    vector<vector<int> > cellEs;	//{{v1,v2},...}
    vector<vector<int> > cellFs;	//{{e1,...},...}
    vector<vector<int> > sorted_verIdOnCellE;	//{{v1,..},...}

    inline float EuclideanDistSquare1Cell2Seg(int c_v, int s_v);
    void reOrderCellEinCellF();

    vector<vector<double> >tetV;
    vector<vector<float> >tetVh;
    vector<vector<int> >tetVmarker;
    vector<int>tetnV;

    vector<vector<int> >tetQ;

    inline double getTetVCoord(int ci, int vi, int dim);
    inline float dotProductSum(double x1, double y1, double z1, double x2, double y2, double z2);

    inline void setPara(string _outDir, float _teta, float _randwb, int _nSmBefLoop, int _nLoop, int _nSmInLoop, int genus);
    inline void setFilename(const char * filename);
    int t_celli;
    float t_teta;
    string t_fname;
    string t_savedir;
    int t_nSmBefLoop;
    int t_nLoop;
    int t_nSmInLoop;
    int t_genus;

    inline bool isVOnLeftSideOfSeg(const Eigen::Vector3f& v, const Eigen::Vector3f& sv1, const Eigen::Vector3f& sv2, const Eigen::Vector3f& n);
};

inline bool Arrangement::isVOnLeftSideOfSeg(const Eigen::Vector3f& v, const Eigen::Vector3f& sv1, const Eigen::Vector3f& sv2, const Eigen::Vector3f& n){
    return (v-sv1).dot(n.cross(sv2-sv1)) > 0;
}

inline bool Arrangement::isSegClosed(int segi){
    return seg2V[segi].front() == seg2V[segi].back();
}
inline bool Arrangement::cycShareSeg(int cyc1, int cyc2){
    return !((cyc2seg_b[cyc1] & cyc2seg_b[cyc2]).isZero());
}
inline float Arrangement::EuclideanDistSquare1Cell2Seg(int c_v, int s_v){
    return pow(Vs[s_v][0] - cellVs[c_v][0], 2)
            + pow(Vs[s_v][1] - cellVs[c_v][1], 2)
            + pow(Vs[s_v][2] - cellVs[c_v][2], 2);
}
inline double Arrangement::getTetVCoord(int ci, int vi, int dim){
    return tetV[ci][vi * 3 + dim];
}
inline float Arrangement::dotProductSum(double u1, double u2, double u3, double v1, double v2, double v3){
    return abs(u1*(v2 - v3) + u2*(v3 - v1) + u3*(v1 - v2));
}

inline void Arrangement::setPara(string _outDir, float _teta, float _randwb, int _nSmBefLoop, int _nLoop, int _nSmInLoop, int genus){
    t_savedir = _outDir;
    t_teta = _teta;
    t_randwb = _randwb;
    t_genus = genus;
    t_nSmBefLoop = _nSmBefLoop;
    t_nLoop = _nLoop;
    t_nSmInLoop = _nSmInLoop;
}

inline void Arrangement::setFilename(const char * _filename){
    string fullname(_filename);

    t_fname = "";
    int pos = fullname.find_last_of(".");
    if (pos != string::npos){
        t_fname = fullname.substr(0, pos);
    }
    pos = t_fname.find_last_of("\\/");
    if (pos != string::npos){
        t_fname.erase(0, pos+1);
    }
}

#endif
