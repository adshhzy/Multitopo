#include <vector>
#include <unordered_map>
#include <fstream>
#include <map>
#include "Cell.h"
#include "tetgen.h"
#include "../DMCore/Arrangement/Arrangement.h"
#include "../Utility/geo_curv.h"
#include "../Utility/geo_sur.h"
#include "../Utility/a_multisur.h"
#include "MMLevelSet.h"
using namespace std;
typedef unordered_map<int, int> indexmap;
typedef vector<int> rgb;
typedef vector<vector<int>> fGroup;
typedef vector<vector<int>> eGroup;
typedef vector<vector<int>> tGroup;
typedef vector<vector<float>> vGroup;

class ParallelArrangement {
public:
    ParallelArrangement() : _nCell(0), _nMat(0) {}
    ~ParallelArrangement() {}


    //connector

    void loadEntireTet(tetgenio &out);
    void loadInfoFromSpaceDivision(Arrangement &ar, int MMprop_f2vmethod, string &outfoldername);

    // Load parallel cross sections (colored images) and the spacings between
    // adjacent cross sections.
    void loadCrossSections(const vector<string>& filenames,
                           vector<float> spacings, float tetVolLimit);
    //  Load parallel cross sections (colored images) and use uniform spacing to
    //  create the cells.
    void loadCrossSections(const vector<string>& filenames, float spacing,
                           float tetVolLimit);

    void GetCell(int celli, vector<vector<float>>& cellVs,
                 vector<vector<int>>& cellTs, vector<int>& cellInitLabel,
                 int& cellNMat);
    void GetCell(int celli, vector<vector<float>>& cellVs,vector<int>& cellVMarkers,
                 vector<vector<int>>& cellTs, vector<int>& cellInitLabel,
                 int& cellNMat);

    void ReportCrossSectionCurve(int celli, const fGroup& activeFs,
                                 const vGroup& crvVs, const eGroup& crvEs);

    void ExportToMathematica(const char* filename);

    void ExportCrossSectionCurves(const char* filename);
    void ExportCrossSectionCurves(int celli, const char* filename);
    void ExportCrossSectionCurvesNonParallel(const char* filename);

    void CombineInterfaces(const vector<Cell>& cells, const vector<int>& topoIDs,
                           vector<vector<float>>& sufVs,
                           vector<vector<int>>& sufFs,
                           vector<vector<int>>& sufFMats,
                           vector<vector<int>>& segEs);

    void GatherAllActiveFAndCrvVE();

    void HandleFaceletLabels(int f2vmethod, vector<int> &OutfacesMat, vector<int> &ReOrientatedFaces);
    void HandleFaceletLabels2(int f2vmethod, vector<int> &OutfacesMat, vector<int> &ReOrientatedFaces);
    void CtrEdgesRetracing();
    void CtrNetWorkReconstruction(vector<int>&ReOrientatedFaces, vector<int> &facesMat, vector<int> &verticesMat);
    void CtrNetWorkReconstruction2(vector<int>&ReOrientatedFaces, vector<int> &facesMat, vector<int> &verticesMat);

    void OutputCellbox2(vector<int> &ReOrientatedFaces,vector<int> &facesMat);

    void MapComponentLoops2CompactLoop(const int celli, CellTopology &TPwithinASingleCell);
    void MapComponentLoops2CompactLoop(const int celli, const vector<vector<int> > &comps, vector<vector<int>> &comps2compactloop);
    void OutputLoops(int loopi, vector<double>&loopVs, vector<uint>&loopEs, int &label);
    int MakeupLoop2SegsLists();
    void MapLoopToSegs(int celli, vector< vector<int> >&_bdLoopE2V);
    void MapActiveFCP(int celli, vector< int >&_ActFV, vector<vector<float> > &_ActFCP, vector<int>&_ActNonE2V, bool isForMainpipe = true);
    void CreateTotalSurf(vector<vector<double>>&TopoVs, vector<vector<uint>>&TopoFs, vector<vector<int>>&TopoFMs, vector<vector<uint> > &TopoCtrs,
                         bool isOutputFiles = true, n_rf::SufStructure *pSuf = NULL);

    void GetCellTopo(int celli, vector<CellTopology>&cellTopo, vector<vector<int>> &label2bdLoops, vector<int> &cellNexplored2Topogroup);
    int MergeTwoCellTopo(int stepi,CellTopology &preCell,CellTopology &curCell,CellTopology &mergeCell);
    void PrintPickInfo(int celli, int topoi);

    vector<vector<int>>_globalEdges;
    vector<vector<int>>_globalEdgesMat;
    vector<int>_tetEdges;
    vector<int>_tetEdgeMarkers;
    vector<int>_tetVMarkers;
    vector<set<int>>_tetVMarkers_overlap;

    vector<rgb> _cols;

    vector<float>_frameVstack;
    vector<int>_framef2Cs;
    vector<vector<float>> _frameVs;
    vector<vector<int>> _frameFs;
    vector<vector<int>> _frameCells;
    vector<vector<float>> _frameF2PlaneParas;
    vector<vector<float>> _Cs2PlaneParas;
    vector<bool>_framefisCs;
    vector<vector<int>>_frameF2Cell;


    vector<vector<float>> _tetVs;
    vector<fGroup> _frameF2tetFs;  // Boundary tet faces on each frame faces.
    vector<tGroup> _cell2tetTs;
    vector<vGroup> _cell2tetVs;  //  not used in the end.
    vector<indexmap> _cell2MapTetL2CellL;//{map }
    vector<indexmap> _cell2MapTetV2CellV;

    vector<int> _constTetV2Label;

    vector<vector<int>> _cell2labels;//{1,3,n_mat}
    vector<vector<int>> _cell2labeling;//{0,1,2,n_vertices, consensuse}
    vector<vector<int>> _cell2tetVIds;

    // For cross section curves.
    vector<vector<int>> _frameCell2csFrameIDs;
    vector<fGroup> _frameF2ActiveFs;
    vector<vGroup> _frameF2CrvVs;
    vector<eGroup> _frameF2CrvEs;
    vector<bool> _isFrameFReported;

    vGroup _allCrvVs;
    eGroup _allCrvEs;
    map<vector<int>, int> _mapActF2CrvVi;

    int _nCell;
    int _nMat;
    int _nCrossSection;
    int _nFac;
    int _nLoops;

    n_rf::Surface mergeCs;
    n_rf::PartialCurve newCtrNet;


    vector<int>mergeCsV2Tv;
    vector< vector< vector<int> > >loop2Segs;
    vector< vector<int> > seg_dualedges2Tv;
    vector< int > ActiveFV2Tv;
    vector< int > ActiveNonEV2Tv;
    vector< float > ActiveFCP;
    vector< int > ActiveFCP_setup;
    map<int,vector<float>>mapActiveNonE2CP;
    vector<vector<int>>FbdCtrV2ActiveF_percell;
    vector<vector<int>>EbdCtrV2ActiveNonE_percell;
    vector<int>numofActiveF_percell;
    vector<int>numofActiveNonE_percell;


    vector< vector<int> >compactloop2Segs;
    vector< vector<int> >loop2compactloopN;
    vector<int> compactloop2label;

    vector<vector<CellTopology>>_cellTopo;
    vector<vector<CellTopology>>_cellTopo_ori;

    vector< vector< vector<int> > > _label2bdLoops;


    vector< vector<int> >  _bdLoops2lable;



    vector< vector<vector<int>> >mergeLoop2Seg;
    vector< vector<vector<int>> >newloop2mergeLoop;
    vector< vector<vector<int>> >oldloop2mergeLoop;
    vector< vector<vector<int>> >label2mergeLoop;
    vector< vector<int> >mergeLoop2label;

    vector< vector<vector<vector<int>>> >label2mergeloopComp;
    vector< vector<vector<int>> >label2mergeLoopCompVE;
    vector< vector<int> >mergeLoop2VE;

    vector<vector<bool>>mergeloopVoid;

    vector<int>pickTopoInd;
    vector<int>maxExistStep;

    vector<vector<int>>_cellNexplored2Topogroup;


    vector<int>DP_label2NComp;
    vector<vector<int>>DP_label2CGenus;
    vector<int>DP_constraintTopos;
    vector<vector<vector<int>>>DP_label2Cloops;

    vector<vector<int>>DP_groupping;
    vector<int>DP_grouppinggenus;
    vector<int>DP_loops2groupping;

    double time_stage1;

    string outfoldername;

public:

    void DynamicProgramming(string protocalname,vector<int>&constraintTopos);
    void DynamicProgramming(vector<int> &pickTopo);
    bool DynamicProgramming();




    void MergeSameTopo(vector<CellTopology>&topos, int stepi);
    void LoopMergingDP();
    int CheckComponentNum(vector<vector<int> > &pickSet, vector<vector<int> > &preSetSeg);
    bool CheckValidation(int stepi,CellTopology &curCell);
    void PrintMergeTopoInfo(int stepi,CellTopology &curCell);
    void PrintBeMergeTopoInfo(int stepi,CellTopology &curCell);

    void Comp2label(vector<vector<int>>&comps, vector<int> &loop2label, vector<int> &relabel);
    void FindLoopDP(int stepi, vector<int> &pickSeg, vector<vector<int> > &reloop);
    void MakeContourfromImage(const vector<string>& filenames, string outfilename,double zdis);

    void OutputNewCtrs(string filename, const vector<double> &vertices, const vector<uint> &edge2vertices);

    void OutputMappingInfo(string filepath);

    void WriteToMappingInfo();

    void OutputDPInfo(string filepath, vector<int> pickTopoIndex, double totalscore, double time, vector<int> isoComp2VEF, vector<int> isoComps2label);
};
