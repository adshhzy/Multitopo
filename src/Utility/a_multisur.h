#ifndef A_MULTISUR_H
#define A_MULTISUR_H




#include "geo_sur.h"
#include "geo_curv.h"
namespace n_rf {

class SufStructure{

public:

    int nCell;
    vector<double>Vs;
    vector<uint>Fs;
    vector<int>FMs;
    vector<int>Ctrs;

    vector<vector<int>>mappingsToGlobal;

    void Construct(const vector<double>&vertices,const vector<unsigned int>&faces2vertices,
                   const vector<int>&facesMat, const vector<int>&CtrEdges,
                   const vector<vector<int>>&in_mappingsToGlobal);

};





class MultiSurface{
private:
    enum{S_NORMAL,S_SURFCUT};
    int Cmode;
private:
    int numofSurf;

    int sparsecoff;
private:
    string modelname;
    string prepath;
    string ext;
public:
    Surface mixSurf;

    Curve crossSection;
    Curve NonmanifoldNet;
public:
    vector<int> surf_group1;
    vector<int> surf_group2;
    vector<uint>curvee2v;

    vector<double>colordegree;
private:
    vector<vector<int> > hidden_group;
    vector<vector<int> > hidden_groupv;

private:
    //vector<Curve> NonmanifoldNet;
    vector<Surface> surfbuffer;
    vector< vector<uint> >inversemap;
    vector<vector<double>>Vpos_specialcellmode;
    vector<vector<uint>>CtrE_specialcellmode;
private:
    int n_curveedges;

    vector<uint> curve_inversemap;


    vector<bool>staticV;

    vector<double>smallVnormal;

    double lscale;
    double pcenter[3];




private:
    vector< vector<double>* > display_vertices;
    vector< vector<double>* >display_normal;
    vector< vector<uint>* >display_edges;


    vector< vector<uint>* >display_field_dot;
    vector< vector<unsigned char>* >display_vcolor;
    vector< vector<uint>* >display_faces;

private:
    CrossSections cutCs;

public:
    void WriteCutCrossSection(string filename);
public:

    MultiSurface(){
        Cmode = S_NORMAL;
    }
    void ActivateSufCutMode(){
        Cmode = S_SURFCUT;
    }
    void CutMixSurfaceByPlane(double *para, int sampleStep = 1, int deleteIso = 0);
    void CutMixSurfaceByPlaneTranslation(double *p_trl,double *para,int sampleStep, int deleteIso);
    void CutMixSurfaceByBatchPlane(vector<vector<double> > &paras);
    void DelLastContour();
    void DelAllContour();



    bool clearup();
    bool ReadFile(string filename);
    int splitContour(string infilename, string outfilename, vector<double> &points, vector<uint> &edges);
    int ReadandStackContour(string infilename, vector<double>&points, vector<uint> &edges);



    bool readSufFile(string filename, bool isreArrMat = false);
    bool readAoffFile(string filename, bool isreArrMat = false);
    bool readCslFile(string prefixname, bool isreArrMat = false);

    bool ReAllocation(bool ismanifold = true, int colormethod = 0, char* inputcolors = NULL, bool isImportContour = true);
    bool HiddenGroup_SpecialAllocation();


    int Initialize(string filename, bool defaultpath = false, bool isAutoSmooth = true, bool isSetnMat = false, int nMat = 0, bool isCtrl = false, double lscale = 0, double *pcenter = NULL, int colormethod = 0, char *inputcolor = NULL);
    int ReadBatchResult(string inFpre, string insurend, string incontourend, int n_m);
    int ReadBatchResult(string inFpre, string insurend, string incontourend, string infaceMat, string inVMat, int n_f, bool ismanifold = true);
    void ImportFromSuf(SufStructure &pSuf, bool isAutoSmooth = false, bool isRescale = false, bool isReAllocation = false);
    int ReadBatchResult_specialCellmode(int codeN, string inFpre, string insurend, string inf2Cells, string infaceMat, string inCtrE);

    int InitializeForCutPlane(string filename, bool issmooth = false);
    int InitializeForCutPlane(string filename,string writefilename,double x,double y,double z);
    void MergeMaterials();
    void ChangeMaterials();

    void BuildDisplay(infoSurfDisp info,bool rebuildCs = true);
    void BuildDisplay_specialCellmode(int celli);

    void ComputeSmallMVnormal();
    void SmoothSurfNet_Laplacian(int iter);
    void SmoothSurfNet_JuFair(int iter);

    void SurfaceFFFairing(double lambda,int iter,vector<bool>* staticv);

    void LiepaRefinement(double alpha);

    int GetNumofSurfaces(){return numofSurf;}

    void GetScaleInfoViaCrossSection(double &lscale, double *pcenter);
    void GetScaleInfoViaCrossSection_csl(double &lscale,double *pcenter);

    void GetColorDegree(vector<double>&out_colordegree);


    bool WriteSuf(string filename);
    bool WriteObj(string filename);

    bool WritePickMatObj(string filename,int pickmat);

public:
    vector<vector<double>*>* getDisplayVertices(){return &display_vertices;}
    vector<vector<double>*>* getDisplayVerticesNormal(){return &display_normal;}
    vector<vector<uint>*>* getDisplayEdges(){return &display_edges;}
    vector<vector<uint>*>* getDisplayFaces(){return &display_faces;}
    vector<vector<uchar>*>* getDisplayColor(){return &display_vcolor;}



};






class MultiCellTopo{

public:
    vector< vector<MultiSurface> > CellTopo;

    int nCell;
    int nMat;
    vector<int>nTopo;
    vector<int>curTopo;

    vector<double> mat_colordegree;

    Curve crossS;

    double lscale, pcenter[3];


public:
    void BuildDisplay(infoSurfDisp info, vector<int> &topo, int pickMat);
    void ReadCellTopo(string nameprefix, vector<int>n_Topos,bool isSmooth = true); 
    void ReadCellTopo_picked(string nameprefix, vector<int>&n_Topos, vector<int>&picked_Topos, bool isSmooth = true);

    void GetAllCellTopos(vector<int>&pickTopos, vector<vector<double>>&TopoVs, vector<vector<uint>>&TopoFs, vector<vector<int>>&TopoFMs, vector<vector<uint> > &TopoCtrs);

    void ReReadSingleCellTopo(string nameprefix, int celli, int c_nTopo);

    void CutSurfaceByPlane(int celli,int Topoi,double *para,vector<double>&outCtrV,vector<uint>&outCtrE,vector<int>&outCtrEMat);

    void GetRescaleInfo(double &outlscale, double *outpcenter);

    void GetColorDegree(vector<double>&out_colordegree);
    double GetLabel2Colordegree(int label);
    void UpdateVerticesPosition(vector<double>&newVpos, vector<int>&pickTopos, vector<vector<int> > &mappingToGlobal, vector<bool> &isChanged);
    void GlobalSmoothing();


    void WriteAllSurface(string nameprefix);



private:
    vector< vector<double>* > display_vertices;
    vector< vector<double>* >display_normal;
    vector< vector<uint>* >display_edges;


    vector< vector<uint>* >display_field_dot;
    vector< vector<unsigned char>* >display_vcolor;
    vector< vector<uint>* >display_faces;


public:
    vector<vector<double>*>* getDisplayVertices(){return &display_vertices;}
    vector<vector<double>*>* getDisplayVerticesNormal(){return &display_normal;}
    vector<vector<uint>*>* getDisplayEdges(){return &display_edges;}
    vector<vector<uint>*>* getDisplayFaces(){return &display_faces;}
    vector<vector<uchar>*>* getDisplayColor(){return &display_vcolor;}


};

















}//n_rf















#endif // A_MULTISUR_H
