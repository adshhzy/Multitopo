#ifndef GEORF_H
#define GEORF_H

#include<vector>
#include<math.h>
#include"my_mesh.h"
#include"InfoStruct.h"
using namespace std;
using namespace MyUtility;
namespace n_rf {


class Curve{
public:

    int n_vertices, n_edges;
    vector<double>vertices;
    vector<uint>edges;



    vector<double>tangent;

    double disscale = 2.;
    double colordegree = -1;

private:

    vector< uint >vertices2edges;
    vector< uint >vertices2vertices;
    vector< unsigned long >vertices2edges_accumulate;
    vector<uchar>vertices2edgesPN;
    vector< uint >vertices2vertices_inverse;


private:
    bool isloop;

private:

    bool isbuild;
    bool isSettangent;
    bool isBuildFrenet;
    double scale;



    vector<double>display_vertices;
    vector<double>display_vnormals;
    vector<uint>display_edges;
    vector<uint>display_faces;
    vector<unsigned char>v_color;


    double avedis;
    static bool isload;
    static Mesh sphere;
    static Mesh cylinder;
    static Mesh cone;
    int n_ver_f1;
    int n_col_f1;
    int n_face_f1;
    int n_vnor_f1;


    vector<double>linespeed;
    vector<double>curvature;
    vector<double>torsion;
    vector<double>Frenetnormal;
    vector<double>Frenetbinormal;
    vector<double>normCurvature;



public:
    vector<double>edge_len;
    vector<double>vertices_len;
    vector<double>inverse_edge_len;

public:
    vector< vector<uint> >seg2ver;
    vector< vector<uint> >loops2seg;
    vector< uint >keyVer;
    vector< vector<uint> >seg2edges;
    vector< int >superedges;
    vector< uint >seg2Mat;

    vector< vector<uint> >cell2Seg;
    vector< uint >seg2Cells;
    vector<vector<int>>SegNSeg;
    int nCell;
    int nSeg;

public:

    void CurveNetAnalayze(vector<uint> &edgesMat, vector<int>&edges2Cells, int _nCell);
    void CurveNetSimplification(int s, vector<double>&out_ver, vector<uint>&out_edges, vector<uint> &outMat, int deleteIso = 0);
    void CurveNormalParalellTransport(vector<double>&oriNor,vector<double>&rmf);



private:
    bool load();

public:
    void reset();
    void setparameters();


    bool ReadCurve(string filename);
    bool ReadCurve(ifstream &fin, int p0, int p1, int p2);
    bool SaveCurve(string filename);
    bool ImportCurve(vector<double>& vertices_in, vector<uint>& edges_in);
    bool ImportCurve(vector<double>& vertices_in, vector<uint>& edges_in,double lscale, double *pcenter);
    bool ImportCurve(string filename);
    bool ValidationCheck();

    bool AddCurve(Curve& a);




    void Rescale(double &lscale, double *pcenter);
    void RescaleUniform();
    void BuildDisplay(bool isrescale = true, bool isuseball=true, bool isshowtangent = false, bool isspecifiedcolor = false);
    void BuildDisplay(vector<double>&vcolordegree, vector<double>&ecolordegree, int method = 0);

    bool EstimateFirstDerivativeO5(vector<double>& invect, vector<double>& derivative, int interval,bool isloop);

    void BuildEdges(bool isseq = true);
    bool SortLoopEdges();
    void BuildTangent();

    void EstimateFrenetFrame();
    void EliminateDuplication(double thres = 1e-5);
    bool isBuildDisplay();
    bool Isloop(){return isloop;}
    void ChangeLoop(bool isloop){this-> isloop = isloop;}

public:
    void BuildDisplay(infoSet info);
    void CurveScalarFunctionSmoothing(vector<double>&scalarF,int maxiter);

public:
    vector<double>* getDisplayVertices(){return &display_vertices;}
    vector<double>* getDisplayVerticesNormal(){return &display_vnormals;}
    vector<uint>* getDisplayEdges(){return &display_edges;}
    vector<uint>* getDisplayFaces(){return &display_faces;}
    vector<unsigned char>* getDisplayColor(){return &v_color;}
public:
    void BuildToyData(bool isloop);
    void BuildToyData(int ind,int numV,double start,double end);


public:
    Curve();

    void vectorize(int first, int second,double *e){
        auto vdim = 3;
        auto p_v1 = &(vertices[first*vdim]);
        auto p_v2 = &(vertices[second*vdim]);
        minusVec(p_v1,p_v2,e,vdim);
    }
    void vectorize(double *v1, int second,double *e){
        auto vdim = 3;
        auto p_v2 = &(vertices[second*vdim]);
        minusVec(v1,p_v2,e,vdim);
    }
    double VerticesDistance(int first, int second)const{
        double dist = 0.0;
        auto p_v1 = &(vertices[first*3]);
        auto p_v2 = &(vertices[second*3]);
        for (int i = 0; i < 3; ++i)
            dist += pow((p_v1[i] - p_v2[i] ), 2);
        return sqrt(dist);
    }
    double EdgeSquareLen(int eind)const{
        double dist = 0.0;
        auto p_ev = edges.data()+eind*2;
        auto p_v1 = &(vertices[p_ev[0]*3]);
        auto p_v2 = &(vertices[p_ev[1]*3]);
        for (int i = 0; i < 3; ++i)
            dist += pow((p_v1[i] - p_v2[i] ), 2);
        return dist;
    }


    inline double* v_begin(int v_ind){ return &(vertices[v_ind * 3]); }
    inline double* v_end(int v_ind){ return &(vertices[(v_ind+1) * 3]); }
    inline uint* ev_begin(int e_ind){ return &(edges[e_ind * 2]); }
    inline uint* ev_end(int e_ind){ return &(edges[e_ind * 2 + 2]); }
    inline double* t_begin(int v_ind){ return &(tangent[v_ind * 3]); }
    inline double* t_end(int v_ind){ return &(tangent[(v_ind+1) * 3]); }



    inline double get_curvature(int v_ind){ return curvature[v_ind]; }
    inline double get_linespeed(int v_ind){ return linespeed[v_ind]; }
    inline double* nor_begin(int v_ind){ return &Frenetnormal[v_ind*3]; }
    inline double* nor_end(int v_ind){ return &Frenetnormal[v_ind*3+3]; }
    inline double* binor_begin(int v_ind){ return &Frenetbinormal[v_ind*3]; }
    inline double* binor_end(int v_ind){ return &Frenetbinormal[v_ind*3+3]; }


    inline double elen_begin(int e_ind){ return (edge_len[e_ind]); }
    inline double einvlen_begin(int e_ind){ return (inverse_edge_len[e_ind]); }

public:
    inline uint* ev_begin(uint e_ind){ return &(edges[e_ind * 2]); }
    inline uint* ev_end(uint e_ind){ return &(edges[e_ind * 2 + 2]); }
    inline uint* ve_begin(uint v_ind){ return &(vertices2edges[vertices2edges_accumulate[v_ind]]); }
    inline uint* ve_end(uint v_ind){ return &(vertices2edges[vertices2edges_accumulate[v_ind+1]]); }
    inline uint ve_num(uint v_ind){ return vertices2edges_accumulate[v_ind+1]-vertices2edges_accumulate[v_ind]; }
    inline uchar* vepn_begin(uint v_ind){ return &(vertices2edgesPN[vertices2edges_accumulate[v_ind]]); }
    inline uchar* vepn_end(uint v_ind){ return &(vertices2edgesPN[vertices2edges_accumulate[v_ind+1]]); }


    inline uint* vv_begin(uint v_ind){ return &(vertices2vertices[vertices2edges_accumulate[v_ind]]); }
    inline uint* vv_end(uint v_ind){ return &(vertices2vertices[vertices2edges_accumulate[v_ind+1]]); }
    inline uint* vvinv_begin(uint v_ind){ return &(vertices2vertices_inverse[vertices2edges_accumulate[v_ind]]); }
    inline uint* vvinv_end(uint v_ind){ return &(vertices2vertices_inverse[vertices2edges_accumulate[v_ind+1]]); }
    inline uint vv_num(uint v_ind){ return vertices2edges_accumulate[v_ind+1]-vertices2edges_accumulate[v_ind]; }
};

class Contour{
public:
    int n_vertices,n_edges,n_materials;
    double plane_para[4];
    vector<double>vertices;
    vector<uint>edges;
    vector<int>m1;
    vector<int>m2;
public:
    double maxx,maxy;

public:
    Contour():n_vertices(0),n_edges(0),n_materials(0){}
    void ReadContour(ifstream &in);
    void splitContour(vector<Contour> &outC);
    void WriteContour(ofstream &out);
    void CalSelfProperty();
    bool ValidationCheck();
    void BuildContourFromImage(vector<vector<int>>&im,double z);
    void CurveSmoothing(int maxiter);
    void ImportContour(double *para, vector<double>&inCtrV, vector<uint>&inCtrE, vector<int>&inCtrEMat);
    void MergeMaterial(vector<int>&mappingMat);
    void RotatingandShifting(double r_angle, vector<double> &shiftbeforeR, vector<double> &raxis, vector<double> &shiftv);
    void AddContour(Contour &mc);
    void InversePlaneNormal();
    void MergeBoundaryBox(vector<double>&inCtrV, vector<uint>&inCtrE, vector<int>&inCtrEMat);
    void OutputContour(double *para, vector<double>&outCtrV, vector<uint>&outCtrE, vector<int>&outCtrEMat);
    void LableTrianglesWithWindingNumber(vector<double> &testFaceCenter, vector<int>&outfaceMat);
    void ComputeMat2VerticesLoop(vector<vector<double>>&verticesSeqs, vector<int>&loops2mats);

    void WriteContourToFCM(ofstream &out);
    void WriteContourToFCMtCtr(ofstream &out);

    inline uint* ev_begin(uint e_ind){ return &(edges[e_ind * 2]); }
    inline uint* ev_end(uint e_ind){ return &(edges[e_ind * 2 + 2]); }
    inline double* v_begin(int v_ind){ return &(vertices[v_ind * 3]); }
    inline double* v_end(int v_ind){ return &(vertices[(v_ind+1) * 3]); }
    double VerticesDistance(int first, int second)const{
        double dist = 0.0;
        auto p_v1 = &(vertices[first*3]);
        auto p_v2 = &(vertices[second*3]);
        for (int i = 0; i < 3; ++i)
            dist += pow((p_v1[i] - p_v2[i] ), 2);
        return sqrt(dist);
    }
};
class CrossSections{
public:
    int n_materials,n_CrossSec,total_nver,total_nedges;
    vector<Contour>CSs;
    vector<int>acc_nver;
public:
    double width,height;
public:
    CrossSections():n_CrossSec(0){}
    CrossSections(vector<Contour>&a):CSs(a),n_CrossSec(a.size()){
        n_materials=-1;
        for(auto &b:CSs)n_materials = max(n_materials,b.n_materials);
    }
    void naiveConstruction();
    void CalSelfProperty();
    void Stackup(vector<double>&points,vector<uint>&edges);

    void MapMaterials();
    int WriteCrossSections(string filename);
    void ReadCrossSectionsfromImages(vector<vector<vector<int>>>&im, double zdis);
    int SpecialMerging(string infilename, vector<double>&points, vector<uint> &edges);
    void ManipulateCrossSections();
    bool ReadCrossSections(string filename);
    void ShiftCrossSections(vector<vector<double>>&trlvec);

    int WriteCrossSectionsToYaron(string filename);
    int WriteCrossSectionsToFCM(string filename, int nMat);
    int WriteCrossSectionsToFCMtCtr(string filename);

    void WriteCrossSectionToObj(string filename);

    void RescaleToUniform();

    bool ReadCrossSections( vector<vector<double>>&planeinfo, vector<vector<double>>&plane2vertices,
                                vector<vector<uint>>&plane2edges, vector<vector<int>>&plane2edgesMat);


    bool ReadCrossSections( vector<double>&vertices,vector<vector<double>>&planeinfo,
                                vector<vector<uint>>&plane2edges, vector<vector<int>>&plane2edgesMat);

    void SyntheticCircles(int ncircle, int nlayers);

    void SyntheticCircles_multilabelled(int nmMat, int nlayers);

    void SyntheticCircles_nonparallel(int ncircles, int nlayers);
};





class PartialCurve{
private:
    bool isEnableFrames = false;

public:
    Curve mixCurve;
    Curve DisplayCurve;
    Curve MergeCurve;
    Curve FrameCurve;

    vector<int> mergingorder;

    vector<vector<bool>>existSeg;
    vector<vector<uint>>disSeg;
    vector<vector<bool>>mergeSeg;
    vector<vector<uint>>dismergeSeg;

    vector<double>frameVs;
    vector<vector<int>>frameF2Es;
    vector<vector<int>>frameCell2Fs;



public:
    vector<double>display_vertices;
    vector<double>display_vnormals;
    vector<uint>display_edges;
    vector<uint>display_faces;
    vector<unsigned char>display_colors;

public:
//    vector<double>* getDisplayVertices(){return DisplayCurve.getDisplayVertices();}
//    vector<double>* getDisplayVerticesNormal(){return DisplayCurve.getDisplayVerticesNormal();}
//    vector<uint>* getDisplayEdges(){return DisplayCurve.getDisplayEdges();}
//    vector<uint>* getDisplayFaces(){return DisplayCurve.getDisplayFaces();}
//    vector<unsigned char>* getDisplayColor(){return DisplayCurve.getDisplayColor();}

    vector<double>* getDisplayVertices(){return &display_vertices;}
    vector<double>* getDisplayVerticesNormal(){return &display_vnormals;}
    vector<uint>* getDisplayEdges(){return &display_edges;}
    vector<uint>* getDisplayFaces(){return &display_faces;}
    vector<unsigned char>* getDisplayColor(){return &display_colors;}

    int Processing();
    void BuildDisplay(bool isrescale = true, bool isuseball=true,bool isShowExist = true, bool isShowMerge = true);
    int ImportCurve(vector<double>& vertices_in, vector<uint>& edges_in, vector<uint> &edgesMat, vector<int>&edges2Cells, int nCells);
    void Relocation(int visitStep);

    void BuildDisplay(int visitStep,bool isrescale = true, bool isuseball=true){
        Relocation(visitStep);
        BuildDisplay( isrescale ,  isuseball);
    }
    void ReadFrameCurves(string prefix, string fvendfix,string ffendfix, string c2fendfix);
};

//class CurveNet{
//private:
//    string R_filename;
//public:
//    enum{SINGLE,SPECIFIED,MULTI};
//    int mode,pickid;
//    Frame singlecurve;
//    vector<Frame>curves;
//    bool isbuild;
//public:

//    vector<double>display_vertices;
//    vector<uint>display_edges;
//    vector<uint>display_faces;
//    vector<unsigned char>v_color;


//    CurveNet();
//public:
//    int ReadCurve(string filename,int minNum);
//    bool SaveCurve(string filename);
//    bool ReadNormalFile(string filename);

//    void ReadMatlabSol(bool isunified, string fin);
//public:
//    bool GenerateToydata(int ind,int numV,double start,double end);
//    bool excuteRMF(bool isrotate);
//    bool isBuildDisplay(){return isbuild;}
//public:
//    void BuildDisplay(infoSet info);
//    bool compute(infoFrame info);
//public:
//    int MultiAlpha(infoFrame info);

//public:
//    void ReSetConstrain();
//    void InitializePicking(double *startp, double* endp);
//    void SetConstrain(int angle);
//    void DeleteConstrain();
//    void ConfirmConstrain(bool isset);
//    bool IsConstrainPicking();
//    void setPickingVetices(int ind);
//public:
//    void InitializeBezier(int numPoints, int numControlPoints, bool isloop);
//    bool HitControlPoints(double *startp, double* endp);
//    bool MoveControlPoints(double *startp, double *endp);
//    bool EndControlProcessing(bool isvalid);
//    bool isOnControl();
//    bool isPicking();

//    void outputAnglediffN(ofstream &fout);
//    void outputAnglediffN(FILE* fout);
//public:
//    vector<double>* getBezierDisplayVertices();
//    vector<uint>* getBezierDisplayEdges();
//    vector<uint>* getBezierDisplayFaces();
//    vector<unsigned char>* getBezierDisplayColor();

//public:
//    vector<double>* getDisplayVertices();
//    vector<double>* getDisplayVerticesNormals();
//    vector<uint>* getDisplayEdges();
//    vector<uint>* getDisplayFaces();
//    vector<unsigned char>* getDisplayColor();



//};





























}//n_rf

#endif // GEORF_H
