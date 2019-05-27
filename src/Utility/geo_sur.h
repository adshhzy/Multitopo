#ifndef GEO_SUR_H
#define GEO_SUR_H


#include"InfoStruct.h"
#include<vector>
#include<map>
#include<math.h>
#include"my_mesh.h"
#include<eigen3/Eigen/Dense>
using namespace std;
using namespace MyUtility;

namespace n_rf {


class Surface: public Mesh{
public:
    bool reinitflag;
private:
    bool isbuilddisp;
    bool isresetColor;
    bool isFaceRenderMode;
    bool ismarkC;

    int nv = 0;
    int nc = 0;
    int ne = 0;
    int nf = 0;
    int nn = 0;
public:
    int colormethod = 0;
    double colordegree;
    vector<double> invert_faces_colordegree;
    vector<double> invert_faces_colordegree_backside;
    vector<double> invert_vertices_colordegree;
    //vector<uchar> invert_vertices_colordegree;

public:
    vector<double>hidden_ctrV;
    vector<uint>hidden_ctrE;
    vector<double>hidden_ctrVnor;
    vector<uint>hidden_ctrELable;
    vector<int>hidden_ctrE2Cells;
    vector<int>hidden_F2Cells;
    vector<int>hidden_F2CSs;
    vector<uint>hidden_ctrActiveE;
    vector<uint>hidden_ctrActiveE2V;
    vector<uint>hidden_ctrActiveF2V;

    vector<int>hidden_ctrVCreator;
    vector<int>hidden_ctrActiveNonE2V;
    map<int,vector<float>>mapActiveF2Vpos;
    map<int,vector<float>>mapActiveNonE2Vpos;
    vector<int>hidden_writeCtrE;
    vector<int>hidden_CtrE2CSs;

    vector<int>hiddenF2Mat;
    vector<int>hiddenV2Mat;
    vector<double>hidden_colordegreeT;
    vector<int>hiddencontainer_ctrE;
private:
    static bool isload;
    static Mesh sphere;

    int sparsecoff;
private:
    string modelname;
    string prepath;
    string ext;
private:
    vector<double> faces_center;
    vector<double> edge_len;
    vector<double> edge_vec;
    vector<double> edge_cot;
    double ave_edge_len;
public:
    vector<int> face_Material1;
    vector<int> face_Material2;
private:

    vector<double>vertices_area;
    vector<double>faces_area;
    vector<double>faces_inverse_area;

public:

    vector<uchar> weighted_color;

    vector<uchar> weighted_fcolor;
    vector<uchar> weighted_fcolorBackside;
    vector<uchar> sparseShowField;
    vector<double> displayInvertfaceNormal;
private:
    vector<double> display_vertices;
    vector<double> display_normal;
    vector<uint> display_edges;


    vector<uint>display_field_dot;
    vector<unsigned char> display_vcolor;
    vector<uint> display_faces;
private:

    vector<double>id_pick_vertices;
    vector<uint>id_pick_faces;
    vector<unsigned char>id_pick_color;
public:



public:
    void clearup();
    void importSurface(vector<double>&vv, vector<uint>&fv, bool isReOrientate=false, bool isBuild = true, bool isRescale = false);
    void MMpropagation(int f2vmethod,vector<uint>&ctredges,vector<int>&ctredgesMat,vector<int>&facesMat,vector<int>&verticesMat);
    void TopoPreservingLabling(vector<int>&facesMat,vector<int>&verticesMat);
    void TopoPreservingLabling2(vector<int>&facesMat,vector<int>&verticesMat);
    void hiddenCtr(vector<int> &verticesLable, bool isOnlyGatherCtrInfo = false);
    void hiddenCtrSmoothing();
    int importSurface(string inFpre, string insurend, string inf2Cs, string inf2Cells, string infaceMat, string inVMat, string inCtrE = string(), bool isrescale = true);
    void importSurface(vector<double>&vv, vector<uint>&fv, vector<int>&vmat, vector<int>&f2Cell, vector<int> &f2CS);
    void importSurface(vector<double>&vv, vector<uint>&fv, vector<int>&vmat, bool isReOrientate=false, bool isRescale = false);

    int importSurface_specialcontainer(int codeN, string inFpre, string insurend, string inf2Cells, string infaceMat, string inCtrE);


    void importSurface(vector<double>&vv, vector<uint>&fv, vector<int>&vmat, vector<int>&fmat, vector<double> &colordegreeT);
    void exportSurface(int csi,vector<double>&vv, vector<uint>&fv, vector<int>&vmat, vector<int>&fmat, vector<double> &colordegreeT);

    void exportSurfaceCell(int celli,vector<double>&vv, vector<uint>&fv, vector<int>&fmat, vector<int>&ctrE,vector<double> &colordegreeT);
private:
    bool saveObjx(string filename);
    void load();
    bool readSufFile(string filename);


    void sparseSampling(int a);


    //bool readEobjfile(string filename);
    //bool readObjxfile(string filename);


    void BuildDisplay(int colormethod, double colordegree, bool isfield, bool isnormal, bool issurf, bool iswire, bool issingularity, bool ismark, int length, int width, int upnormal);
    void BuildFacesCenter();
    void ComputeEdgeLength();
    void ComputeArea();

    void WeightColor(const double thres,const vector<double>&weights,vector<uchar>&out_color);
    void meshScalarfunctionSmoothing(vector<double>& scalarF,int maxiter);
    void testSlicer(int thres);
private:
    void BuildPickID();
public:
    vector<double>* getPickIDVertices(){return &id_pick_vertices;}
    vector<uint>* getPickIDFaces(){return &id_pick_faces;}
    vector<unsigned char>* getPickIDColor(){return &id_pick_color;}
public:
    int PickFaceViaColor(unsigned char* pcolor);



public:

    void SetDisplayTransparency(uchar alpha);
    void BuildDisplay(infoSurfDisp info);
    int Initialize(string filename);
    bool ReadFile(string filename);
    bool SaveInterface(string filename);

public:
    Surface();
public:
    vector<double>* getDisplayVertices(){return &display_vertices;}
    vector<double>* getDisplayVerticesNormal(){return &display_normal;}
    vector<uint>* getDisplayEdges(){return &display_edges;}
    vector<uint>* getDisplayFaces(){return &display_faces;}
    vector<unsigned char>* getDisplayColor(){return &display_vcolor;}


    bool isBuildDisp(){return isbuilddisp;}

    void GetPerFaceColorDegree(vector<double>&facedegree);
    void GetPerFaceColorDegreeBackside(vector<double>&facedegree);

    void SetDisplayNormal(vector<bool>&isinvertFnormal);

    void GetPerVertexColorDegree(vector<double>&verticesdegree);
public:



    inline double* fcent_begin(int f_ind){ return &(faces_center[f_ind * 3]); }
    inline double* fcent_end(int f_ind){ return &(faces_center[(f_ind+1) * 3]); }



    inline double elen_begin(int e_ind){ return (edge_len[e_ind]); }
    inline double *evec_begin(int e_ind){return &(edge_vec[e_ind*3]);}

    inline double ecot_begin(int e_ind){ return (edge_cot[e_ind]); }
    inline double* ecotp_begin(int e_ind){ return &(edge_cot[e_ind]); }


    inline double v_area(int v_ind){return vertices_area[v_ind];}
    inline double f_area(int v_ind){return faces_area[v_ind];}





public:

    void ShowOptInfo(int isacc);

public:

    double computeRotationAlongAxis(const double angle, const double *norAxis, const double *vec, double *vecout);

public:

    void ReOrientateConvexShape();

    void CutSurfaceByPlanePointDistance(double *para, vector<double>&outV2planeDist);
    void CutSurfaceByPlane(double *para,vector<double>&outCtrV,vector<uint>&outCtrE,vector<int>&outCtrEMat);
    void CutSurfaceByBatchPlane(vector<vector<double>>paras, vector<double>&out_CtrV, vector<vector<int>>&out_CtrE, vector<vector<int>>&out_CtrEMat, vector<double> &vnormals);

    void CutSurfaceByPlaneTranslation(double* p_trl, double *para, vector<double>&outCtrV, vector<uint>&outCtrE, vector<int>&outCtrEMat);

public:
    void SpecialRouteforFCM();

};









}// n_rf




#endif // GEO_SUR_H
