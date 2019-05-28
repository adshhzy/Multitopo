#include "a_multisur.h"
#include <stdio.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<limits>
#include<cstdlib>
#include <functional>
#include<time.h>
#include<string>
#include <sstream>
#include <iterator>
#include<set>
#include<map>
#include "tetgen.h"
#include "readers.h"
#include "../MMCore/ParallelArrangement.h"
#include "../DMCore/LuCtr2Suf/JuFair/meshFair.h"
#include "../MMCore/UnionFind.h"
#include<random>


double randomdouble();

extern ParallelArrangement pArr;
namespace n_rf {

void SufStructure::Construct(const vector<double>&vertices,const vector<unsigned int>&faces2vertices,
                             const vector<int>&facesMat, const vector<int>&CtrEdges,
                             const vector<vector<int>>&in_mappingsToGlobal){

    Vs = vertices;
    Fs = faces2vertices;
    FMs = facesMat;
    Ctrs = CtrEdges;
    mappingsToGlobal = in_mappingsToGlobal;
    nCell = mappingsToGlobal.size();
}


int MultiSurface::Initialize(string filename, bool defaultpath, bool isAutoSmooth,bool isSetnMat, int nMat, bool isCtrl, double lscale, double *pcenter,int colormethod, char* inputcolor){
    clearup();
    bool rere = false;
    if(defaultpath){

        ReadFile("/Users/Research/Geometry/MM/OutputSuf.suf");
        //ReadFile("/Users/Research/Geometry/MM/result/chickenheart/suf_0_5.suf");

        rere = true;
    }
    else{
        rere = ReadFile(filename);
    }
    if(!rere)return -1;

    mixSurf.setparameters();

    mixSurf.BuildNeighborTable_nonmanifold();

    if(ext ==".csl"){

        crossSection.Rescale(this->lscale,this->pcenter);
        //cout<<"crossSection: "<<this->lscale<<endl;

        //mixSurf.GetRescaleInfo(this->lscale,this->pcenter);
        //cout<<"mixSurf: "<<this->lscale<<endl;

        mixSurf.ReScale(this->lscale,this->pcenter);

    }
    else if(ext!=".suf")mixSurf.ReScale_uniform(1.0);

    if(ext==".suf"){


        if(!isCtrl){

            GetScaleInfoViaCrossSection(this->lscale,this->pcenter);

            mixSurf.ReScale(this->lscale,this->pcenter);
            //mixSurf.ReScale_uniform(1.0);
        }
        else mixSurf.ReScale(lscale,pcenter);



        if(isAutoSmooth){
            //ComputeSmallMVnormal();
            //mixSurf.Fairing(true,false,false,50,0.5,0.1,&staticV);
            ComputeSmallMVnormal();
            SurfaceFFFairing(0.5,20,&staticV);
        }


        if(isSetnMat)numofSurf = nMat;



    }

    if(ext ==".csl")ReAllocation(false,colormethod,inputcolor,false);
    else ReAllocation(false,colormethod,inputcolor);
    ComputeSmallMVnormal();

    cout<<"number of Materia: "<<numofSurf<<endl;

    return numofSurf;

}

int MultiSurface::InitializeForCutPlane(string filename,bool issmooth){
    clearup();
    bool rere = false;

    //rere = readSufFile(filename,true);
    rere = ReadFile(filename);
    if(!rere)return -1;

    {
        mixSurf.setparameters();



        MergeMaterials();

        mixSurf.BuildNeighborTable_nonmanifold();

        //ChangeMaterials();
        if(issmooth){
            ComputeSmallMVnormal();
            mixSurf.Fairing(true,false,false,20,0.5,0.1,&staticV);
            ComputeSmallMVnormal();
            SurfaceFFFairing(0.5,100,&staticV);
        }
        //ComputeSmallMVnormal();



        //mixSurf.ReScale(1.0,0.5,0.5);
        mixSurf.ReScale_uniform(1.0);


        ReAllocation();
    }


    cout<<"number of Materia: "<<numofSurf<<endl;

    return numofSurf;


}



int MultiSurface::InitializeForCutPlane(string filename,string writefilename,double x,double y,double z){
    clearup();
    bool rere = false;

    //rere = readSufFile(filename,true);
    rere = ReadFile(filename);
    if(!rere)return -1;
    {
        mixSurf.setparameters();

        mixSurf.BuildNeighborTable_nonmanifold();



        mixSurf.ReScale(x,y,z);

        ReAllocation();
    }

    writeObjFile(writefilename,mixSurf.vertices,mixSurf.faces2vertices);

    cout<<"number of Materia: "<<numofSurf<<endl;

    return numofSurf;


}


void MultiSurface::ImportFromSuf(SufStructure &pSuf,bool isAutoSmooth,bool isRescale, bool isReAllocation){



    cout<<"ImportFromSuf"<<endl;

    clearup();


    mixSurf.vertices=pSuf.Vs;

    mixSurf.faces2vertices=pSuf.Fs;

    mixSurf.setparameters();


    surf_group1.resize(mixSurf.n_faces);
    surf_group2.resize(mixSurf.n_faces);
    for(int i=0;i<mixSurf.n_faces;++i){
        surf_group1[i] = pSuf.FMs[i*2];
        surf_group2[i] = pSuf.FMs[i*2+1];

    }
    curvee2v.clear();
    for(auto a:pSuf.Ctrs)curvee2v.push_back(a);

    n_curveedges = curvee2v.size()/2;


    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;


    set<int>matSet;
    for(auto a:surf_group1)matSet.insert(a);
    for(auto a:surf_group2)matSet.insert(a);

    if(1){
        map<int,int>mappM;
        int newMInd = 0;
        for(auto a:matSet)mappM[a] = newMInd++;
        for(auto &a:surf_group1)a = mappM[a];
        for(auto &a:surf_group2)a = mappM[a];
        numofSurf = matSet.size();

    }else {
        numofSurf = *max_element(matSet.begin(),matSet.end())+1;
    }


    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;


    if(isRescale){
        GetScaleInfoViaCrossSection(this->lscale,this->pcenter);
        mixSurf.ReScale(this->lscale,this->pcenter);
    }
    //cout<<"numofSurf: "<<numofSurf<<endl;

    if(isAutoSmooth){
        //ComputeSmallMVnormal();
        //mixSurf.Fairing(true,false,false,50,0.5,0.1,&staticV);
        mixSurf.BuildNeighborTable_nonmanifold();
        ComputeSmallMVnormal();
        SurfaceFFFairing(0.5,100,&staticV);
        //SurfaceFFFairing(0.5,100,&staticV);
    }

    ComputeSmallMVnormal();

    if(isReAllocation)ReAllocation();


}
bool readVector_uint(string filename,vector<uint>&vvv){
    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }
    int nnum = 0;
    reader>>nnum;
    nnum*=2;
    vvv.resize(nnum);
    for(int i=0;i<nnum;++i)reader>>vvv[i];
    reader.close();
    return true;
}
int MultiSurface::ReadBatchResult(string inFpre, string insurend, string incontourend, int n_m){

    numofSurf = n_m;

    surf_group1.clear();
    surf_group2.clear();
    mixSurf.clearup();
    curvee2v.clear();
    for(int i=1;i<n_m;++i){
        vector<uint>vvv;
        readVector_uint(inFpre + to_string(i) + incontourend,vvv);
        for(auto &a:vvv)a+=mixSurf.n_vertices;
        for(auto a:vvv)curvee2v.push_back(a);


        Mesh a;
        a.readfile(inFpre + to_string(i) + insurend );
        //a.BuildNeighborTable();
        //a.ReOrientFaces();
        mixSurf.addMesh(a);


        int be = surf_group1.size();
        int nf = a.n_faces;
        surf_group1.resize(be+nf);
        for(int j=be;j<surf_group1.size();++j)surf_group1[j] = i;


    }
    surf_group2.resize(surf_group1.size(),0);
    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;

    mixSurf.ReScale_uniform(1.0);
    mixSurf.BuildNeighborTable_nonmanifold();
    cout<<mixSurf.n_faces<<endl;

    //curvee2v.resize(10);
    //for(int i=0;i<curvee2v.size();++i)curvee2v[i]=i;

    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;
    //for(auto a:curvee2v)cout<<a<<' ';cout<<endl;


    ReAllocation();

    return n_m;

}

int MultiSurface::ReadBatchResult(string inFpre, string insurend, string incontourend, string infaceMat, string inVMat, int n_f,bool ismanifold){

    numofSurf = n_f;

    surf_group1.clear();
    surf_group2.clear();
    mixSurf.clearup();
    curvee2v.clear();

    hidden_group.resize(n_f);
    hidden_groupv.resize(n_f);
    for(int i=0;i<n_f;++i){
        vector<uint>vvv;
        readVector_uint(inFpre + to_string(i) + incontourend,vvv);
        for(auto &a:vvv)a+=mixSurf.n_vertices;
        for(auto a:vvv)curvee2v.push_back(a);


        Mesh a;
        a.readfile(inFpre + to_string(i) + insurend );
        //a.BuildNeighborTable();
        //a.ReOrientFaces();
        mixSurf.addMesh(a);


        int be = surf_group1.size();
        int nf = a.n_faces;
        surf_group1.resize(be+nf);
        for(int j=be;j<surf_group1.size();++j)surf_group1[j] = i;


        readVecFile(inFpre + to_string(i) + infaceMat,hidden_group[i]);
        readVecFile(inFpre + to_string(i) + inVMat,hidden_groupv[i]);



    }
    surf_group2.resize(surf_group1.size(),0);

    mixSurf.ReScale_uniform(1.0);
    mixSurf.BuildNeighborTable_nonmanifold();
    cout<<mixSurf.n_faces<<endl;

    //curvee2v.resize(10);
    //for(int i=0;i<curvee2v.size();++i)curvee2v[i]=i;

    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;
    //for(auto a:curvee2v)cout<<a<<' ';cout<<endl;



    ReAllocation();
    HiddenGroup_SpecialAllocation();




    return n_f;
}

int MultiSurface::ReadBatchResult_specialCellmode(int codeN, string inFpre, string insurend, string inf2Cells, string infaceMat, string inCtrE){



    codeN = 200000;
    mixSurf.readfile(inFpre + to_string(codeN) + insurend );
    mixSurf.ReScale_uniform(1.0);

    mixSurf.BuildUp(false);
    //hiddenF2Mat.clear();hiddenV2Mat.clear();hidden_F2CSs.clear();hidden_F2Cells.clear();hiddencontainer_ctrE.clear();
    vector<int>hiddenF2Mat,verticesMat,faces2Cs,hidden_F2Cells,hiddencontainer_ctrE;
    readVecFile(inFpre + to_string(codeN) + infaceMat,hiddenF2Mat);
    //readVecFile(inFpre + to_string(codeN) + inVMat,hiddenV2Mat);
    //readVecFile(inFpre + to_string(codeN) + inf2Cs,hidden_F2CSs);
    readVecFile(inFpre + to_string(codeN) + inf2Cells,hidden_F2Cells);
    if(!inCtrE.empty())readContourEdgeTxtFile(inFpre + to_string(codeN) + inCtrE,hiddencontainer_ctrE);

    surf_group1 = hiddenF2Mat;surf_group2 = hiddenF2Mat;

    int nMat = *max_element(hiddenF2Mat.begin(),hiddenF2Mat.end())+1;
    colordegree.resize(nMat);
    for(int i=0;i<nMat;++i){
        colordegree[i] = 340*double(i)/double(nMat);
    }
    colordegree[0] = -1;






    numofSurf = *max_element(hidden_F2Cells.begin(),hidden_F2Cells.end())+1;
    vector< vector<uint> >mm_faces(numofSurf);
    vector< vector<uint> >mm_facesmapping(numofSurf);
    vector< vector<uint> >inver_faces_mat(numofSurf);
    for(int i =0;i<mixSurf.n_faces;i++){

        int m1 = hidden_F2Cells[i*2];
        int m2 = hidden_F2Cells[i*2 +1 ];
        auto p_fv = mixSurf.fv_begin(i);


        //if(i>8093)for(int j=0;j<9;++j)cout<<p_fv[j]<<endl;
        if(m1!=-1)for(int j=0;j<3;++j){mm_faces[m1].push_back(p_fv[j]);}



        //if(i>8093)cout<<i<<' '<<m1<<' '<<m2<<endl;
        if(m2!=-1)for(int j=0;j<3;++j){mm_faces[m2].push_back(p_fv[j]);}
        //if(m2!=-1){mm_faces[m2].push_back(p_fv[2]);mm_faces[m2].push_back(p_fv[1]);mm_faces[m2].push_back(p_fv[0]);}

        if(m1!=-1)mm_facesmapping[m1].push_back(i);
        if(m2!=-1)mm_facesmapping[m2].push_back(i);

        if(1){
            if(m1!=-1)inver_faces_mat[m1].push_back(hiddenF2Mat[i]);
            if(m2!=-1)inver_faces_mat[m2].push_back(hiddenF2Mat[i]);
        }

    }

    cout<<"ReadBatchResult_specialCellmode debug 1"<<endl;
    vector< int > newverticeInd (mixSurf.n_vertices,-1);
    inversemap.clear();surfbuffer.clear();Vpos_specialcellmode.clear();CtrE_specialcellmode.clear();
    inversemap.resize(numofSurf);
    surfbuffer.resize(numofSurf);
    Vpos_specialcellmode.resize(numofSurf+1);CtrE_specialcellmode.resize(numofSurf+1);

    int numofPickup = 0;
    vector<uint>&e2v = CtrE_specialcellmode[0];
    vector<double>&vpos = Vpos_specialcellmode[0];

    for(auto &a:newverticeInd)a = -1;
    for(auto a:hiddencontainer_ctrE)newverticeInd[a] = 0;
    for(int j=0;j< newverticeInd.size();++j)if(newverticeInd[j]==0){
        newverticeInd[j] = numofPickup++;
        auto p_v = mixSurf.v_begin(j);
        for(int k=0;k<3;++k)vpos.push_back(p_v[k]);
    }
    for(auto a:hiddencontainer_ctrE)e2v.push_back(newverticeInd[a]);

    for(int i =0;i<numofSurf;i++){

        auto &pinversemap = inversemap[i];
        auto &psurfbuffer = surfbuffer[i];
        auto &pmm_faces = mm_faces[i];
        for(auto &a:newverticeInd)a = -1;
        pinversemap.clear();
        int numofPickup =0;

        for(auto a:pmm_faces)newverticeInd[a] = 0;
        for(auto &a : newverticeInd)if(a==0){
            a = numofPickup;
            ++numofPickup;
        }
        pinversemap.resize(numofPickup);
        for(int j=0;j<newverticeInd.size();++j)if(newverticeInd[j]!=-1){
            pinversemap[newverticeInd[j]] = j;
        }

        for(auto &a:pmm_faces)a = newverticeInd[a];
        psurfbuffer.faces2vertices = pmm_faces;

        psurfbuffer.vertices.resize(numofPickup*3);
        for(int j = 0;j<numofPickup;++j){
            copyVec(mixSurf.v_begin(pinversemap[j]),psurfbuffer.v_begin(j));
        }



        numofPickup = 0;
        vector<uint>&e2v = CtrE_specialcellmode[i+1];
        vector<double>&vpos = Vpos_specialcellmode[i+1];
        for(int j=0,endj = hiddencontainer_ctrE.size()/2;j<endj;++j){
            int a1 = hiddencontainer_ctrE[j*2];
            int a2 = hiddencontainer_ctrE[j*2+1];

            if(newverticeInd[a1]>=0 && newverticeInd[a2]>=0){
                e2v.push_back(a1);
                e2v.push_back(a2);
            }
        }

        for(auto &a:newverticeInd)a = -1;
        for(auto a:e2v)newverticeInd[a] = 0;
        for(int j=0;j< newverticeInd.size();++j)if(newverticeInd[j]!=-1){
            newverticeInd[j] = numofPickup++;
            auto p_v = mixSurf.v_begin(j);
            for(int k=0;k<3;++k)vpos.push_back(p_v[k]);
        }
        for(auto &a:e2v)a = newverticeInd[a];

    }

    for(int i=0;i<numofSurf;++i){
        auto &a = surfbuffer[i];

        a.setparameters();
        a.ReOrientFaces();

        a.colordegree = colordegree[i];
        a.colormethod = 2;

        auto &c_inver_faces_mat = inver_faces_mat[i];
        vector<double>inver_faces_colordegree(c_inver_faces_mat.size());
        for(int j=0;j<c_inver_faces_mat.size();++j)inver_faces_colordegree[j] = colordegree[c_inver_faces_mat[j]];
        a.GetPerFaceColorDegree(inver_faces_colordegree);

        a.BuildUp(false);
        a.reinitflag = true;

        //cout<<a.weighted_fcolor.size()<<endl;
        //a.isbuilddisp = false;
    }

    for(int i =0;i<numofSurf;i++){

        auto &pinversemap = inversemap[i];
        auto &psurfbuffer = surfbuffer[i];
        auto &pfacesmapping = mm_facesmapping[i];


        for(int j=0;j<pfacesmapping.size();++j){
            auto p_fvd = mixSurf.fv_begin(pfacesmapping[j]);
            auto p_fvs = psurfbuffer.fv_begin(j);
            for(int k=0;k<3;++k)p_fvd[k] = pinversemap[p_fvs[k]];
        }
    }


    ComputeSmallMVnormal();
    mixSurf.reinitflag = true;
    mixSurf.colordegree = colordegree[0];
    mixSurf.colormethod = 3;
    vector<double>faces_colordegree(surf_group1.size());
    //for(int j=0;j<surf_group1.size();++j)faces_colordegree[j] = colordegree[min(surf_group1[j],surf_group2[j])==0?max(surf_group1[j],surf_group2[j]):min(surf_group1[j],surf_group2[j])];
    for(int j=0;j<surf_group1.size();++j)faces_colordegree[j] = colordegree[hiddenF2Mat[j]];
    mixSurf.GetPerFaceColorDegree(faces_colordegree);
    for(int j=0;j<surf_group2.size();++j)faces_colordegree[j] = colordegree[hiddenF2Mat[j]];
    mixSurf.GetPerFaceColorDegreeBackside(faces_colordegree);


    return numofSurf;




}










void MultiSurface::SmoothSurfNet_Laplacian(int iter){
    mixSurf.SurfaceLaplacianFairing(0.5,iter,&staticV);

    ReAllocation();
}
void MultiSurface::SmoothSurfNet_JuFair(int iter){

    SurfaceFFFairing(0.5,iter,&staticV);

    ReAllocation();
}

void MultiSurface::SurfaceFFFairing(double lambda,int iter,vector<bool>* staticv){

    vector<bool>isInverseNor(mixSurf.n_faces);
    for(int i=0;i<mixSurf.n_faces;++i){
        if(surf_group1[i]>surf_group2[i])isInverseNor[i] = false;
        else isInverseNor[i] = true;
    }

    mixSurf.SurfaceFFFairing(lambda,iter,&isInverseNor,&staticV);


}

void MultiSurface::LiepaRefinement(double alpha){

    LuMesh m_meshFair;
    vector<float> m_verLu(mixSurf.vertices.size());
    intvector m_faceLu(mixSurf.faces2vertices.size());
    intvector m_faceMat(surf_group1.size()*2);
    intvector m_ctrmedgeLu(curvee2v.size());
    bool m_doSwap = true;

    auto tmpcurv = curvee2v;

    auto &mverOri = mixSurf.vertices;
    auto &mfaceOri = mixSurf.faces2vertices;
    for(int i=0;i<m_verLu.size();++i)m_verLu[i] = mverOri[i];
    for(int i=0;i<m_faceLu.size();++i)m_faceLu[i] = mfaceOri[i];
    for(int i=0;i<m_ctrmedgeLu.size();++i)m_ctrmedgeLu[i] = curvee2v[i];
    for(int i=0;i<surf_group1.size();++i){m_faceMat[i*2] = surf_group1[i];m_faceMat[i*2+1] = surf_group2[i];}

    //m_ctrmedgeLu.clear();
    m_meshFair.InputData(m_verLu,m_faceLu,m_faceMat,m_ctrmedgeLu, m_doSwap);

    cout<<"Liepa refineMent Start!"<<endl;

    m_meshFair.LiepaRefine(alpha);

    cout<<"Liepa refineMent finish!"<<endl;


    m_meshFair.OutputData(m_verLu,m_faceLu,m_faceMat);

    clearup();

    mverOri.resize(m_verLu.size());
    mfaceOri.resize(m_faceLu.size());
    for(int i=0;i<m_verLu.size();++i)mverOri[i] = m_verLu[i];
    for(int i=0;i<m_faceLu.size();++i)mfaceOri[i] = m_faceLu[i];

    surf_group1.resize(m_faceMat.size()/2);surf_group2.resize(m_faceMat.size()/2);
    for(int i=0;i<surf_group1.size();++i){ surf_group1[i] = m_faceMat[i*2];surf_group2[i] = m_faceMat[i*2+1];}






    //writeObjFile(string("../liepatest"),mverOri,mfaceOri);
    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;

    mixSurf.setparameters();
    mixSurf.BuildNeighborTable_nonmanifold();

    cout<<mixSurf.n_faces<<' '<<surf_group1.size()<<endl;

    curvee2v = tmpcurv;
    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;



    ReAllocation();




}
bool MultiSurface::clearup(){

    mixSurf.clearup();

    surf_group1.clear();
    surf_group2.clear();


    surfbuffer.clear();
    inversemap.clear();


    curvee2v.clear();
    curve_inversemap.clear();

    crossSection.reset();
    NonmanifoldNet.reset();

    staticV.clear();

    Vpos_specialcellmode.clear();
    CtrE_specialcellmode.clear();




    display_vertices.clear();
    display_normal.clear();
    display_edges.clear();


    display_field_dot.clear();
    display_vcolor.clear();
    display_faces.clear();



}


bool MultiSurface::ReAllocation(bool ismanifold, int colormethod, char *inputcolors,bool isImportContour){


    //cout<<"ReAllocation"<<endl;
    vector< int > newverticeInd (mixSurf.n_vertices,-1);
    for(auto a:curvee2v)newverticeInd[a] = 0;

    int numofPickup = 0;
    for(auto &a : newverticeInd)if(a==0){
        a = numofPickup;
        ++numofPickup;
    }
    curve_inversemap.resize(numofPickup);
    for(int j=0;j<newverticeInd.size();++j)if(newverticeInd[j]!=-1){
        curve_inversemap[newverticeInd[j]] = j;
    }

    vector<uint>e2v;
    for(auto &a:curvee2v)e2v.push_back(newverticeInd[a]);

    vector<double>veveve(numofPickup*3);
    for(int j = 0;j<numofPickup;++j){
        copyVec(mixSurf.v_begin(curve_inversemap[j]),veveve.data()+3*j);
    }

    if(isImportContour)crossSection.ImportCurve(veveve,e2v);

    int nonvind = 0;
    veveve.clear();
    vector<int>nonv(mixSurf.n_vertices,-1);
    vector<uint>nonedges;
    if(0){
        for(int i=0;i<mixSurf.n_vertices;++i)if(mixSurf.vertices_non_manifold[i]){
            nonv[i] = nonvind++;
            auto p_v = mixSurf.v_begin(i);
            for(int j=0;j<3;++j)veveve.push_back(p_v[j]);
        }
        for(int i=0;i<mixSurf.n_edges;++i)if(mixSurf.edge_non_manifold[i]){
            auto p_ev = mixSurf.ev_begin(i);
            nonedges.push_back(nonv[p_ev[0]]);nonedges.push_back(nonv[p_ev[1]]);
        }

        NonmanifoldNet.ImportCurve(veveve,nonedges);
    }




    //cout<<"ReAllocation "<<mixSurf.n_faces<<endl;
    //for(auto a: mixSurf.faces2vertices)cout<<a<<' ';cout<<endl;
    //cout<<mixSurf.faces2vertices.size()/3-mixSurf.n_faces<<endl;

    vector< vector<uint> >mm_faces(numofSurf);
    vector< vector<uint> >inver_faces_mat(numofSurf);
    for(int i =0;i<mixSurf.n_faces;i++){

        int m1 = surf_group1[i];
        int m2 = surf_group2[i];
        auto p_fv = mixSurf.fv_begin(i);


        //if(i>8093)for(int j=0;j<9;++j)cout<<p_fv[j]<<endl;
        for(int j=0;j<3;++j){mm_faces[m1].push_back(p_fv[j]);}



        //if(i>8093)cout<<i<<' '<<m1<<' '<<m2<<endl;
        //for(int j=0;j<3;++j){mm_faces[m2].push_back(p_fv[j]);}
        mm_faces[m2].push_back(p_fv[2]);mm_faces[m2].push_back(p_fv[1]);mm_faces[m2].push_back(p_fv[0]);

        if(1){
            inver_faces_mat[m1].push_back(m2);
            inver_faces_mat[m2].push_back(m1);
        }else{
            inver_faces_mat[m1].push_back(m1==0?m2:(m2==0?m1:m2));
            inver_faces_mat[m2].push_back(m2==0?m1:(m1==0?m2:m1));
        }

    }



    //cout<<"reallocation debug 1"<<endl;
    inversemap.resize(numofSurf);
    surfbuffer.resize(numofSurf);
    for(int i =0;i<numofSurf;i++){

        auto &pinversemap = inversemap[i];
        auto &psurfbuffer = surfbuffer[i];
        auto &pmm_faces = mm_faces[i];
        for(auto &a:newverticeInd)a = -1;
        pinversemap.clear();
        int numofPickup =0;

        for(auto a:pmm_faces)newverticeInd[a] = 0;
        for(auto &a : newverticeInd)if(a==0){
            a = numofPickup;
            ++numofPickup;
        }
        pinversemap.resize(numofPickup);
        for(int j=0;j<newverticeInd.size();++j)if(newverticeInd[j]!=-1){
            pinversemap[newverticeInd[j]] = j;
        }

        for(auto &a:pmm_faces)a = newverticeInd[a];
        psurfbuffer.faces2vertices = pmm_faces;

        psurfbuffer.vertices.resize(numofPickup*3);
        for(int j = 0;j<numofPickup;++j){
            copyVec(mixSurf.v_begin(pinversemap[j]),psurfbuffer.v_begin(j));
        }
    }

    colordegree.resize(numofSurf);
    //colormethod = 0;
    if(colormethod ==0){
        for(int i=0;i<numofSurf;++i){
            colordegree[i] = 340*double(i)/double(numofSurf);
        }
    }if(colormethod ==1){
        for(int i=0;i<numofSurf;++i){
            colordegree[i] = 340*double(i)/double(numofSurf-1);
        }
    }
    colordegree[0] = -1;
   // colordegree[1] = 174;


    for(int i=0;i<numofSurf;++i){
        auto &a = surfbuffer[i];




        a.setparameters();

        a.colordegree = colordegree[i];
        a.colormethod = 1;

        auto &c_inver_faces_mat = inver_faces_mat[i];
        vector<double>inver_faces_colordegree(c_inver_faces_mat.size());
        for(int j=0;j<c_inver_faces_mat.size();++j)inver_faces_colordegree[j] = colordegree[c_inver_faces_mat[j]];
        a.GetPerFaceColorDegree(inver_faces_colordegree);

        if(1){a.BuildUp(false);}
        else{
            a.ComputeFaceNormal(false);
            //a.GetRescaleInfo(lscale, pcenter);
        }
        a.reinitflag = true;

        //cout<<a.weighted_fcolor.size()<<endl;
        //a.isbuilddisp = false;
    }
    ComputeSmallMVnormal();
    mixSurf.reinitflag = true;
    mixSurf.colordegree = colordegree[0];
    mixSurf.colormethod = 3;
    vector<double>faces_colordegree(surf_group1.size());
    //for(int j=0;j<surf_group1.size();++j)faces_colordegree[j] = colordegree[min(surf_group1[j],surf_group2[j])==0?max(surf_group1[j],surf_group2[j]):min(surf_group1[j],surf_group2[j])];
    for(int j=0;j<surf_group1.size();++j)faces_colordegree[j] = colordegree[min(surf_group1[j],surf_group2[j])];
    mixSurf.GetPerFaceColorDegree(faces_colordegree);
    //mixSurf.GetPerFaceColorDegreeBackside(faces_colordegree);


    for(int j=0;j<surf_group2.size();++j)faces_colordegree[j] = colordegree[max(surf_group1[j],surf_group2[j])];
    mixSurf.GetPerFaceColorDegreeBackside(faces_colordegree);
    //mixSurf.GetPerFaceColorDegree(faces_colordegree);




    //for(int i=0;i<2;++i)surfbuffer[2].ReOrientFaces();
    //surfbuffer[2].SaveInterface(string("ads"));

    //    for(int i=0;i<numofSurf;++i){
    //        cout<<surfbuffer[i].weighted_fcolor.size()<<endl;
    //    }
    //cout<<"reallocation debug end "<<surfbuffer.size()<<endl;

}


bool MultiSurface::HiddenGroup_SpecialAllocation(){


    cout<<"HiddenGroup_SpecialAllocation"<<endl;

    vector<uint>maxM(numofSurf,0);
    for(int i=1;i<numofSurf;++i)maxM[i] = *max_element(hidden_group[i].begin(),hidden_group[i].end());
    int n_m = *max_element(maxM.begin(),maxM.end());
    vector<double>colordegree(n_m);
    for(int i=0;i<n_m;++i){
        colordegree[i] = 360*double(i)/double(n_m+1)+90;
    }
    for(int i=1;i<numofSurf;++i){
        auto &a = surfbuffer[i];

        a.setparameters();

        a.colordegree = 50;
        a.reinitflag = true;


        vector<double>inver_faces_colordegree(a.n_faces);
        for(int j=0;j<inver_faces_colordegree.size();++j)inver_faces_colordegree[j] = colordegree[hidden_group[i][j]];
        a.GetPerFaceColorDegree(inver_faces_colordegree);

        vector<double>inver_vertices_colordegree(a.n_vertices);
        for(int j=0;j<inver_vertices_colordegree.size();++j)inver_vertices_colordegree[j] = colordegree[hidden_groupv[i][j]];
        a.GetPerVertexColorDegree(inver_vertices_colordegree);

        a.hiddenCtr(hidden_groupv[i]);




        //cout<<a.weighted_fcolor.size()<<endl;



        //a.isbuilddisp = false;
    }

}


bool MultiSurface::ReadFile(string filename){



    string filepath;
    string modname;
    string extname;

    SplitFileName(filename,filepath,modname,extname);
    cout<<modname<<' '<<extname<<' '<<filepath<<endl;

    bool issuc = false;


    //issuc = ReadFile(filename);

    if(extname==".suf"){issuc = readSufFile(filename);}
    else if(extname==".aoff"){issuc = readAoffFile(filename,true);}
    else if(extname==".csl"){issuc = readCslFile(filepath+modname,true);}
    else {
        //        if(extname==".obj" || extname==".off"){
        //            surfbuffer.resize(1);
        //            issuc = surfbuffer[0].ReadFile(filename);

        //        }
        issuc = mixSurf.readfile(filename);
        surf_group1.clear();
        surf_group2.clear();
        surf_group1.resize(mixSurf.n_faces,1);
        surf_group2.resize(mixSurf.n_faces,0);
        mixSurf.face_Material1 = surf_group1;
        mixSurf.face_Material2 = surf_group2;

        numofSurf=2;
    }

    if(issuc){
        modelname = modname;
        prepath = filepath;
        ext=extname;
    }

    return issuc;
}


bool MultiSurface::readSufFile(string filename,bool isreArrMat){

    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Suf file " << filename << endl;
        return false;
    }
    mixSurf.clearup();

    cout<<"Reading Suf File"<<endl;



    reader>>mixSurf.n_vertices;
    reader>>mixSurf.n_faces;

    mixSurf.vertices.resize(mixSurf.n_vertices*3);

    mixSurf.faces2vertices.clear();

    surf_group1.clear();
    surf_group2.clear();

    for(int i =0;i<mixSurf.vertices.size();i++){
        reader>>mixSurf.vertices[i];
    }

    int ivalue;
    for(int i =0;i<mixSurf.n_faces;i++){
        for(int j =0;j<3;++j){reader>>ivalue;mixSurf.faces2vertices.push_back(ivalue);}
        reader>>ivalue;surf_group1.push_back(ivalue);
        reader>>ivalue;surf_group2.push_back(ivalue);
    }

    reader>>n_curveedges;
    curvee2v.resize(n_curveedges*2);
    for(int i =0;i<curvee2v.size();i++){
        reader>>curvee2v[i];
    }

    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;


    set<int>matSet;
    for(auto a:surf_group1)matSet.insert(a);
    for(auto a:surf_group2)matSet.insert(a);

    if(isreArrMat){
        map<int,int>mappM;
        int newMInd = 0;
        for(auto a:matSet)mappM[a] = newMInd++;
        for(auto &a:surf_group1)a = mappM[a];
        for(auto &a:surf_group2)a = mappM[a];
        numofSurf = matSet.size();

    }else {
        numofSurf = *max_element(matSet.begin(),matSet.end())+1;
    }


    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;












    cout<<"numofSurf: "<<numofSurf<<endl;
    return true;






}

bool MultiSurface::WriteSuf(string filename){

    vector<int >fmat;
    for(int i=0;i<mixSurf.n_faces;++i){
        fmat.push_back(surf_group1[i]);fmat.push_back(surf_group2[i]);
    }

    vector<int>ctre;
    for(auto a:curvee2v)ctre.push_back(a);

    return writeSufFile(filename,mixSurf.vertices,mixSurf.faces2vertices,fmat,ctre);




}

bool MultiSurface::WriteObj(string filename){




    vector<uint>faces2vertices;
    for(int i=0;i<mixSurf.n_faces;++i){
        auto p_fv = mixSurf.fv_begin(i);
        if(surf_group1[i]<surf_group2[i]){
            for(int j=0;j<3;++j)faces2vertices.push_back(p_fv[j]);

        }else{
            for(int j=2;j>=0;--j)faces2vertices.push_back(p_fv[j]);

        }



    }

    return writeObjFile(filename,mixSurf.vertices,faces2vertices);

}

bool MultiSurface::WritePickMatObj(string filename,int pickmat){



//    vector<double>newvertices;
//    vector<uint>faces2vertices;
//    vector<int>pickvertices(mixSurf.n_vertices,-1);
//    for(int i=0;i<mixSurf.n_faces;++i){

//        if(surf_group1[i]==pickmat || surf_group2[i] == pickmat){
//            auto p_fv = mixSurf.fv_begin(i);
//            for(int j=0;j<3;++j)pickvertices[p_fv[j]] = 0;
//            if(surf_group1[i]<surf_group2[i]){
//                for(int j=0;j<3;++j)faces2vertices.push_back(p_fv[j]);
//            }else{
//                for(int j=2;j>=0;--j)faces2vertices.push_back(p_fv[j]);
//            }
//        }

//    }
//    for(int i=0;i<mixSurf.n_vertices;++i){


//    }

    Surface &picksurf = surfbuffer[pickmat];

    return writeObjFile(filename,picksurf.vertices,picksurf.faces2vertices);

}

bool MultiSurface::readAoffFile(string filename,bool isreArrMat){

    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Aoff file " << filename << endl;
        return false;
    }
    mixSurf.clearup();

    cout<<"Reading Aoff File"<<endl;


    string titl;
    int tmp;


    reader>>titl;

    reader>>mixSurf.n_vertices;
    reader>>mixSurf.n_faces;

    reader>>tmp;
    mixSurf.vertices.resize(mixSurf.n_vertices*3);

    mixSurf.faces2vertices.clear();

    surf_group1.clear();
    surf_group2.clear();

    for(int i =0;i<mixSurf.vertices.size();i++){
        reader>>mixSurf.vertices[i];
    }

    int ivalue;
    for(int i =0;i<mixSurf.n_faces;i++){
        reader>>tmp;
        for(int j =0;j<3;++j){reader>>ivalue;mixSurf.faces2vertices.push_back(ivalue);}

    }
    reader>>titl;
    for(int i =0;i<mixSurf.n_faces;i++){
        reader>>ivalue;surf_group1.push_back(ivalue);
        reader>>ivalue;surf_group2.push_back(ivalue);
    }

    curvee2v.clear();


    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;


    set<int>matSet;
    for(auto a:surf_group1)matSet.insert(a);
    for(auto a:surf_group2)matSet.insert(a);

    if(isreArrMat){
        map<int,int>mappM;
        int newMInd = 0;
        for(auto a:matSet)mappM[a] = newMInd++;
        for(auto &a:surf_group1)a = mappM[a];
        for(auto &a:surf_group2)a = mappM[a];
        numofSurf = matSet.size();

    }else {
        numofSurf = *max_element(matSet.begin(),matSet.end())+1;
    }


    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;










    cout<<"numofSurf: "<<numofSurf<<endl;
    return true;






}



bool MultiSurface::readCslFile(string prefixname,bool isreArrMat){


    string cslfilename = prefixname + string(".csl");
    ifstream reader(cslfilename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Csl file " << cslfilename << endl;
        return false;
    }
    mixSurf.clearup();

    cout<<"Reading Aoff File"<<endl;


    string header;
    reader>>header;
    int nP, nMat;
    reader>>nP; reader>>nMat;
    vector<vector<double>>planeinfo(nP);
    vector<int>plane2nV(nP);
    vector<int>plane2nE(nP);
    vector<vector<double>>plane2vertices(nP);
    vector<vector<unsigned int>>plane2edges(nP);
    vector<vector<int>>plane2edgesMat(nP);


    for(int i=0;i<nP;++i){
        int val;
        int ncomp;
        reader>>val;
        reader>>plane2nV[i];
        reader>>ncomp;
        auto &p_planeinfo = planeinfo[i];
        p_planeinfo.resize(4);
        for(int j=0;j<4;++j){
            reader>>p_planeinfo[j];
        }
        //        reader>>plane2nV[i];
        //        reader>>plane2nE[i];
        auto &p_plane2vertices = plane2vertices[i];
        p_plane2vertices.resize(plane2nV[i]*3);
        for(int nv = plane2nV[i]*3,j=0;j<nv;++j){
            reader>>p_plane2vertices[j];
        }
        reader>>plane2nE[i];
        reader>>val;
        auto &p_plane2edges = plane2edges[i];
        auto &p_plane2edgesMat = plane2edgesMat[i];
        //        p_plane2edges.resize(plane2nE[i]*2);
        //        p_plane2edgesMat.resize(plane2nE[i]*2);
        int val1,val2;
        reader>>val1;
        for(int j=1;j<plane2nE[i];++j){
            reader>>val2;
            p_plane2edges.push_back(val1);p_plane2edges.push_back(val2);
            p_plane2edgesMat.push_back(0);p_plane2edgesMat.push_back(1);

            val1 = val2;
        }

        p_plane2edges.push_back(val1);p_plane2edges.push_back(p_plane2edges[0]);
        p_plane2edgesMat.push_back(0);p_plane2edgesMat.push_back(1);

    }



    /*******************************************************/
    /*******************************************************/


    CrossSections CS;
    CS.ReadCrossSections(planeinfo,plane2vertices,plane2edges,plane2edgesMat);
    vector<double>vvvv;
    vector<uint>e2v;
    CS.Stackup(vvvv,e2v);
    crossSection.ImportCurve(vvvv,e2v);




    mixSurf.readOfffile(prefixname + string(".off"));


    curvee2v.clear();
    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;


    surf_group1.clear();
    surf_group2.clear();
    surf_group1.resize(mixSurf.n_faces,1);
    surf_group2.resize(mixSurf.n_faces,0);
    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;

    numofSurf=2;

    cout<<"numofSurf: "<<numofSurf<<endl;
    return true;

}


void MultiSurface::GetScaleInfoViaCrossSection(double &lscale,double *pcenter){
    for(int i=0;i<3;++i)pcenter[i] = 0;
    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    int nSv = 0;
    for(int i=0;i< mixSurf.n_vertices;++i)if(staticV[i]){
        auto point = mixSurf.v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);

        for(int j=0;j<3;++j)pcenter[j] += point[j];
        ++nSv;
    }
    for(int j=0;j<3;++j)pcenter[j] /=nSv;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    lscale = 1/largestdis;

    cout<<pcenter[0]<<"  "<<pcenter[1]<<"  "<<pcenter[2]<<"   "<<largestdis<<endl;

}

void MultiSurface::GetScaleInfoViaCrossSection_csl(double &lscale,double *pcenter){
    for(int i=0;i<3;++i)pcenter[i] = 0;
    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    int nSv = 0;
    for(int i=0;i< crossSection.n_vertices;++i){
        auto point = crossSection.v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);

        for(int j=0;j<3;++j)pcenter[j] += point[j];
        ++nSv;
    }
    for(int j=0;j<3;++j)pcenter[j] /=nSv;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    lscale = 1/largestdis;

    cout<<pcenter[0]<<"  "<<pcenter[1]<<"  "<<pcenter[2]<<"   "<<largestdis<<endl;

}



void MultiSurface::GetColorDegree(vector<double>&out_colordegree){


    out_colordegree = colordegree;

}



void MultiSurface::WriteCutCrossSection(string filename){

    if(Cmode!=S_SURFCUT)return;

    cout<<"write!!!!"<<endl;

    //cutCs.WriteCrossSectionsToYaron(string("../rootinfo/"));

    //cutCs.ManipulateCrossSections();
    //cutCs.WriteCrossSections("../ConvertFolder/root_");

    //cutCs.WriteCrossSections("../ConvertFolder/silicium_");

    //cutCs.WriteCrossSectionsToFCM("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/cutmode/tctrloop.txt",numofSurf);

    //cutCs.WriteCrossSectionsToFCMtCtr("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/cutmode/tctr.txt");

    //cutCs.WriteCrossSectionsToFCM("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/mousebraincutforsurface/tctrloop.txt",numofSurf);

    //cutCs.WriteCrossSectionsToFCMtCtr("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/mousebraincutforsurface/tctr.txt");

    //cutCs.WriteCrossSectionsToFCM("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/newliver/tctrloop.txt",numofSurf);

    //cutCs.WriteCrossSectionsToFCMtCtr("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/newliver/tctr.txt");

    //cutCs.WriteCrossSectionToObj(filename);


    //cutCs.WriteCrossSections("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/mousebraincut/mousebraincut");
    cutCs.WriteCrossSections("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/newliver/newliver");

    //cutCs.WriteCrossSectionToObj("/Users/Research/Geometry/RBF/data/DataSet/surfacecut/torus.obj");
}
void MultiSurface::DelLastContour(){
    if(Cmode!=S_SURFCUT)return;
    auto &CSs = cutCs.CSs;
    if(CSs.size()>0){
        CSs.resize(CSs.size()-1);
        vector<double>outCtrV;vector<uint>outCtrE;vector<int>outCtrEMat;
        cutCs.Stackup(outCtrV,outCtrE);



        cout<<"CSs.size(): "<<CSs.size()<<endl;
        //for(int i=0;i<4;++i)cout<<para[i]<<' ';cout<<endl;
        crossSection.ImportCurve(outCtrV,outCtrE);

        crossSection.BuildDisplay(false,true);
    }


}
void MultiSurface::DelAllContour(){
    if(Cmode!=S_SURFCUT)return;
    auto &CSs = cutCs.CSs;
    if(CSs.size()>0){
        CSs.resize(0);
        vector<double>outCtrV;vector<uint>outCtrE;vector<int>outCtrEMat;
        cutCs.Stackup(outCtrV,outCtrE);
        cout<<"CSs.size(): "<<CSs.size()<<endl;
        //for(int i=0;i<4;++i)cout<<para[i]<<' ';cout<<endl;
        crossSection.ImportCurve(outCtrV,outCtrE);

        crossSection.BuildDisplay(false,true);
    }


}
void MultiSurface::CutMixSurfaceByPlane(double *para,int sampleStep, int deleteIso){

    if(Cmode!=S_SURFCUT)return;
    auto &CSs = cutCs.CSs;
    CSs.resize(CSs.size()+1);
    Contour &ccCs = CSs[CSs.size()-1];
    vector<double>outCtrV;vector<uint>outCtrE;vector<int>outCtrEMat;
    mixSurf.CutSurfaceByPlane(para,outCtrV,outCtrE,outCtrEMat);

    if(sampleStep>1){
        int n_edges = outCtrE.size()/2;

        vector<uint>edgesMat;
        vector<int>edges2Cell(n_edges*2,0);
        for(int i=0;i<n_edges;++i)edges2Cell[i*2]=1;
        for(auto a:outCtrEMat)edgesMat.push_back(a);
        Curve simCurve;
        simCurve.ImportCurve(outCtrV,outCtrE);
        simCurve.CurveNetAnalayze(edgesMat,edges2Cell,2);
        simCurve.CurveNetSimplification(sampleStep,outCtrV,outCtrE,edgesMat,deleteIso);
        outCtrEMat.clear();
        for(auto a:edgesMat)outCtrEMat.push_back(a);
    }

    ccCs.ImportContour(para,outCtrV,outCtrE,outCtrEMat);

    ccCs.CurveSmoothing(200);
    ccCs.ValidationCheck();

    //crossSection.ImportCurve(ccCs.vertices,ccCs.edges);


    cutCs.Stackup(outCtrV,outCtrE);


    cout<<"CSs.size(): "<<CSs.size()<<endl;
    for(int i=0;i<4;++i)cout<<para[i]<<' ';cout<<endl;
    crossSection.ImportCurve(outCtrV,outCtrE);

    crossSection.BuildDisplay(false,true);
}

void MultiSurface::CutMixSurfaceByPlaneTranslation(double *p_trl,double *para,int sampleStep, int deleteIso){

    if(Cmode!=S_SURFCUT)return;
    auto &CSs = cutCs.CSs;
    CSs.resize(CSs.size()+1);
    Contour &ccCs = CSs[CSs.size()-1];
    vector<double>outCtrV;vector<uint>outCtrE;vector<int>outCtrEMat;
    mixSurf.CutSurfaceByPlaneTranslation(p_trl,para,outCtrV,outCtrE,outCtrEMat);

    if(sampleStep>1){
        int n_edges = outCtrE.size()/2;

        vector<uint>edgesMat;
        vector<int>edges2Cell(n_edges*2,0);
        for(int i=0;i<n_edges;++i)edges2Cell[i*2]=1;
        for(auto a:outCtrEMat)edgesMat.push_back(a);
        Curve simCurve;
        simCurve.ImportCurve(outCtrV,outCtrE);
        simCurve.CurveNetAnalayze(edgesMat,edges2Cell,2);
        simCurve.CurveNetSimplification(sampleStep,outCtrV,outCtrE,edgesMat,deleteIso);
        outCtrEMat.clear();
        for(auto a:edgesMat)outCtrEMat.push_back(a);
    }

    ccCs.ImportContour(para,outCtrV,outCtrE,outCtrEMat);

    ccCs.CurveSmoothing(500);
    ccCs.ValidationCheck();

    //crossSection.ImportCurve(ccCs.vertices,ccCs.edges);


    cutCs.Stackup(outCtrV,outCtrE);


    cout<<"CSs.size(): "<<CSs.size()<<endl;
    for(int i=0;i<4;++i)cout<<para[i]<<' ';cout<<endl;
    crossSection.ImportCurve(outCtrV,outCtrE);

    crossSection.BuildDisplay(false,true);
}


void MultiSurface::CutMixSurfaceByBatchPlane(vector<vector<double>>&paras){

    if(Cmode!=S_SURFCUT)return;

    vector<double>outCtrV;vector<vector<int>>outCtrE;vector<vector<int>>outCtrEMat;
    vector<double>outCtrVnor;
    mixSurf.CutSurfaceByBatchPlane(paras,outCtrV,outCtrE,outCtrEMat,outCtrVnor);


    vector<uint>out_CtrE;
    vector<int>out_CtrEMat;
    for(auto &a:outCtrE)for(auto b:a)out_CtrE.push_back(b);
    for(auto &a:outCtrEMat)for(auto b:a)out_CtrEMat.push_back(b);
    vector<double>rmf;
    if(true){
        int n_edges = out_CtrE.size()/2;

        vector<uint>edgesMat;
        vector<int>edges2Cell(n_edges*2,0);
        for(int i=0;i<n_edges;++i)edges2Cell[i*2]=1;
        for(auto a:out_CtrEMat)edgesMat.push_back(a);
        Curve simCurve;
        simCurve.ImportCurve(outCtrV,out_CtrE);
        simCurve.CurveNetAnalayze(edgesMat,edges2Cell,2);
        simCurve.CurveNormalParalellTransport(outCtrVnor,rmf);
        simCurve.CurveNetSimplification(1,outCtrV,out_CtrE,edgesMat,0);
        //outCtrEMat.clear();
        //for(auto a:edgesMat)out_CtrEMat.push_back(a);
    }


    string filenamessss("../output/ferretCut/ferret_6cut");
    writeCurNetFile(filenamessss + string("_exact"),outCtrV,outCtrE,outCtrEMat,paras,outCtrVnor);
    writeCurNetFile(filenamessss + string("_rmf"),outCtrV,outCtrE,outCtrEMat,paras,rmf);



    crossSection.ImportCurve(outCtrV,out_CtrE);

    crossSection.disscale = 0.3;
    //crossSection.tangent = outCtrVnor;
    //crossSection.tangent.clear();
    crossSection.BuildDisplay(false,true,true);
    cout<<"batch cut with normal"<<endl;
}


void MultiSurface::BuildDisplay(infoSurfDisp info, bool rebuildCs){

    //mixSurf.BuildDisplay(info);
    //for(auto &a:surfbuffer)cout<<a.weighted_fcolor.size()<<endl;
    for(auto &a:surfbuffer)a.BuildDisplay(info);
    if(rebuildCs){crossSection.BuildDisplay(false,true,true);}


    mixSurf.BuildDisplay(info);
    //cout<<"aaaaaaaa"<<endl;
    NonmanifoldNet.BuildDisplay(false,true,true);
    display_vertices.resize(surfbuffer.size()+3);
    display_edges.resize(surfbuffer.size()+3);
    display_faces.resize(surfbuffer.size()+3);
    display_normal.resize(surfbuffer.size()+3);
    display_vcolor.resize(surfbuffer.size()+3);

    for(int i=0;i<surfbuffer.size();++i){
        //if(i!=1)continue;
        display_vertices[i+1] = surfbuffer[i].getDisplayVertices();
        display_edges[i+1] = surfbuffer[i].getDisplayEdges();
        display_faces[i+1] = surfbuffer[i].getDisplayFaces();
        display_normal[i+1] = surfbuffer[i].getDisplayVerticesNormal();
        display_vcolor[i+1] = surfbuffer[i].getDisplayColor();
    }
    //    for(int i=0;i<surfbuffer.size();++i){
    //        //if(i!=1)continue;
    //        display_vertices[i+1] = mixSurf.getDisplayVertices();
    //        display_edges[i+1] = mixSurf.getDisplayEdges();
    //        display_faces[i+1] = mixSurf.getDisplayFaces();
    //        display_normal[i+1] = mixSurf.getDisplayVerticesNormal();
    //        display_vcolor[i+1] = mixSurf.getDisplayColor();
    //    }


    //    display_vertices[0] = mixSurf.getDisplayVertices();
    //    display_edges[0] = mixSurf.getDisplayEdges();
    //    display_faces[0] = mixSurf.getDisplayFaces();
    //    display_normal[0] = mixSurf.getDisplayVerticesNormal();
    //    display_vcolor[0] = mixSurf.getDisplayColor();

    display_vertices[0] = crossSection.getDisplayVertices();
    display_edges[0] = crossSection.getDisplayEdges();
    display_faces[0] = crossSection.getDisplayFaces();
    display_normal[0] = crossSection.getDisplayVerticesNormal();
    display_vcolor[0] = crossSection.getDisplayColor();

    int nonind = surfbuffer.size()+1;
    display_vertices[nonind] = NonmanifoldNet.getDisplayVertices();
    display_edges[nonind] = NonmanifoldNet.getDisplayEdges();
    display_faces[nonind] = NonmanifoldNet.getDisplayFaces();
    display_normal[nonind] = NonmanifoldNet.getDisplayVerticesNormal();
    display_vcolor[nonind] = NonmanifoldNet.getDisplayColor();

    int mixind = surfbuffer.size()+2;
    display_vertices[mixind] = mixSurf.getDisplayVertices();
    display_edges[mixind] = mixSurf.getDisplayEdges();
    display_faces[mixind] = mixSurf.getDisplayFaces();
    display_normal[mixind] = mixSurf.getDisplayVerticesNormal();
    display_vcolor[mixind] = mixSurf.getDisplayColor();

    //    display_vertices[mixind] = NonmanifoldNet.getDisplayVertices();
    //    display_edges[mixind] = NonmanifoldNet.getDisplayEdges();
    //    display_faces[mixind] = NonmanifoldNet.getDisplayFaces();
    //    display_normal[mixind] = NonmanifoldNet.getDisplayVerticesNormal();
    //    display_vcolor[mixind] = NonmanifoldNet.getDisplayColor();

    //    mixSurf.getDisplayVertices();
    //    mixSurf.getDisplayEdges();
    //    mixSurf.getDisplayFaces();
    //    mixSurf.getDisplayVerticesNormal();
    //    mixSurf.getDisplayColor();

}
void MultiSurface::BuildDisplay_specialCellmode(int celli){


    //cout<<celli<<endl;

    crossSection.ImportCurve(Vpos_specialcellmode[celli+1],CtrE_specialcellmode[celli+1]);

    crossSection.BuildDisplay(false,true);
    display_vertices[0] = crossSection.getDisplayVertices();
    display_edges[0] = crossSection.getDisplayEdges();
    display_faces[0] = crossSection.getDisplayFaces();
    display_normal[0] = crossSection.getDisplayVerticesNormal();
    display_vcolor[0] = crossSection.getDisplayColor();

}

void MultiSurface::ComputeSmallMVnormal(){

    mixSurf.ComputeFaceNormal(false);
    vector<bool>pisinverN(mixSurf.n_faces);
    for(int i=0;i<mixSurf.n_faces;++i)
        if(surf_group1[i]>surf_group2[i])pisinverN[i]=true;
        else pisinverN[i]=false;
    mixSurf.SetDisplayNormal(pisinverN);
    //mixSurf.ComputeFaceNormal(true,&pisinverN);
    //mixSurf.ComputeFaceNormal(true,NULL);



}
int MultiSurface::ReadandStackContour(string infilename, vector<double>&points, vector<uint> &edges){
    ifstream reader(infilename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Contour file " << infilename << endl;
        return -1;
    }
    int n_c;
    reader>>n_c;

    vector<Contour>mainC(n_c);

    for(int i=0;i<n_c;++i)mainC[i].ReadContour(reader);

    int n_materials = mainC[0].n_materials;




    CrossSections CCC(mainC);
    CCC.Stackup(points,edges);
    return 0;

}

void MultiSurface::MergeMaterials(){



    if(0){
        vector<set<int>>contactness(numofSurf);
        for(int i=0;i<mixSurf.n_faces;++i){
            contactness[surf_group1[i]].insert(surf_group2[i]);
            contactness[surf_group2[i]].insert(surf_group1[i]);
        }

        for(int i=0;i<mixSurf.n_faces;++i){
            cout<<i<<": ";
            for(auto a:contactness[i])cout<<a<<' ';cout<<endl;
        }

        exit(11);
    }
    vector<int>megermap(numofSurf,0);

    if(0){
        for(int i=0;i<megermap.size();++i)megermap[i] = i;



        megermap[1] = 0;megermap[2] = 0;//megermap[5] = 4;
        //    megermap[1] = 0;
        //    megermap[2] = 0;megermap[9] = 0;megermap[12] = 0;

        //    megermap[10] = 7;megermap[17] = 7;megermap[7] = 7;megermap[15] = 7;megermap[16] = 7;megermap[8] = 7;megermap[11] = 7;
        //    megermap[3] = 7;

        //    megermap[4] = 7;megermap[5] = 7;megermap[13] = 7;megermap[14] = 7;


        for(auto &a:surf_group1)a = megermap[a];
        for(auto &a:surf_group2)a = megermap[a];

        vector<bool>restFace(mixSurf.n_faces,false);

        for(int i=0;i<mixSurf.n_faces;++i){
            if(surf_group1[i]!=surf_group2[i])restFace[i]=true;
        }

        vector<int>newVer(mixSurf.n_vertices,-1);
        for(int i=0;i<mixSurf.n_faces;++i)if(restFace[i]){
            auto p_fv= mixSurf.fv_begin(i);
            for(int j=0;j<3;++j)newVer[p_fv[j]] = 0;

        }
        int neInd = 0;
        for(auto &a:newVer)if(a==0)a=neInd++;

        vector<double>verticesPos;
        vector<uint>faces2vertices;

        for(int i=0;i<mixSurf.n_vertices;++i)if(newVer[i]!=-1){
            auto p_v = mixSurf.v_begin(i);
            for(int j=0;j<3;++j)verticesPos.push_back(p_v[j]);
        }
        for(int i=0;i<mixSurf.n_faces;++i)if(restFace[i]){
            auto p_fv= mixSurf.fv_begin(i);
            for(int j=0;j<3;++j)faces2vertices.push_back(newVer[p_fv[j]]);
        }

        vector<int>surf_groupN1,surf_groupN2;
        for(int i=0;i<mixSurf.n_faces;++i)if(restFace[i]){
            surf_groupN1.push_back(surf_group1[i]);
            surf_groupN2.push_back(surf_group2[i]);
        }


        mixSurf.faces2vertices = faces2vertices;
        mixSurf.vertices = verticesPos;
        mixSurf.setparameters();
        surf_group1 = surf_groupN1;
        surf_group2 = surf_groupN2;
    }

    if(1){
        set<int>matSet;
        for(auto a:surf_group1)matSet.insert(a);
        for(auto a:surf_group2)matSet.insert(a);

        map<int,int>mappM;
        int newMInd = 0;
        for(auto a:matSet)mappM[a] = newMInd++;
        for(auto &a:surf_group1)a = mappM[a];
        for(auto &a:surf_group2)a = mappM[a];
        numofSurf = matSet.size();

    }

    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;








}

void MultiSurface::ChangeMaterials(){

    vector<bool>pickF(mixSurf.n_faces,false);
    vector<int>pickFs;
    for(int i=0;i<pickF.size();++i){
        if(surf_group1[i]==0 || surf_group2[i]==0 )pickF[i] = true;
    }
    for(int i=0;i<pickF.size();++i)if(pickF[i]){
        pickFs.push_back(i);
    }
    UnionFind unifind;
    unifind.SetElements(pickFs);
    for(int i=0;i<mixSurf.n_faces;++i)if(pickF[i]){
        auto p_ff = mixSurf.ffnon_begin(i);
        auto p_ffn = mixSurf.ffnon_end(i);
        for(;p_ff!=p_ffn;++p_ff)if(pickF[*p_ff])unifind.Union(*p_ff,i);
    }

    vector<vector<int>>comp;
    unifind.ExtractComponents(comp);
    vector<int>* pc;
    if(comp[0].size()<comp[1].size())pc = &(comp[0]);else pc = &(comp[1]);
    for(int i=0;i<pc->size();++i){
        if(surf_group1[pc->at(i)]==0)surf_group1[pc->at(i)]  = numofSurf;
        else surf_group2[pc->at(i)]  = numofSurf;
    }
    vector<int >fmat;
    for(int i=0;i<mixSurf.n_faces;++i){
        fmat.push_back(surf_group1[i]);fmat.push_back(surf_group2[i]);
    }
    numofSurf++;
    vector<int>ctre;

    writeSufFile("../meshEx/MeshesSuf/hepatic.suf",mixSurf.vertices,mixSurf.faces2vertices,fmat,ctre);


}


int MultiSurface::splitContour(string infilename, string outfilename, vector<double>&points, vector<uint> &edges){

    ifstream reader(infilename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Contour file " << infilename << endl;
        return -1;
    }
    int n_c;
    reader>>n_c;

    vector<Contour>mainC(n_c);

    for(int i=0;i<n_c;++i)mainC[i].ReadContour(reader);

    int n_materials = mainC[0].n_materials;

    vector< vector<Contour> >mainSC(n_c,vector<Contour>(n_materials));


    CrossSections CCC(mainC);
    CCC.Stackup(points,edges);
    //CCC.MapMaterials();
    //CCC.WriteCrossSections(outfilename);
    //CCC.naiveConstruction();

    for(int i=0;i<n_c;++i)mainC[i].splitContour(mainSC[i]);

    vector<int>m_n_c(n_materials,0);

    for(int i=0;i<n_materials;++i)for(int j=0;j<n_c;++j)if(mainSC[j][i].n_vertices!=0)m_n_c[i]++;

    for(int i=0;i<n_materials;++i){
        string outCname = outfilename + string("_m") + to_string(i) + string(".contour");
        ofstream outer(outCname.data(), ofstream::out);
        if (!outer.good()) {
            cout << "Can not open the Contour file " << outCname << endl;
            return -1;
        }
        outer<<m_n_c[i]<<endl;
        for(int j=0;j<n_c;++j)mainSC[j][i].WriteContour(outer);
        outer.close();
    }

    return n_materials;

}
void Contour::ReadContour(ifstream &in){

    for(int i=0;i<4;++i)in>>plane_para[i];
    in>>n_vertices;
    in>>n_edges;

    vertices.resize(n_vertices*3);
    edges.resize(n_edges*2);
    m1.resize(n_edges);
    m2.resize(n_edges);

    for(int i=0;i<vertices.size();++i)in>>vertices[i];
    for(int i=0;i<n_edges;++i){
        in>>edges[i*2];
        in>>edges[i*2+1];
        in>>m1[i];
        in>>m2[i];
    }

    set<int>mmm;
    for(auto a:m1)mmm.insert(a);
    for(auto a:m2)mmm.insert(a);

    if(0)n_materials = mmm.size();
    else n_materials = *max_element(mmm.begin(),mmm.end())+1;
    //for(auto )
    //cin>>n_materials;

}
void Contour::splitContour(vector<Contour> &outC){

    if(n_vertices==0 || n_edges==0 || n_materials<=0)return;

    outC.resize(n_materials);
    for(int i=0;i<n_materials;++i){

        vector<int>vpick(n_vertices,-2);
        vector<bool>epick(n_edges,false);
        auto &pvertices = outC[i].vertices;
        auto &pedges = outC[i].edges;
        auto &pm1 = outC[i].m1;
        auto &pm2 = outC[i].m2;

        for(int j = 0; j<n_edges;++j){
            if(m1[j]==i || m2[j] == i)epick[j] = true;
            else epick[j] = false;
        }

        for(int j = 0; j<n_edges;++j){
            if(epick[j]){
                auto p_ev = ev_begin(j);
                vpick[p_ev[0]] = vpick[p_ev[1]] = -1;
            }
        }

        int vind = 0;
        for(int j = 0; j<n_vertices;++j){
            if(vpick[j]==-1){vpick[j] = vind;++vind;}
        }

        int eind = 0;
        for(int j = 0; j<n_edges;++j){
            if(m1[j]==i || m2[j] == i){epick[j] = true;++eind;}
            else epick[j] = false;
        }

        outC[i].n_vertices = vind;
        outC[i].n_edges = eind;
        outC[i].n_materials = 2;
        pvertices.resize(0);
        pedges.resize(0);

        pm1.resize(0);
        pm2.resize(0);

        for(int j=0;j<4;++j)outC[i].plane_para[j] = plane_para[j];
        for(int j = 0; j<n_edges;++j)if(epick[j]){
            auto p_ev = ev_begin(j);
            pedges.push_back(vpick[p_ev[0]]);
            pedges.push_back(vpick[p_ev[1]]);

            if(i==0){
                pm1.push_back(m1[j]==0?0:1);
                pm2.push_back(m2[j]==0?0:1);
            }else{
                pm1.push_back(m1[j]==i?1:0);
                pm2.push_back(m2[j]==i?1:0);
            }
        }


        for(int j = 0; j<n_vertices;++j){
            if(vpick[j]>-1){
                auto p_v = v_begin(j);
                for(int k=0;k<3;++k)pvertices.push_back(p_v[k]);


            }
        }

    }


}
void Contour::WriteContour(ofstream &out){

    for(int i=0;i<4;++i)out<<plane_para[i]<<' ';out<<endl;
    out<<n_vertices<<' ';
    out<<n_edges<<endl;


    bool direction = 1;

    for(int i=0;i<n_vertices;++i){
        int ind = i*3;
        for(int j=0;j<3;++j)out<<vertices[ind+j]<<' ';
        out<<endl;
    }
    for(int i=0;i<n_edges;++i){
        out<<edges[i*2]<<' ';
        out<<edges[i*2+1]<<' ';
        if(direction){
            out<<m1[i]<<' ';
            out<<m2[i]<<endl;
        }else{
            out<<m2[i]<<' ';
            out<<m1[i]<<endl;
        }
    }

}



void Contour::CalSelfProperty(){
    maxx = maxy = -100000;
    for(int i=0;i<vertices.size();i+=3){
        maxx = max(maxx,vertices[i]);
    }
    for(int i=1;i<vertices.size();i+=3){
        maxy = max(maxy,vertices[i]);
    }
    n_edges = edges.size()/2;
    n_vertices = vertices.size()/3;
}
bool Contour::ValidationCheck(){
    for(auto a:edges)assert(a<n_vertices);
    double baa=0;
    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        if((baa = VerticesDistance(p_ev[0],p_ev[1]))<1e-6){
            cout<<"Contour degenerate edge: "<<baa<<' '<<p_ev[0]<<' '<<p_ev[1]<<endl;
        }
    }
    cout<<"Contour validation check"<<endl;

}

void Contour::BuildContourFromImage(vector<vector<int>>&im, double z){

    if(im.size()==0)return;

    plane_para[0] = 0;plane_para[1] = 0;plane_para[2] = 1;plane_para[3] = -z;

    int height = im.size();
    int width = im[0].size();

    int bw = width*2-1;
    int bh = height*2-1;

    vector<vector<int>>fbuf(bh);
    for(int i=0;i<bh;++i){
        fbuf[i].resize(bw,-1);
    }


    for(int i = 0; i<height-1;++i ){
        auto &buf1 = fbuf[2*i];
        auto &im1 = im[i];
        for(int j = 0; j<width-1;++j ){
            if(im1[j]!=im1[j+1])buf1[j*2+1] = 0;
        }
    }
    for(int i = 0; i<height-1;++i ){
        auto &buf1 = fbuf[2*i+1];
        auto &im1 = im[i];
        auto &im2 = im[i+1];
        for(int j = 0; j<width-1;++j ){
            if(im1[j]!=im2[j])buf1[j*2] = 0;
        }
    }
    for(int i = 0; i<height-1;++i ){
        auto &buf1 = fbuf[2*i];
        auto &buf2 = fbuf[2*i+1];
        auto &buf3 = fbuf[2*i+2];
        for(int j = 1; j<bw;j+=2 ){
            if(buf1[j]==0 || buf3[j]==0)buf2[j] = 0;
        }
        for(int j = 0; j<bw-2;j+=2 ){
            if(buf2[j]==0 || buf2[j+2]==0)buf2[j+1] = 0;
        }
    }
    int vind = 0;
    vertices.clear();
    for(int i=1;i<bh;i+=2){
        auto &buf1 = fbuf[i];
        for(int j = 1; j<bw;j+=2 )if(buf1[j]==0){
            vertices.push_back(j);
            vertices.push_back(i);
            vertices.push_back(z);
            buf1[j] = vind++;
        }
    }

    edges.clear();m1.clear();m2.clear();
    for(int i = 0; i<height-1;++i ){
        auto &buf1 = fbuf[2*i+1];
        auto &im1 = im[i];
        auto &im2 = im[i+1];
        for(int j = 1; j<width-1;++j ){
            if(im1[j]!=im2[j]){
                edges.push_back(buf1[j*2-1]);
                edges.push_back(buf1[j*2+1]);
                m1.push_back(im2[j]);
                m2.push_back(im1[j]);

            }
        }
    }
    for(int i = 1; i<height-1;++i ){
        auto &buf1 = fbuf[2*i-1];
        auto &buf2 = fbuf[2*i+1];
        auto &im1 = im[i];
        for(int j = 0; j<width-1;++j ){
            if(im1[j]!=im1[j+1]){
                edges.push_back(buf1[j*2+1]);
                edges.push_back(buf2[j*2+1]);
                m1.push_back(im1[j]);
                m2.push_back(im1[j+1]);

            }
        }
    }


    n_vertices = vertices.size()/3;
    n_edges = edges.size()/2;
    if(1){
        vector<uint>edgesMat;
        vector<int>edges2Cell(n_edges*2,0);
        for(int i=0;i<n_edges;++i)edges2Cell[i*2]=1;
        for(int i=0;i<n_edges;++i){
            edgesMat.push_back(m1[i]);
            edgesMat.push_back(m2[i]);
        }
        Curve simCurve;
        simCurve.ImportCurve(vertices,edges);
        simCurve.CurveNetAnalayze(edgesMat,edges2Cell,2);
        simCurve.CurveNetSimplification(6,vertices,edges,edgesMat);//6
        n_vertices = vertices.size()/3;
        n_edges = edges.size()/2;
        m1.resize(n_edges);m2.resize(n_edges);
        for(int i=0;i<n_edges;++i){
            m1[i] = edgesMat[2*i];
            m2[i] = edgesMat[2*i+1];
        }
    }
    CurveSmoothing(10800);

}

void Contour::ImportContour(double *para,vector<double>&inCtrV,vector<uint>&inCtrE,vector<int>&inCtrEMat){
    for(int i=0;i<4;++i)plane_para[i] = para[i];
    vertices=inCtrV;
    edges = inCtrE;
    n_vertices = vertices.size()/3;
    n_edges = edges.size()/2;
    m1.resize(n_edges);m2.resize(n_edges);
    for(int i=0;i<n_edges;++i){
        m1[i] = inCtrEMat[2*i];
        m2[i] = inCtrEMat[2*i+1];
    }

}

void Contour::OutputContour(double *para, vector<double>&outCtrV, vector<uint>&outCtrE, vector<int>&outCtrEMat){
    for(int i=0;i<4;++i)para[i] = plane_para[i];
    outCtrV = vertices;
    outCtrE = edges;
    outCtrEMat.clear();
    for(int i=0;i<n_edges;++i){
        outCtrEMat.push_back(m1[i]);
        outCtrEMat.push_back(m2[i]);
    }


}


void Contour::CurveSmoothing(int maxiter){


    vector<vector<int>>vnn(n_vertices);


    for(int i=0;i<n_edges;++i){
        int v1 = edges[i*2],v2 = edges[i*2+1];
        vnn[v1].push_back(v2);
        vnn[v2].push_back(v1);
    }



    vector<double>buffer1(vertices.size());
    vector<double>buffer2(vertices.size());
    vector<double>*pre = &buffer1,*cur = &buffer2;
    *pre = vertices;

    double kaipa = 0.1,lamda = 0.5;
    auto gama = 1 / (kaipa - 1 / lamda);

    auto oneiter = [&vnn,this](double coef,vector<double>*pre,vector<double>*cur){
        double ocoef = 1- coef;
        vector<double>aaa(3);
        for(int i =0;i<n_vertices;++i){
            for(int j=0;j<3;++j)aaa[j] = 0;
            auto &vv = vnn[i];
            if(vnn.size()!=2){
                for(int j=0;j<3;++j)cur->at(i*3+j) = pre->at(i*3+j);
                continue;
            }
            for(int k=0;k<vv.size();++k){
                for(int j=0;j<3;++j){aaa[j]+=pre->at((vv[k]*3)+j);}
            }
            for(int j=0;j<3;++j)aaa[j]/=(vv.size());
            for(int j=0;j<3;++j)cur->at(i*3+j) = ocoef * pre->at(i*3+j) + coef * aaa[j];
        }

    };

    for(int iter=0;iter<maxiter;++iter){

        oneiter(lamda,pre,cur);
        swap(pre,cur);
        oneiter(gama,pre,cur);
        swap(pre,cur);

    }
    vertices=(*pre);




}

void Contour::MergeBoundaryBox(vector<double>&inCtrV, vector<uint>&inCtrE, vector<int>&inCtrEMat){


    if(inCtrV.size()==0 || inCtrE.size()==0)return;
    vector<vector<uint>>v2v(n_vertices);
    vector<vector<uint>>v2e(n_vertices);

    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        v2v[p_ev[0]].push_back(p_ev[1]);
        v2v[p_ev[1]].push_back(p_ev[0]);

        v2e[p_ev[0]].push_back(i);
        v2e[p_ev[1]].push_back(i);
    }

    int nv = inCtrV.size()/3;
    int ne = inCtrE.size()/2;


    vector<vector<uint>>in_v2v(nv);
    vector<vector<uint>>in_v2e(nv);

    for(int i=0;i<ne;++i){
        auto p_ev = inCtrE.data()+i*2;
        in_v2v[p_ev[0]].push_back(p_ev[1]);
        in_v2v[p_ev[1]].push_back(p_ev[0]);

        in_v2e[p_ev[0]].push_back(i);
        in_v2e[p_ev[1]].push_back(i);
    }

    for(auto &a:in_v2v)assert(a.size()==2);






    auto closestPoint = [this,&nv](double *pv, vector<double>&inV){


        vector<double>dist(nv);
        auto p_dv = inV.data();
        for(int i=0;i<nv;++i){
            dist[i] = vecSquareDist(pv,p_dv+i*3);
        }

        return min_element(dist.begin(),dist.end())-dist.begin();

    };



    vector<int>newAddInd(nv,-1);
    vector<int>attachV;
    vector<int>attachinV;
    vector<bool>modifiedV(nv,false);
    vector<bool>modifiedE(ne,false);
    vector<uint>addE;
    vector<int>addEMat;
    for(int i=0;i<n_vertices;++i){
        assert(v2v[i].size()>0);
        if(v2v[i].size()>1)continue;

        int cInd = closestPoint(v_begin(i),inCtrV);

        assert(newAddInd[cInd]==-1);//extreme cases may appear
        modifiedV[cInd] = true;
        newAddInd[cInd] = i;
        attachV.push_back(i);
        attachinV.push_back(cInd);

    }

    cout<<"attachV: "<<attachV.size()<<endl;

    int outsideM = numeric_limits<int>::max();

    if(attachV.size()==0){
        addE = inCtrE;
        addEMat = inCtrEMat;
        for(auto &a:addEMat)if(a==1)a = outsideM;
    }else {


        for(int i=0;i<attachV.size();++i){

            int oriV = attachV[i],inV = attachinV[i];
            int ml,mr;

            int e = v2e[oriV][0];
            if(ev_begin(e)[1]==oriV){ml = m1[e];mr = m2[e];}
            else {ml = m2[e];mr = m1[e];}

            auto &p_ve = in_v2e[inV];
            auto &p_vv = in_v2v[inV];
            for(int j=0;j<2;++j)if(!modifiedV[p_vv[j]]){
                int ine = p_ve[j];
                int nextv = p_vv[j];
                int currv = inV;
                int inml,inmr;
                auto p_ev = inCtrE.data()+ine*2;
                auto p_em = inCtrEMat.data()+ine*2;
                if(p_ev[0]==inV){inml = p_em[0];inmr = p_em[1];}
                else {inml = p_em[1];inmr = p_em[0];}
                bool isleftturn = (inml == 0);

                while(true){
                    addE.push_back(currv);addE.push_back(nextv);
                    if(isleftturn){addEMat.push_back(ml);addEMat.push_back(outsideM);}
                    else {addEMat.push_back(outsideM);addEMat.push_back(mr);}
                    if(modifiedV[nextv])break;
                    modifiedV[nextv] = true;
                    int tmpv = nextv;
                    auto &pvvnext = in_v2v[nextv];
                    nextv = pvvnext[0]==currv?pvvnext[1]:pvvnext[0];
                    currv = tmpv;

                }

            }
        }

    }
    int addInd = n_vertices;
    for(int i=0;i<nv;++i)if(newAddInd[i]==-1){
        newAddInd[i] = addInd++;
        auto p_v = inCtrV.data()+i*3;
        for(int j=0;j<3;++j)vertices.push_back(p_v[j]);
    }
    for(auto &a:addE)a = newAddInd[a];
    edges.insert(edges.end(),addE.begin(),addE.end());
    for(int i=0;i<addEMat.size()/2;++i){
        m1.push_back(addEMat[i*2]);
        m2.push_back(addEMat[i*2+1]);
    }
    CalSelfProperty();
    cout<<"aa"<<endl;


}
void Contour::LableTrianglesWithWindingNumber(vector<double>&testFaceCenter, vector<int>&outfaceMat){


    if(n_vertices == 0 || n_edges == 0){
        outfaceMat.clear();
        outfaceMat.resize(testFaceCenter.size()/3,numeric_limits<int>::max());
        return;
    }
    vector<vector<int>>edgesloops;
    vector<vector<int>>verticesloops;
    vector<int>loops2mat;

    vector<vector<int>>lable2edges;

    vector<vector<uint>>v2v(n_vertices);
    vector<vector<uint>>v2e(n_vertices);
    vector<vector<uint>>e2e(n_edges);

    //map<pair<int,int>,int>edgesmap;
    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        v2v[p_ev[0]].push_back(p_ev[1]);
        v2v[p_ev[1]].push_back(p_ev[0]);

        v2e[p_ev[0]].push_back(i);
        v2e[p_ev[1]].push_back(i);

        //edgesmap[make_pair(p_ev[0],p_ev[1])] = i;
        //edgesmap[make_pair(p_ev[1],p_ev[0])] = i;
    }
    for(auto &v2ee : v2e){
        for(auto a:v2ee)for(auto b:v2ee)if(a!=b)e2e[a].push_back(b);
    }

    set<int>labelset;vector<int>labelvec;
    map<int,int>labelmap;
    for(int i=0;i<n_edges;++i){
        labelset.insert(m1[i]);
        labelset.insert(m2[i]);
    }
    for(auto a:labelset)labelvec.push_back(a);
    for(int i=0;i<labelvec.size();++i)labelmap[labelvec[i]] = i;
    lable2edges.resize(labelvec.size());
    for(int i=0;i<n_edges;++i){
        lable2edges[labelmap[m1[i]]].push_back(i);
        lable2edges[labelmap[m2[i]]].push_back(i);
    }


    for(int i=0;i<lable2edges.size();++i){

        int mat = i;
        vector<int>&ele  = lable2edges[i];

        UnionFind unifind;
        unifind.SetElements(ele);
        vector<bool>chosenEdges( n_edges, false);

        for(auto b:ele)chosenEdges[b] = true;

        for(auto b:ele)for(auto c:e2e[b])if(chosenEdges[c])unifind.Union(b,c);

        vector<vector<int>>tmp_edgesloops;
        unifind.ExtractComponents(tmp_edgesloops);
        for(auto &b:tmp_edgesloops){
            edgesloops.push_back(b);
            loops2mat.push_back(mat);
        }

    }

    verticesloops.resize(loops2mat.size());
    vector<int>countVn(n_vertices,0);
    for(int i=0;i<loops2mat.size();++i){
        vector<int>&v2l  = verticesloops[i];
        vector<int>&ele  = edgesloops[i];
        vector<bool>chosenEdges( n_edges, false);
        vector<bool>chosenVs( n_vertices, false);
        for(auto b:ele)chosenEdges[b] = true;
        for(auto b:ele)for(int j=0;j<2;++j)chosenVs[edges[b*2+j]] = true;

        int pickE = 0;
        int curV,preV;
        if(m2[ele[pickE]] == labelvec[loops2mat[i]]){curV = edges[ele[pickE]*2+1];preV = edges[ele[pickE]*2];}
        else {curV = edges[ele[pickE]*2];preV = edges[ele[pickE]*2+1];}

        int ne = ele.size();
        for(int j=0;j<ne;++j){
            countVn[curV]++;
            v2l.push_back(curV);
            for(auto b:v2v[curV])if(chosenVs[b] && b!=preV){
                preV = curV;
                curV = b;
                break;
            }
        }
        //        ele.clear();
        //        int curE = ele[0],preE = ele[0];
        //        for(int i=0;i<ne;++i){
        //            ele.push_back(curE);
        //            for(auto b:e2e[curE])if(chosenEdges[b] && b!=preE){
        //                preE = curE;
        //                curE = b;
        //                break;
        //            }
        //        }

    }
    //for(int i=0;i<n_vertices;++i)assert(countVn[i]==v2e[i].size());

    //    for(int i=0;i<verticesloops.size();++i){
    //        auto &a=verticesloops[i];
    //        int nv = a.size();
    //        int curV = a[1],preV = a[0];
    //        int gmat = labelvec[loops2mat[i]];
    //        int mat,matinvese;
    //        int kkk = 0;
    //        for(int j=1;j<nv;++j){
    //            int e = edgesmap[make_pair(curV,preV)];

    //            if(edges[e*2]==preV){mat = m2[e];matinvese = m1[e];}
    //            else {mat = m1[e];matinvese = m2[e];}
    //            //assert(gmat==mat);
    //            if(gmat!=mat && gmat==matinvese){
    //                ++kkk;
    //            }
    //            preV = curV;
    //            curV = a[j+1];

    //        }
    //        curV = a[0],preV = a[nv-1];
    //        int e = edgesmap[make_pair(curV,preV)];
    //        if(edges[e*2]==preV){mat = m2[e];matinvese = m1[e];}
    //        else {mat = m1[e];matinvese = m2[e];}
    //        //assert(gmat==mat);
    //        if(gmat!=mat && gmat==matinvese){
    //            ++kkk;
    //        }
    //        cout<<"kkk: "<<kkk<<endl;

    //    }

    //unordered_map<int,double>windingnumber;
    auto calculateWindingnumber = [this](double *p_vpos,vector<int>&vseq){
        int nv = vseq.size();
        vector<double>rays(nv*3);
        auto p_rayd = rays.data();
        for(int i=0;i<nv;++i){
            auto p_rays = p_rayd+i*3;
            minusVec(v_begin(vseq[i]),p_vpos,p_rays);
            normalize(p_rays);
        }
        vector<double>angles(nv);
        double cs[3];
        double isreverse = -1;

        for(int i=0;i<nv-1;++i){
            angles[i] = angleNor(p_rayd+i*3,p_rayd+i*3+3);
            cross(plane_para,p_rayd+i*3,cs);
            if(dot(cs,p_rayd+i*3+3)*isreverse<0)angles[i] = -angles[i];
        }
        angles[nv-1] = angleNor(p_rayd+(nv-1)*3,p_rayd);
        cross(plane_para,p_rayd+(nv-1)*3,cs);
        if(dot(cs,p_rayd)*isreverse<0)angles[nv-1] = -angles[nv-1];

        double windingnumber = 0;
        for(auto a:angles)windingnumber+=a;
        windingnumber/=(2*my_PI);
        return windingnumber;


    };
    auto calculateWindingnumberBatch = [&calculateWindingnumber,this,&loops2mat,&verticesloops](double *p_vpos,vector<double>&outwindingnumber){

        outwindingnumber.resize(loops2mat.size());
        for(int i=0;i<loops2mat.size();++i){
            outwindingnumber[i] = calculateWindingnumber(p_vpos,verticesloops[i]);
        }

    };

    int ntcv = testFaceCenter.size()/3;
    outfaceMat.resize(ntcv);
    for(int i=0;i<ntcv;++i){
        auto p_tcv = testFaceCenter.data()+i*3;
        vector<double>outwindingnumber;
        calculateWindingnumberBatch(p_tcv,outwindingnumber);
        //cout<<"ori: ";for(auto a:outwindingnumber)cout<<a<<' ';cout<<endl;
        vector<double>windingnumberEveryMat(labelvec.size(),0);for(int j=0;j<labelvec.size();++j)windingnumberEveryMat[j] = 0;
        for(int j=0;j<loops2mat.size();++j)windingnumberEveryMat[loops2mat[j]] += outwindingnumber[j];


        windingnumberEveryMat[labelvec.size()-1]+=1;
        //for(auto a:windingnumberEveryMat)cout<<a<<' ';cout<<endl;

        for(auto &a:windingnumberEveryMat)a = fabs(1-a);

        outfaceMat[i] = labelvec[min_element(windingnumberEveryMat.begin(),windingnumberEveryMat.end())-windingnumberEveryMat.begin()];

    }

    cout<<endl;






}



void Contour::ComputeMat2VerticesLoop(vector<vector<double>>&verticesSeqs, vector<int>&loops2mats){


    if(n_vertices == 0 || n_edges == 0){
        verticesSeqs.clear();
        loops2mats.clear();
        return;
    }
    vector<vector<int>>edgesloops;
    vector<vector<int>>verticesloops;
    vector<int>loops2mat;

    vector<vector<int>>lable2edges;

    vector<vector<uint>>v2v(n_vertices);
    vector<vector<uint>>v2e(n_vertices);
    vector<vector<uint>>e2e(n_edges);

    //map<pair<int,int>,int>edgesmap;
    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        v2v[p_ev[0]].push_back(p_ev[1]);
        v2v[p_ev[1]].push_back(p_ev[0]);

        v2e[p_ev[0]].push_back(i);
        v2e[p_ev[1]].push_back(i);

        //edgesmap[make_pair(p_ev[0],p_ev[1])] = i;
        //edgesmap[make_pair(p_ev[1],p_ev[0])] = i;
    }
    for(auto &v2ee : v2e){
        for(auto a:v2ee)for(auto b:v2ee)if(a!=b)e2e[a].push_back(b);
    }

    set<int>labelset;vector<int>labelvec;
    map<int,int>labelmap;
    for(int i=0;i<n_edges;++i){
        labelset.insert(m1[i]);
        labelset.insert(m2[i]);
    }
    for(auto a:labelset)labelvec.push_back(a);
    for(int i=0;i<labelvec.size();++i)labelmap[labelvec[i]] = i;
    lable2edges.resize(labelvec.size());
    for(int i=0;i<n_edges;++i){
        lable2edges[labelmap[m1[i]]].push_back(i);
        lable2edges[labelmap[m2[i]]].push_back(i);
    }


    for(int i=0;i<lable2edges.size();++i){

        int mat = i;
        vector<int>&ele  = lable2edges[i];

        UnionFind unifind;
        unifind.SetElements(ele);
        vector<bool>chosenEdges( n_edges, false);

        for(auto b:ele)chosenEdges[b] = true;

        for(auto b:ele)for(auto c:e2e[b])if(chosenEdges[c])unifind.Union(b,c);

        vector<vector<int>>tmp_edgesloops;
        unifind.ExtractComponents(tmp_edgesloops);
        for(auto &b:tmp_edgesloops){
            edgesloops.push_back(b);
            loops2mat.push_back(mat);
        }

    }

    verticesloops.resize(loops2mat.size());
    vector<int>countVn(n_vertices,0);
    for(int i=0;i<loops2mat.size();++i){
        vector<int>&v2l  = verticesloops[i];
        vector<int>&ele  = edgesloops[i];
        vector<bool>chosenEdges( n_edges, false);
        vector<bool>chosenVs( n_vertices, false);
        for(auto b:ele)chosenEdges[b] = true;
        for(auto b:ele)for(int j=0;j<2;++j)chosenVs[edges[b*2+j]] = true;

        int pickE = 0;
        int curV,preV;
        if(m2[ele[pickE]] == labelvec[loops2mat[i]]){curV = edges[ele[pickE]*2+1];preV = edges[ele[pickE]*2];}
        else {curV = edges[ele[pickE]*2];preV = edges[ele[pickE]*2+1];}

        int ne = ele.size();
        for(int j=0;j<ne;++j){
            countVn[curV]++;
            v2l.push_back(curV);
            for(auto b:v2v[curV])if(chosenVs[b] && b!=preV){
                preV = curV;
                curV = b;
                break;
            }
        }
    }


    loops2mats = loops2mat;
    for(auto &a:loops2mats)a = labelvec[a];
    verticesSeqs.clear();
    verticesSeqs.resize(loops2mats.size());
    for(int i=0;i<loops2mats.size();++i){
        auto &vsq = verticesSeqs[i];
        vsq.clear();
        auto &vlp = verticesloops[i];
        for(int j=0;j<vlp.size();++j){
            auto p_v = v_begin(vlp[j]);
            for(int k=0;k<3;++k)vsq.push_back(p_v[k]);
        }
    }


}



void Contour::MergeMaterial(vector<int>&mappingMat){

    if(1){
        cout<<n_materials<<endl;
        vector<set<int>>contactness(n_materials);
        for(int i=0;i<n_edges;++i){
            contactness[m1[i]].insert(m2[i]);
            contactness[m2[i]].insert(m1[i]);
        }

        for(int i=0;i<n_materials;++i){
            cout<<i<<": ";
            for(auto a:contactness[i])cout<<a<<' ';cout<<endl;
        }

        //exit(11);
    }

    for(auto &a:m1)a = mappingMat[a];
    for(auto &a:m2)a = mappingMat[a];

    vector<bool>restEdge(n_edges,false);

    for(int i=0;i<n_edges;++i){
        if(m1[i]!=m2[i])restEdge[i]=true;
    }

    vector<int>newVer(n_vertices,-1);
    for(int i=0;i<n_edges;++i)if(restEdge[i]){
        auto p_ev= ev_begin(i);
        for(int j=0;j<2;++j)newVer[p_ev[j]] = 0;

    }
    int neInd = 0;
    for(auto &a:newVer)if(a==0)a=neInd++;

    vector<double>verticesPos;
    vector<uint>edges2vertices;

    for(int i=0;i<n_vertices;++i)if(newVer[i]!=-1){
        auto p_v = v_begin(i);
        for(int j=0;j<3;++j)verticesPos.push_back(p_v[j]);
    }
    for(int i=0;i<n_edges;++i)if(restEdge[i]){
        auto p_ev= ev_begin(i);
        for(int j=0;j<2;++j)edges2vertices.push_back(newVer[p_ev[j]]);
    }

    vector<int>_groupN1,_groupN2;
    for(int i=0;i<n_edges;++i)if(restEdge[i]){
        _groupN1.push_back(m1[i]);
        _groupN2.push_back(m2[i]);
    }


    edges = edges2vertices;
    vertices = verticesPos;

    m1 = _groupN1;
    m2 = _groupN2;

    n_vertices = vertices.size()/3;
    n_edges = edges.size()/2;

    CurveSmoothing(100);



}
void Contour::AddContour(Contour &mc){

    int offset = vertices.size()/3;
    for(auto a:mc.vertices)vertices.push_back(a);
    for(auto a:mc.edges)edges.push_back(a+offset);
    for(auto a:mc.m1)m1.push_back(a);
    for(auto a:mc.m2)m2.push_back(a);
    CalSelfProperty();

}
void Contour::InversePlaneNormal(){

    for(int i=0;i<4;++i)plane_para[i] = -plane_para[i];
    swap(m1,m2);

}
void Contour::RotatingandShifting(double r_angle, vector<double>&shiftbeforeR,vector<double>&raxis, vector<double>&shiftv){

    vector<Eigen::Vector3d>ori(n_vertices);
    vector<Eigen::Vector3d>tranformed(n_vertices);

    for(int j = 0;j<n_vertices;++j){
        auto p_spv = v_begin(j);
        for(int k=0;k<3;++k)ori[j](k) = p_spv[k] + shiftbeforeR[k];
    }

    Eigen::Vector3d ra;
    for(int k=0;k<3;++k)ra(k) = raxis[k];


    Eigen::AngleAxisd t(r_angle,ra);


    for(int j = 0;j<n_vertices;++j){
        tranformed[j] = t*ori[j];
    }
    for(int j = 0;j<n_vertices;++j){
        auto p_spv = v_begin(j);
        for(int k=0;k<3;++k)p_spv[k] = tranformed[j](k) - shiftbeforeR[k] + shiftv[k];
    }
    Eigen::Vector3d vnormal;
    for(int i=0;i<3;++i)vnormal(i) = plane_para[i];
    vnormal = t*vnormal;
    for(int i=0;i<3;++i)plane_para[i] = vnormal(i);


    //plane_para[3]-=dot(shiftv.data(),plane_para);
    plane_para[3]=-dot(vertices.data(),plane_para);


}

int CrossSections::SpecialMerging(string infilename, vector<double>&points, vector<uint> &edges){
    ifstream reader(infilename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Contour file " << infilename << endl;
        return -1;
    }
    int n_c;
    reader>>n_c;

    CSs.resize(n_c);


    for(int i=0;i<n_c;++i)CSs[i].ReadContour(reader);

    vector<vector<int>>mappingMat(n_c);
    if(0){
        for(int i=0;i<n_c;++i){
            mappingMat[i].resize(18,0);
            //for(int j=10;j<17;++j)mappingMat[i][j] = 1;
            for(int j=0;j<18;++j)mappingMat[i][j] = j;
        }
        {
            for(int i=4;i<18;++i)mappingMat[0][i] = 3;
            mappingMat[0][15] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[1][i] = 3;
            mappingMat[1][8] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[2][i] = 3;
            mappingMat[2][8] = 8;
            mappingMat[2][7] = 7;
        }

        {
            //for(int i=4;i<18;++i)mappingMat[2][i] = 3;
            mappingMat[3][8] = 0;
            mappingMat[3][10] = 0;
            mappingMat[3][11] = 0;
            mappingMat[3][16] = 6;
        }
        {

            mappingMat[4][12] = 11;
            mappingMat[4][17] = 11;
            mappingMat[4][2] = 6;
            mappingMat[4][10] = 1;

        }
    }

    if(1){
        for(int i=0;i<n_c;++i){
            mappingMat[i].resize(18,0);
            //for(int j=10;j<17;++j)mappingMat[i][j] = 1;
            for(int j=0;j<18;++j)mappingMat[i][j] = j;
        }
        {
            for(int i=4;i<18;++i)mappingMat[0][i] = 1;
            mappingMat[0][13] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[1][i] = 3;
            mappingMat[1][15] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[2][i] = 3;
            mappingMat[2][8] = 8;
        }
        {
            for(int i=4;i<18;++i)mappingMat[3][i] = 3;
            mappingMat[3][8] = 8;
            mappingMat[3][7] = 7;
        }

        {
            //for(int i=4;i<18;++i)mappingMat[2][i] = 3;
            mappingMat[4][8] = 0;
            mappingMat[4][10] = 0;
            mappingMat[4][11] = 0;
            mappingMat[4][16] = 6;
        }
        {

            mappingMat[5][12] = 11;
            mappingMat[5][17] = 11;
            mappingMat[5][2] = 6;
            mappingMat[5][10] = 1;

        }
        {

            mappingMat[6][12] = 11;
            mappingMat[6][17] = 11;
            mappingMat[6][2] = 6;
            mappingMat[6][10] = 1;

        }
    }

    int testCC = 6;
    for(int i=0;i<n_c;++i){
        //if(i!=testCC)continue;
        CSs[i].MergeMaterial(mappingMat[i]);
    }
    int n_materials = CSs[0].n_materials;

    if(0){
        points = CSs[testCC].vertices;
        edges = CSs[testCC].edges;
    }else Stackup(points,edges);
    return 0;


}
void CrossSections::ReadCrossSectionsfromImages(vector<vector<vector<int>>>&im,double zdis){

    n_CrossSec = im.size();

    CSs.resize(n_CrossSec);
    for(int i=0;i<n_CrossSec;++i){

        CSs[i].BuildContourFromImage(im[i],(i+1)*zdis);
    }




}






bool CrossSections::ReadCrossSections(string filename){
    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        return false;
    }

    fin>>n_CrossSec;

    //n_CrossSec = 5;


    CSs.resize(n_CrossSec);
    for(int i=0;i<n_CrossSec;++i){
        CSs[i].ReadContour(fin);
    }
    CalSelfProperty();

    return true;
}
bool CrossSections::ReadCrossSections(vector<vector<double>>&planeinfo, vector<vector<double>>&plane2vertices,
                                      vector<vector<uint> > &plane2edges, vector<vector<int>>&plane2edgesMat){



    n_CrossSec = planeinfo.size();

    CSs.resize(n_CrossSec);
    for(int i=0;i<n_CrossSec;++i){
        CSs[i].ImportContour(planeinfo[i].data(),plane2vertices[i],plane2edges[i],plane2edgesMat[i]);
    }
    CalSelfProperty();

    return true;


}

bool CrossSections::ReadCrossSections( vector<double>&vertices,vector<vector<double>>&planeinfo,
                            vector<vector<uint>>&plane2edges, vector<vector<int>>&plane2edgesMat){

//    vector<int>v_inplane(n_v,-1);
//    int ind = 0;
//    vector<double>Vpos3D;
//    vector<uint>E2V;

//    for(auto a:E[i])v_inplane[a] = 0;
//    for(int j=0;j<n_v;++j)if(v_inplane[j]!=-1){
//        v_inplane[j] = ind++;
//        auto p_v = V.data()+j*3;
//        for(int k=0;k<3;++k)Vpos3D.push_back(p_v[k]);
//    }

//    for(auto a:E[i])E2V.push_back(v_inplane[a]);

    int n_v =    vertices.size()/3;

    n_CrossSec = planeinfo.size();
    CSs.clear();
    CSs.resize(n_CrossSec);
    for(int i=0;i<n_CrossSec;++i){


        vector<int>v_inplane(n_v,-1);

        int ind = 0;
        vector<double>Vp;
        vector<uint>E2V;

        for(auto a:plane2edges[i])v_inplane[a] = 0;
        for(int j=0;j<n_v;++j)if(v_inplane[j]!=-1){
                v_inplane[j] = ind++;
                auto p_v = vertices.data()+j*3;
                for(int k=0;k<3;++k)Vp.push_back(p_v[k]);
        }
        for(auto a:plane2edges[i])E2V.push_back(v_inplane[a]);
        CSs[i].ImportContour(planeinfo[i].data(),Vp,E2V,plane2edgesMat[i]);

    }
    CalSelfProperty();

    return true;



}






void CrossSections::ShiftCrossSections(vector<vector<double>>&trlvec){
    vector<double>rota({0,0,1});
    vector<double>shiftbeforeR({0,0,0});

    for(int i=0;i<CSs.size();++i){

        CSs[i].RotatingandShifting(0,shiftbeforeR,rota,trlvec[i]);
    }


}
void CrossSections::ManipulateCrossSections(){

    //RotatingandShifting

    vector<double>rota({0,0,1});
    vector<double>shiftbeforeR({-500,-500,0});
    vector<double>shiftv({0,0,0});

    if(0){
        srand(0);
        vector<double>shiftbeforeR({0,0,0});
        vector<double>rota({0,1,0});
        vector<Contour>tmpfake = CSs;
        CSs.clear();
        CSs.resize(9);
        CSs[0] = tmpfake[3];
        shiftv[1] = -0.56;
        CSs[0].RotatingandShifting(0,shiftbeforeR,rota,shiftv);

        for(int i=1;i<8;++i)CSs[i] = tmpfake[i-1];
        CSs[8] = tmpfake[3];
        shiftv[1] = 0.56;
        CSs[8].RotatingandShifting(0,shiftbeforeR,rota,shiftv);

        n_CrossSec = 9;
        for(int nre = 0;nre<4;++nre){
            n_CrossSec +=8;
            CSs.resize(n_CrossSec);
            for(int i=n_CrossSec-8;i<n_CrossSec;++i){
                vector<double>shiftbeforeR({0,0,0});
                vector<double>rota({0,1,0});
                CSs[i] = CSs[i-8*(nre+1)];
                //shiftv[2] = 0.01*(randomdouble()-0.5);shiftv[0] = 0.01*(randomdouble()-0.5);

                shiftv[1] = 1.16*(nre+1) + 0.005*(randomdouble());
                CSs[i].RotatingandShifting(0,shiftbeforeR,rota,shiftv);
                // for(int i=0;i<CSs[0].vertices.size();++i)cout<<CSs[0].vertices[i]<<" ";cout<<endl;
            }
        }

    }

    if(1){
        vector<double>shiftbeforeR({-0,-0,0});
        vector<double>shiftv11({0,0,0});
        for(int i=0;i<CSs.size();++i){
            shiftv11[1] = i*170 + 10*(randomdouble());
            shiftv11[0] =  2.5-5*(randomdouble());
            CSs[i].RotatingandShifting(0,shiftbeforeR,rota,shiftv11);
        }

    }

    if(0){
        //n_CrossSec = 2;
        //CSs.resize(n_CrossSec);
        {
            vector<double>shiftbeforeR({-500,-500,0});
            vector<double>rota({0,1,0});
            //CSs[1] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[1].RotatingandShifting(90/180.*3.1415,shiftbeforeR,rota,shiftv);
            // for(int i=0;i<CSs[0].vertices.size();++i)cout<<CSs[0].vertices[i]<<" ";cout<<endl;
        }

    }
    if(0){
        n_CrossSec = 5;
        CSs.resize(n_CrossSec);
        {
            CSs[2] = CSs[1];
            shiftv[2] = 100;shiftv[0] = 0;shiftv[1] = 0;
            CSs[2].RotatingandShifting(20/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        {
            CSs[3] = CSs[1];
            shiftv[2] = 200;shiftv[0] = 30;shiftv[1] = -60;
            CSs[3].RotatingandShifting(20/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        {
            CSs[4] = CSs[0];
            shiftv[2] = 400;shiftv[0] = 60;shiftv[1] = -30;
            CSs[4].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv);
        }
    }

    if(0){
        n_CrossSec = 5;
        CSs.resize(n_CrossSec);
        {
            CSs[1] = CSs[0];
            shiftv[2] = 100;shiftv[0] = 75;shiftv[1] = 70;
            CSs[1].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv);
        }
        {
            CSs[2] = CSs[0];
            shiftv[2] = 200;shiftv[0] = 00;shiftv[1] = 00;
            CSs[2].RotatingandShifting(20/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        {
            CSs[3] = CSs[0];
            shiftv[2] = 300;shiftv[0] = -80;shiftv[1] = -90;
            CSs[3].RotatingandShifting(-20/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        {
            CSs[4] = CSs[0];
            shiftv[2] = 400;shiftv[0] = -10;shiftv[1] = -10;
            CSs[4].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv);
        }

        //        {
        //            CSs[5] = CSs[0];
        //            shiftv[2] = 500;shiftv[0] = 0;shiftv[1] = -0;
        //            CSs[5].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv);
        //        }
    }

    if(0){
        n_CrossSec = 4;
        CSs.resize(n_CrossSec);
        vector<double>shiftbeforeRBB({(-800-1000)/2,-00,-00});
        vector<double>rotaRR({0,1,0});
        vector<double>shiftv11({0,0,-100});
        CSs[0].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv11);

        {
            Contour tmpcc = CSs[0];
            shiftv[2] = 0;shiftv[0] =800;shiftv[1] = 0;
            tmpcc.RotatingandShifting(-180/180.*3.1415,shiftbeforeR,rotaRR,shiftv);
            tmpcc.InversePlaneNormal();
            CSs[0].AddContour(tmpcc);
        }

        if(n_CrossSec>1){
            CSs[1] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[1].RotatingandShifting(90/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }

        if(n_CrossSec>2){
            CSs[2] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[2].RotatingandShifting(45/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }

        if(n_CrossSec>3){
            CSs[3] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[3].RotatingandShifting(135/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }


    }
    if(0){
        n_CrossSec = 3;
        CSs.resize(n_CrossSec);
        int middistance = 1500;
        vector<double>shiftbeforeR({-500,-500,0});
        vector<double>shiftbeforeRBB({0,-00,-00});
        shiftbeforeRBB[0] = (-middistance-1000)/2;
        vector<double>rotaRR({0,1,0});
        vector<double>shiftv11({0,0,-200});
        CSs[0].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv11);

        {
            Contour tmpcc = CSs[0];
            shiftv[2] = 0;shiftv[0] =middistance;shiftv[1] = 0;
            tmpcc.RotatingandShifting(-180/180.*3.1415,shiftbeforeR,rotaRR,shiftv);
            tmpcc.InversePlaneNormal();
            CSs[0].AddContour(tmpcc);
        }

        if(n_CrossSec>1){
            CSs[1] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[1].RotatingandShifting(60/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }

        if(n_CrossSec>2){
            CSs[2] = CSs[0];
            shiftv[2] = 50;shiftv[0] = 50;shiftv[1] = 0;
            CSs[2].RotatingandShifting(120/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }

        if(n_CrossSec>3){
            CSs[3] = CSs[0];
            shiftv[2] = 0;shiftv[0] = 0;shiftv[1] = 0;
            CSs[3].RotatingandShifting(135/180.*3.1415,shiftbeforeRBB,rotaRR,shiftv);
        }


    }

    if(0){
        vector<double>shiftbeforeR({-0,-0,0});
        vector<double>shiftv11({10,0,0});
        for(int i=0;i<CSs.size();++i){
            shiftv11[0] = i*20;
            CSs[i].RotatingandShifting(0,shiftbeforeR,rota,shiftv11);
        }

    }
}


void CrossSections::SyntheticCircles(int ncircle, int nlayers){

    double r_smallcircle = 0.6;//1.2
    double r_bigcircle = 4.0;
    double ratio_coef = 0.5;


    int nv_percicle = 20;


    double step_bigcircle = 2*my_PI / ncircle;
    double step_small_circle = 2*my_PI / nv_percicle;

    CSs.clear();
    CSs.resize(1);
    Contour& onec = CSs[0];

    vector<double>v_smallcircle;
    vector<uint>e_smallcircle;
    vector<int>eMat1_smallcircle(nv_percicle,0),eMat2_smallcircle(nv_percicle,1),eMat3_smallcircle(nv_percicle,2),eMat4_smallcircle(nv_percicle,3);

    for(int i=0;i<nv_percicle;++i){
        v_smallcircle.push_back(r_smallcircle * cos(i*step_small_circle));
        v_smallcircle.push_back(r_smallcircle * sin(i*step_small_circle));
        v_smallcircle.push_back(0);
    }


    for(int i=0;i<nv_percicle-1;++i){
        e_smallcircle.push_back(i);e_smallcircle.push_back(i+1);
    }
    e_smallcircle.push_back(nv_percicle-1);e_smallcircle.push_back(0);


    onec.plane_para[0] = onec.plane_para[1] = onec.plane_para[3] = 0;
    onec.plane_para[2] = -1;

    auto &pv = onec.vertices;
    auto &pe = onec.edges;
    auto &pemat1 =  onec.m1;auto &pemat2 =  onec.m2;
    for(int i=0;i<ncircle;++i){

        double center[3];
        center[0] = r_bigcircle* cos(i*step_bigcircle);
        center[1] = r_bigcircle* sin(i*step_bigcircle);
        center[2] = 0;
        uint offset = pv.size()/3;
        for(int j=0;j<nv_percicle;++j){
            auto p_vsmall = v_smallcircle.data()+j*3;
            for(int k=0;k<3;++k)pv.push_back(p_vsmall[k] + center[k]);
        }
        for(auto e:e_smallcircle)pe.push_back(e+offset);
        pemat1.insert(pemat1.end(),eMat1_smallcircle.begin(),eMat1_smallcircle.end());
        if(1)pemat2.insert(pemat2.end(),eMat2_smallcircle.begin(),eMat2_smallcircle.end());
        if(0){
            if(i<4)pemat2.insert(pemat2.end(),eMat2_smallcircle.begin(),eMat2_smallcircle.end());
            else if(i<8)pemat2.insert(pemat2.end(),eMat3_smallcircle.begin(),eMat3_smallcircle.end());
            else if(i<12)pemat2.insert(pemat2.end(),eMat4_smallcircle.begin(),eMat4_smallcircle.end());
        }
        if(0){
            if(i<6)pemat2.insert(pemat2.end(),eMat2_smallcircle.begin(),eMat2_smallcircle.end());
            else if(i<12)pemat2.insert(pemat2.end(),eMat3_smallcircle.begin(),eMat3_smallcircle.end());
        }
    }

    if(0){
        uint offset = pv.size()/3;
        for(int j=0;j<nv_percicle;++j){
            auto p_vsmall = v_smallcircle.data()+j*3;
            for(int k=0;k<3;++k)pv.push_back(p_vsmall[k]);
        }
        for(auto e:e_smallcircle)pe.push_back(e+offset);
        pemat1.insert(pemat1.end(),eMat1_smallcircle.begin(),eMat1_smallcircle.end());
        pemat2.insert(pemat2.end(),eMat2_smallcircle.begin(),eMat2_smallcircle.end());
    }



    onec.CalSelfProperty();


    n_CrossSec = nlayers;
    vector<double>rota({0,0,1});
    vector<double>shiftbeforeR({0,0,0});

    vector<vector<double>>shiftv(nlayers,vector<double>(3));



    srand(120);
    vector<double>randomoscillation_height(nlayers,0);
    vector<double>randomoscillation_angle(nlayers,0);
    if(1)for(int i=0;i<nlayers;++i){
        randomoscillation_height[i] = ratio_coef*1.4*(randomdouble()-0.5);
        randomoscillation_angle[i] = ratio_coef*10*(randomdouble()-0.5);
    }


    CSs.resize(n_CrossSec);
    for(int i=1;i<nlayers;++i)
    {
        CSs[i] = CSs[0];
        shiftv[i][2] = i*2.5*ratio_coef + randomoscillation_height[i];

        shiftv[i][0] = 0;shiftv[i][1] = 0;

        CSs[i].RotatingandShifting(randomoscillation_angle[i]/180.*3.1415,shiftbeforeR,rota,shiftv[i]);
    }


    CalSelfProperty();


}



void CrossSections::SyntheticCircles_multilabelled(int nmMat, int nlayers){

    double r_smallcircle = 0.9;//0.6
    double r_bigcircle = 4.0;
    double ratio_coef = 1;


    int nv_percicle = 20;


    int ncircle = 6;
    double step_bigcircle = 2*my_PI / ncircle;
    double step_small_circle = 2*my_PI / nv_percicle;

    CSs.clear();
    CSs.resize(1);
    Contour& onec = CSs[0];

    vector<double>v_smallcircle;
    vector<uint>e_smallcircle;
    vector<vector<int>>eMat_smallcircle(nmMat+1);
    for(int i=0;i<nmMat+1;++i){

        eMat_smallcircle[i].clear();
        eMat_smallcircle[i].resize(nv_percicle,i);
    }
    //vector<int>eMat1_smallcircle(nv_percicle,0),eMat2_smallcircle(nv_percicle,1),eMat3_smallcircle(nv_percicle,2),eMat4_smallcircle(nv_percicle,3);

    for(int i=0;i<nv_percicle;++i){
        v_smallcircle.push_back(r_smallcircle * cos(i*step_small_circle));
        v_smallcircle.push_back(r_smallcircle * sin(i*step_small_circle));
        v_smallcircle.push_back(0);
    }


    for(int i=0;i<nv_percicle-1;++i){
        e_smallcircle.push_back(i);e_smallcircle.push_back(i+1);
    }
    e_smallcircle.push_back(nv_percicle-1);e_smallcircle.push_back(0);


    onec.plane_para[0] = onec.plane_para[1] = onec.plane_para[3] = 0;
    onec.plane_para[2] = -1;

    auto &pv = onec.vertices;
    auto &pe = onec.edges;
    auto &pemat1 =  onec.m1;auto &pemat2 =  onec.m2;

    double biasAngle = 65;
    vector<double>anglelist({0,biasAngle,180-biasAngle,180,180+biasAngle,360-biasAngle});
    vector<double>radiulist({1.,0.8,0.8,1.,0.8,0.8});
    for(auto &a:anglelist)a=a/180*3.14;

    for(int iMat = 0;iMat<nmMat;++iMat){
        double shiftcenter[3] = {0,0,0};
        //shiftcenter[0] = r_bigcircle * iMat *4./3.;
        shiftcenter[1] = r_bigcircle * iMat * 2.2;


        for(int i=0;i<ncircle;++i){


            double center[3];
            center[0] = r_bigcircle* radiulist[i]* cos(anglelist[i]) + shiftcenter[0];
            center[1] = r_bigcircle* radiulist[i]* sin(anglelist[i]) + shiftcenter[1];
            center[2] = 0 + shiftcenter[2];
            uint offset = pv.size()/3;
            for(int j=0;j<nv_percicle;++j){
                auto p_vsmall = v_smallcircle.data()+j*3;
                for(int k=0;k<3;++k)pv.push_back(p_vsmall[k] + center[k]);
            }
            for(auto e:e_smallcircle)pe.push_back(e+offset);
            pemat1.insert(pemat1.end(),eMat_smallcircle[0].begin(),eMat_smallcircle[0].end());
            pemat2.insert(pemat2.end(),eMat_smallcircle[iMat+1].begin(),eMat_smallcircle[iMat+1].end());

        }
    }




    onec.CalSelfProperty();


    n_CrossSec = nlayers;
    vector<double>rota({0,0,1});
    vector<double>shiftbeforeR({0,0,0});

    vector<vector<double>>shiftv(nlayers,vector<double>(3));



    srand(120);
    vector<double>randomoscillation_height(nlayers,0);
    vector<double>randomoscillation_x(nlayers,0);
    vector<double>randomoscillation_y(nlayers,0);
    if(1)for(int i=0;i<nlayers;++i){
        randomoscillation_height[i] = ratio_coef*0.6*(randomdouble()-0.5);
        //randomoscillation_angle[i] = ratio_coef*10*(randomdouble()-0.5);

        randomoscillation_x[i] = ratio_coef*0.3*(randomdouble()-0.5);
        randomoscillation_y[i] = ratio_coef*0.3*(randomdouble()-0.5);
    }


    CSs.resize(n_CrossSec);
    for(int i=1;i<nlayers;++i)
    {
        CSs[i] = CSs[0];
        shiftv[i][2] = i*2.0*ratio_coef + randomoscillation_height[i];

        shiftv[i][0] = randomoscillation_x[i];shiftv[i][1] = randomoscillation_y[i];

        CSs[i].RotatingandShifting(0/180.*3.1415,shiftbeforeR,rota,shiftv[i]);
    }


    CalSelfProperty();


}

void CrossSections::SyntheticCircles_nonparallel(int ncircles, int nlayers){

    double r_smallcircle = 1.0;//0.6
    double r_bigcircle = 4.0;
    double ratio_coef = 1;

    double center_step = 3.5;

    int nv_percicle = 20;


    int nmMat = 1;
    double step_bigcircle = 2*my_PI / ncircles;
    double step_small_circle = 2*my_PI / nv_percicle;

    CSs.clear();


    vector<double>v_smallcircle;
    vector<double>v_smallcircle_vertical;
    vector<uint>e_smallcircle;
    vector<vector<int>>eMat_smallcircle(nmMat+1);
    for(int i=0;i<nmMat+1;++i){

        eMat_smallcircle[i].clear();
        eMat_smallcircle[i].resize(nv_percicle,i);
    }
    //vector<int>eMat1_smallcircle(nv_percicle,0),eMat2_smallcircle(nv_percicle,1),eMat3_smallcircle(nv_percicle,2),eMat4_smallcircle(nv_percicle,3);

    for(int i=0;i<nv_percicle;++i){
        v_smallcircle.push_back(r_smallcircle * cos(i*step_small_circle));
        v_smallcircle.push_back(r_smallcircle * sin(i*step_small_circle));
        v_smallcircle.push_back(0);
    }
    for(int i=0;i<nv_percicle;++i){
        v_smallcircle_vertical.push_back(0);
        v_smallcircle_vertical.push_back(r_smallcircle * cos(i*step_small_circle));
        v_smallcircle_vertical.push_back(r_smallcircle * sin(i*step_small_circle));

    }


    for(int i=0;i<nv_percicle-1;++i){
        e_smallcircle.push_back(i);e_smallcircle.push_back(i+1);
    }
    e_smallcircle.push_back(nv_percicle-1);e_smallcircle.push_back(0);






    srand(120);
    vector<double>randomoscillation_height(nlayers,0);
    vector<double>randomoscillation_x(nlayers,0);
    vector<double>randomoscillation_y(nlayers,0);
    for(int i=0;i<nlayers;++i){
        randomoscillation_height[i] = ratio_coef*0.6*(randomdouble()-0.5);
        randomoscillation_x[i] = ratio_coef*0.3*(randomdouble()-0.5);
        randomoscillation_y[i] = ratio_coef*0.3*(randomdouble()-0.5);
    }
    for(int iplane = 0;iplane<nlayers;++iplane){
        double h = iplane*center_step + 0.8*(randomdouble()-0.5);
        CSs.resize(CSs.size()+1);
        Contour& onec = CSs[CSs.size()-1];
        onec.plane_para[0] = 0;onec.plane_para[1] = 0;onec.plane_para[3] = -h;
        onec.plane_para[2] = -1;

        auto &pv = onec.vertices;
        auto &pe = onec.edges;
        auto &pemat1 =  onec.m1;auto &pemat2 =  onec.m2;



        for(int i=0;i<nlayers;++i){
            for(int j=0;j<nlayers;++j){

                double center[3] = {i*center_step,j*center_step,h};
                center[0] += 0.25*(randomdouble()-0.5);
                center[1] += 0.5*(randomdouble()-0.5);

                uint offset = pv.size()/3;
                for(int j=0;j<nv_percicle;++j){
                    auto p_vsmall = v_smallcircle.data()+j*3;
                    for(int k=0;k<3;++k)pv.push_back(p_vsmall[k] + center[k]);
                }
                for(auto e:e_smallcircle)pe.push_back(e+offset);
                pemat1.insert(pemat1.end(),eMat_smallcircle[0].begin(),eMat_smallcircle[0].end());
                pemat2.insert(pemat2.end(),eMat_smallcircle[1].begin(),eMat_smallcircle[1].end());


            }
        }

        onec.CalSelfProperty();
    }


    for(int iplane = 0;iplane<nlayers-1;++iplane){
        double h = (iplane+0.5)*center_step + 0.2*(randomdouble()-0.5);
        CSs.resize(CSs.size()+1);
        Contour& onec = CSs[CSs.size()-1];
        onec.plane_para[0] = -1;onec.plane_para[1] = 0;onec.plane_para[3] = -h;
        onec.plane_para[2] = 0;

        auto &pv = onec.vertices;
        auto &pe = onec.edges;
        auto &pemat1 =  onec.m1;auto &pemat2 =  onec.m2;



        for(int i=0;i<nlayers-1;++i){
            for(int j=0;j<nlayers-1;++j){

                double center[3] = {h,(i+0.5)*center_step,(j+0.5)*center_step};
                center[1] += 0.5*(randomdouble()-0.5);
                center[2] += 0.5*(randomdouble()-0.5);

                uint offset = pv.size()/3;
                for(int j=0;j<nv_percicle;++j){
                    auto p_vsmall = v_smallcircle_vertical.data()+j*3;
                    for(int k=0;k<3;++k)pv.push_back(p_vsmall[k] + center[k]);
                }
                for(auto e:e_smallcircle)pe.push_back(e+offset);
                pemat1.insert(pemat1.end(),eMat_smallcircle[0].begin(),eMat_smallcircle[0].end());
                pemat2.insert(pemat2.end(),eMat_smallcircle[1].begin(),eMat_smallcircle[1].end());


            }
        }

        onec.CalSelfProperty();
    }






    n_CrossSec = nlayers;


    CalSelfProperty();


}

void CrossSections::CalSelfProperty(){
    for(auto &a:CSs)a.CalSelfProperty();
    width = height = -1;

    n_CrossSec = CSs.size();
    for(int i=0;i<n_CrossSec;i++)width = max(width,CSs[i].maxx);
    for(int i=0;i<n_CrossSec;i++)height = max(height,CSs[i].maxy);
    width *= 2;height *= 2;
    total_nver = 0;
    total_nedges = 0;
    for(auto &a:CSs)total_nedges+=a.n_edges;
    acc_nver.resize(n_CrossSec+1,0);
    for(int i=0;i<n_CrossSec;i++)acc_nver[i+1] = CSs[i].n_vertices+acc_nver[i];
    total_nver = acc_nver[n_CrossSec];

    set<int>matts;
    for(auto &a:CSs)for(auto b:a.m1)matts.insert(b);
    for(auto &a:CSs)for(auto b:a.m2)matts.insert(b);
    n_materials = matts.size();


}
void CrossSections::Stackup(vector<double>&points, vector<uint> &edges){

    CalSelfProperty();

    points.clear();

    if(1){
        for(int i=0;i<n_CrossSec;i++)for(int j=0;j<CSs[i].n_vertices*3;++j)points.push_back( CSs[i].vertices[j]);
        edges.clear();
        for(int i=0;i<n_CrossSec;i++)for(int j=0;j<CSs[i].n_edges*2;++j){
            edges.push_back( CSs[i].edges[j]+acc_nver[i]);
        }
    }else{
        vector<int>pickCs({0,2,3,5,6});
        for(auto i:pickCs)for(int j=0;j<CSs[i].n_vertices*3;++j)points.push_back( CSs[i].vertices[j]);
        edges.clear();
        for(auto i:pickCs)for(int j=0;j<CSs[i].n_edges*2;++j){
            edges.push_back( CSs[i].edges[j]+acc_nver[i]);
        }
    }

}

void ReadCrossSection(string filename){





}

void CrossSections::MapMaterials(){

    map<int,int>mapping;
    mapping[11]=0;
    mapping[15] = 2;
    mapping[17] = 1;

    for(auto &a:CSs){
        for(auto &b:a.m1)b = mapping[b];
        for(auto &b:a.m2)b = mapping[b];
    }


}
int CrossSections::WriteCrossSections(string filename){
    string outCname = filename + string("mapping.contour");
    ofstream outer(outCname.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }

    outer<<CSs.size()<<endl;
    for(auto &b:CSs)b.WriteContour(outer);
    outer.close();
    return 0;
    for(int i=0;i<CSs.size();++i){
        string outCname = filename + string("mapping")+to_string(i)+ string(".contour");
        ofstream outer(outCname.data(), ofstream::out);
        if (!outer.good()) {
            cout << "Can not open the Contour file " << outCname << endl;
            return -1;
        }

        outer<<1<<endl;
        CSs[i].WriteContour(outer);
        outer.close();
    }

    return 0;

}


int CrossSections::WriteCrossSectionsToYaron(string filepath){

    Curve AnalyzeCurve;

    CalSelfProperty();

    vector<double>points;
    vector<uint>edges;
    vector<uint>edgeMat;
    vector<int>edgeCells;
    vector<uint>acc_nver_(n_CrossSec+1,0);
    for(int i=0;i<n_CrossSec;i++)acc_nver_[i+1] = CSs[i].n_vertices+acc_nver_[i];
    points.clear();

    for(int i=0;i<n_CrossSec;i++)for(auto a:CSs[i].vertices)points.push_back(a);
    edges.clear();
    for(int i=0;i<n_CrossSec;i++)for(auto a:CSs[i].edges)edges.push_back(a+acc_nver_[i]);

    for(int i=0;i<n_CrossSec;i++)for(int j=0;j<CSs[i].n_edges;++j){
        edgeCells.push_back(i+1);edgeCells.push_back(0);
        edgeMat.push_back(CSs[i].m1[j]);edgeMat.push_back(CSs[i].m2[j]);
    }
    AnalyzeCurve.ImportCurve(points,edges);
    AnalyzeCurve.CurveNetAnalayze(edgeMat,edgeCells,n_CrossSec+1);


    string filename = filepath + "ctr";
    writeObjFile_line(filename,points,edges);

    auto& seg2edges = AnalyzeCurve.seg2edges;

    filename = filepath + "loop2edges";
    writeVVecFile(filename,seg2edges);


    vector<vector<uint>>plane2loops(n_CrossSec);
    for(int planei=0;planei<n_CrossSec;++planei)plane2loops[planei] = AnalyzeCurve.cell2Seg[planei+1];
    filename = filepath + "plane2loops";
    writeVVecFile(filename,plane2loops);





}

void Contour::WriteContourToFCM(ofstream &out){


    vector<vector<double>>verticesSeqs;
    vector<int>loops2mats;

    ComputeMat2VerticesLoop(verticesSeqs,loops2mats);

    //    out<<"{0,0,0}";
    //    return;
    // {{v, e},{v, e}}
    //loop1 = {loop1, eloop1, 1};
    //{loop1, loop2, loop3}
    //
    out<<"{";
    for(int i=0;i<loops2mats.size();++i){
        out<<"{";
        //v
        //out<<"{{0,0,0}},";
        out<<"{";
        auto p_v = verticesSeqs[i].data();
        int nv = verticesSeqs[i].size()/3;
        int j=0;
        for(j=0;j<nv-1;++j){
            out<<"{";
            for(int k=0;k<2;++k)out<<p_v[j*3+k]<<", ";
            out<<p_v[j*3+2];
            out<<"}, ";
        }
        out<<"{";
        for(int k=0;k<2;++k)out<<p_v[j*3+k]<<", ";
        out<<p_v[j*3+2];
        out<<"}";
        out<<"}, ";

        //e
        //out<<"{{0,1}},";

        out<<"{";
        for(int k=1;k<nv;++k){
            out<<"{";
            out<<k<<", "<<k+1;
            out<<"}, ";
        }
        out<<"{";out<<nv<<", "<<0;out<<"}";
        out<<"},";

        out<<loops2mats[i]+1;

        if(i!=loops2mats.size()-1)out<<"}, ";
        else out<<"}";
    }
    out<<"}";


}

int CrossSections::WriteCrossSectionsToFCM(string filename,int nMat){

    string outCname = filename;
    ofstream out(outCname.data(), ofstream::out);
    if (!out.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }
    cout << "write:  " << outCname << endl;

    out<<std::setprecision(12);
    out<<std::fixed;
    out<<"{";
    //loops
    out<<"{";
    for(int i=0;i<CSs.size();++i){
        CSs[i].WriteContourToFCM(out);
        if(i!=CSs.size()-1)out<<",";

    }

    //out<<"},";
    out<<"},";

    //plane paras
    out<<"{";
    for(int i=0;i<CSs.size();++i){
        out<<"{";
        out<<"{0,0,0},";
        //          auto p_cv = CSs[i].v_begin(0);
        //          out<<"{";
        //          out<<p_cv[0]<<","<<p_cv[1]<<","<<p_cv[2];
        //          out<<"},";
        out<<"{";
        for(int j=0;j<3;++j)out<<CSs[i].plane_para[j]<<",";
        out<<CSs[i].plane_para[3];
        out<<"},";
        vector<double>tmp(3,1);
        for(int j=0;j<3;++j)if(fabs(CSs[i].plane_para[j])>1e-1){
            int v1,v2;
            if(j==0){v1=1;v2=2;}else if(j==1){v1=0;v2=2;}else if(j==2){v1=0;v2=1;}
            tmp[j] = (0-CSs[i].plane_para[v1]*tmp[v1]-CSs[i].plane_para[v2]*tmp[v2])/CSs[i].plane_para[j];
            break;
        }
        normalize(tmp.data());
        out<<"{"<<tmp[0]<<","<<tmp[1]<<","<<tmp[2]<<"}";

        if(i!=CSs.size()-1) out<<"},";
        else out<<"}";
    }
    out<<"},";

    out<<nMat<<"}";
    out.close();


    cout<<"write fin"<<endl;
    return 0;

}

void Contour::WriteContourToFCMtCtr(ofstream &out){





    int nV = vertices.size()/3;
    out<<"{";
    out<<"{";

    //v
    out<<"{";
    for(int i=0;i<nV;++i){
        auto p_v = v_begin(i);
        out<<"{";
        out<<p_v[0]<<","<<p_v[1]<<","<<p_v[2];
        if(i!=nV-1)out<<"},";else out<<"}";
    }
    out<<"},";

    //e
    int ne = edges.size()/2;
    out<<"{";
    for(int i=0;i<ne;++i){
        auto p_e = ev_begin(i);
        out<<"{";
        out<<p_e[0]+1<<","<<p_e[1]+1;
        if(i!=ne-1)out<<"},";else out<<"}";
    }
    out<<"}";


    out<<"}";
    out<<"}";

}

int CrossSections::WriteCrossSectionsToFCMtCtr(string filename){

    string outCname = filename;
    ofstream out(outCname.data(), ofstream::out);
    if (!out.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }
    cout << "write:  " << outCname << endl;
    out<<std::setprecision(12);
    out<<std::fixed;
    out<<"{";
    //loops
    out<<"{";
    for(int i=0;i<CSs.size();++i){
        CSs[i].WriteContourToFCMtCtr(out);
        if(i!=CSs.size()-1)out<<",";

    }

    //out<<"},";
    out<<"},";

    //plane paras
    out<<"{";
    for(int i=0;i<CSs.size();++i){
        out<<"{";
        //          auto p_cv = CSs[i].v_begin(0);
        //          out<<"{";
        //          out<<p_cv[0]<<","<<p_cv[1]<<","<<p_cv[2];
        //          out<<"},";
        out<<"{0,0,0},";

        out<<"{";
        for(int j=0;j<3;++j)out<<CSs[i].plane_para[j]<<",";
        out<<CSs[i].plane_para[3];
        out<<"},";
        vector<double>tmp(3,1);
        for(int j=0;j<3;++j)if(fabs(CSs[i].plane_para[j])>1e-1){
            int v1,v2;
            if(j==0){v1=1;v2=2;}else if(j==1){v1=0;v2=2;}else if(j==2){v1=0;v2=1;}
            tmp[j] = (0-CSs[i].plane_para[v1]*tmp[v1]-CSs[i].plane_para[v2]*tmp[v2])/CSs[i].plane_para[j];
            break;
        }
        normalize(tmp.data());
        out<<"{"<<tmp[0]<<","<<tmp[1]<<","<<tmp[2]<<"}";

        if(i!=CSs.size()-1) out<<"},";
        else out<<"}";
    }
    out<<"}";

    out<<"}";

    out.close();


    cout<<"write fin"<<endl;
    return 0;

}

void CrossSections::WriteCrossSectionToObj(string filename){


    ofstream out(filename.data(), ofstream::out);
    if (!out.good()) {
        cout << "Can not open the Contour file " << filename << endl;
        return;
    }
    cout << "write:  " << filename << endl;

    vector<double>Vs;
    vector<uint>Es;
    Stackup(Vs,Es);


    int nv = Vs.size()/3;
    int ne = Es.size()/2;


    out<<std::fixed;
    for(int i=0;i<nv;++i){
        out<<"v ";
        int ind = i*3;
        out<<Vs[ind]<<" "<<Vs[ind+1]<<" "<<Vs[ind+2]<<endl;
    }
    for(int i=0;i<ne;++i){
        out<<"l ";
        int ind = i*2;
        out<<Es[ind]+1<<" "<<Es[ind+1]+1<<endl;
    }

    out.close();



}

void CrossSections::naiveConstruction(){

    cout<<"naiveConstruction"<<endl;
    //--------------------- Construct the frame ---------------------
    // <1> frameVs.
    // Two addtional layers are added on the top and bottom. The spaceing of each
    // addtional layer to it adjacent layer is set to be the minimum spacing among
    // all the spacings.
    // |spacings| = _nCrossSection-1; |extSpacings| = _nCrossSection+2 = nLayer.
    vector<float> spacings(n_CrossSec-1, 100.0f);

    vector<vector<double> > _frameVs;
    vector<vector<int> > _frameFs;
    vector<vector<int>> _frameCells;
    vector<vector<double> > _frameF2PlaneParas;


    CalSelfProperty();


    float bdSpacing = *min_element(spacings.begin(), spacings.end());
    vector<float> extSpacings = spacings;
    // For the second layer (first true layer).
    extSpacings.insert(extSpacings.begin(), bdSpacing);
    // For the first layer.
    extSpacings.insert(extSpacings.begin(), 0.0);
    // For the last layer.
    extSpacings.push_back(bdSpacing);
    int nLayer = n_CrossSec + 2;
    int _nCell = n_CrossSec + 1;
    for (int i = 1; i < extSpacings.size(); ++i) {
        extSpacings[i] += extSpacings[i - 1];
    }
    _frameVs.reserve(nLayer * 4);
    for (int i = 0; i < nLayer; ++i) {
        _frameVs.push_back({0, 0, extSpacings[i]});
        _frameVs.push_back({1.0 * width, 0., extSpacings[i]});
        _frameVs.push_back({1.0 * width, 1.0 * height, extSpacings[i]});
        _frameVs.push_back({0, 1.0 * height, extSpacings[i]});
    }
    // <2> frameFs.
    // frameFs = shared Fs and side Fs.
    // Shared Fs.
    for(auto a:extSpacings)cout<<a<<' ';cout<<endl;
    for (int i = 0; i < nLayer; ++i) {
        _frameFs.push_back({i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 3});
    }
    // Side Fs.
    for (int i = 0; i < _nCell; ++i) {
        int pre = i * 4;
        _frameFs.push_back({pre + 0, pre + 1, pre + 5, pre + 4});
        _frameFs.push_back({pre + 1, pre + 2, pre + 6, pre + 5});
        _frameFs.push_back({pre + 2, pre + 3, pre + 7, pre + 6});
        _frameFs.push_back({pre + 3, pre + 0, pre + 4, pre + 7});
    }
    // <3> frameCells.
    for (int i = 0; i < _nCell; ++i) {
        int pre = i * 4 + nLayer;
        _frameCells.push_back({i, i + 1, pre, pre + 1, pre + 2, pre + 3});
    }
    // <4> frameF2PlaneParas.
    // Four parameters used to define each plane {a,b,c,d} in the form of
    // ax+by+cz+d=0.
    for (int i = 0; i < nLayer; ++i) {
        _frameF2PlaneParas.push_back({0, 0, 1, -1.0f * extSpacings[i]});
    }

    //return;

    // ----------------------- Generate tets ------------------------
    // Find the centroid of each cell for classifying tets in to cells.
    vector<vector<double>> cellCenter(_nCell, vector<double>(3));
    for (int i = 0; i < _nCell; ++i) {
        auto& res = cellCenter[i];
        for (int p = i * 4; p < i * 4 + 8; ++p) {
            std::transform(res.begin(), res.end(), _frameVs[p].begin(), res.begin(),
                           std::plus<double>());
        }
        for(auto &a:res)a/=8;
        //Utility::MultiplyVectorsByValue(0.125f, res);
    }
    // Prepare the data structure for tetgen.
    tetgenio in, out;
    tetgenio::facet* f;
    tetgenio::polygon* p;
    // PLC: pointlist.
    cout<<"PLC: pointlist. "<<_frameVs.size()<<' '<<total_nver<<endl;
    in.numberofpoints = _frameVs.size()+total_nver;
    //in.numberofpoints = _frameVs.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    in.pointmarkerlist = new int[in.numberofpoints];
    for (int i = 0; i < _frameVs.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            in.pointlist[i * 3 + j] = _frameVs[i][j];
        }
        in.pointmarkerlist[i] = -1;
    }
    for(int i=0;i<n_CrossSec;++i){
        int offset = (_frameVs.size()+acc_nver[i])*3;
        for(int j=0;j<CSs[i].vertices.size();++j)in.pointlist[offset+j] = CSs[i].vertices[j];
        offset = offset/3;
        for(int j=0;j<CSs[i].n_vertices;++j)in.pointmarkerlist[offset+j] = -10-i;
    }


    // PLC: facetlist.
    cout<<"PLC: facetlist."<<endl;
    in.numberoffacets = _frameFs.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (int i = 0; i < _frameFs.size(); ++i) {
        f = &in.facetlist[i];
        f->numberofpolygons = 1;
        if(i<1 || i>n_CrossSec)f->numberofpolygons = 1;
        else f->numberofpolygons=1+CSs[i-1].n_edges;

        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = _frameFs[i].size();
        p->vertexlist = new int[p->numberofvertices];
        for (int j = 0; j < _frameFs[i].size(); ++j) {
            p->vertexlist[j] = _frameFs[i][j];
        }
        in.facetmarkerlist[i] = i;

        if(1)if(!(i<1 || i>n_CrossSec)){
            int offset = _frameVs.size()+acc_nver[i-1];
            for(int j=0;j<CSs[i-1].n_edges;++j){
                p = &f->polygonlist[j+1];
                p->numberofvertices = 2;p->vertexlist = new int[p->numberofvertices];
                auto p_ev = CSs[i-1].ev_begin(j);
                for (int k = 0; k < 2; ++k)p->vertexlist[k] = p_ev[k]+ offset;
            }
        }
    }

    // PLC: region.
    cout<<"PLC: region."<<endl;
    in.numberofregions = _frameCells.size();
    in.regionlist = new REAL[in.numberofregions * 5];
    for (int i = 0; i < _frameCells.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            in.regionlist[i * 5 + j] = cellCenter[i][j];
        }
        // Region attribute (marker).
        // Starting from -1. {-1, -2, -3, ...}.
        in.regionlist[i * 5 + 3] = -i - 1;
        // Region volume constraint (not used if specify overall constraint after -a
        // switch).
        in.regionlist[i * 5 + 4] = 100000;
    }


    char tetin[] = "../tetgen1.4.3/c++tetin";
    // in.save_nodes(tetin);
    // in.save_poly(tetin);

    // Tetrahedralize the PLC.
    // Switches are chosen to read PLC (p), do quality mesh generation (q) with a
    // specified quality bound (1.414), apply a maximum volume constraint
    // (a10000) and omit terminal output except errors (Q).
    // -Qpq1.4a10000
    float tetVolLimit = 5000.0f;
    string switchesStr = "Qpq1.414Aa" + to_string(tetVolLimit);
    // char switches[switchesS.size()] = switchesS.c_str();
    // char switches[] = "Qpq1.414Aa10000";
    std::vector<char> tmp(switchesStr.begin(), switchesStr.end());
    tmp.push_back('\0');
    char* switches = &tmp[0];
    // char switches[] = "Qpq1.414Aa10";
    tetrahedralize(switches, &in, &out);


    char tetout[] = "/Users/Research/Geometry/MM/ConvertFolder/tetout";
    out.save_nodes(tetout);
    out.save_elements(tetout);
    out.save_faces(tetout);



}







/**********************************************************************/



void MultiCellTopo::BuildDisplay(infoSurfDisp info,vector<int>&topo,int pickMat){
    //mixSurf.BuildDisplay(info);
    //for(auto &a:surfbuffer)cout<<a.weighted_fcolor.size()<<endl;
    //for(auto &a:CellTopo)for(auto&b : a)b.BuildDisplay(info,false);

    for(int i=0;i<CellTopo.size();++i)if(nTopo[i]>0)for(auto&b : CellTopo[i])b.BuildDisplay(info,false);



    display_vertices.resize(CellTopo.size()*2+1);
    display_edges.resize(CellTopo.size()*2+1);
    display_faces.resize(CellTopo.size()*2+1);
    display_normal.resize(CellTopo.size()*2+1);
    display_vcolor.resize(CellTopo.size()*2+1);

    int pickMatInd = pickMat+1;
    if(pickMat==-1)pickMatInd = nMat + 2;
    //cout<<pickMat<<endl;
    //cout<<pickMatInd<<endl;
    //pickMatInd = 1;
    for(int i=0;i<CellTopo.size();++i)if(nTopo[i]>0){
        //if(i!=1)continue;
        int pickTopo = topo[i];
        display_vertices[i+1] = (CellTopo[i][pickTopo].getDisplayVertices())->at(pickMatInd);
        display_edges[i+1] = (CellTopo[i][pickTopo].getDisplayEdges())->at(pickMatInd);
        display_faces[i+1] = (CellTopo[i][pickTopo].getDisplayFaces())->at(pickMatInd);
        display_normal[i+1] =( CellTopo[i][pickTopo].getDisplayVerticesNormal())->at(pickMatInd);
        display_vcolor[i+1] = (CellTopo[i][pickTopo].getDisplayColor())->at(pickMatInd);
        //cout<<"iii: "<<i<<' '<<(CellTopo[i][pickTopo].getDisplayVertices())->size()<<endl;
    }

    for(int i=0;i<CellTopo.size();++i)if(nTopo[i]>0){
        //if(i!=1)continue;
        int pickTopo = topo[i];
        int nnind = (CellTopo[i][pickTopo].getDisplayVertices())->size()-2;
        display_vertices[i+1+CellTopo.size()] = (CellTopo[i][pickTopo].getDisplayVertices())->at(nnind);
        display_edges[i+1+CellTopo.size()] = (CellTopo[i][pickTopo].getDisplayEdges())->at(nnind);
        display_faces[i+1+CellTopo.size()] = (CellTopo[i][pickTopo].getDisplayFaces())->at(nnind);
        display_normal[i+1+CellTopo.size()] =( CellTopo[i][pickTopo].getDisplayVerticesNormal())->at(nnind);
        display_vcolor[i+1+CellTopo.size()] = (CellTopo[i][pickTopo].getDisplayColor())->at(nnind);
        //cout<<"iii: "<<i<<' '<<(CellTopo[i][pickTopo].getDisplayVertices())->size()<<endl;
    }




    //    display_vertices[0] = mixSurf.getDisplayVertices();
    //    display_edges[0] = mixSurf.getDisplayEdges();
    //    display_faces[0] = mixSurf.getDisplayFaces();
    //    display_normal[0] = mixSurf.getDisplayVerticesNormal();
    //    display_vcolor[0] = mixSurf.getDisplayColor();

    //    display_vertices[0] = CellTopo[0][0][0];
    //    display_edges[0] = crossSection.getDisplayEdges();
    //    display_faces[0] = crossSection.getDisplayFaces();
    //    display_normal[0] = crossSection.getDisplayVerticesNormal();
    //    display_vcolor[0] = crossSection.getDisplayColor();

    crossS.BuildDisplay(false,true);
    //    display_vertices[0] = (CellTopo[0][0].getDisplayVertices())->at(0);
    //    display_edges[0] = (CellTopo[0][0].getDisplayEdges())->at(0);
    //    display_faces[0] = (CellTopo[0][0].getDisplayFaces())->at(0);
    //    display_normal[0] =( CellTopo[0][0].getDisplayVerticesNormal())->at(0);
    //    display_vcolor[0] = (CellTopo[0][0].getDisplayColor())->at(0);

    display_vertices[0] = crossS.getDisplayVertices();
    display_edges[0] = crossS.getDisplayEdges();
    display_faces[0] = crossS.getDisplayFaces();
    display_normal[0] = crossS.getDisplayVerticesNormal();
    display_vcolor[0] = crossS.getDisplayColor();



}
void MultiCellTopo::ReadCellTopo(string nameprefix, vector<int>n_Topos, bool isSmooth){

    CellTopo.clear();
    CellTopo.resize(n_Topos.size());
    string ed(".suf");

    Mesh tmp;
    MultiSurface a;

    for(int i=0;i<n_Topos.size();++i)if(n_Topos[i]>0){

        string s_cell = string("_")+to_string(i);
        string s_topo = string("_")+to_string(0);
        nameprefix+s_cell+s_topo+ed;

        a.readSufFile(nameprefix+s_cell+s_topo+ed);
        tmp.addMesh(a.mixSurf);
        nMat = max(nMat,a.GetNumofSurfaces());
    }

    tmp.GetRescaleInfo(lscale, pcenter);


    //cout<<lscale<<endl;
    //cout<<pcenter[0]<<' '<<pcenter[1]<<' '<<pcenter[2]<<endl;



    nTopo = n_Topos;
    for(int i=0;i<n_Topos.size();++i){
        CellTopo[i].resize(n_Topos[i]);
        string s_cell = string("_")+to_string(i);
        for(int j=0;j<n_Topos[i];++j){
            string s_topo = string("_")+to_string(j);
            CellTopo[i][j].Initialize(nameprefix+s_cell+s_topo+ed,false,false,true,nMat,true,lscale,pcenter);
        }
    }


    for(auto &a:CellTopo)if(a.size()!=0){
        a[0].GetColorDegree(mat_colordegree);
    }
    //CellTopo[0][0].GetColorDegree(mat_colordegree);

    crossS.reset();


    for(int i=0;i<n_Topos.size();++i)if(n_Topos[i]>0){

        crossS.AddCurve(CellTopo[i][0].crossSection);
    }
    crossS.setparameters();
    crossS.BuildEdges(false);


    if(!isSmooth)return;

    GlobalSmoothing();


    //for(auto a:TopoCurInd)cout<<a<<' ';cout<<endl;
    //for(auto &b:pSuf.mappingsToGlobal){for(auto a:b)cout<<a<<' ';cout<<endl;}


}


void MultiCellTopo::ReadCellTopo_picked(string nameprefix, vector<int>&n_Topos, vector<int>&picked_Topos, bool isSmooth){

    CellTopo.clear();
    CellTopo.resize(n_Topos.size());
    string ed(".suf");

    Mesh tmp;
    MultiSurface a;

    for(int i=0;i<n_Topos.size();++i)if(n_Topos[i]>0){

        string s_cell = string("_")+to_string(i);
        string s_topo = string("_")+to_string(0);
        nameprefix+s_cell+s_topo+ed;

        a.readSufFile(nameprefix+s_cell+s_topo+ed);
        tmp.addMesh(a.mixSurf);
        nMat = max(nMat,a.GetNumofSurfaces());
    }

    tmp.GetRescaleInfo(lscale, pcenter);


    //cout<<lscale<<endl;
    //cout<<pcenter[0]<<' '<<pcenter[1]<<' '<<pcenter[2]<<endl;



    nTopo = n_Topos;
    for(int i=0;i<n_Topos.size();++i){
        CellTopo[i].resize(n_Topos[i]);
        string s_cell = string("_")+to_string(i);

        assert(picked_Topos[i]<n_Topos[i]);

        string s_topo = string("_")+to_string(picked_Topos[i]);
        CellTopo[i][picked_Topos[i]].Initialize(nameprefix+s_cell+s_topo+ed,false,false,true,nMat,true,lscale,pcenter);
    }


    for(int i=0;i<n_Topos.size();++i)if(CellTopo[i].size()!=0){
        CellTopo[i][picked_Topos[i]].GetColorDegree(mat_colordegree);
    }
    //CellTopo[0][0].GetColorDegree(mat_colordegree);

    crossS.reset();


    for(int i=0;i<n_Topos.size();++i)if(n_Topos[i]>0){

        crossS.AddCurve(CellTopo[i][picked_Topos[i]].crossSection);
    }
    crossS.setparameters();
    crossS.BuildEdges(false);


    if(!isSmooth)return;

    GlobalSmoothing();


    //for(auto a:TopoCurInd)cout<<a<<' ';cout<<endl;
    //for(auto &b:pSuf.mappingsToGlobal){for(auto a:b)cout<<a<<' ';cout<<endl;}


}

void MultiCellTopo::WriteAllSurface(string nameprefix){


    for(int i=0;i<nTopo.size();++i){

        string s_cell = string("_")+to_string(i);
        for(int j=0;j<nTopo[i];++j){
            string s_topo = string("_")+to_string(j);
            CellTopo[i][j].WriteSuf(nameprefix+s_cell+s_topo);
        }
    }

}

void MultiCellTopo::ReReadSingleCellTopo(string nameprefix, int celli, int c_nTopo){


    string ed(".suf");

    nTopo[celli] = c_nTopo;
    CellTopo[celli].clear();
    CellTopo[celli].resize(c_nTopo);
    string s_cell = string("_")+to_string(celli);
    for(int j=0;j<c_nTopo;++j){
        string s_topo = string("_")+to_string(j);
        CellTopo[celli][j].Initialize(nameprefix+s_cell+s_topo+ed,false,false,true,nMat,true,lscale,pcenter);
    }

    SufStructure pSuf;
    MultiSurface multiSuf;
    vector<vector<double>>TopoVs; vector<vector<uint>>TopoFs;
    vector<vector<int> > TopoFMs;vector<vector<uint>>TopoCtrs;
    vector<int>TopoCurInd(nTopo.size(),-1);
    for(int j=0;j<nTopo.size();++j)if(TopoCurInd[j]+1<nTopo[j]){TopoCurInd[j]+=1;}
    vector<bool>hasChange(nTopo.size(),false);hasChange[celli] = true;
    for(int i=0;i<c_nTopo;++i){
        TopoCurInd[celli]=i;
        GetAllCellTopos(TopoCurInd, TopoVs, TopoFs,TopoFMs,TopoCtrs);
        pArr.CreateTotalSurf(TopoVs, TopoFs,TopoFMs,TopoCtrs,false,&pSuf);
        multiSuf.ImportFromSuf(pSuf,true);
        UpdateVerticesPosition(multiSuf.mixSurf.vertices,TopoCurInd,pSuf.mappingsToGlobal,hasChange);
    }


}

void MultiCellTopo::UpdateVerticesPosition(vector<double>&newVpos,vector<int>&pickTopos,vector<vector<int>>&mappingToGlobal,vector<bool>&isChanged){

    for(int i=0;i<pickTopos.size();++i)if(pickTopos[i]>=0)if(isChanged[i]){

        vector<double>&toVs = CellTopo[i][pickTopos[i]].mixSurf.vertices;
        auto &mappingsTog = mappingToGlobal[i];
        auto p_newVd = newVpos.data();
        //toVs.clear();toVs.resize(mappingsTog.size()*3,0);
        auto p_vd = toVs.data();
        for(int j=0;j<mappingsTog.size();++j){

            copyVec(p_newVd+mappingsTog[j]*3,p_vd+j*3);
        }
        CellTopo[i][pickTopos[i]].ComputeSmallMVnormal();
        CellTopo[i][pickTopos[i]].ReAllocation(false);
        //cout<<"CACACA"<<endl;
    }


}
void MultiCellTopo::GlobalSmoothing(){



    SufStructure pSuf;
    MultiSurface multiSuf;

    vector<vector<double>>TopoVs; vector<vector<uint>>TopoFs;
    vector<vector<int> > TopoFMs;vector<vector<uint>>TopoCtrs;

    vector<int>TopoCurInd(nTopo.size(),-1);

    int maxNCellTopo = *max_element(nTopo.begin(),nTopo.end());
    for(int i=0;i<maxNCellTopo;++i){
        vector<bool>hasChange(nTopo.size(),false);
        for(int j=0;j<nTopo.size();++j)if(TopoCurInd[j]+1<nTopo[j]){TopoCurInd[j]+=1;hasChange[j] = true;}
        GetAllCellTopos(TopoCurInd, TopoVs, TopoFs,TopoFMs,TopoCtrs);
        pArr.CreateTotalSurf(TopoVs, TopoFs,TopoFMs,TopoCtrs,false,&pSuf);
        multiSuf.ImportFromSuf(pSuf,true);
        UpdateVerticesPosition(multiSuf.mixSurf.vertices,TopoCurInd,pSuf.mappingsToGlobal,hasChange);
    }

}

void MultiCellTopo::GetAllCellTopos(vector<int>&pickTopos, vector<vector<double>>&TopoVs, vector<vector<uint>>&TopoFs, vector<vector<int> > &TopoFMs,vector<vector<uint>>&TopoCtrs){

    if(pickTopos.size()!=CellTopo.size()){
        cout<<"please input a complete choice of all cells"<<endl;
        return;
    }
    TopoVs.clear();
    TopoFs.clear();
    TopoCtrs.clear();
    TopoFMs.clear();


    TopoVs.resize(CellTopo.size());
    TopoFs.resize(CellTopo.size());
    TopoFMs.resize(CellTopo.size());
    TopoCtrs.resize(CellTopo.size());

    for(int i=0;i<CellTopo.size();i++)if(pickTopos[i]>=0){

        TopoVs[i] = CellTopo[i][pickTopos[i]].mixSurf.vertices;
        TopoFs[i] = CellTopo[i][pickTopos[i]].mixSurf.faces2vertices;
        int nf = CellTopo[i][pickTopos[i]].mixSurf.n_faces;
        auto &sm1 =  CellTopo[i][pickTopos[i]].surf_group1;
        auto &sm2 =  CellTopo[i][pickTopos[i]].surf_group2;
        for(int j=0;j<nf;++j){
            TopoFMs[i].push_back(sm1[j]);
            TopoFMs[i].push_back(sm2[j]);
        }

        TopoCtrs[i] = CellTopo[i][pickTopos[i]].curvee2v;


    }






}


void MultiCellTopo::CutSurfaceByPlane(int celli,int Topoi,double *para,vector<double>&outCtrV,vector<uint>&outCtrE,vector<int>&outCtrEMat){

    CellTopo[celli][Topoi].mixSurf.CutSurfaceByPlane(para,outCtrV,outCtrE,outCtrEMat);

}

void MultiCellTopo::GetRescaleInfo(double &outlscale, double *outpcenter){

    outlscale = lscale;
    for(int i=0;i<3;++i)outpcenter[i] = pcenter[i];


}

void MultiCellTopo::GetColorDegree(vector<double>&out_colordegree){


    out_colordegree = mat_colordegree;

}
double MultiCellTopo::GetLabel2Colordegree(int label){

    return mat_colordegree[label];
}
}//n_rf
