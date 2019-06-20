#include "geo_sur.h"
#include <stdio.h>
#include<iostream>
#include<fstream>
#include<limits>
#include<cstdlib>
#include <functional>
#include<time.h>
#include<list>
#include <sstream>
#include <eigen3/Eigen/Geometry>
#include "readers.h"
#include<set>
#include "gurobi_c++.h"

#include"../MMCore/UnionFind.h"
namespace n_rf {
bool Surface::isload = false;
Mesh Surface::sphere;

double randf() {
    return ((double)(rand()/(double)RAND_MAX)-0.5)*2.;
}
Surface::Surface():Mesh(){
    isbuilddisp = false;

    isresetColor = true;
    reinitflag = true;

    isFaceRenderMode = false;
    ismarkC = false;

    load();

}

void Surface::load(){
    if(!isload){
        //sphere.createToy(3);
        //sphere.BuildNeighborTable();
        //sphere.ComputeFaceNormal(true);
        isload = true;
    }
}


void Surface::BuildFacesCenter(){
    if(n_vertices==0)return;
    if(n_faces==0)return;

    faces_center.resize(n_faces*3);
    for(uint i = 0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        _TriangleMidpoint(v_begin(p_fv[0]),v_begin(p_fv[1]),v_begin(p_fv[2]),fcent_begin(i));
    }
}
void Surface::ComputeArea(){
    if(n_vertices==0)return;
    if(n_faces==0)return;
    vertices_area.resize(n_vertices);
    faces_area.resize(n_faces);
    faces_inverse_area.resize(n_faces);
    const double scale_up = 1.;
    for(uint i = 0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        faces_area[i]=TriangleArea(p_fv[0],p_fv[1],p_fv[2]);
    }
    for(int i=0;i<n_vertices;++i){
        double b = 0;
        for(auto p_vf = vf_begin(i);p_vf != vf_end(i);++p_vf)b+=faces_area[*p_vf];
        //cout<<b<<' ';cout<<endl;
    }
    double scaleup = 1e-5;
    for(uint i = 0;i<n_faces;++i){
        faces_inverse_area[i]=scaleup/faces_area[i];
    }


    for(uint i = 0;i<n_vertices;++i){
        double a=0;
        for(auto p_vf = vf_begin(i);p_vf != vf_end(i);++p_vf){
            a+=faces_area[*p_vf];
        }
        vertices_area[i]=a*scale_up/3;
    }


}
void Surface::ComputeEdgeLength(){
    if(n_vertices==0)return;
    if(n_faces==0)return;
    if(n_edges==0)return;

    edge_len.resize(n_edges);
    edge_vec.resize(n_edges*3);
    for(uint i = 0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        vectorize(p_ev[0],p_ev[1],evec_begin(i));
        edge_len[i]=normVec(evec_begin(i));

    }

    edge_cot.resize(n_edges);
    for(uint i = 0;i<n_edges;++i){
        edge_cot[i] = EdgeCottan(i);
    }

    double ave_edge_lentmp = 0.;
    ave_edge_lentmp = 0.;
    for(auto a:edge_len)ave_edge_lentmp+=a;
    ave_edge_lentmp/=n_edges;
    ave_edge_lentmp*=2.;

    ave_edge_len = 0.;
    double numofE = 0;
    for(auto a:edge_len)if(a<ave_edge_lentmp){ave_edge_len+=a;numofE++;}
    ave_edge_len/=numofE;
    //cout<<"numofE: "<<(numofE==n_edges)<<' '<<ave_edge_len<<' '<<ave_edge_lentmp<<endl;



}



//bool Surface::readObjxfile(string filename){
//    ifstream reader(filename.data(), ifstream::in);
//    if (!reader.good()) {
//        cout << "Can not open the Objx file " << filename << endl;
//        return false;
//    }
//    clearup();
//    auto readVertices = [this](stringstream &ss){
//        double dvalue;
//        for(int i=0;i<3;++i){ss>>dvalue;vertices.push_back(dvalue);}
//    };
//    auto readFace = [this](stringstream &ss){
//        int ivalue;
//        for(int i=0;i<3;++i){
//            ss>>ivalue;faces2vertices.push_back(ivalue-1);
//        }
//    };
//    auto readVerticesField = [this](stringstream &ss){
//        double dvalue;
//        for(int i=0;i<3;++i){ss>>dvalue;vertices_field.push_back(dvalue);}
//    };


//    string oneline;

//    cout<<"reading: "<<filename<<endl;

//    while( getline( reader, oneline ))
//    {
//        stringstream ss( oneline );
//        string token;

//        ss >> token;

//        if( token == "v"  ) { readVertices( ss ); continue; } // vertex
//        if( token == "vt" ) {  continue; } // texture coordinate
//        if( token == "vn" ) {  continue; } // vertex normal
//        if( token == "vf" ) { readVerticesField( ss ); continue; } // tangent vector
//        if( token == "f"  ) { readFace( ss ); continue; } // face
//        if( token[0] == '#' ) continue; // comment
//        if( token == "o" ) continue; // object name
//        if( token == "g" ) continue; // group name
//        if( token == "s" ) continue; // smoothing group
//        if( token == "mtllib" ) continue; // material library
//        if( token == "usemtl" ) continue; // material
//        if( token == "k" ) continue; // field degree
//        if( token == "fs" ) continue; // field singularity
//        if( token == "" ) continue; // empty string
//        if( token == "vp" ) continue;// principal field, ignore here

//        cerr << "Error: does not appear to be a valid Wavefront OBJ file!" << endl;
//        cerr << "(Offending line: " << oneline << ")" << endl;
//        return false;
//    }
//    //cout<<"nfvec "<<nfvec<<endl;


//    reader.close();

//    cout<<vertices_field.size()<<' '<<vertices.size()<<' '<<faces2vertices.size()<<endl;
//    setparameters();
//    cout<<n_vertices<<' '<<n_faces<<endl;

//    auto p_fv = fv_begin(0);
//    cout<<p_fv[0]<<' '<<p_fv[1]<<' '<<p_fv[2]<<endl;

//    isfacesfield = false;
//    return true;

//}

void Surface::clearup(){

    reset();



    faces_center.clear();
    edge_len.clear();
    edge_vec.clear();
    edge_cot.clear();

    vertices_area.clear();
    faces_area.clear();
    faces_inverse_area.clear();
    face_Material1.clear();
    face_Material2.clear();

    weighted_color.clear();
    weighted_fcolor.clear();

    display_vertices.clear();
    display_normal.clear();
    display_edges.clear();

    display_vcolor.clear();
    display_faces.clear();

    invert_faces_colordegree.clear();
    invert_vertices_colordegree.clear();

    isbuilddisp = false;
    //cout<<"clearup!!!!!!!!"<<endl;


}



double Surface::computeRotationAlongAxis(const double angle, const double *norAxis, const double *vec, double *vecout){

    Eigen::Vector3d axis,invec;
    for(int i=0;i<3;++i)axis(i)= norAxis[i];
    for(int i=0;i<3;++i)invec(i)= vec[i];
    normalize(axis.data());

    Eigen::AngleAxisd tb(angle,axis);
    Eigen::Vector3d outvec = tb*invec;
    for(int i=0;i<3;++i)vecout[i]= outvec(i);


}



void Surface::BuildPickID(){
    if(n_vertices==0)return;
    if(n_faces==0)return;

    id_pick_faces.resize(  n_faces*3);
    for(uint i=0;i<id_pick_faces.size();++i){
        id_pick_faces[i] = i;
    }


    id_pick_vertices.resize(n_faces*3*3);
    for(uint i=0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        auto ind = i*3*3;

        for(uint j=0;j<3;++j){
            auto p_v = v_begin(p_fv[j]);
            for(int k=0;k<3;++k)id_pick_vertices[ind+j*3+k] = p_v[k];
        }

    }


    id_pick_color.resize(n_faces*3*3*4);
    unsigned char cc[3];

    for(uint i=0;i<n_faces;++i){

        cc[0] = ( i & 0xFF );
        cc[1] = ( i & 0xFF00 ) >> 8;
        cc[2] = ( i & 0xFF0000 ) >> 16;
        cc[3] = ( i & 0xFF000000 ) >> 24;
        auto ind = i*3*4;

        for(uint j=0;j<3;++j){
            id_pick_color[ind+j*4]=cc[0];
            id_pick_color[ind+j*4+1]=cc[1];
            id_pick_color[ind+j*4+2]=cc[2];
            id_pick_color[ind+j*4+3]=cc[3];
        }


    }

    cout<<"BuildPickID!!!!"<<endl;

}





int Surface::PickFaceViaColor(unsigned char* pcolor){


    int pickIndex = pcolor[0] | (pcolor[1]<<8) | (pcolor[2]<<16) | (pcolor[4]<<24);



    cout<<"PickFaceViaColor: "<<pickIndex<<endl;
    if(pickIndex>=n_faces)pickIndex=-1;
    return pickIndex;

}


bool Surface::readSufFile(string filename){

    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Suf file " << filename << endl;
        return false;
    }
    reset();

    cout<<"Reading Suf File"<<endl;



    reader>>n_vertices;
    reader>>n_faces;

    vertices.resize(n_vertices*3);
    face_Material1.clear();
    face_Material2.clear();
    faces2vertices.clear();

    for(int i =0;i<vertices.size();i++){
        reader>>vertices[i];
    }

    int ivalue;
    for(int i =0;i<n_faces;i++){
        for(int j =0;j<3;++j){reader>>ivalue;faces2vertices.push_back(ivalue);}
        reader>>ivalue;face_Material1.push_back(ivalue);
        reader>>ivalue;face_Material2.push_back(ivalue);
    }

    vector<bool>mmmmm(32,false);
    for(int i =0;i<n_faces;i++){
        mmmmm[face_Material1[i]] = true;
    }
    int numofM = 0;
    for(int i =0;i<32;i++){
        if(!mmmmm[i]){numofM = i;break;}
    }

    vector< vector<uint> >mm_faces(numofM);
    for(int i =0;i<n_faces;i++){
        int m1 = face_Material1[i];
        int m2 = face_Material2[i];
        auto p_fv = fv_begin(i);
        for(int j=0;j<3;++j){mm_faces[m1].push_back(p_fv[j]);}
        for(int j=0;j<3;++j){mm_faces[m2].push_back(p_fv[j]);}


    }

    faces2vertices = mm_faces[0];
    n_faces = faces2vertices.size()/3;
    //cout<<n_vertices<<n_faces<<endl;




}


void Surface::importSurface(vector<double>&vv, vector<uint>&fv, bool isReOrientate, bool isBuild, bool isRescale){

    reset();

    cout<<"importSurfaces"<<endl;
    vertices = vv;
    faces2vertices = fv;

    setparameters();
    if(isRescale)ReScale_uniform();
    if(!isBuild)return;

    BuildNeighborTable();
    if(isReOrientate)ReOrientFaces();
    ComputeEdgefromFace();
    ComputeFaceNormal();
    return;

}

void Surface::MMpropagation(int f2vmethod, vector<uint>&ctredges, vector<int>&ctredgesMat, vector<int>&facesMat, vector<int>&verticesMat){


    if(0){
        facesMat.clear();facesMat.resize(n_faces,1);
        verticesMat.clear();verticesMat.resize(n_vertices,1);
        return;
    }


    int n_ctrE = ctredges.size()/2;
    int n_ctrEE = ctredgesMat.size()/2;
    assert(n_ctrE==n_ctrEE);
    cout<<n_ctrE<<endl;
    vector<uint>ctrmapEInd(n_ctrE,-1);
    for(int i=0;i<n_ctrE;++i){
        uint v1 = ctredges[i*2],v2 = ctredges[i*2+1];
        int vvnum = vv_num(v1);
        auto p_vv = vv_begin(v1);
        for(int j=0;j<vvnum;++j)if(p_vv[j]==v2){
            ctrmapEInd[i] = ve_begin(v1)[j];
            break;
        }
        assert(ctrmapEInd[i]!=-1);
    }



    //cout<<"ctrmapEInd: ";for(auto a:ctrmapEInd)cout<<a<<' ';cout<<endl;

    auto ctrEMat = ctredgesMat;
    facesMat.clear();facesMat.resize(n_faces,-2);

    int n_filled = 0;
    double vecE[3],vecF[3],vecC[3];
    for(int i=0;i<n_ctrE;++i){
        uint eind = ctrmapEInd[i];
        auto p_ef = ef_begin(eind);
        auto p_ev = &(ctredges[i*2]);
        uint v1 = p_ev[0],v2 = p_ev[1];
        auto p_fv = fv_begin(p_ef[0]);
        uint v3 = 999999;
        for(int k=0;k<3;++k)if(p_fv[k]!=v1 && p_fv[k]!=v2){
            v3 = p_fv[k];break;
        }
        vectorize(v2,v1,vecE);
        vectorize(v3,v1,vecF);
        cross(vnor_begin(v1),vecE,vecC);
        if(dot(vecC,vecF)<0)swap(ctrEMat[i*2],ctrEMat[i*2+1]);
    }
    for(int i=0;i<n_ctrE;++i){
        uint eind = ctrmapEInd[i];
        uint efnum = ef_num(eind);
        auto p_ef = ef_begin(eind);
        for(uint j=0;j<efnum;++j){
            auto &fm = facesMat[p_ef[j]];
            if(fm==-2){
                fm = ctrEMat[i*2+j];
            }else{
                //if material didn't match, fail!
                assert(fm==ctrEMat[i*2+j]);
                if(fm!=ctrEMat[i*2+j])
                    cout<<i<<endl;
            }
        }
    }



    list<int>faceBuff;
    vector<bool>isVisited(n_faces,false);
    vector<bool>isInitBound(n_faces,false);
    for(auto a:facesMat)if(a!=-2)++n_filled;
    for(int i=0;i<n_faces;++i)if(facesMat[i]!=-2)isVisited[i] = isInitBound[i] = true;
    assert(n_filled>0);
    while(true){
        int seed = -1;
        for(int i=0;i<n_faces;++i)if(facesMat[i]==-2 || !isVisited[i]){
            seed=i;break;
        }
        if(seed == -1)break;
        vector<int>componentBuf;
        componentBuf.push_back(seed);
        //cout<<"seed: "<<seed<<endl;
        faceBuff.push_back(seed);
        int curM = -1;
        while(!faceBuff.empty()){
            int curfa = faceBuff.front();
            faceBuff.pop_front();
            //cout<<curfa<<' '<<isVisited[curfa]<<' '<<n_faces-n_filled<<endl;
            if(isVisited[curfa])continue;
            isVisited[curfa] = true;
            componentBuf.push_back(curfa);
            ++n_filled;
            //cout<<"curfa "<<curfa<<' '<<curM<<' '<<faceBuff.empty()<<endl;

            for(auto p_ff = ff_begin(curfa);p_ff!= ff_end(curfa);++p_ff){
                if(isInitBound[*p_ff]){
                    if(curM!=-1);//assert(facesMat[*p_ff] == curM);
                    else curM = facesMat[*p_ff];
                }
                if(isVisited[*p_ff])continue;
                faceBuff.push_back(*p_ff);
            }
        }
        assert(curM!=-1);
        for(auto a:componentBuf)facesMat[a] = curM;

    }

    for(int i=0;i<n_faces;++i)if(!isInitBound[i])
        for(auto p_ff = ff_begin(i);p_ff!= ff_end(i);++p_ff){
            assert(facesMat[*p_ff]==facesMat[i]);
        }

    //cout<<"facesMat: ";for(auto a:facesMat)cout<<a<<' ';cout<<endl;

    if(f2vmethod==3){

        verticesMat.resize(n_vertices);
        vector<int>orderMapping({6,5,4,3,12,11,0});
        unordered_map<int,int>inversemapping;
        for(int i=0;i<orderMapping.size();++i)inversemapping[orderMapping[i]] = i;
        for(int i=0;i<n_vertices;++i){
            int maxM = -1;
            for(auto p_vf = vf_begin(i);p_vf!=vf_end(i);++p_vf)maxM = max(maxM,orderMapping[facesMat[*p_vf]]);
            verticesMat[i] = inversemapping[maxM];
        }

    }
    else if(f2vmethod ==2){
        TopoPreservingLabling(facesMat,verticesMat);


    }if(f2vmethod == 4){
        TopoPreservingLabling2(facesMat,verticesMat);


    }else{
        verticesMat.resize(n_vertices);
        for(int i=0;i<n_vertices;++i){
            if(f2vmethod == 0){
                int maxM = -1;
                for(auto p_vf = vf_begin(i);p_vf!=vf_end(i);++p_vf)maxM = max(maxM,facesMat[*p_vf]);
                verticesMat[i] = maxM;
            }else if(f2vmethod == 1){

                int minM = 9999;
                for(auto p_vf = vf_begin(i);p_vf!=vf_end(i);++p_vf)minM = min(minM,facesMat[*p_vf]);
                verticesMat[i] = minM;
            }

        }
    }


}
void Surface::TopoPreservingLabling2(vector<int>&facesMat,vector<int>&verticesMat){

    BuildUp(false);
    set<int>lableset;
    for(auto a:facesMat)lableset.insert(a);



    UnionFind unifind;
    unordered_map<int,int>aWeight;
    vector<int>labelv;
    for(auto a:lableset){

        vector<int>epool,vpool;
        int stateO,stateA,stateB;

        for(int  i=0;i<n_edges;++i)if(efnon_num(i)==2){
            auto p_ef = efnon_begin(i);
            if( (facesMat[p_ef[0]]!=a && facesMat[p_ef[1]]==a) || (facesMat[p_ef[0]]==a && facesMat[p_ef[1]]!=a))epool.push_back(i);
        }
        vector<bool>vpick(n_vertices,false);
        for(auto b:epool){
            auto p_ev = ev_begin(b);
            for(int i=0;i<2;++i)vpick[p_ev[i]] = true;
        }
        for(int  i=0;i<n_vertices;++i)if(vpick[i])vpool.push_back(i);
        unifind.SetElements(vpool);
        for(auto b:epool){
            auto p_ev = ev_begin(b);
            unifind.Union(p_ev[0],p_ev[1]);
        }
        stateO = unifind.GetNumOfComponents();



        vector<int>vpoolsA(n_vertices,0);
        vector<int>vpoolsB(n_vertices,1);
        for(int i=0;i<n_faces;++i){
            bool fm = facesMat[i] == a;
            auto fv = fv_begin(i);
            if(fm)for(int j=0;j<3;++j)vpoolsA[fv[j]] = 1;
            if(!fm)for(int j=0;j<3;++j)vpoolsB[fv[j]] = 0;
        }

        hiddenCtr(vpoolsA,true);
        if(hidden_ctrV.size()==0)stateA=-1;
        else{

            vector<int>ele(hidden_ctrV.size()/3);
            for(int i =0;i<ele.size();++i)ele[i]=i;
            unifind.SetElements(ele);
            for(int i=0;i<hidden_ctrE.size()/2;++i)unifind.Union(hidden_ctrE[i*2],hidden_ctrE[i*2+1]);
            stateA = unifind.GetNumOfComponents();
        }

        hiddenCtr(vpoolsB,true);
        if(hidden_ctrV.size()==0)stateB=-1;
        else{
            vector<int>ele(hidden_ctrV.size()/3);
            for(int i =0;i<ele.size();++i)ele[i]=i;
            unifind.SetElements(ele);
            for(int i=0;i<hidden_ctrE.size()/2;++i)unifind.Union(hidden_ctrE[i*2],hidden_ctrE[i*2+1]);
            stateB = unifind.GetNumOfComponents();

            ele.clear();
            for(int i=0;i<n_vertices;++i)if(vpoolsB[i]==1)ele.push_back(i);
            unifind.SetElements(ele);
            for(int  i=0;i<n_edges;++i){
                auto p_ev = ev_begin(i);
                if(vpoolsB[p_ev[0]]==1 && vpoolsB[p_ev[1]]==1)unifind.Union(p_ev[0],p_ev[1]);
            }
            vector<vector<int>>comptmp;
            unifind.ExtractComponents(comptmp);
            for(auto &c:comptmp)if(c.size()<3){stateB=-1;break;}
        }

        if(stateA==stateB)aWeight[a] = 0;
        else {
            if(stateO==stateA)aWeight[a] = 1;
            else aWeight[a] = 0;
        }
        if(a==2)aWeight[a] = 022;
        //if(a==3)aWeight[a] = 026;
        //else if(a==6)aWeight[aWeight.size()-1] = 0;
        labelv.push_back(a);

    }

    auto comparator = [&aWeight](int a, int b) { return aWeight[a] < aWeight[b]; };
    sort(labelv.begin(), labelv.end(), comparator);

    unordered_map<int,int>indexxx;
    for(int i=0;i<labelv.size();++i)indexxx[labelv[i]] = i;

    vector<int>priorityW(n_vertices,-1);

    verticesMat.resize(n_vertices);
    for(int i=0;i<n_faces;++i){
        auto fmW = indexxx[facesMat[i]];
        auto fv = fv_begin(i);
        for(int j=0;j<3;++j)if(fmW>priorityW[fv[j]]){
            verticesMat[fv[j]] = facesMat[i];
            priorityW[fv[j]] = fmW;
        }

    }






}
void Surface::TopoPreservingLabling(vector<int>&facesMat,vector<int>&verticesMat){

    BuildUp(false);
    set<int>labelset;
    for(auto a:facesMat)labelset.insert(a);



    UnionFind unifind;
    unordered_map<int,int>aWeight;
    vector<int>labelv;
    for(auto a:labelset){

        vector<int>epool,vpool;
        int stateO,stateA,stateB;

        for(int  i=0;i<n_edges;++i)if(efnon_num(i)==2){
            auto p_ef = efnon_begin(i);
            if( (facesMat[p_ef[0]]!=a && facesMat[p_ef[1]]==a) || (facesMat[p_ef[0]]==a && facesMat[p_ef[1]]!=a))epool.push_back(i);
        }
        vector<bool>vpick(n_vertices,false);
        for(auto b:epool){
            auto p_ev = ev_begin(b);
            for(int i=0;i<2;++i)vpick[p_ev[i]] = true;
        }
        for(int  i=0;i<n_vertices;++i)if(vpick[i])vpool.push_back(i);
        unifind.SetElements(vpool);
        for(auto b:epool){
            auto p_ev = ev_begin(b);
            unifind.Union(p_ev[0],p_ev[1]);
        }
        stateO = unifind.GetNumOfComponents();



        vector<int>vpoolsA(n_vertices,0);
        vector<int>vpoolsB(n_vertices,1);
        for(int i=0;i<n_faces;++i){
            bool fm = facesMat[i] == a;
            auto fv = fv_begin(i);
            if(fm)for(int j=0;j<3;++j)vpoolsA[fv[j]] = 1;
            if(!fm)for(int j=0;j<3;++j)vpoolsB[fv[j]] = 0;
        }

        hiddenCtr(vpoolsA,true);
        if(hidden_ctrV.size()==0)stateA=-1;
        else{

            vector<int>ele(hidden_ctrV.size()/3);
            for(int i =0;i<ele.size();++i)ele[i]=i;
            unifind.SetElements(ele);
            for(int i=0;i<hidden_ctrE.size()/2;++i)unifind.Union(hidden_ctrE[i*2],hidden_ctrE[i*2+1]);
            stateA = unifind.GetNumOfComponents();
        }

        hiddenCtr(vpoolsB,true);
        if(hidden_ctrV.size()==0)stateB=-1;
        else{
            vector<int>ele(hidden_ctrV.size()/3);
            for(int i =0;i<ele.size();++i)ele[i]=i;
            unifind.SetElements(ele);
            for(int i=0;i<hidden_ctrE.size()/2;++i)unifind.Union(hidden_ctrE[i*2],hidden_ctrE[i*2+1]);
            stateB = unifind.GetNumOfComponents();

            ele.clear();
            for(int i=0;i<n_vertices;++i)if(vpoolsB[i]==1)ele.push_back(i);
            unifind.SetElements(ele);
            for(int  i=0;i<n_edges;++i){
                auto p_ev = ev_begin(i);
                if(vpoolsB[p_ev[0]]==1 && vpoolsB[p_ev[1]]==1)unifind.Union(p_ev[0],p_ev[1]);
            }
            vector<vector<int>>comptmp;
            unifind.ExtractComponents(comptmp);
            for(auto &c:comptmp)if(c.size()<3){stateB=-1;break;}
        }

        if(stateA==stateB)aWeight[a] = 0;
        else {
            if(stateO==stateA)aWeight[a] = 1;
            else aWeight[a] = 0;
        }

        labelv.push_back(a);


    }

    //    for(auto a:lableset)cout<<a<<' ';cout<<endl;
    //    for(auto a:aWeight)cout<<a<<' ';cout<<endl;



    auto comparator = [&aWeight](int a, int b) { return aWeight[a] < aWeight[b]; };
    sort(labelv.begin(), labelv.end(), comparator);

    unordered_map<int,int>indexxx;
    for(int i=0;i<labelv.size();++i)indexxx[labelv[i]] = i;

    vector<int>priorityW(n_vertices,-1);

    verticesMat.resize(n_vertices);
    for(int i=0;i<n_faces;++i){
        auto fmW = indexxx[facesMat[i]];
        auto fv = fv_begin(i);
        for(int j=0;j<3;++j)if(fmW>priorityW[fv[j]]){
            verticesMat[fv[j]] = facesMat[i];
            priorityW[fv[j]] = fmW;
        }

    }






}

void Surface::hiddenCtr(vector<int>&verticesLable,bool isOnlyGatherCtrInfo){

    if(faces_center.size()!=n_vertices*3)BuildFacesCenter();
    assert(verticesLable.size()==n_vertices);

    auto maxn = *max_element(verticesLable.begin(),verticesLable.end());
    set<int>nLable;
    for(auto a:verticesLable)nLable.insert(a);
    if(maxn==0 || nLable.size()<=1){
        hidden_ctrE.clear();
        hidden_ctrV.clear();
        hidden_ctrVnor.clear();
        return;
    }

    hidden_ctrE.clear();
    hidden_ctrV.clear();
    hidden_ctrVnor.clear();
    hidden_ctrActiveE.clear();
    vector<int>activeF(n_faces,-1);
    vector<int>activeNonE(n_edges,-1);
    hidden_ctrActiveF2V.clear();
    hidden_ctrVCreator.clear();
    hidden_ctrActiveNonE2V.clear();
    int aVind = 0,aEind = 0;
    for(int i=0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        if(verticesLable[p_fv[0]]!=verticesLable[p_fv[1]] || verticesLable[p_fv[1]]!=verticesLable[p_fv[2]] || verticesLable[p_fv[0]]!=verticesLable[p_fv[2]] ){
            activeF[i]=aVind;++aVind;
            for(int j=0;j<3;++j)hidden_ctrV.push_back(faces_center[i*3+j]);
            for(int j=0;j<3;++j)hidden_ctrVnor.push_back(faces_normal[i*3+j]);
            for(int j=0;j<3;++j)hidden_ctrActiveF2V.push_back(p_fv[j]);
            hidden_ctrVCreator.push_back(i+1);
        }
    }


    vector<int>edgeCenterInd(n_edges,-1);
    vector<uint>dualEdges;
    vector<int>referFaces;


    double ecent[3];
    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        if(verticesLable[p_ev[0]]!=verticesLable[p_ev[1]]){
            auto p_ef = efnon_begin(i);
            int efnum = efnon_num(i);
            if(efnum==2){
                if(activeF[p_ef[0]]!=-1 && activeF[p_ef[1]]!=-1){
                    hidden_ctrE.push_back(activeF[p_ef[0]]);hidden_ctrE.push_back(activeF[p_ef[1]]);
                    //hidden_ctrActiveE.push_back(i);
                    dualEdges.push_back(i);referFaces.push_back(p_ef[0]);
                }

            }else{
                edgeCenterInd[i] = aVind;++aVind;
                VerticesMidpoint(p_ev[0],p_ev[1],ecent);
                for(int j=0;j<3;++j)hidden_ctrV.push_back(ecent[j]);
                for(int j=0;j<3;++j)hidden_ctrVnor.push_back(faces_normal[3*p_ef[0]+j]);
                hidden_ctrVCreator.push_back(-i);
                for(int j=0;j<efnum;++j)if(activeF[p_ef[j]]!=-1){
                    hidden_ctrE.push_back(activeF[p_ef[j]]);hidden_ctrE.push_back(edgeCenterInd[i]);
                    //hidden_ctrActiveE.push_back(i);
                    dualEdges.push_back(i);referFaces.push_back(p_ef[j]);
                }

                activeNonE[i]=aEind++;
                hidden_ctrActiveNonE2V.push_back(p_ev[0]);hidden_ctrActiveNonE2V.push_back(p_ev[1]);
            }
        }
    }

    int n_ctrE = dualEdges.size();
    double evec[3],dualevec[3],crosVec[3];
    auto p_ctrv = hidden_ctrV.data();
    for(int i=0;i<n_ctrE;++i){
        auto p_ev = &(hidden_ctrE[i*2]);
        auto p_dualev = ev_begin(dualEdges[i]);
        minusVec(p_ctrv+p_ev[1]*3,p_ctrv+3*p_ev[0],evec);
        minusVec(v_begin(p_dualev[1]),v_begin(p_dualev[0]),dualevec);
        cross(evec,fnor_begin(referFaces[i]),crosVec);
        if(dot(crosVec,dualevec)>0){
            hidden_ctrELable.push_back(verticesLable[p_dualev[0]]);
            hidden_ctrELable.push_back(verticesLable[p_dualev[1]]);
        }else{
            hidden_ctrELable.push_back(verticesLable[p_dualev[1]]);
            hidden_ctrELable.push_back(verticesLable[p_dualev[0]]);
        }

    }

    if(hidden_F2Cells.size()!=0){
        hidden_ctrE2Cells.clear();
        for(int i=0;i<n_ctrE;++i){
            auto re = referFaces[i];
            hidden_ctrE2Cells.push_back(hidden_F2Cells[re*2]);
            hidden_ctrE2Cells.push_back(hidden_F2Cells[re*2+1]);
        }

    }

    hidden_ctrActiveE = dualEdges;
    hidden_ctrActiveE2V.clear();
    for(auto a:hidden_ctrActiveE){
        hidden_ctrActiveE2V.push_back(edges[a*2]);
        hidden_ctrActiveE2V.push_back(edges[a*2+1]);
    }

    if(isOnlyGatherCtrInfo)return;
    hiddenCtrSmoothing();
    mapActiveF2Vpos.clear();
    mapActiveNonE2Vpos.clear();

    for(int i=0;i<hidden_ctrVCreator.size();++i){
        vector<float>vpos(3);
        auto p_v = hidden_ctrV.data()+i*3;
        for(int j=0;j<3;++j)vpos[j] = p_v[j];
        if(hidden_ctrVCreator[i]>0){
            mapActiveF2Vpos[activeF[hidden_ctrVCreator[i]-1]] = vpos;
        }else{
            mapActiveNonE2Vpos[activeNonE[-hidden_ctrVCreator[i]]] = vpos;
        }
    }

    int aFind = aVind - aEind;
    vector<int>iiiiFInd(aFind),iiiiEInd(aEind);
    for(int i=0;i<aFind;++i)iiiiFInd[i] = i;
    for(int i=0;i<aEind;++i)iiiiEInd[i] = aFind + i;
    hidden_writeCtrE.clear();

    for(auto a:hidden_ctrE)if(hidden_ctrVCreator[a]>0){
        hidden_writeCtrE.push_back(iiiiFInd[activeF[hidden_ctrVCreator[a]-1]]);
    }else{
        hidden_writeCtrE.push_back(iiiiEInd[activeNonE[-hidden_ctrVCreator[a]]]);
    }

    hidden_CtrE2CSs.resize(n_ctrE);
    if(hidden_F2CSs.size()!=n_faces){
        cout<<"hiddenCtr: No hidden_F2CSs!"<<endl;
        hidden_F2CSs.clear();hidden_F2CSs.resize(n_faces,0);
    }
    for(int i=0;i<n_ctrE;++i){

        hidden_CtrE2CSs[i] = hidden_F2CSs[referFaces[i]];
    }
    //cout<<aFind<<' '<<aEind<<endl;
    //exit(0);



}

void combinationHelper(int s, int e, int k, vector<vector<int>>& res,
                       vector<int>& buf) {
    if (k <= 0) {
        res.push_back(buf);
        return;
    }
    for (int i = s; i <= e; ++i) {
        buf.push_back(i);
        combinationHelper(i + 1, e, k - 1, res, buf);
        buf.pop_back();
    }
}

void combinations(int n, int k, vector<vector<int>>& res) {
    res.clear();
    vector<int> buf;
    combinationHelper(0, n - 1, min(n, k), res, buf);
}
void Surface::hiddenCtrSmoothing(){

    cout<<"hiddenCtrSmoothing"<<endl;
    int n_hiddenedges = hidden_ctrE.size()/2;
    int n_hiddenvertices = hidden_ctrV.size()/3;
    vector<vector<int>>hidden_vertices2vertices(n_hiddenvertices);
    for(int i=0;i<n_hiddenedges;++i){
        auto p_ev = hidden_ctrE.data()+i*2;
        hidden_vertices2vertices[p_ev[0]].push_back(p_ev[1]);
        hidden_vertices2vertices[p_ev[1]].push_back(p_ev[0]);
    }
    int n_afv = 0,n_aev = 0;
    vector<int>varInd(n_hiddenvertices);
    for(int i=0;i<hidden_ctrVCreator.size();++i){
        if(hidden_ctrVCreator[i] >0)varInd[i] = n_afv++;
        else varInd[i] = n_aev++;
    }





    try {
        // Solve for the new locations of each orgCrvVs and store them in crvVs.
        // One curve vertex v can be represented as a linear interpolation of the 3
        // vertices of the triangle, in which v lies.
        // xi = a*v1 + b*v2 + c*v3, where a+b+c=1, a,b,c>0.
        // So the variable we are going to solve are those a,b,c for each crvV.
        // obj:
        //      min (smoothTerm + lambda*intactTerm)
        //      intactTerm = Sum_i( ||xi - triCentroidi||^2 )
        //      smoothTerm = Sum_{i,k,j}( ||xk - (xi - xj)/2||^2 )
        // s.t.:
        //      ai+bi+ci = 1
        //      ai, bi, ci >=0 (Actually we want > 0 to avoid duplicated points, but
        //      it is not possible in opt. intactTerm plays some role of keeping
        //      this from happening.)
        //
        // double lambda = 1.0;
        double lambda = 0.0;
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);
        // Quiet output from Gurobi.
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        vector<GRBVar> A(n_afv);
        vector<GRBVar> B(n_afv);
        vector<GRBVar> C(n_afv);
        vector<GRBVar> E(n_aev);

        vector<GRBVar> X(n_hiddenvertices);
        vector<GRBVar> Y(n_hiddenvertices);
        vector<GRBVar> Z(n_hiddenvertices);



        // Create variables.
        for (int i = 0; i < n_afv; ++i) {
            A[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            B[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            C[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
        }
        for (int i = 0; i < n_hiddenvertices; ++i) {
            X[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            Y[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            Z[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }
        for (int i = 0; i < n_aev; ++i) {
            E[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
        }
        model.update();

        // Set objective, obj = smoothTerm + lambda*intactTerm.
        GRBQuadExpr obj = 0;
        GRBQuadExpr intactTerm = 0;
        for (int i = 0; i < n_hiddenvertices; ++i) {
            auto p_centroid = hidden_ctrV.data()+3*i;
            intactTerm += (X[i] - p_centroid[0]) * (X[i] - p_centroid[0]) +
                    (Y[i] - p_centroid[1]) * (Y[i] - p_centroid[1]) +
                    (Z[i] - p_centroid[2]) * (Z[i] - p_centroid[2]);
        }
        GRBQuadExpr smoothTerm = 0;
        GRBQuadExpr elasticTerm = 0;
        vector<vector<int>> res;
        int nonCounter = 0;
        for (int j = 0; j < n_hiddenvertices; ++j) {
            auto &v2v = hidden_vertices2vertices[j];
            int nv2v = hidden_vertices2vertices[j].size();
            double ww = 1./double(nv2v);
            double ww2 = hidden_ctrVCreator[j]>0 ? 1. : 1.;
            GRBLinExpr tmp;
            if(nv2v<3){
                tmp = X[j];
                for(auto a:v2v)tmp -= ww*X[a];
                smoothTerm += ww2 * tmp * tmp;
                tmp = Y[j];
                for(auto a:v2v)tmp -= ww*Y[a];
                smoothTerm += ww2 * tmp * tmp;
                tmp = Z[j];
                for(auto a:v2v)tmp -= ww*Z[a];
                smoothTerm += ww2 * tmp * tmp;
            }else{
                combinations(nv2v,2,res);
                for(auto &cbs:res){
                    int v0 = v2v[cbs[0]],v1 = v2v[cbs[1]];
                    if(hidden_ctrVCreator[v0]<1 || hidden_ctrVCreator[v1]<1)continue;
                    if(hidden_F2CSs[hidden_ctrVCreator[v0]-1]!=hidden_F2CSs[hidden_ctrVCreator[v1]-1])continue;
                    tmp -= X[j] - 0.5*X[v0] - 0.5*X[v1];
                    smoothTerm += ww2 * tmp * tmp;
                    tmp -= Y[j] - 0.5*Y[v0] - 0.5*Y[v1];
                    smoothTerm += ww2 * tmp * tmp;
                    tmp -= Z[j] - 0.5*Z[v0] - 0.5*Z[v1];
                    smoothTerm += ww2 * tmp * tmp;
                    ++nonCounter;
                }
            }
            if(0){
                for(auto a:v2v){
                    tmp = X[j] - X[a];
                    elasticTerm -= tmp * tmp;
                    tmp = Y[j] - Y[a];
                    elasticTerm -= tmp * tmp;
                    tmp = Z[j] - Z[a];
                    elasticTerm -= tmp * tmp;
                }
            }
        }


        //cout<<nonCounter<<endl;exit(0);


        obj = smoothTerm + lambda * intactTerm ;
        model.setObjective(obj, GRB_MINIMIZE);

        // Add linear constraint.
        for (int i = 0; i < n_hiddenvertices; ++i){
            if(hidden_ctrVCreator[i]>0){
                const int activeF = hidden_ctrVCreator[i]-1;
                auto p_fv = fv_begin(activeF);
                auto v0 = v_begin(p_fv[0]);
                auto v1 = v_begin(p_fv[1]);
                auto v2 = v_begin(p_fv[2]);
                int var_i = varInd[i];
                model.addConstr(A[var_i] + B[var_i] + C[var_i] == 1.0);
                model.addConstr(X[i] - v0[0] * A[var_i] - v1[0] * B[var_i] - v2[0] * C[var_i] == 0.0);
                model.addConstr(Y[i] - v0[1] * A[var_i] - v1[1] * B[var_i] - v2[1] * C[var_i] == 0.0);
                model.addConstr(Z[i] - v0[2] * A[var_i] - v1[2] * B[var_i] - v2[2] * C[var_i] == 0.0);
            }else {
                const int activeE = -hidden_ctrVCreator[i];
                auto p_ev = ev_begin(activeE);
                auto v0 = v_begin(p_ev[0]);
                auto v1 = v_begin(p_ev[1]);
                int var_i = varInd[i];
                model.addConstr(X[i] - v0[0] * E[var_i] - v1[0] * (1-E[var_i]) == 0.0);
                model.addConstr(Y[i] - v0[1] * E[var_i] - v1[1] * (1-E[var_i]) == 0.0);
                model.addConstr(Z[i] - v0[2] * E[var_i] - v1[2] * (1-E[var_i]) == 0.0);

            }
        }

        // Optimize model.
        model.optimize();

        // Extract Results.

        cout << "Obj: " << model.get(GRB_IntAttr_Status) << endl;
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            //cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
            for (int i = 0; i < n_hiddenvertices; ++i) {
                auto p_vp = hidden_ctrV.data()+i*3;

                p_vp[0] = X[i].get(GRB_DoubleAttr_X);
                p_vp[1] = Y[i].get(GRB_DoubleAttr_X);
                p_vp[2] = Z[i].get(GRB_DoubleAttr_X);
            }
        }
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }


}


void Surface::SpecialRouteforFCM(){


    string filename("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/toGurobiV.txt");
    ifstream reader(filename.data(), ofstream::in);
    if (!reader.good()) {
        cout << "Can not open file: " << filename << endl;
    }else {
        cout << "Reading: "<<filename<<endl;
    }


    int tmp;
    int nV;
    reader>>nV;
    vector<double>vpos(nV*3);
    for(auto &a:vpos)reader>>a;
    reader.close();

    reader.open("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/toGurobiE.txt");
    int nP;

    reader>>nP;
    vector<vector<int>>plane2edges(nP);
    for(int i=0;i<nP;i++){
        reader>>tmp;
        plane2edges[i].resize(tmp*2);
        for(auto &a:plane2edges[i]){
            reader>>a;--a;
        }
    }
    reader.close();

    reader.open("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/toGurobiL.txt");
    int nL;

    reader>>nL;
    vector<vector<int>>lines(nL);
    for(int i=0;i<nL;i++){
        reader>>tmp;
        lines[i].resize(tmp);
        for(auto &a:lines[i]){
            reader>>a;--a;
        }
    }
    reader.close();

    reader.open("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/toGurobivertValsOld.txt");
    int nMat;

    reader>>nV;reader>>nP;reader>>nMat;
    vector<vector<double>>oldval(nV);
    for(int i=0;i<nV;i++){
        oldval[i].resize(nP*nMat);
        for(auto &a:oldval[i]){
            reader>>a;
        }
    }
    reader.close();


    reader.open("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/toGurobivertValsNew.txt");

    reader>>nV;reader>>nP;reader>>nMat;
    vector<vector<double>>newval(nV);
    for(int i=0;i<nV;i++){
        newval[i].resize(nP*nMat);
        for(auto &a:newval[i]){
            reader>>a;
        }
    }
    reader.close();


    //nP=1;
    vector<vector<int>>plane2vers(nP);
    vector<vector<int>>plane2verInv(nP);
    for(int i=0;i<nP;++i){
        auto &Vflag = plane2verInv[i];
        Vflag.clear();Vflag.resize(nV,-1);
        for(auto a:plane2edges[i])Vflag[a] = 0;
        auto& vers = plane2vers[i];
        vers.clear();
        int vind = 0;
        for(int j=0;j<nV;++j)if(Vflag[j]!=-1){vers.push_back(j);Vflag[j]=vind++;}
    }





    vector<int>nConsv(nP);
    vector<int>nUnConsv(nP);
    vector<vector<double>>pEn(nP);
    vector<vector<double>>pTran(nP);

    vector<int>nConsVertices(nP);
    vector<int>nUnConsVertices(nP);


    for(int i=0;i<nP;++i){

        FILE* fp = NULL;
        string finame = "/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/newH"+to_string(i+1)+".dat";
        fp = fopen(finame.data(),"r");
        if(!fp)cout<<"open fail"<<endl;


        fread(&(nConsv[i]),sizeof(int),1,fp);
        pEn[i].resize(nConsv[i]*nConsv[i]);
        fread(pEn[i].data(),sizeof(double),nConsv[i]*nConsv[i],fp);


        fread(&(nUnConsv[i]),sizeof(int),1,fp);
        pTran[i].resize(nUnConsv[i]*nConsv[i]);
        fread(pTran[i].data(),sizeof(double),nUnConsv[i]*nConsv[i],fp);

        fclose(fp);

        nConsVertices[i] = nConsv[i]/nMat;
        nUnConsVertices[i] = nUnConsv[i]/nMat;
    }

    double threshold = 1e-3;
    for(int i=0;i<nP;++i){

        int count=0;
        for(auto a:pEn[i])if(a<threshold)++count;
        cout<<double(count)/pEn[i].size()<<' ';

        count=0;
        for(auto a:pTran[i])if(a<threshold)++count;
        cout<<double(count)/pTran[i].size()<<endl;
    }


    int nTotalCons = 0;
    vector<int>mapCons2Binary(nV,-1);
    for(int i=0;i<nL;++i){
        for(auto a:lines[i])mapCons2Binary[a] = 0;
    }
    for(auto &a:mapCons2Binary)if(a!=-1)a=nTotalCons++;


    vector<vector<int>>mapplaneV2globalV(nP);
    for(int i=0;i<nP;++i){
        vector<bool>ispick(nV,false);
        auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
        mapplaneV2globalV_local.clear();
        for(auto a:plane2edges[i])ispick[a] = true;
        for(int j=0;j<nV;++j)if(ispick[j])mapplaneV2globalV_local.push_back(j);

    }

    vector<vector<double>>planev2oldval(nP);
    for(int i=0;i<nP;++i){
        auto &pvov = planev2oldval[i];
        auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
        int npv = mapplaneV2globalV_local.size();
        pvov.resize(npv*nMat);
        for(int j=0;j<npv;++j){
            auto pv = pvov.data()+j*nMat;
            auto pvo = oldval[mapplaneV2globalV_local[j]].data()+i*nMat;
            for(int k=0;k<nMat;++k){
                pv[k] = pvo[k];
            }
        }
    }

    vector<vector<double>>planev2newval(nP);
    for(int i=0;i<nP;++i){
        auto &pvov = planev2newval[i];
        auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
        int npv = mapplaneV2globalV_local.size();
        pvov.resize(npv*nMat);
        for(int j=0;j<npv;++j){
            auto pv = pvov.data()+j*nMat;
            auto pvo = newval[mapplaneV2globalV_local[j]].data()+i*nMat;
            for(int k=0;k<nMat;++k){
                pv[k] = pvo[k];
            }
        }
    }

    vector<int>binaryVarsInit(nTotalCons,-1);
    for(int i=0;i<nP;++i){
        auto &pvov = planev2newval[i];
        auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
        int npv = mapplaneV2globalV_local.size();
        pvov.resize(npv*nMat);
        for(int j=0;j<npv;++j){
            auto pv = pvov.begin()+j*nMat;
            int vbinary = mapCons2Binary[mapplaneV2globalV_local[j]];
            if(vbinary==-1)continue;
            int maxind = max_element(pv,pv+nMat)-pv;
            //if(binaryVarsInit[vbinary]!=-1)assert(binaryVarsInit[vbinary]==maxind);
            binaryVarsInit[vbinary] = maxind;


        }
    }

    cout<<"read end"<<endl;


    vector<vector<double>>solval = newval;
    
    /*************************************************************/
    /*************************************************************/


    try{

        cout << "optimize 0" << endl;
        GRBEnv env = GRBEnv();
        cout << "optimize 0" << endl;
        GRBModel model = GRBModel(env);
        cout << "optimize 0" << endl;
        // Quiet output from Gurobi.
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);

        vector<vector<GRBVar>> conVars(nP);
        vector<GRBVar>binaryVars(nTotalCons*nMat);



        cout << "optimize 0" << endl;

        for(int i=0;i<nP;++i){
            auto &conVar = conVars[i];
            conVar.resize(nConsv[i]);
            auto &pvov = planev2newval[i];
            auto &pvovold = planev2oldval[i];

            for(int j=0;j<nConsv[i];++j){
                conVar[j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, pvov[j]-pvovold[j], GRB_CONTINUOUS);
            }
        }
        //        for(int j=0,k = nTotalCons*nMat;j<k;++j){
        //            binaryVars[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        //        }
        for(int j=0;j<nTotalCons;++j)for(int k = 0;k<nMat;++k){
            if(k==binaryVarsInit[j])binaryVars[j*nMat+k] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY);
            else binaryVars[j*nMat+k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }


        cout << "optimize 1" << endl;
        GRBQuadExpr obj = 0;


        for(int k=0;k<nP;++k){

            auto &en = pEn[k];
            auto &variables = conVars[k];
            int nCon = nConsv[k];
            for(int i=0;i<nCon;++i){
                for(int j=0;j<nCon;++j){
                    double a = en[i*nCon+j];
                    obj+=a*variables[i]*variables[j];
                }
            }
            //            double a = 0;
            //            for(int i=0;i<nCon;++i){
            //                for(int j=0;j<nCon;++j){
            //                     a+=fabs(en[i*nCon+j]-en[j*nCon+i]);

            //                }
            //            }

            //            cout<<a<<endl;
        }


        cout << "optimize 2" << endl;
        model.setObjective(obj, GRB_MINIMIZE);
        cout << "optimize 3" << endl;

        //        for(int k=0;k<nP;++k){
        //            auto &variables = conVars[k];
        //            int nCon = nConsv[k];
        //            int nConvertice = nConsVertices[k];
        //            auto &mapplaneV2globalV_local = mapplaneV2globalV[k];
        //            for(int v = 0;v<nConvertice;++v){
        //                int bias = nMat*v;
        //                int binarybias = nMat*mapCons2Binary[mapplaneV2globalV_local[v]];
        //                for(int i=0;i<nMat;++i)for(int j=i+1;j<nMat;++j){
        //                    //model.addQConstr((binaryVars[binarybias+i]-binaryVars[binarybias+j])*(variables[bias+i]-variables[bias+j])>= 0);
        //                    //model.add((binaryVars[binarybias+i]-binaryVars[binarybias+j])*(variables[bias+i]-variables[bias+j])>= 0);
        //                }
        //            }
        //        }
        if(0){
            for(int k=0;k<nP;++k){
                auto &variables = conVars[k];
                int nCon = nConsv[k];
                int nConvertice = nConsVertices[k];
                auto &mapplaneV2globalV_local = mapplaneV2globalV[k];
                auto &oldvalue = planev2oldval[k];
                for(int v = 0;v<nConvertice;++v){
                    int bias = nMat*v;
                    int binarybias = nMat*mapCons2Binary[mapplaneV2globalV_local[v]];
                    auto olv = oldvalue.data()+v*nMat;
                    for(int i=0;i<nMat;++i){
                        for(int j=0;j<nMat;++j)if(j!=i){
                            //GRBLinExpr a = variables[bias+i]-variables[bias+j];
                            model.addGenConstrIndicator(binaryVars[binarybias+i],1,
                                    (variables[bias+i]+olv[i])-(variables[bias+j]+olv[j])>=0.);
                        }
                    }
                }
            }
        }else{


            for(int k=0;k<nP;++k){
                auto &variables = conVars[k];
                int nCon = nConsv[k];
                int nConvertice = nConsVertices[k];
                auto &oldvalue = planev2oldval[k];
                auto &mapplaneV2globalV_local = mapplaneV2globalV[k];
                for(int v = 0;v<nConvertice;++v){
                    int bias = nMat*v;
                    int binarybias = nMat*mapCons2Binary[mapplaneV2globalV_local[v]];
                    auto olv = oldvalue.data()+v*nMat;
                    int picki=binaryVarsInit[mapCons2Binary[mapplaneV2globalV_local[v]]];
                    for(int j=0;j<nMat;++j)if(j!=picki){
                        model.addConstr((variables[bias+picki]+olv[picki])-(variables[bias+j]+olv[j])>=0.);
                    }

                }
            }

        }
        for(int i=0;i<nTotalCons;++i){
            GRBLinExpr a = 0;
            int binarybias = nMat*i;
            for(int j=0;j<nMat;++j)a+=binaryVars[binarybias+j];
            model.addConstr(a==1);

        }


        cout << "optimize begin" << endl;
        model.optimize();
        cout << "optimize end" << endl;

        // Extract Results.

        cout << "Obj: " << model.get(GRB_IntAttr_Status) << endl;
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
            //            for(int i=0;i<nTotalCons;++i){
            //                int binarybias = nMat*i;
            //                for(int j=0;j<nMat;++j)cout<<binaryVars[binarybias+j].get(GRB_DoubleAttr_X) << " ";
            //                cout<<binaryVarsInit[i]<<' ';cout<<endl;
            //            }
            //                        for(int i=0;i<nP;++i){
            //                            auto &conVar = conVars[i];
            //                            conVar.resize(nConsv[i]);

            //                            auto &pvov = planev2newval[i];
            //                            auto &pvovold = planev2oldval[i];
            //                            for(int j=0;j<nConsv[i];++j)if(fabs((conVar[j].get(GRB_DoubleAttr_X)-(pvov[j]-pvovold[j]))/(pvov[j]-pvovold[j]))>-1000){
            //                                cout<<(conVar[j].get(GRB_DoubleAttr_X)-(pvov[j]-pvovold[j]))/(pvov[j]-pvovold[j])<<" "<<conVar[j].get(GRB_DoubleAttr_X)<<' '<<(pvov[j]-pvovold[j])<<endl;
            //                            }
            //                            cout<<endl;
            //                        }

        }

        vector<vector<double>>planev2solveval(nP);

        for(int i=0;i<nP;++i){
            auto &pvov = planev2solveval[i];
            auto &pvovold = planev2oldval[i];
            auto &pvovnew = planev2newval[i];
            auto &mapplaneV2globalV_local = mapplaneV2globalV[i];
            int npv = mapplaneV2globalV_local.size();
            pvov.resize(npv*nMat);
            auto &conVar = conVars[i];
            for(int j=0;j<nConsVertices[i];++j){
                auto psolv = conVar.data()+j*nMat;
                auto pv = pvov.data()+j*nMat;

                auto pvo = solval[mapplaneV2globalV_local[j]].data()+i*nMat;
                for(int k=0;k<nMat;++k){
                    pv[k] = pvo[k] = psolv[k].get(GRB_DoubleAttr_X);
                }
            }

            auto &pTra = pTran[i];
            int w = nUnConsv[i],h = nConsv[i];
            vector<double>matrixMul(w);
            for(int j=0;j<w;++j){
                auto pv = pvov.data();
                double aaa = 0;
                for(int k=0;k<h;++k)aaa+=pTra[k*w+j]*pv[k];
                matrixMul[j] = -aaa;
            }
            int offset = nConsVertices[i];
            for(int j=0;j<nUnConsVertices[i];++j){
                auto pvnew = pvovnew.data()+(j+offset)*nMat;

                auto pvold = pvovold.data()+(j+offset)*nMat;
                auto pv = pvov.data()+(j+offset)*nMat;
                auto pvo = solval[mapplaneV2globalV_local[j+offset]].data()+i*nMat;
                auto pmat = matrixMul.data()+j*nMat;
                for(int k=0;k<nMat;++k){
                    pv[k] = pvo[k] = pmat[k]+pvold[k];
                    //cout<<pvnew[k] - pv[k]<<' ';
                }
            }

            for(int j=0;j<nConsVertices[i];++j){
                auto psolv = conVar.data()+j*nMat;
                auto pv = pvov.data()+j*nMat;
                auto pvold = pvovold.data()+j*nMat;

                auto pvo = solval[mapplaneV2globalV_local[j]].data()+i*nMat;
                for(int k=0;k<nMat;++k){
                    pv[k] = pvo[k] = psolv[k].get(GRB_DoubleAttr_X)+pvold[k];
                }
            }


        }







    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }




    string nfilename("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Gurobisol.txt");
    ofstream writer(nfilename.data(), ofstream::out);
    if (!writer.good()) {
        cout << "Can not open file: " << nfilename << endl;
    }else {
        cout << "writing: "<<nfilename<<endl;
    }

    writer<< std::fixed;
    writer<<"{";
    for(int i =0;i<nV;++i){
        auto &pval = solval[i];
        writer<<"{";
        for(int j=0;j<nP;++j){
            auto pv = pval.data()+nMat*j;
            bool isNull = true;
            for(int k=0;k<nMat;++k)isNull  = isNull&&(pv[k]==-1);
            if(isNull)writer<<"Null";
            else{
                writer<<"{";

                for(int k=0;k<nMat;++k){
                    writer<<pv[k];
                    if(k!=nMat-1)writer<<",";
                }
                writer<<"}";
            }
            if(j!=nP-1)writer<<",";
        }

        writer<<"}";
        if(i!=nV-1)writer<<",";
    }
    writer<<"}";


}


void Surface::exportSurface(int csi,vector<double>&vv, vector<uint>&fv, vector<int>&vmat, vector<int>&fmat,vector<double> &colordegreeT){

    if(csi==-1){
        vv = vertices;
        fv = faces2vertices;
        fmat = hiddenF2Mat;
        vmat = hiddenV2Mat;
        colordegreeT = hidden_colordegreeT;
        return;

    }
    vv.clear();fv.clear();vmat.clear();fmat.clear();
    colordegreeT = hidden_colordegreeT;

    vector<int>pickV(n_vertices,-1);
    for(int i=0;i<n_faces;++i)if(hidden_F2CSs[i] == csi){
        auto p_fv = fv_begin(i);
        for(int j=0;j<3;++j){
            fv.push_back(p_fv[j]);
            pickV[p_fv[j]] = 0;
        }
        fmat.push_back(hiddenF2Mat[i]);
    }
    int ind = 0;
    for(int i=0;i<n_vertices;++i)if(pickV[i] == 0){
        pickV[i] = ind++;
        auto p_v = v_begin(i);
        for(int j=0;j<3;++j)vv.push_back(p_v[j]);
        vmat.push_back(hiddenV2Mat[i]);
    }
    for(auto &a:fv)a = pickV[a];


}

void Surface::exportSurfaceCell(int celli, vector<double>&vv, vector<uint>&fv, vector<int>&fmat, vector<int> &ctrE, vector<double> &colordegreeT){
    //    if(celli==-1){
    //        vv = vertices;
    //        fv = faces2vertices;
    //        fmat = hiddenF2Mat;
    //        colordegreeT = hidden_colordegreeT;
    //        return;

    //    }
    vv.clear();fv.clear();fmat.clear();ctrE.clear();
    colordegreeT = hidden_colordegreeT;

    vector<int>pickV(n_vertices,-1);
    for(int i=0;i<n_faces;++i)if(hidden_F2Cells[i*2] == celli || hidden_F2Cells[i*2 +1] == celli ){
        auto p_fv = fv_begin(i);
        for(int j=0;j<3;++j){
            fv.push_back(p_fv[j]);
            pickV[p_fv[j]] = 0;
        }
        fmat.push_back(hiddenF2Mat[i]);
    }
    int ind = 0;
    for(int i=0;i<n_vertices;++i)if(pickV[i] == 0){
        pickV[i] = ind++;
        auto p_v = v_begin(i);
        for(int j=0;j<3;++j)vv.push_back(p_v[j]);
        //vmat.push_back(hiddenV2Mat[i]);
    }
    for(auto &a:fv)a = pickV[a];

    for(int i=0,endi = hiddencontainer_ctrE.size()/2;i<endi;++i){
        int a1 = pickV[hiddencontainer_ctrE[i*2]];
        int a2 = pickV[hiddencontainer_ctrE[i*2+1]];
        if(a1>=0 && a2>=0){
            ctrE.push_back(a1);
            ctrE.push_back(a2);
        }

    }



}

void Surface::importSurface(vector<double>&vv, vector<uint>&fv, vector<int>&verticesMat, vector<int>&facesMat,vector<double>&colordegreeT){

    reset();
    vertices = vv;
    faces2vertices = fv;
    hiddenF2Mat = facesMat;
    hiddenV2Mat = verticesMat;


    BuildUp(false);

    int nMat = *max_element(hiddenF2Mat.begin(),hiddenF2Mat.end())+1;
    hidden_colordegreeT = colordegreeT;

    hidden_F2CSs.clear();
    hidden_F2CSs.resize(n_faces,0);


    colordegree = 50;
    reinitflag = true;
    colormethod = 2;


    vector<double>inver_faces_colordegree(n_faces);
    for(int j=0;j<inver_faces_colordegree.size();++j)inver_faces_colordegree[j] = hidden_colordegreeT[hiddenF2Mat[j]];
    GetPerFaceColorDegree(inver_faces_colordegree);

    vector<double>inver_vertices_colordegree(n_vertices);
    for(int j=0;j<inver_vertices_colordegree.size();++j)inver_vertices_colordegree[j] = hidden_colordegreeT[hiddenV2Mat[j]];
    GetPerVertexColorDegree(inver_vertices_colordegree);

    hiddenCtr(hiddenV2Mat);
    //hiddenCtrSmoothing();



}

int Surface::importSurface_specialcontainer(int codeN,string inFpre, string insurend, string inf2Cells, string infaceMat,string inCtrE){

    codeN = 200000;
    readfile(inFpre + to_string(codeN) + insurend );
    ReScale_uniform(1.0);

    BuildUp(false);
    //Fairing(true,false,false,20,0.5,0.1);
    hiddenF2Mat.clear();hiddenV2Mat.clear();hidden_F2CSs.clear();hidden_F2Cells.clear();hiddencontainer_ctrE.clear();
    //vector<int>facesMat,verticesMat,faces2Cs,faces2Cells,Ctredges;
    readVecFile(inFpre + to_string(codeN) + infaceMat,hiddenF2Mat);
    //readVecFile(inFpre + to_string(codeN) + inVMat,hiddenV2Mat);
    //readVecFile(inFpre + to_string(codeN) + inf2Cs,hidden_F2CSs);
    readVecFile(inFpre + to_string(codeN) + inf2Cells,hidden_F2Cells);
    if(!inCtrE.empty())readContourEdgeTxtFile(inFpre + to_string(codeN) + inCtrE,hiddencontainer_ctrE);

    int nMat = *max_element(hiddenF2Mat.begin(),hiddenF2Mat.end())+1;
    hidden_colordegreeT.resize(nMat);
    for(int i=0;i<nMat;++i){
        hidden_colordegreeT[i] = 340*double(i)/double(nMat-1);
    }
    hidden_colordegreeT[0] = -1;

    return *max_element(hidden_F2Cells.begin(),hidden_F2Cells.end())+1;


}

int Surface::importSurface(string inFpre, string insurend, string inf2Cs, string inf2Cells, string infaceMat, string inVMat, string inCtrE, bool isrescale){



    int codeN = 100000;
    readfile(inFpre + to_string(codeN) + insurend );
    if(isrescale)ReScale_uniform(1.0);

    BuildUp(false);
    //Fairing(true,false,false,20,0.5,0.1);
    hiddenF2Mat.clear();hiddenV2Mat.clear();hidden_F2CSs.clear();hidden_F2Cells.clear();
    vector<int>facesMat,verticesMat,faces2Cs,faces2Cells,Ctredges;
    readVecFile(inFpre + to_string(codeN) + infaceMat,hiddenF2Mat);
    readVecFile(inFpre + to_string(codeN) + inVMat,hiddenV2Mat);
    readVecFile(inFpre + to_string(codeN) + inf2Cs,hidden_F2CSs);
    readVecFile(inFpre + to_string(codeN) + inf2Cells,hidden_F2Cells);
    if(!inCtrE.empty())readContourEdgeTxtFile(inFpre + to_string(codeN) + inCtrE,Ctredges);


    int nMat = *max_element(hiddenF2Mat.begin(),hiddenF2Mat.end())+1;
    hidden_colordegreeT.resize(nMat);
    for(int i=0;i<nMat;++i){
        hidden_colordegreeT[i] = 340*double(i)/double(nMat-1);
    }

    int nCs = *max_element(hidden_F2CSs.begin(),hidden_F2CSs.end())+1;
    vector<double>colordegreeCs(nCs);
    for(int i=0;i<nCs;++i){
        colordegreeCs[i] = 360*double(i)/double(nCs+1)+90;
    }


    colordegree = 50;
    reinitflag = true;
    colormethod = 2;


    vector<double>inver_faces_colordegree(n_faces);
    if(1)for(int j=0;j<inver_faces_colordegree.size();++j)inver_faces_colordegree[j] = hidden_colordegreeT[hiddenF2Mat[j]];
    else for(int j=0;j<inver_faces_colordegree.size();++j)inver_faces_colordegree[j] = colordegreeCs[hidden_F2CSs[j]];
    GetPerFaceColorDegree(inver_faces_colordegree);

    vector<double>inver_vertices_colordegree(n_vertices);
    for(int j=0;j<inver_vertices_colordegree.size();++j)inver_vertices_colordegree[j] = hidden_colordegreeT[hiddenV2Mat[j]];
    GetPerVertexColorDegree(inver_vertices_colordegree);

    hiddenCtr(hiddenV2Mat);
    //hiddenCtrSmoothing();



    vector<int>fakefacesMat;
    for(auto a:hiddenF2Mat){
        fakefacesMat.push_back(a);
        fakefacesMat.push_back(nMat);
    }


    writeSufFile(inFpre+string("space_arrangement"),vertices,faces2vertices,fakefacesMat,Ctredges);

    return nCs;

}

void Surface::importSurface(vector<double>&vv, vector<uint>&fv, vector<int>&vmat, bool isReOrientate, bool isRescale){

    vertices = vv;
    faces2vertices = fv;
    if(isReOrientate){
        BuildUp(true);
        ReOrientFaces();
    }
    BuildUp(false);
    if(isRescale)ReScale_uniform(1.0);



    hidden_F2Cells.resize(n_faces*2,0);
    for(int i=0;i<n_faces;++i)hidden_F2Cells[i*2+1] = 1;
    hidden_F2CSs.resize(n_faces,0);

    int nMat = *max_element(vmat.begin(),vmat.end())+1;
    vector<double>colordegreeT(nMat);
    for(int i=0;i<nMat;++i){
        colordegreeT[i] = 360*double(i)/double(nMat+1)+90;
    }

    int nCs = 1;
    vector<double>colordegreeCs(nCs);
    for(int i=0;i<nCs;++i){
        colordegreeCs[i] = 360*double(i)/double(nCs+1)+90;
    }


    colordegree = 50;
    reinitflag = true;


    vector<double>inver_faces_colordegree(n_faces,colordegree);
    //if(1)for(int j=0;j<inver_faces_colordegree.size();++j)inver_faces_colordegree[j] = colordegreeT[hidden_F2CSs[j]];
    //else for(int j=0;j<inver_faces_colordegree.size();++j)inver_faces_colordegree[j] = colordegreeCs[hidden_F2CSs[j]];
    GetPerFaceColorDegree(inver_faces_colordegree);

    vector<double>inver_vertices_colordegree(n_vertices);
    for(int j=0;j<inver_vertices_colordegree.size();++j)inver_vertices_colordegree[j] = colordegreeT[vmat[j]];
    GetPerVertexColorDegree(inver_vertices_colordegree);

    //hiddenCtr(vmat);
    face_Material1.clear();
    face_Material2.clear();

    face_Material1.resize(n_faces,0);face_Material2.resize(n_faces,0);


}


void Surface::importSurface(vector<double>&vv, vector<uint>&fv, vector<int>&vmat, vector<int> &f2Cell, vector<int> &f2CS){

    vertices = vv;
    faces2vertices = fv;

    setparameters();
    BuildUp(false);

    hidden_F2Cells = f2Cell;
    hidden_F2CSs = f2CS;

    hiddenCtr(vmat);

}

void Surface::CutSurfaceByPlanePointDistance(double *para,vector<double>&outV2planeDist){

    outV2planeDist.resize(n_vertices);
    for(int i=0;i<n_vertices;++i){
        auto p_v = v_begin(i);
        double &a = outV2planeDist[i];
        for(int j=0;j<3;++j)a+=p_v[j]*para[j];
        a+=para[3];
    }

}
void Surface::CutSurfaceByPlane(double *para,vector<double>&outCtrV,vector<uint>&outCtrE,vector<int>&outCtrEMat){

    vector<bool>vonPlaneSide(n_vertices);

    vector<int>activeE(n_edges,-1);
    //vector<bool>activeF(n_faces,false);
    outCtrV.clear();outCtrE.clear();outCtrEMat.clear();
    for(int i=0;i<n_vertices;++i){
        auto p_v = v_begin(i);
        double a = 0;
        for(int j=0;j<3;++j)a+=p_v[j]*para[j];
        a+=para[3];
        vonPlaneSide[i] = a>0;

    }

    int ctrvInd = 0;
    double intersetPoint[3];
    auto intersetLinewithPlane = [](double *pv1,double *pv2,double *plane_para,double *intsetV){
        double dirc[3],tmp[3];
        minusVec(pv2,pv1,dirc);
        normalize(dirc);
        double s = -(dot(pv1,plane_para)+plane_para[3])/dot(dirc,plane_para);
        product(s,dirc,intsetV);
        add(intsetV,pv1,intsetV);
    };
    for(int i= 0;i<n_edges;++i){
        auto p_ev =ev_begin(i);
        if(vonPlaneSide[p_ev[0]]!=vonPlaneSide[p_ev[1]]){
            activeE[i] = ctrvInd++;
            intersetLinewithPlane(v_begin(p_ev[0]),v_begin(p_ev[1]),para,intersetPoint);
            for(int j = 0;j<3;++j)outCtrV.push_back(intersetPoint[j]);
        }

    }

    auto p_vdata = outCtrV.data();
    for(int i=0;i<n_faces;++i){
        auto p_fe = fe_begin(i);
        vector<int>pickEss;
        for(int j=0;j<3;++j)if(activeE[p_fe[j]]!=-1)pickEss.push_back(activeE[p_fe[j]]);

        if(pickEss.size()==0)continue;
        assert(pickEss.size()==2);

        for(auto e:pickEss){
            outCtrE.push_back(e);
        }

        auto p_v1 = p_vdata+pickEss[0]*3;
        auto p_v2 = p_vdata+pickEss[1]*3;

        double vec[3],cr[3];
        minusVec(p_v2,p_v1,vec);
        cross(para,vec,cr);
        if(dot(cr,fnor_begin(i))>0){
            outCtrEMat.push_back(face_Material1[i]);
            outCtrEMat.push_back(face_Material2[i]);
        }else{
            outCtrEMat.push_back(face_Material2[i]);
            outCtrEMat.push_back(face_Material1[i]);
        }
    }

}

void Surface::CutSurfaceByPlaneTranslation(double *p_trl, double *para, vector<double>&outCtrV, vector<uint>&outCtrE, vector<int>&outCtrEMat){

    vector<bool>vonPlaneSide(n_vertices);

    vector<int>activeE(n_edges,-1);
    //vector<bool>activeF(n_faces,false);
    outCtrV.clear();outCtrE.clear();outCtrEMat.clear();

    for(int i=0;i<n_vertices;++i){
        auto p_v = v_begin(i);
        for(int j=0;j<3;++j) p_v[j]+=p_trl[j];
    }

    for(int i=0;i<n_vertices;++i){
        auto p_v = v_begin(i);
        double a = 0;
        for(int j=0;j<3;++j)a+=p_v[j]*para[j];
        a+=para[3];
        vonPlaneSide[i] = a>0;

    }

    int ctrvInd = 0;
    double intersetPoint[3];
    auto intersetLinewithPlane = [](double *pv1,double *pv2,double *plane_para,double *intsetV){
        double dirc[3],tmp[3];
        minusVec(pv2,pv1,dirc);
        normalize(dirc);
        double s = -(dot(pv1,plane_para)+plane_para[3])/dot(dirc,plane_para);
        product(s,dirc,intsetV);
        add(intsetV,pv1,intsetV);
    };
    for(int i= 0;i<n_edges;++i){
        auto p_ev =ev_begin(i);
        if(vonPlaneSide[p_ev[0]]!=vonPlaneSide[p_ev[1]]){
            activeE[i] = ctrvInd++;
            intersetLinewithPlane(v_begin(p_ev[0]),v_begin(p_ev[1]),para,intersetPoint);
            for(int j = 0;j<3;++j)outCtrV.push_back(intersetPoint[j]);
        }

    }

    auto p_vdata = outCtrV.data();
    for(int i=0;i<n_faces;++i){
        auto p_fe = fe_begin(i);
        vector<int>pickEss;
        for(int j=0;j<3;++j)if(activeE[p_fe[j]]!=-1)pickEss.push_back(activeE[p_fe[j]]);

        if(pickEss.size()==0)continue;
        assert(pickEss.size()==2);

        for(auto e:pickEss){
            outCtrE.push_back(e);
        }

        auto p_v1 = p_vdata+pickEss[0]*3;
        auto p_v2 = p_vdata+pickEss[1]*3;

        double vec[3],cr[3];
        minusVec(p_v2,p_v1,vec);
        cross(para,vec,cr);
        if(dot(cr,fnor_begin(i))>0){
            outCtrEMat.push_back(face_Material1[i]);
            outCtrEMat.push_back(face_Material2[i]);
        }else{
            outCtrEMat.push_back(face_Material2[i]);
            outCtrEMat.push_back(face_Material1[i]);
        }
    }

    for(int i=0;i<n_vertices;++i){
        auto p_v = v_begin(i);
        for(int j=0;j<3;++j) p_v[j]-=p_trl[j];
    }

}

void Surface::CutSurfaceByBatchPlane(vector<vector<double>>paras,vector<double>&out_CtrV,vector<vector<int>>&out_CtrE,vector<vector<int>>&out_CtrEMat,
                                     vector<double>&vnormals){



    out_CtrV.clear();out_CtrE.clear();out_CtrEMat.clear();vnormals.clear();
    out_CtrE.resize(paras.size());out_CtrEMat.resize(paras.size());
    auto intersetLinewithPlane = [](double *pv1,double *pv2,double *plane_para,double *intsetV){
        double dirc[3];
        minusVec(pv2,pv1,dirc);
        normalize(dirc);
        double s = -(dot(pv1,plane_para)+plane_para[3])/dot(dirc,plane_para);
        product(s,dirc,intsetV);
        add(intsetV,pv1,intsetV);
    };
    int ctrvInd = 0;
    double intersetPoint[3];
    vector<vector<int>>ActiveF2E(n_faces);
    for(int cp=0;cp<paras.size();++cp){
        auto &outCtrE = out_CtrE[cp];
        auto &outCtrEMat = out_CtrEMat[cp];
        double *para = paras[cp].data();
        vector<bool>vonPlaneSide(n_vertices);
        vector<int>activeE(n_edges,-1);
        //vector<bool>activeF(n_faces,false);


        for(int i=0;i<n_vertices;++i){
            auto p_v = v_begin(i);
            double a = 0;
            for(int j=0;j<3;++j)a+=p_v[j]*para[j];
            a+=para[3];
            vonPlaneSide[i] = a>0;

        }




        for(int i= 0;i<n_edges;++i){
            auto p_ev =ev_begin(i);
            if(vonPlaneSide[p_ev[0]]!=vonPlaneSide[p_ev[1]]){
                activeE[i] = ctrvInd++;
                intersetLinewithPlane(v_begin(p_ev[0]),v_begin(p_ev[1]),para,intersetPoint);
                for(int j = 0;j<3;++j)out_CtrV.push_back(intersetPoint[j]);

                auto p_ef = efnon_begin(i);
                add(fnor_begin(p_ef[0]),fnor_begin(p_ef[1]),intersetPoint);
                for(int j = 0;j<3;++j)vnormals.push_back(intersetPoint[j]);

            }

        }

        auto p_vdata = out_CtrV.data();
        for(int i=0;i<n_faces;++i){
            auto p_fe = fe_begin(i);
            vector<int>pickEss;
            for(int j=0;j<3;++j)if(activeE[p_fe[j]]!=-1)pickEss.push_back(activeE[p_fe[j]]);

            if(pickEss.size()==0)continue;
            assert(pickEss.size()==2);
            ActiveF2E[i].push_back(cp);
            ActiveF2E[i].push_back(outCtrE.size()/2);
            for(auto e:pickEss){
                outCtrE.push_back(e);
            }


            auto p_v1 = p_vdata+pickEss[0]*3;
            auto p_v2 = p_vdata+pickEss[1]*3;

            double vec[3],cr[3];
            minusVec(p_v2,p_v1,vec);
            cross(para,vec,cr);
            if(dot(cr,fnor_begin(i))>0){
                outCtrEMat.push_back(face_Material1[i]);
                outCtrEMat.push_back(face_Material2[i]);
            }else{
                outCtrEMat.push_back(face_Material2[i]);
                outCtrEMat.push_back(face_Material1[i]);
            }
        }


    }

    double tri2para[4],intersection[3];
    auto p_v=out_CtrV.data();
    int numIn = 0,numAc = 0;
    for(int i=0;i<n_faces;++i)if(ActiveF2E[i].size()>2){
        //assert(ActiveF2E[i].size()==4);
        ++numAc;
        auto pd = ActiveF2E[i];

        pointNor2para4(v_begin(fv_begin(i)[0]),fnor_begin(i),tri2para);

        vector<double *>pev(ActiveF2E[i].size());
        for(int j=0;j<ActiveF2E[i].size()/2;++j){
            int cp = pd[j*2],ce = 2*pd[j*2+1];
            pev[j*2] = p_v+3*out_CtrE[cp][ce];
            pev[j*2+1] = p_v+3*out_CtrE[cp][ce+1];
        }

        int ll[2];bool findone = false;
        for(int j =0;j<ActiveF2E[i].size()/2;++j)
            for(int k=j+1;k<ActiveF2E[i].size()/2;++k)
                if(planeSegIntersectionTest(pev[j*2],pev[j*2+1],pev[k*2],pev[k*2+1])){
                    assert(findone == false);
                    ll[0] = j;ll[1] = k;
                    findone = true;
                }
        if(!findone)continue;


        threeplaneIntersection(paras[pd[ll[0]*2]].data(),paras[pd[ll[1]*2]].data(),tri2para,intersection);
        for(int j = 0;j<3;++j)out_CtrV.push_back(intersection[j]);

        auto p_fn = fnor_begin(i);
        for(int j = 0;j<3;++j)vnormals.push_back(p_fn[j]);

        for(int j=0;j<2;++j){
            int cp = pd[ll[j]*2],ce = 2*pd[ll[j]*2+1];

            int &s2 = out_CtrE[cp][ce+1];
            out_CtrE[cp].push_back(ctrvInd);
            out_CtrE[cp].push_back(s2);
            out_CtrE[cp][ce+1] = ctrvInd;
            out_CtrEMat[cp].push_back(out_CtrEMat[cp][ce]);
            out_CtrEMat[cp].push_back(out_CtrEMat[cp][ce+1]);

        }

        //        if(ctrvInd==1801){
        //            for(int j=0;j<ActiveF2E[i].size()/2;++j){
        //                int cp = pd[j*2],ce = 2*pd[j*2+1];
        //                cout<<out_CtrE[cp][ce]<<' '<<out_CtrE[cp][ce+1]<<endl;
        //            }
        //            cout<<"paius"<<endl;
        //        }

        //        for(int j=0;j<ActiveF2E[i].size()/2;++j){
        //            int cp = pd[j*2],ce = 2*pd[j*2+1];

        //            int &s2 = out_CtrE[cp][ce+1];
        //            out_CtrE[cp].push_back(ctrvInd);
        //            out_CtrE[cp].push_back(s2);
        //            out_CtrE[cp][ce+1] = ctrvInd;
        //            out_CtrEMat[cp].push_back(out_CtrEMat[cp][ce]);
        //            out_CtrEMat[cp].push_back(out_CtrEMat[cp][ce+1]);

        //        }

        //        if(ctrvInd==1801){
        //            for(int j=0;j<ActiveF2E[i].size()/2;++j){
        //                int cp = pd[j*2],ce = 2*pd[j*2+1];
        //                cout<<out_CtrE[cp][ce]<<' '<<out_CtrE[cp][ce+1]<<endl;
        //            }
        //            cout<<"paius"<<endl;
        //        }
        ++ctrvInd;
        ++numIn;

    }

    cout<<"numIn: "<<numIn<<endl;
    cout<<"numAc: "<<numAc<<endl;





}

void Surface::ReOrientateConvexShape(){

    BuildUp(true);
    ReOrientFaces();
    double center[3] = {0,0,0};

    for(int i=0;i<n_vertices;++i){
        auto p_v = v_begin(i);
        for(int j=0;j<3;++j)center[j]+=p_v[j];
    }

    for(int j=0;j<3;++j)center[j]/=n_vertices;
    double f0[3] = {0,0,0};

    auto p_fv = fv_begin(0);
    for(int i=0;i<3;++i){
        auto p_v = v_begin(p_fv[i]);
        for(int j=0;j<3;++j)f0[j]+=p_v[j];
    }
    for(int j=0;j<3;++j)f0[j]/=3;

    double fcc[3];
    minusVec(f0,center,fcc);

    if(dot(fnor_begin(0),fcc)<0)InverseFaces();

    face_Material1.clear();
    face_Material2.clear();
    face_Material1.resize(n_faces,1);
    face_Material1.resize(n_faces,0);


    BuildUp(false);
}


}//n_rf
