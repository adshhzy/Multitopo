#include "my_mesh.h"
//#include <hash_map>
#include <list>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include<iostream>
#include<fstream>
#include <sstream>
#include <functional>
#include <assert.h>
#include <eigen3/Eigen/Geometry>
#include<random>


double randomdouble() {return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);}

void randomRotmatrix(vector<vector<double>>&outmat){

    double u1 = randomdouble();
    double u2 = randomdouble();
    double u3 = randomdouble();

    double theta = acos(2*u1-1);
    double phi = 2*my_PI*u2;
    double delta = 2*my_PI*u3;

    Eigen::AngleAxisd rollAngle(theta, Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd yawAngle(phi, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd pitchAngle(delta, Eigen::Vector3d::UnitX());

    Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;

    Eigen::Matrix3d rotationMatrix = q.matrix();

    outmat.clear();
    outmat.resize(3,vector<double>(4,0));
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)outmat[i][j] = rotationMatrix(i,j);



}


bool Mesh::readfile(string filename){

    string filepath;
    string modname;
    string extname;

    SplitFileName(filename,filepath,modname,extname);

    if(extname == ".obj" )return readObjfile(filename);
    else if(extname == ".off" )return readOfffile(filename);

    return false;
}

bool Mesh::createToy(int toyind){

    auto RenderCone = [this](int n){

        n_vertices = 3*n+1;
        vertices.clear();
        vertices.reserve(n_vertices*3);
        double thetastep = 2*my_PI/n;
        double r = 1;
        for(int i = 0;i<n;++i){
            vertices.push_back(0);
            vertices.push_back(0);
            vertices.push_back(1);
        }
        for(int i = 0;i<2*n;++i){
            vertices.push_back(r*cos(thetastep*i));
            vertices.push_back(r*sin(thetastep*i));
            vertices.push_back(0);
        }
        for(int i=0;i<3;++i)vertices.push_back(0);

        n_faces = 2*n;
        faces2vertices.clear();
        faces2vertices.reserve(n_faces*3);
        for(int i = 0;i<n-1;++i){
            faces2vertices.push_back(i);
            faces2vertices.push_back(n+i);
            faces2vertices.push_back(n+1+i);
        }
        faces2vertices.push_back(n-1);
        faces2vertices.push_back(2*n-1);
        faces2vertices.push_back(n);

        for(int i = 0;i<n-1;++i){
            faces2vertices.push_back(3*n);
            faces2vertices.push_back(2*n+i);
            faces2vertices.push_back(2*n+i+1);
        }
        faces2vertices.push_back(3*n);
        faces2vertices.push_back(3*n-1);
        faces2vertices.push_back(2*n);


    };

    auto RenderCyclinderShell = [this](int n){

        n_vertices = 2*n;
        vertices.clear();
        vertices.reserve(n_vertices*3);
        double thetastep = 2*my_PI/n;
        double r = 1;
        for(int i = 0;i<n;++i){
            vertices.push_back(r*cos(thetastep*i));
            vertices.push_back(r*sin(thetastep*i));
            vertices.push_back(0);
        }
        for(int i = 0;i<n;++i){
            vertices.push_back(r*cos(thetastep*i));
            vertices.push_back(r*sin(thetastep*i));
            vertices.push_back(1);
        }


        n_faces = 2*n;
        faces2vertices.clear();
        faces2vertices.reserve(n_faces*3);
        for(int i = 0;i<n-1;++i){
            faces2vertices.push_back(i+1);
            faces2vertices.push_back(n+1+i);
            faces2vertices.push_back(n+i);
            faces2vertices.push_back(i);
            faces2vertices.push_back(i+1);
            faces2vertices.push_back(n+i);
        }
        faces2vertices.push_back(0);
        faces2vertices.push_back(n);
        faces2vertices.push_back(2*n-1);
        faces2vertices.push_back(n-1);
        faces2vertices.push_back(0);
        faces2vertices.push_back(2*n-1);



    };


    auto RenderSphere = [this](){
        double ppp[] = {0,0,-1,
                        0.7236,-0.52572,-0.447215,
                        -0.276385,-0.85064,-0.447215,
                        -0.894425,0,-0.447215,
                        -0.276385,0.85064,-0.447215,
                        0.7236,0.52572,-0.447215,
                        0.276385,-0.85064,0.447215,
                        -0.7236,-0.52572,0.447215,
                        -0.7236,0.52572,0.447215,
                        0.276385,0.85064,0.447215,
                        0.894425,0,0.447215,
                        0,0,1,
                        0.425323,-0.309011,-0.850654,
                        -0.162456,-0.499995,-0.850654,
                        0.262869,-0.809012,-0.525738,
                        0.425323,0.309011,-0.850654,
                        0.850648,0,-0.525736,
                        -0.52573,0,-0.850652,
                        -0.688189,-0.499997,-0.525736,
                        -0.162456,0.499995,-0.850654,
                        -0.688189,0.499997,-0.525736,
                        0.262869,0.809012,-0.525738,
                        0.951058,0.309013,0,
                        0.951058,-0.309013,0,
                        0.587786,-0.809017,0,
                        0,-1,0,
                        -0.587786,-0.809017,0,
                        -0.951058,-0.309013,0,
                        -0.951058,0.309013,0,
                        -0.587786,0.809017,0,
                        0,1,0,
                        0.587786,0.809017,0,
                        0.688189,-0.499997,0.525736,
                        -0.262869,-0.809012,0.525738,
                        -0.850648,0,0.525736,
                        -0.262869,0.809012,0.525738,
                        0.688189,0.499997,0.525736,
                        0.52573,0,0.850652,
                        0.162456,-0.499995,0.850654,
                        -0.425323,-0.309011,0.850654,
                        -0.425323,0.309011,0.850654,
                        0.162456,0.499995,0.850654};
        uint fff[] = {1,14,12,
                      13,12,14,
                      14,2,13,
                      12,13,0,
                      12,16,1,
                      16,12,15,
                      15,5,16,
                      15,12,0,
                      2,18,13,
                      17,13,18,
                      18,3,17,
                      13,17,0,
                      3,20,17,
                      19,17,20,
                      20,4,19,
                      17,19,0,
                      4,21,19,
                      15,19,21,
                      21,5,15,
                      19,15,0,
                      16,23,1,
                      23,16,22,
                      22,10,23,
                      5,22,16,
                      14,25,2,
                      25,14,24,
                      24,6,25,
                      1,24,14,
                      18,27,3,
                      27,18,26,
                      26,7,27,
                      2,26,18,
                      20,29,4,
                      29,20,28,
                      28,8,29,
                      3,28,20,
                      21,31,5,
                      31,21,30,
                      30,9,31,
                      4,30,21,
                      10,32,23,
                      24,23,32,
                      32,6,24,
                      23,24,1,
                      6,33,25,
                      26,25,33,
                      33,7,26,
                      25,26,2,
                      7,34,27,
                      28,27,34,
                      34,8,28,
                      27,28,3,
                      8,35,29,
                      30,29,35,
                      35,9,30,
                      29,30,4,
                      9,36,31,
                      22,31,36,
                      36,10,22,
                      31,22,5,
                      32,38,6,
                      38,32,37,
                      37,11,38,
                      10,37,32,
                      33,39,7,
                      39,33,38,
                      38,11,39,
                      6,38,33,
                      34,40,8,
                      40,34,39,
                      39,11,40,
                      7,39,34,
                      35,41,9,
                      41,35,40,
                      40,11,41,
                      8,40,35,
                      36,37,10,
                      37,36,41,
                      41,11,37,
                      9,41,36};
        n_vertices = 42;
        vertices.resize(n_vertices*3);
        for(int i=0;i<vertices.size();++i)vertices[i] = ppp[i];
        n_faces = 80;
        faces2vertices.resize(n_faces*3);
        for(int i=0;i<faces2vertices.size();++i)faces2vertices[i] = fff[i];

    };

    cout<<"create"<<endl;
    switch(toyind){
    case 1:RenderCone(10);break;
    case 2:RenderCyclinderShell(8);break;
    case 3:RenderSphere();break;
    default:cout<<"wrong toy index!"<<endl;break;



    };
    cout<<"create done"<<endl;
    setparameters();
    cout<<"create done2"<<endl;
    return true;


}

void Mesh::BuildUp(bool ismanifold_in){
    ismanifold = ismanifold_in;

    setparameters();
    if(ismanifold){
        BuildNeighborTable();
        ComputeEdgefromFace();
    }else BuildNeighborTable_nonmanifold();

    ComputeFaceNormal();
    GetRescaleInfo(lscale, pcenter);
}


void Mesh::ReScale_uniform(double lscale,bool isMoveCent){
    setparameters();
    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    double centers[3] = {0,0,0};
    for (int i = 0; i<n_vertices; i++){
        auto point = v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);
        centers[0]+=point[0];centers[1]+=point[1];centers[2]+=point[2];
        //if(centers[0]!=centers[0]){cout<<centers[0]<<"  "<<point[0]<<" "<<i/dim<<endl;}
    }
    centers[0]/=n_vertices;centers[1]/=n_vertices;centers[2]/=n_vertices;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2 / lscale;
    //cout<<centers[0]<<"  "<<centers[1]<<"  "<<centers[2]<<"   "<<largestdis<<endl;

    largestdis = 1/largestdis;
    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-centers[i])*largestdis;
    if(!isMoveCent){
        for(int i=0;i<3;++i)centers[i]*=largestdis;
        for (int j = 0; j<n_vertices*3; j+=3)
            for(int i = 0; i<3; i++)
                vertices[j+i]=(vertices[j+i]+centers[i]);
    }

}
void Mesh::GetRescaleInfo(double &lscale, double *pcenter){

    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    double centers[3] = {0,0,0};
    for (int i = 0; i<n_vertices; i++){
        auto point = v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);
        centers[0]+=point[0];centers[1]+=point[1];centers[2]+=point[2];
        //if(centers[0]!=centers[0]){cout<<centers[0]<<"  "<<point[0]<<" "<<i/dim<<endl;}
    }
    centers[0]/=n_vertices;centers[1]/=n_vertices;centers[2]/=n_vertices;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    //cout<<centers[0]<<"  "<<centers[1]<<"  "<<centers[2]<<"   "<<largestdis<<endl;

    lscale = 1/largestdis;
    for(int i=0;i<3;++i)pcenter[i] = centers[i];


}
void Mesh::ReScale(double lscale,double *pcenter){

    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-pcenter[i])*lscale;


}

void Mesh::ReScale_uniform_directratio(double ratio){
    setparameters();
    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    double centers[3] = {0,0,0};
    for (int i = 0; i<n_vertices; i++){
        auto point = v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);
        centers[0]+=point[0];centers[1]+=point[1];centers[2]+=point[2];
        //if(centers[0]!=centers[0]){cout<<centers[0]<<"  "<<point[0]<<" "<<i/dim<<endl;}
    }
    centers[0]/=n_vertices;centers[1]/=n_vertices;centers[2]/=n_vertices;

    //cout<<centers[0]<<"  "<<centers[1]<<"  "<<centers[2]<<"   "<<largestdis<<endl;


    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-centers[i])*ratio;

}

void Mesh::ReScale(double xscale,double yscale,double zscale){
    setparameters();
    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    double centers[3] = {0,0,0};
    for (int i = 0; i<n_vertices; i++){
        auto point = v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);
        centers[0]+=point[0];centers[1]+=point[1];centers[2]+=point[2];
        //if(centers[0]!=centers[0]){cout<<centers[0]<<"  "<<point[0]<<" "<<i/dim<<endl;}
    }
    centers[0]/=n_vertices;centers[1]/=n_vertices;centers[2]/=n_vertices;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    double ratiodis[3];
    ratiodis[0] = xscale/largestdis;ratiodis[1] = yscale/largestdis;ratiodis[2] = zscale/largestdis;
    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-centers[i])*ratiodis[i];

}

void Mesh::MoveBottom(bool xx,bool yy,bool zz){
    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    double centers[3] = {0,0,0};
    for (int i = 0; i<n_vertices; i++){
        auto point = v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);
        centers[0]+=point[0];centers[1]+=point[1];centers[2]+=point[2];
        //if(centers[0]!=centers[0]){cout<<centers[0]<<"  "<<point[0]<<" "<<i/dim<<endl;}
    }

    if(xx)for (int i = 0; i<n_vertices; i++)v_begin(i)[0] -=xmin;
    if(yy)for (int i = 0; i<n_vertices; i++)v_begin(i)[1] -=ymin;
    if(zz)for (int i = 0; i<n_vertices; i++)v_begin(i)[2] -=zmin;


}


void Mesh::CheckValidation(){


    for(auto a:faces2vertices)if(a>=n_vertices){
        cout<<"Wrong Surfaces: "<<n_vertices<<' '<<a<<endl;
        assert(a<n_vertices);
    }




}


void Mesh::BuildNeighborTable(){

    
    setparameters();
    CheckValidation();
    //cout<<"build: "<<n_faces<<' '<<faces.size()<<endl;
    int facedim = dimtable[2];

    vertices2faces.clear();
    vertices2faces_accumulate.clear();
    vertices2faces_accumulate.resize(n_vertices,0);
    for (uint i = 0; i < n_faces; i++){
        auto p_f = &(faces2vertices[facedim * i]);
        for (uint j = 0; j < facedim; ++j){
            vertices2faces_accumulate[p_f[j]]++;
        }
    }
    for (uint i = 1; i < n_vertices; i++){
        vertices2faces_accumulate[i]+=vertices2faces_accumulate[i-1];
    }
    vertices2faces_accumulate.insert(vertices2faces_accumulate.begin(),0);
    vertices2faces.resize(vertices2faces_accumulate[n_vertices]);

    vector<uchar>cur_num(n_vertices,0);
    
    for (uint i = 0; i < n_faces; i++){
        auto p_f = &(faces2vertices[facedim * i]);
        for (uint j = 0; j < facedim; ++j){
            auto v = p_f[j];
            vertices2faces[vertices2faces_accumulate[v]+cur_num[v]] = i;
            cur_num[v]++;
        }
    }

    vertices2faces_inverse.resize(vertices2faces.size());
    for(int i=0;i<n_vertices;++i){
        auto p_vfinv = vfinv_begin(i);
        auto p_vfend = vf_end(i);
        for(auto p_vf = vf_begin(i);p_vf!=p_vfend;++p_vf,++p_vfinv){
            auto p_fv = fv_begin(*p_vf);
            for(int j=0;j<3;++j)if(p_fv[j]==i){
                (*p_vfinv)=j;break;
            }
        }
    }

    
    faces2faces.clear();
    faces2faces.resize(n_faces * f2f_step, UINT_MAX);
    for(int i = 0;i<faces2faces.size();i+=f2f_step)faces2faces[i]=0;
    faces_in_manifold.clear();
    faces_in_manifold.resize(n_faces, false);
    vector<bool>neighborface(n_faces, false);
    vector<uint>used_indices;
    
    for (uint i = 0; i < n_faces; i++){
        used_indices.clear();
        auto p_f2f = &(faces2faces[f2f_step * i]);
        auto p_f = &(faces2vertices[facedim * i]);
        for (uint j = 0; j < facedim; ++j){
            for(auto p_vf = vf_begin(p_f[j]);p_vf!=vf_end(p_f[j]);++p_vf)
            {
                if (*p_vf != i){
                    if (neighborface[*p_vf]){
                        if (p_f2f[0] == facedim){ cout << "un-manifold!" << endl;faces_in_manifold[i] = true; break; }
                        else{ ++p_f2f[0]; p_f2f[p_f2f[0]] = *p_vf; }
                    }
                    else { neighborface[*p_vf] = true; used_indices.push_back(*p_vf); }
                }
            }
        }
        for (auto used : used_indices)neighborface[used] = false;
    }
    




    faces_flinge.clear();
    faces_flinge.resize(n_faces, false);
    uint flinge_face = 0, un_manifold = 0;
    for (uint i = 0; i < n_faces; i++)if (faces2faces[i *f2f_step] != facedim){ flinge_face++; faces_flinge[i] = true; }
    for (auto in : faces_in_manifold)if (in)un_manifold++;

    //    vector<uint> vertices2facesTemp(vertices2faces.size());
    //    for(int i=0;i<n_vertices;++i){
    //        auto p_vfb = vf_begin(i);
    //        auto p_vfe = vf_end(i);
    //        auto p_vdbtemp =  &(vertices2facesTemp[vertices2faces_accumulate[v_ind]]);
    //  }





    cout <<"total faces: "<<n_faces<< " un_manifold: " << un_manifold << " flinge_face: " << flinge_face << endl;
    isbuildneighbortable = true;
}



void Mesh::BuildNeighborTable_nonmanifold(){

    if(n_vertices==0)return;
    if(n_faces==0)return;


    CheckValidation();
    vertices2faces_accumulate.resize(n_vertices+1,0);
    vector<uint>vertices2facesnum(n_vertices,0);

    for(uint i=0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        for(int j=0;j<3;++j){
            vertices2facesnum[p_fv[j]]++;
        }
    }

    vertices2faces_accumulate[0]=0;
    for(uint i=0;i<n_vertices;++i){
        vertices2faces_accumulate[i+1] = vertices2faces_accumulate[i]+vertices2facesnum[i];
    }

    vector<uint>curind(n_vertices,0);
    vertices2faces.resize(vertices2faces_accumulate[n_vertices]);
    vertices2faces_inverse.resize(vertices2faces_accumulate[n_vertices]);
    for(uint i=0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        for(int j=0;j<3;++j){
            auto ind = p_fv[j];
            vf_begin(ind)[curind[ind]] = i;
            vfinv_begin(ind)[curind[ind]] = j;
            ++curind[ind];
        }
    }


    vector<bool>ispick(n_vertices,false);
    vector<uint>curbuffer;
    vertices2vertices.clear();
    vertices2vertices.reserve(n_vertices*5);
    vertices2edges_accumulate.resize(n_vertices+1);
    vertices2edges_accumulate[0]=0;
    for(uint i=0;i<n_vertices;++i){
        curbuffer.clear();
        ispick[i]=true;
        auto vfnum = vf_num(i);
        auto p_vf=vf_begin(i);
        for(int j=0;j<vfnum;++j){
            auto find = p_vf[j];
            auto p_fv = fv_begin(find);
            for(int k=0;k<3;++k){
                if(!ispick[p_fv[k]]){
                    curbuffer.push_back(p_fv[k]);
                    ispick[p_fv[k]] = true;
                }
            }
        }

        vertices2edges_accumulate[i+1] = vertices2edges_accumulate[i]+curbuffer.size();
        for(auto a:curbuffer)vertices2vertices.push_back(a);
        for(auto a:curbuffer)ispick[a]=false;
        ispick[i]=false;

    }

    vertices2vertices_inverse.resize(vertices2edges_accumulate[n_vertices]);

    for(int i=0;i<n_vertices;++i){
        auto p_vvinv = vvinv_begin(i);
        auto p_vvend = vv_end(i);
        for(auto p_vv = vv_begin(i);p_vv!=p_vvend;++p_vv,++p_vvinv){
            auto p_vvi = vv_begin(*p_vv);
            auto num = vv_num(*p_vv);
            for(int j=0;j<num;++j)if(p_vvi[j]==i){
                (*p_vvinv)=j;break;
            }
        }
    }

    vertices2edges.resize(vertices2edges_accumulate[n_vertices]);
    vertices2edgesPN.resize(vertices2edges_accumulate[n_vertices]);
    for(auto &a:vertices2edges)a=UINTFLAG;
    edges.clear();
    edges.reserve(n_vertices*2);
    uint eind = 0;
    for(int i=0;i<n_vertices;++i){
        auto p_vvinv = vvinv_begin(i);
        auto p_vvend = vv_end(i);
        auto p_ve = ve_begin(i);
        auto p_vepn = vepn_begin(i);
        for(auto p_vv = vv_begin(i);p_vv!=p_vvend;++p_vv,++p_vvinv,++p_ve,++p_vepn){
            if((*p_ve)!=UINTFLAG){continue;}
            auto p_vei = ve_begin(*p_vv);
            auto p_vepni = vepn_begin(*p_vv);
            p_vei[*p_vvinv] = (*p_ve) = eind;
            ++eind;
            edges.push_back(i);edges.push_back(*p_vv);
            (*p_vepn) = true;(*p_vepni) = false;
        }
    }
    n_edges = eind;

    edge_non_manifold.clear();
    edge_non_manifold.resize(n_edges,false);
    efnum_nonmanifold.clear();
    efnum_nonmanifold.resize(n_edges,0);
    vertices_non_manifold.clear();
    vertices_non_manifold.resize(n_vertices,false);

    vector<vector<uint> >edges2faces_vector(n_edges);
    edges2faces_nonmanifold.clear();
    edges2faces_accumulate_nonmanifold.clear();
    edges2faces_accumulate_nonmanifold.resize(n_edges+1,0);
    faces.resize(n_faces*3);
    vector<int>curfacesE(n_faces,0);
    uint num_none = 0;

    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        for(auto p_vf = vf_begin(p_ev[0]);p_vf != vf_end(p_ev[0]);++p_vf){
            auto p_fv = fv_begin(*p_vf);
            for(int k=0;k<3;++k){
                if(p_fv[k]==p_ev[1]){
                    ++efnum_nonmanifold[i];
                    edges2faces_vector[i].push_back(*p_vf);
                    faces[3*(*p_vf)+curfacesE[*p_vf]]=i;++curfacesE[*p_vf];
                    break;
                }
            }
        }
        if(efnum_nonmanifold[i]>2){
            //cout<<efnum_nonmanifold[i]<<' ';
            edge_non_manifold[i]=true;num_none++;
            vertices_non_manifold[p_ev[0]] = vertices_non_manifold[p_ev[1]] = true;
        }
    }
    for(int i=0;i<n_edges;++i)edges2faces_accumulate_nonmanifold[i+1] = edges2faces_accumulate_nonmanifold[i]+efnum_nonmanifold[i];
    edges2faces_nonmanifold.reserve(edges2faces_accumulate_nonmanifold[n_edges]);
    for(auto &a:edges2faces_vector)for(auto &b:a)edges2faces_nonmanifold.push_back(b);
    //cout<<endl;
    uint num_nonv = 0;
    for(auto a:vertices_non_manifold)if(a)++num_nonv;
    //cout<<"number of edges: "<<n_edges<<endl;
    //cout<<"number of non-manifold edges: "<<num_none<<endl;
    //cout<<"number of non-manifold vertices: "<<num_nonv<<endl;
    isbuildneighbortable = true;



    vertices2edges_accumulate_nonmanifold.clear();
    vertices2edges_accumulate_nonmanifold.resize(n_vertices+1,0);

    for(int i=0;i<n_vertices;++i){
        auto p_ve = ve_begin(i);
        auto venum = vv_num(i);
        if(vertices_non_manifold[i]){
            for(int j=0;j<venum;++j){
                if(edge_non_manifold[p_ve[j]])++(vertices2edges_accumulate_nonmanifold[i+1]);
            }
        }else{
            vertices2edges_accumulate_nonmanifold[i+1] = venum;
        }

    }
    for(uint i=0;i<n_vertices;++i)vertices2edges_accumulate_nonmanifold[i+1] += vertices2edges_accumulate_nonmanifold[i];
    vertices2vertices_nonmanifold.resize(vertices2edges_accumulate_nonmanifold[n_vertices]);
    for(int i=0;i<n_vertices;++i){
        auto p_ve = ve_begin(i);
        auto p_vv = vv_begin(i);
        auto venum = ve_num(i);
        auto p_vvnon = vvnon_begin(i);

        if(vertices_non_manifold[i]){
            int numnon=0;
            for(int j=0;j<venum;++j){
                if(edge_non_manifold[p_ve[j]]){
                    p_vvnon[numnon]=p_vv[j];
                    ++numnon;
                }
            }
            //cout<<numnon<<' ';
        }else{
            for(int j=0;j<venum;++j)p_vvnon[j]=p_vv[j];

        }
    }


    //    vector< uint >faces2faces_nonmanifold;
    //    vector<  unsigned long >faces2faces_accumulate_nonmanifold;

    //    vector< uint >edges2faces_nonmanifold;
    //    vector<  unsigned long >edges2faces_accumulate_nonmanifold;

    faces2faces_nonmanifold.clear();
    faces2faces_accumulate_nonmanifold.clear();
    faces2faces_accumulate_nonmanifold.resize(n_faces+1,0);
    vector<int> faces2facesnum(n_faces,0);


    for(int i=0;i<n_faces;++i){
        for(auto p_fe = fe_begin(i);p_fe!=fe_end(i);++p_fe){
            faces2facesnum[i]+=efnon_num(*p_fe)-1;
        }
    }
    for(int i=0;i<n_faces;++i)faces2faces_accumulate_nonmanifold[i+1]=faces2faces_accumulate_nonmanifold[i]+ faces2facesnum[i];
    faces2facesnum.clear();
    faces2facesnum.resize(n_faces,0);

    faces2faces_nonmanifold.resize(faces2faces_accumulate_nonmanifold[n_faces]);

    for(int i=0;i<n_faces;++i){
        auto p_ff = ffnon_begin(i);
        int ind = 0;
        for(auto p_fe = fe_begin(i);p_fe!=fe_end(i);++p_fe){
            uint* p_ef = efnon_begin(*p_fe);
            auto efnum = efnon_num(*p_fe);
            for(int j=0;j<efnum;++j){
                if(p_ef[j]!=i)p_ff[ind++] = p_ef[j];
            }
        }
        if( ind!=ffnon_num(i) ){cout<<"err nonmanifold neighborhood"<<endl;exit(-1);}
    }


}


void Mesh::ComputeDefectAngle(){
    if (!isbuildEdges){ ComputeEdgefromFace(false); }
    vertices_defect.resize(n_vertices,-1.0);

    double e1[3],e2[3],*p_e;
    bool ise2 = false;
    for(uint i = 0;i<n_vertices;++i){
        if(isflingevertives(i))continue;
        vertices_defect[i]=0;
        auto p_v = v_begin(i);
        //int a =0;
        for(auto p_vf = vf_begin(i);p_vf!=vf_end(i);++p_vf){
            p_e = e1;
            for(auto p_fv = fv_begin(*p_vf);p_fv!=fv_end(*p_vf);++p_fv){
                if((*p_fv)==i)continue;
                minusVec(v_begin(*p_fv),p_v,p_e);
                p_e = e2;
            }
            normalize(e1);normalize(e2);
            vertices_defect[i]+=angleNor(e1,e2);
            //++a;
        }
        //cout<<a<<' ';
    }
    //cout<<endl;
    //for(uint i = 0;i<n_vertices;++i)vertices_defect[i] = 2*PI - vertices_defect[i];


}
void Mesh::DeleteInmanifoldedFaces(){
    cout<<"not implemented!"<<endl;return;
    if (!isbuildneighbortable){ BuildNeighborTable(); }

    //vector<face2>
    for (int i = 0; i < n_faces; ++i){




    }


}

inline bool is_sequenced(uint first, uint second){
    if (first < second){
        if (first == 1 && second == 2)return true;
        else if (first == 2 && second == 3)return true;
        else if (first == 1 && second == 3)return false;
    }
    else {
        if(first == 2 && second == 1)return false;
        else if (first == 3 && second == 2)return false;
        else if (first == 3 && second == 1)return true;
    }
    return true;
}
inline void swap_index(uint* first, uint* second){
    auto v_first = *first;
    *first = *second;
    *second = v_first;
}
void Mesh::ReOrientFaces(){
    if(n_faces == 0)return;
    if (!isbuildneighbortable){ BuildNeighborTable(); }


    //vector<bool>isset = faces_in_manifold;
    vector<bool>isset(n_faces,false);
    uint seed = 0;

    vector<int>used_indices(n_vertices, 0);
    list<uint>waitinglist;
    waitinglist.push_back(seed);
    vector< pair<int, int> >order;
    int set_count = 0;
    int ccc = 0;
    while(true){

        while (!waitinglist.empty()){
            const int findex = waitinglist.front();
            waitinglist.pop_front();
            if(isset[findex])continue;
            isset[findex] = true;
            set_count++;
            auto set_p_fv = fv_begin(findex);
            for (auto p_fv = set_p_fv; p_fv != fv_end(findex); ++p_fv)used_indices[*p_fv] = p_fv - set_p_fv + 1;
            for (auto p_ff = ff_begin(findex); p_ff != ff_end(findex); ++p_ff){
                if (isset[*p_ff])continue;
                order.clear();
                auto np_fvb = fv_begin(*p_ff);
                for (auto np_fv = np_fvb; np_fv != fv_end(*p_ff); ++np_fv){
                    if (used_indices[*np_fv]!=0)order.push_back(make_pair(used_indices[*np_fv], np_fv - np_fvb + 1));

                }
                if (order.size() != 2){ cout << "error in reorintation" << endl; isset[*p_ff] = true; continue; }
                if (is_sequenced(order[0].first, order[1].first) == is_sequenced(order[0].second, order[1].second)){
                    swap(*np_fvb, *(np_fvb + 1));++ccc;
                }
                //isset[*p_ff] = true;
                waitinglist.push_front(*p_ff);
            }
            for (auto p_fv = set_p_fv; p_fv != fv_end(findex); ++p_fv)used_indices[*p_fv] = 0;
        }

        //if(set_count==n_faces)break;
        for(int i = 0; i <n_faces;++i)if(!isset[i])waitinglist.push_back(i);
        //cout<<"orientation"<<endl;
        if(waitinglist.empty())break;
    }


    for(auto a:isset)if(!a)cout<<"orientation:adsadasdasdasdasdasdas"<<endl;


}
void Mesh::InverseFaces(){
    for(int i = 0; i<n_faces;++i){
        auto p_fv = fv_begin(i);
        swap_index(p_fv,p_fv+2);
    }
    for(auto& p_vn : vertices_normal)p_vn *= -1;
    for(auto& p_fn : faces_normal)p_fn *= -1;
}

void Mesh::ComputeFaceNormal(bool isComputeVerticesNormal, vector<bool> *p_isinver){

    if (dimtable[0] != 3 || dimtable[2] != 3){ cout << "Not triangles! Unable to compute normals of faces!" << endl; return; }
    faces_normal.clear();
    faces_normal.resize(3 * n_faces);
    double e1[3], e2[3];
    for (uint i = 0; i < faces_normal.size(); i += 3){
        auto p_fv = &(faces2vertices[i]);
        auto p_n = &(faces_normal[i]);
        auto v1 = &(vertices[3 * p_fv[0]]);
        auto v2 = &(vertices[3 * p_fv[1]]);
        auto v3 = &(vertices[3 * p_fv[2]]);
        for (uint j = 0; j < 3; ++j){
            e1[j] = v2[j] - v1[j];
            e2[j] = v3[j] - v1[j];
        }
        cross(e1,e2,p_n);
        normalize(p_n);
        if(0)if(p_n[0]!=p_n[0]){
            cross(e1,e2,p_n);
            cout<<"Nan face normal: ";for(int j=0;j<3;++j)cout<<p_n[j]<<' ';cout<<endl;
            for(int j=0;j<3;++j)cout<<p_fv[j]<<' ';cout<<endl;
            for(int j=0;j<3;++j)cout<<v1[j]<<' ';cout<<endl;
            for(int j=0;j<3;++j)cout<<v2[j]<<' ';cout<<endl;
            for(int j=0;j<3;++j)cout<<v3[j]<<' ';cout<<endl;
        }

    }
    if(p_isinver!=NULL){
        for(int i=0;i<n_faces;++i)if(p_isinver->at(i))MyUtility::inversevec(fnor_begin(i),fnor_begin(i));
    }
    if(!isComputeVerticesNormal)return;
    if (!isbuildneighbortable){ BuildNeighborTable(); }
    vertices_normal.clear();
    vertices_normal.resize(3 * n_vertices, 0);
    for (uint i = 0; i < n_vertices; i ++){
        auto p_vn = &(vertices_normal[i * 3]);
        for (auto p_vf = vf_begin(i); p_vf != vf_end(i); ++p_vf){
            auto p_fn = get_face_normal(*p_vf);
            for (uint j = 0; j < 3; ++j)p_vn[j] += p_fn[j];
        }
        normalize(p_vn);
        if(p_vn[0]!=p_vn[0]){
            cout<<"Nan Vertices normal: "<<vf_num(i)<<endl;
            if(0){
                //cout<<"Nan Vertices normal: "<<vf_num(i)<<endl;
                p_vn[0]=0;p_vn[1]=0;p_vn[2]=0;
                for (auto p_vf = vf_begin(i); p_vf != vf_end(i); ++p_vf){
                    auto p_fn = get_face_normal(*p_vf);
                    for (uint j = 0; j < 3; ++j)p_vn[j] += p_fn[j];
                }
                for(int j=0;j<3;++j)cout<<p_vn[j]<<' ';cout<<endl;
            }

            p_vn[0]=0;p_vn[1]=0;p_vn[2]=1;
        }
    }

}


void Mesh::SpcecialNormalComputation(vector<bool>*p_isinver,vector<double>&out_vnor,vector<double>&out_mag){


    if (dimtable[0] != 3 || dimtable[2] != 3){ cout << "Not triangles! Unable to compute normals of faces!" << endl; return; }

    vector<double>tmp_faces_normal(3 * n_faces);
    double e1[3], e2[3];
    for (uint i = 0; i < tmp_faces_normal.size(); i += 3){
        auto p_v = &(faces2vertices[i]);
        auto p_n = &(tmp_faces_normal[i]);
        auto v1 = &(vertices[3 * p_v[0]]);
        auto v2 = &(vertices[3 * p_v[1]]);
        auto v3 = &(vertices[3 * p_v[2]]);
        for (uint j = 0; j < 3; ++j){
            e1[j] = v2[j] - v1[j];
            e2[j] = v3[j] - v1[j];
        }
        cross(e1,e2,p_n);
    }
    auto p_fnd = tmp_faces_normal.data();
    if(p_isinver!=NULL){

        for(int i=0;i<n_faces;++i)if(p_isinver->at(i))MyUtility::inversevec(p_fnd+i*3,p_fnd+i*3);
    }
    if (!isbuildneighbortable){ BuildNeighborTable(); }
    out_vnor.clear();
    out_vnor.resize(3 * n_vertices, 0);
    out_mag.clear();
    out_mag.resize(n_vertices);
    for (uint i = 0; i < n_vertices; i ++){
        auto p_vn = &(out_vnor[i * 3]);
        for (auto p_vf = vf_begin(i); p_vf != vf_end(i); ++p_vf){
            auto p_fn = p_fnd+3*(*p_vf);
            for (uint j = 0; j < 3; ++j)p_vn[j] += p_fn[j];
        }

        if(p_vn[0]==p_vn[0])out_mag[i] = sqrt(normVec(p_vn));

        //normalize(p_vn);
        //        if(p_vn[0]!=p_vn[0]){
        //            cout<<"Nan Vertices normal: "<<vf_num(i)<<endl;
        //            if(0){
        //                //cout<<"Nan Vertices normal: "<<vf_num(i)<<endl;
        //                p_vn[0]=0;p_vn[1]=0;p_vn[2]=0;
        //                for (auto p_vf = vf_begin(i); p_vf != vf_end(i); ++p_vf){
        //                    auto p_fn = get_face_normal(*p_vf);
        //                    for (uint j = 0; j < 3; ++j)p_vn[j] += p_fn[j];
        //                }
        //                for(int j=0;j<3;++j)cout<<p_vn[j]<<' ';cout<<endl;
        //            }

        //            p_vn[0]=0;p_vn[1]=0;p_vn[2]=1;
        //        }
    }






}
double Mesh::ComputeTotalArea() {
    double area = 0.0f;
    double e1[3],e2[3],a_unit[3];
    setparameters();

    for(int i = 0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        auto v1 = v_begin(p_fv[0]);
        auto v2 = v_begin(p_fv[1]);
        auto v3 = v_begin(p_fv[2]);
        for (uint j = 0; j < 3; ++j){
            e1[j] = v2[j] - v1[j];
            e2[j] = v3[j] - v2[j];
        }
        cross(e1,e2,a_unit);
        area+=0.5*sqrt(dot(a_unit,a_unit));
    }
    return area;
}
void Mesh::ComputeEdgefromFace(bool isRecomputeNeighborhood){
    
    if(isRecomputeNeighborhood)BuildNeighborTable();
    
    if (!isbuildneighbortable){ BuildNeighborTable(); }
    

    edges2faces.clear();
    edges.clear();

    faces.clear();
    faces.resize(n_faces*dimtable[2],UINT_MAX);
    vector<int>currentedgesnum(n_faces,0);
    int facedim = dimtable[2],edges_index = 0;
    vector<bool>visitedfaces(n_faces,false);
    vector<bool>useVertices(n_vertices,false);

    for(uint i = 0; i<n_faces;++i){
        visitedfaces[i]=true;
        for(auto p_fv = fv_begin(i);p_fv!=fv_end(i);++p_fv)useVertices[*p_fv]=true;
        for(auto p_ff = ff_begin(i);p_ff!=ff_end(i);++p_ff){
            if(visitedfaces[*p_ff])continue;
            for(auto p_fv = fv_begin(*p_ff);p_fv!=fv_end(*p_ff);++p_fv)if(useVertices[*p_fv]){
                edges.push_back(*p_fv);
            }
            if(currentedgesnum[i]<facedim)faces[i*facedim+currentedgesnum[i]] = edges_index;
            ++currentedgesnum[i];
            if(currentedgesnum[*p_ff]<facedim)faces[(*p_ff)*facedim+currentedgesnum[*p_ff]] = edges_index;
            ++currentedgesnum[*p_ff];
            edges2faces.push_back(i);
            edges2faces.push_back(*p_ff);
            ++edges_index;
        }

        for(auto p_fv = fv_begin(i);p_fv!=fv_end(i);++p_fv)useVertices[*p_fv]=false;
    }

    vector<uint>temp;
    vector<uint>temp2;
    map<pair<uint,uint>,bool>e_flag;
    uint flingeedge = 0;


    for(uint i = 0; i<n_faces;++i){
        if(faces2faces[i * f2f_step]!=facedim){
            temp.clear();
            e_flag.clear();
            for(auto p_fv = fv_begin(i);p_fv!=fv_end(i);++p_fv){
                useVertices[*p_fv]=true;
                temp.push_back(*p_fv);
            }
            sort(temp.begin(),temp.end());
            e_flag[make_pair(temp[0],temp[1])]=false;
            e_flag[make_pair(temp[1],temp[2])]=false;
            e_flag[make_pair(temp[0],temp[2])]=false;
            for(auto p_ff = ff_begin(i);p_ff!=ff_end(i);++p_ff){
                temp2.clear();
                for(auto p_fv = fv_begin(*p_ff);p_fv!=fv_end(*p_ff);++p_fv)if(useVertices[*p_fv]){
                    temp2.push_back(*p_fv);
                }

                if(temp2[0]>temp2[1])swap(temp2[0],temp2[1]);
                e_flag[make_pair(temp2[0],temp2[1])]=true;
            }
            for(auto iter = e_flag.begin();iter != e_flag.end();++iter){
                if(!iter->second){
                    edges.push_back(iter->first.first);edges.push_back(iter->first.second);
                    edges2faces.push_back(i);edges2faces.push_back(UINT_MAX);
                    if(currentedgesnum[i]<facedim){faces[i*facedim+currentedgesnum[i]] = edges_index;/*cout<<"hit"<<endl;*/}
                    ++currentedgesnum[i];
                    ++flingeedge;
                    ++edges_index;
                }

            }
            for(auto p_fv = fv_begin(i);p_fv!=fv_end(i);++p_fv)useVertices[*p_fv]=false;

        }

    }


    int numless = 0,numgreater = 0;
    for(auto n:currentedgesnum){
        if(n<3)numless++;
        if(n>3)numgreater++;
    }

    cout<<"less: "<<numless<<' '<<"greater: "<<numgreater<<endl;
    cout<<"new_edges: "<<edges.size()/2<<" flingeedge: "<<flingeedge<<endl;
    setparameters();

    vertices_flinge.clear();
    vertices_flinge.resize(n_vertices,false);
    for(uint i =0;i<n_edges;++i){
        if(isflingeedge(i)){
            for(auto p_ev = ev_begin(i);p_ev!=ev_end(i);++p_ev)vertices_flinge[*p_ev]=true;
        }
    }

    uint tmpedge[3];
    for(uint i=0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        auto p_fe = fe_begin(i);
        for(int j=0;j<3;++j){
            for(int k=0;k<3;++k){
                auto p_ev = ev_begin(p_fe[k]);
                if(p_ev[0]!=p_fv[j] && p_ev[1]!=p_fv[j]){
                    tmpedge[j]=(p_fe[k]);
                }
            }
        }
        for(int j=0;j<3;++j)p_fe[j] = tmpedge[j];
    }




    faces2edgesPN.clear();
    faces2edgesPN.resize(n_faces*dimtable[2],false);
    for(uint i=0;i<n_faces;++i){
        auto p_fepn = fepn_begin(i);
        for(auto p_fe = fe_begin(i);p_fe!=fe_end(i);++p_fe,++p_fepn){
            if((*ef_begin(*p_fe))==i)(*p_fepn)=true;
            else (*p_fepn)=false;

        }
    }


    vertices2edges.clear();
    vertices2edges_accumulate.clear();
    vertices2edges_accumulate.resize(n_vertices,0);

    for (uint i = 0; i < n_edges; i++){
        auto p_ev = ev_begin(i);
        for (uint j = 0; j < 2; ++j){
            ++vertices2edges_accumulate[p_ev[j]];

        }
    }

    //for(auto &a:vertices2edges_accumulate)--a;
    for (uint i = 1; i < n_vertices; i++){
        vertices2edges_accumulate[i]+=vertices2edges_accumulate[i-1];
    }
    vertices2edges_accumulate.insert(vertices2edges_accumulate.begin(),0);
    vertices2edges.resize(vertices2edges_accumulate[n_vertices]);
    vector<uchar>cur_num(n_vertices,0);

    for (uint i = 0; i < n_edges; i++){
        auto p_ev = ev_begin(i);
        for (uint j = 0; j < 2; ++j){
            auto v = p_ev[j];
            vertices2edges[vertices2edges_accumulate[v]+cur_num[v]] = i;
            cur_num[v]++;
        }
    }

    vertices2edgesPN.clear();
    vertices2edgesPN.resize(vertices2edges.size(),false);
    for(uint i=0;i<n_vertices;++i){
        auto p_vepn = vepn_begin(i);
        for(auto p_ve = ve_begin(i);p_ve!=ve_end(i);++p_ve,++p_vepn){
            if((*ev_begin(*p_ve))==i)(*p_vepn)=true;
            else (*p_vepn)=false;

        }
    }


    vertices2vertices.clear();
    vertices2vertices.resize(vertices2edges.size());
    for(int i=0;i<n_vertices;++i){
        auto p_vepn = vepn_begin(i);
        auto p_vv = vv_begin(i);
        for(auto p_ve = ve_begin(i);p_ve!=ve_end(i);++p_ve,++p_vepn,++p_vv)
            if(*p_vepn)(*p_vv) = ev_begin(*p_ve)[1];
            else (*p_vv) = ev_begin(*p_ve)[0];
    }

    vertices2vertices_inverse.clear();
    vertices2vertices_inverse.resize(vertices2edges.size());
    for(int i=0;i<n_vertices;++i){
        auto p_vvinv = vvinv_begin(i);
        auto p_vvend = vv_end(i);
        for(auto p_vv = vv_begin(i);p_vv!=p_vvend;++p_vv,++p_vvinv){
            auto p_vvi = vv_begin(*p_vv);
            auto num = vv_num(*p_vv);
            for(int j=0;j<num;++j)if(p_vvi[j]==i){
                (*p_vvinv)=j;break;
            }
        }
    }
    //for(auto p_vvinv = vvinv_begin(3);p_vvinv!=vvinv_end(3);++p_vvinv)cout<<*p_vvinv<<' ';cout<<endl;

    //for(auto a:vertices2edges_accumulate)cout<<a<<' ';cout<<endl;
    isbuildEdges=true;



}
void Mesh::ComputeEdgefromFaceWithDuplication(){
    edges.clear();
    for(uint i = 0; i<n_faces;++i){
        for(auto p_fv = fv_begin(i);p_fv!=fv_end(i)-1;++p_fv)
        {
            edges.push_back(*p_fv);
            edges.push_back(*(p_fv+1));
        }
        edges.push_back(*(fv_end(i)-1));
        edges.push_back(*(fv_begin(i)));
    }
}

bool Mesh::RayHitTriangle(double* raystart, double* rayend, uint f_ind, double *hitpoint) {
    double ray_d[3],rb2v[3],Intersect[3],e1[3],e2[3],w[3];
    for(int i = 0;i<3;++i)ray_d[i] = rayend[i]-raystart[i];

    auto f_n = get_face_normal(f_ind);
    auto f_v0 = get_vertice(fv_begin(f_ind)[0]);
    auto f_v1 = get_vertice(fv_begin(f_ind)[1]);
    auto f_v2 = get_vertice(fv_begin(f_ind)[2]);
    for(int j = 0;j<3;++j)rb2v[j] = raystart[j]-f_v0[j];

    double a = -dot(f_n,rb2v);
    double b = dot(f_n,ray_d);
    if(fabs(b) < 0.00001){
        if(a == 0){
            if(hitpoint!=NULL)for(int i=0;i<3;++i)hitpoint[i] = (f_v0[i]+f_v1[i]+f_v2[i])/3.;
            return true;
        }
        else return false;
    }
    double r = a / b;
    for(int i = 0;i<3;++i)Intersect[i] = raystart[i] + r * ray_d[i];
    for(int i = 0;i<3;++i)e1[i] = f_v1[i] - f_v0[i];
    for(int i = 0;i<3;++i)e2[i] = f_v2[i] - f_v1[i];

    double    uu, uv, vv, wu, wv, D;
    uu = dot(e1,e1);
    uv = dot(e1,e2);
    vv = dot(e2,e2);
    for(int i = 0;i<3;++i)w[i] = Intersect[i] - f_v0[i];
    wu = dot(w,e1);
    wv = dot(w,e2);
    D = uv * uv - uu * vv;

    double s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)         // I is outside T
        return false;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
        return false;

    if(hitpoint!=NULL)for(int i=0;i<3;++i)hitpoint[i] = Intersect[i];
    return true;

}
bool Mesh::RayHitTarget(double* raystart,double* rayend){

    for(int i = 0; i < n_faces; ++i){
        if(RayHitTriangle(raystart,rayend,i))return true;
    }
    return false;

}
void Mesh::addCellComplex(CellComplex& cc_in){
    for (auto& v : cc_in.vertices)vertices.push_back(v);
    for (auto& e : cc_in.edges)edges.push_back(e + n_vertices);
    for (auto& f : cc_in.faces2vertices)faces2vertices.push_back(f + n_vertices);
    setparameters();
}

void Mesh::addMesh(Mesh& mesh_in){
    if(faces_normal.size()!=0&&mesh_in.faces_normal.size()!=0){
        for (auto& fn : mesh_in.faces_normal)faces_normal.push_back(fn);
    }
    if(vertices_normal.size()!=0&&mesh_in.vertices_normal.size()!=0){
        for (auto& vn : mesh_in.vertices_normal)vertices_normal.push_back(vn);
    }
    this->addCellComplex(mesh_in);
}
void Mesh::TriangleCircumcenter(uint f_ind,double *cc,double &r){
    double *pv[3];
    auto p_fv = fv_begin(f_ind);
    for(int i = 0;i<3;++i){
        pv[i] = get_vertice(p_fv[i]);
    }
    _TriangleCircumcenter(pv[0],pv[1],pv[2],cc,r);
}

void Mesh::TriangleMidpoint(uint f_ind,double *cc){
    double *pv[3];
    auto p_fv = fv_begin(f_ind);
    for(int i = 0;i<3;++i){
        pv[i] = get_vertice(p_fv[i]);
    }
    _TriangleMidpoint(pv[0],pv[1],pv[2],cc,3);

}

double Mesh::point2TriSquareDistance(const double *P,uint find,double *cp){
    auto p_fv = fv_begin(find);
    const double *v1 = v_begin(p_fv[0]);
    const double *v2 = v_begin(p_fv[1]);
    const double *v3 = v_begin(p_fv[2]);
    return MyUtility::point2TriSquareDistance(P,v1,v2,v3,cp);


}

void smoothing(vector<double>& points, vector< uint >& neighborhood, vector< unsigned long>&v2v_acc,double coef,vector<bool>*staticv = NULL){
    auto oldpoints = points;
    int num = v2v_acc.size() - 1;
    uint p_dim = points.size() / num;
    double o_coef = 1 - coef;
    bool isstatic = (staticv !=NULL);
    for (int i = 0; i < num; ++i){
        if(isstatic)if(staticv->at(i))continue;
        auto n_n = v2v_acc[i + 1] - v2v_acc[i];
        if(n_n==0)continue;
        vector<double>midpoint(p_dim, 0);
        uint* p_n = &(neighborhood[v2v_acc[i]]);
        for (int k = 0; k < n_n;++k){
            auto begin = p_n[k]*p_dim;
            for (int j = 0; j < p_dim; ++j){
                midpoint[j] += oldpoints[begin + j];
            }
        }
        for (int j = 0; j < p_dim; ++j){
            midpoint[j] /= n_n;
            points[i * p_dim + j] = o_coef * oldpoints[i * p_dim + j] + coef * midpoint[j];
        }
    }
}
void BuildRingNeighborhood(int rings,vector< uint >& neighborhood, vector< unsigned long>&v2v_acc,
                           vector<uint>&ringneighbor,vector< unsigned long>&ringv2v_acc){

    auto vv_begin = [&neighborhood,&v2v_acc](uint v_ind){ return &(neighborhood[v2v_acc[v_ind]]); };
    auto vv_end = [&neighborhood,&v2v_acc](uint v_ind){ return &(neighborhood[v2v_acc[v_ind+1]]); };
    int n_vertices = v2v_acc.size()-1;
    cout<<"BuildRingNeighborhood "<<n_vertices<<endl;
    ringv2v_acc.resize(n_vertices,0);
    ringneighbor.clear();
    vector<bool>pickrings(n_vertices,false);
    list<uint>vbuffer1,vbuffer2;
    list<uint>*pbuffer_pre = &vbuffer1,*pbuffer_next = &vbuffer2;
    vector<uint>pickbuffer;
    for (uint i = 0; i < n_vertices; i++){
        pbuffer_pre->clear();
        pbuffer_next->clear();
        pbuffer_pre->push_front(i);
        pickbuffer.clear();
        pickrings[i] = true;

        for(int r=0;r<rings;++r){
            while(!pbuffer_pre->empty()){
                auto vind = pbuffer_pre->front();
                pbuffer_pre->pop_front();
                pickbuffer.push_back(vind);
                auto vvend = vv_end(vind);
                for(auto p_vv = vv_begin(vind);p_vv!=vvend;++p_vv){
                    if(pickrings[*p_vv])continue;
                    pbuffer_next->push_back(*p_vv);pickrings[*p_vv]=true;
                    pickbuffer.push_back(*p_vv);
                }
            }
            swap(pbuffer_next,pbuffer_pre);
        }

        for(auto a:pickbuffer)pickrings[a] = false;
        pickbuffer.erase(pickbuffer.begin());
        //cout<<pickbuffer.size()<<endl;
        ringv2v_acc[i+1] = pickbuffer.size();


        for(auto a:pickbuffer)ringneighbor.push_back(a);
    }
    for(int i=0;i<n_vertices;++i)ringv2v_acc[i+1] = ringv2v_acc[i] + ringv2v_acc[i+1];

}
void Mesh::Fairing(bool usingedges, bool iskeepboundary, bool isRecomputeNeighborhood, int round, double lamda, double kaipa,vector<bool>* staticv){
    if(isRecomputeNeighborhood)BuildNeighborTable();
    //if (!isbuildneighbortable)BuildNeighborTable();
    auto* link = &edges;
    int linkdim = 0;
    if (usingedges){
        linkdim = dimtable[1];
    }
    else{
        link = &faces2vertices;
        linkdim = dimtable[2];
    }
    auto& nlink = *link;

    bool ueb = usingedges&&iskeepboundary;
    vector<bool>staticvertices;
    //vector<bool>* staticv = NULL;
    if(ueb){
        staticvertices.resize(n_vertices,false);
        for (int i = 0; i < nlink.size(); i += linkdim){
            if(isflingeedge(i/linkdim)){
                auto p_v = &(nlink[i]);
                for (int j = 0; j < linkdim; ++j)staticvertices[p_v[j]]=true;
            }
        }
        staticv = &staticvertices;
    }

    auto gama = 1 / (kaipa - 1 / lamda);
    // gama = 0;
    //    for (int i = 0; i < round; ++i){
    //        smoothing(vertices, vertices2vertices, vertices2edges_accumulate,lamda,staticv);
    //        smoothing(vertices, neighborhood, v2v_acc,gama,staticv);
    //    }
    //    for (int i = 0; i < round; ++i){
    //        for(int j=0;j<5;++j)smoothing(vertices, vertices2vertices_nonmanifold, vertices2edges_accumulate_nonmanifold,lamda,staticv);
    //        for(int j=0;j<4;++j)smoothing(vertices, vertices2vertices_nonmanifold, vertices2edges_accumulate_nonmanifold,gama,staticv);
    //    }
    //vector<uint>ringneighbor;vector< unsigned long>ringv2v_acc;
    //BuildRingNeighborhood(1,vertices2vertices,vertices2edges_accumulate, ringneighbor,ringv2v_acc);

    for (int i = 0; i < round; ++i){
        //smoothing(vertices, ringneighbor, ringv2v_acc,lamda,staticv);
        //smoothing(vertices, ringneighbor, ringv2v_acc,gama,staticv);
        smoothing(vertices, vertices2vertices_nonmanifold, vertices2edges_accumulate_nonmanifold,lamda,staticv);
        smoothing(vertices, vertices2vertices_nonmanifold, vertices2edges_accumulate_nonmanifold,gama,staticv);
    }
    //cout<<"ringneighbor "<<ringneighbor.size()<<' '<<vertices2vertices.size()<<endl;
}

void computeCenterDiff(const vector<double>& points,const vector< uint >& neighborhood, const vector< unsigned long>&v2v_acc,
                       vector<double>&dis2center){
    int numP = v2v_acc.size()-1;
    uint p_dim = points.size() / numP;
    dis2center.clear();
    dis2center.resize(points.size(),0);

    for (int i = 0; i < numP; ++i){
        auto p_n = &(neighborhood[v2v_acc[i]]);
        auto n_n = v2v_acc[i + 1] - v2v_acc[i];
        auto p_d = dis2center.data()+i*3;
        for (int k = 0; k < n_n; ++k){
            auto begin = p_n[k]*p_dim;
            for (int j = 0; j < p_dim; ++j){
                p_d[j] += points[begin + j];
            }
        }
        for (int j = 0; j < p_dim; ++j){
            p_d[j] /= n_n;
        }
        for (int j = 0; j < p_dim; ++j){
            p_d[j] -= points[i*3 + j];
        }
    }

}
void SurfaceFFFairing_OneIter(vector<double>& points, vector<double>&pointnor, vector<double>&norsqrtmag,
                              vector< uint >& v2v, vector< unsigned long>&v2v_acc,vector<bool>&point_nonmanifold,
                              vector<double>&w2list,double coef,vector<bool>*staticv = NULL){
    //auto oldpoints = points;
    int numV = v2v_acc.size() - 1;

    uint p_dim = points.size() / numV;
    //cout<<p_dim<<endl;
    //double o_coef = 1 - coef;
    //bool isstatic = (staticv !=NULL);
    vector<double>firstdiff;
    vector<double>seconddiff;
    vector<double>tangentdiff(numV*p_dim,0);
    vector<double>dis(numV*p_dim,0);
    computeCenterDiff(points,v2v,v2v_acc,firstdiff);
    //computeCenterDiff(firstdiff,v2v,v2v_acc,seconddiff);

    auto p_v = points.data();
    auto p_first = firstdiff.data();
    auto p_nor = pointnor.data();
    auto p_tang = tangentdiff.data();


    for(int i=0;i<numV;++i){
        auto ind = i*p_dim;
        double sqrtm = norsqrtmag[i];
        double mag = sqrtm*sqrtm;
        if(!point_nonmanifold[i]){
            dis[ind] = dot(p_first+ind,p_nor+ind)/(sqrtm*mag);
            product(1/sqrtm,p_nor+ind,p_nor+ind);
            weightedAddVec(1.,-dis[ind],p_first+ind,p_nor+ind,p_tang+ind);
        }else{
            copyVec(p_first+ind,dis.data()+ind);
        }
    }

    computeCenterDiff(dis,v2v,v2v_acc,seconddiff);
    auto p_sec = seconddiff.data();
    //auto p_oldv = oldpoints.data();
    for(int i=0;i<numV;++i){
        if(staticv!=NULL && staticv->at(i))continue;
        auto ind = i*p_dim;
        auto p_cv = p_v+ind;
        if(!point_nonmanifold[i]){
            weightedAddVec(-coef*p_sec[ind],p_nor+ind,p_cv);
            weightedAddVec(coef,p_tang+ind,p_cv);
        }else{
            weightedAddVec(-w2list[i],p_sec+ind,p_cv);
        }
    }

}
void SurfaceFFFairing_OneIter(vector<double>& points, vector<double>&pointnor,
                              vector< uint >& v2v, vector< unsigned long>&v2v_acc,vector<bool>&point_nonmanifold,
                              vector<double>&w2list,double coef,vector<bool>*staticv = NULL){
    //auto oldpoints = points;
    int numV = v2v_acc.size() - 1;

    uint p_dim = points.size() / numV;
    //cout<<p_dim<<endl;
    //double o_coef = 1 - coef;
    //bool isstatic = (staticv !=NULL);
    vector<double>firstdiff;
    vector<double>seconddiff;
    vector<double>tangentdiff(numV*p_dim,0);
    vector<double>dis(numV*p_dim,0);
    computeCenterDiff(points,v2v,v2v_acc,firstdiff);
    //computeCenterDiff(firstdiff,v2v,v2v_acc,seconddiff);

    auto p_v = points.data();
    auto p_first = firstdiff.data();
    auto p_nor = pointnor.data();
    auto p_tang = tangentdiff.data();


    for(int i=0;i<numV;++i){
        auto ind = i*p_dim;

        if(!point_nonmanifold[i]){
            dis[ind] = dot(p_first+ind,p_nor+ind);
            weightedAddVec(1.,-dis[ind],p_first+ind,p_nor+ind,p_tang+ind);
        }else{
            copyVec(p_first+ind,dis.data()+ind);
        }
    }

    computeCenterDiff(dis,v2v,v2v_acc,seconddiff);
    auto p_sec = seconddiff.data();
    //auto p_oldv = oldpoints.data();
    for(int i=0;i<numV;++i){
        if(staticv!=NULL && staticv->at(i))continue;
        auto ind = i*p_dim;
        auto p_cv = p_v+ind;
        if(!point_nonmanifold[i]){
            weightedAddVec(-coef*p_sec[ind],p_nor+ind,p_cv);
            weightedAddVec(coef,p_tang+ind,p_cv);
        }else{
            weightedAddVec(-w2list[i],p_sec+ind,p_cv);
        }
    }

}
void Mesh::SurfaceFFFairing(double lambda,int iter,vector<bool>* inverNor,vector<bool>* staticv){
    if (!isbuildneighbortable)BuildNeighborTable();

    cout<<"Ju Smoothing"<<endl;

    vector<double>w2list(n_vertices,0);
    vector<double>valence(n_vertices);
    double maxValence = 0.4;
    for(int i=0;i<n_vertices;++i){
        double a = vvnon_num(i);
        if(a!=0)valence[i] = 1/a;else valence[i] = 0;
    }
    for(int i=0;i<n_vertices;++i){
        int a = vvnon_num(i);
        auto p_vv = vvnon_begin(i);
        for(int j=0;j<a;++j){
            w2list[i]+=valence[p_vv[j]];
        }
        w2list[i]*=valence[i];
        w2list[i] += 1;
        w2list[i] = 1/w2list[i];
        if( w2list[i] > maxValence)w2list[i] = maxValence;
    }



    vector<double>tmpvnormal,tmpvnorsqrtMag;
    for(int i=0;i<iter;++i){
        if(1){
            SpcecialNormalComputation(inverNor,tmpvnormal,tmpvnorsqrtMag);

            SurfaceFFFairing_OneIter(vertices,tmpvnormal,tmpvnorsqrtMag,
                                     vertices2vertices_nonmanifold,vertices2edges_accumulate_nonmanifold,vertices_non_manifold,
                                     w2list,lambda,staticv);
        }else{
            ComputeFaceNormal(true,inverNor);
            SurfaceFFFairing_OneIter(vertices,vertices_normal,
                                     vertices2vertices_nonmanifold,vertices2edges_accumulate_nonmanifold,vertices_non_manifold,
                                     w2list,lambda,staticv);
        }
        //cout<<"ads"<<endl;
    }

    ComputeFaceNormal(true);

}

void Mesh::SurfaceLaplacianFairing(double lambda,int iter,vector<bool>* staticv){
    if (!isbuildneighbortable)BuildNeighborTable();

    for(int i=0;i<iter;++i){
        smoothing(vertices, vertices2vertices_nonmanifold, vertices2edges_accumulate_nonmanifold,lambda,staticv);
    }

}



bool Mesh::readObjfile(string filename){

    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the OBJ file " << filename << endl;
        return false;
    }

    reset();

    auto readVertices = [this](stringstream &ss){
        double dvalue;
        for(int i=0;i<3;++i){ss>>dvalue;vertices.push_back(dvalue);}
    };
    auto readFace = [this](stringstream &ss){
        int ivalue;
        string component;
        for(int i=0;i<3;++i){
            ss>>component;
            stringstream tt(component);
            tt>>ivalue;faces2vertices.push_back(ivalue-1);
        }
    };
    //    auto readVerticesField = [this](stringstream &ss){
    //        double dvalue;
    //        for(int i=0;i<3;++i){ss>>dvalue;vertices_field.push_back(dvalue);}
    //    };


    string oneline;

    cout<<"reading: "<<filename<<endl;

    while( getline( reader, oneline ))
    {
        stringstream ss( oneline );
        string token;

        ss >> token;

        if( token == "v"  ) { readVertices( ss ); continue; } // vertex
        if( token == "vt" ) {  continue; } // texture coordinate
        if( token == "vn" ) {  continue; } // vertex normal
        if( token == "vf" ) { /*readVerticesField( ss );*/ continue; } // tangent vector
        if( token == "f"  ) { readFace( ss ); continue; } // face
        if( token[0] == '#' ) continue; // comment
        if( token == "o" ) continue; // object name
        if( token == "g" ) continue; // group name
        if( token == "s" ) continue; // smoothing group
        if( token == "mtllib" ) continue; // material library
        if( token == "usemtl" ) continue; // material
        if( token == "k" ) continue; // field degree
        if( token == "fs" ) continue; // field singularity
        if( token == "" ) continue; // empty string
        if( token == "c" ) continue;


        cerr << "Error: does not appear to be a valid Wavefront OBJ file!" << endl;
        cerr << "(Offending line: " << oneline << ")" << endl;
        return false;
    }
    //cout<<"nfvec "<<nfvec<<endl;


    reader.close();

    setparameters();
    isbuildneighbortable=false;
    return true;

}
bool Mesh::readOfffile(string filename){
    ifstream reader(filename.data(), ofstream::in);
    if (!reader.good()) {
        cout << "Can not open the OFF file " << filename << endl;
        return false;
    }else {
        cout << "Reading: "<<filename<<endl;
    }

    string ss;

    int ivalue,ibuf[20];
    reset();

    reader>>ss;


    //cout<<ss<<endl;
    if(ss!="OFF"){
        cout << "Not OFF file: " << filename << endl;
        return false;
    }


    reader>>n_vertices;
    reader>>n_faces;
    reader>>n_edges;

    cout<<n_vertices<<' '<< n_faces<<' '<<n_edges<<endl;
    vertices.resize(n_vertices*3);
    //faces2vertices.resize(n_faces*3);

    for(int i =0;i<vertices.size();i++){
        reader>>vertices[i];
    }
    for(int i =0;i<n_faces;i++){
        reader>>ivalue;
        if(ivalue==3){for(int j =0;j<3;++j){reader>>ivalue;faces2vertices.push_back(ivalue);}}
        else if(ivalue==4){
            for(int j =0;j<4;++j)reader>>ibuf[j];
            for(int j =0;j<3;++j)faces2vertices.push_back(ibuf[j]);
            for(int j =1;j<4;++j)faces2vertices.push_back(ibuf[j]);
        }
    }
    //for(auto& f:faces)f-=1;
    //for(auto& f:faces2vertices)cout<<f<<' ';

    cout<<"Read Finished!"<<endl;
    reader.close();
    setparameters();
    isbuildneighbortable=false;
    isbuildEdges=false;

    return true;


}
bool Mesh::saveObjFile(string filename){
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open the OBJ file " << filename << endl;
        return false;
    }

    for(int i=0;i<n_vertices;++i){
        auto p_v = v_begin(i);
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }

    for(int i=0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        outer << "f " << p_fv[0]+1<< " "<< p_fv[1]+1 << " "<< p_fv[2]+1 << endl;
    }

    outer.close();
    cout<<"saving finish: "<<filename<<endl;

    return true;
}

bool Mesh::savePLY2File(string filename){
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }

    outer<<n_vertices<<endl;
    outer<<n_faces<<endl;
    for(int i=0;i<n_vertices;++i){
        auto p_v = v_begin(i);
        outer <<  p_v[0] << endl<< p_v[1] << endl<< p_v[2] << endl;
    }

    for(int i=0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        outer << '3' << endl<< p_fv[0]<< endl<< p_fv[1] <<  endl<< p_fv[2] << endl;
    }

    outer.close();
    cout<<"saving finish: "<<filename<<endl;

    return true;


}

/*********************************************************************************************************/
double Mesh::ComputeDihedralAngle(uint e_ind){
    if(isflingeedge(e_ind))return 0;
    auto p_ef = ef_begin(e_ind);
    auto p_fn1 = fnor_begin(p_ef[0]);
    auto p_fn2 = fnor_begin(p_ef[1]);

    auto p_ev = ev_begin(e_ind);
    double evec[3];
    minusVec(v_begin(p_ev[0]),v_begin(p_ev[1]),evec);
    normalize(evec);

    double cn1n2[3];
    cross(p_fn1,p_fn2,cn1n2);


    return atan2( dot(evec,cn1n2), dot(p_fn1,p_fn2) );


}
double Mesh::EdgeCottan(uint e_ind){
    auto p_ev = ev_begin(e_ind);
    uint va = p_ev[0],vb = p_ev[1],vc;
    double w[3],u[3],v[3];
    const double thr = 1e-7;
    auto p_ef = ef_begin(e_ind);
    auto p_fv = fv_begin(p_ef[0]);
    for(int i=0;i<3;++i)if(p_fv[i]!=va &&p_fv[i]!=vb){vc=p_fv[i];break;}

    minusVec(v_begin(va),v_begin(vc),w);
    minusVec(v_begin(vb),v_begin(vc),u);


    cross( w, u,v );
    double cotTheta1 = dot( w, u ) / normVec(v);
    if( cotTheta1 < thr )cotTheta1 = thr;

    if(ef_num(e_ind)==1)return cotTheta1;

    p_fv = fv_begin(p_ef[1]);
    for(int i=0;i<3;++i)if(p_fv[i]!=va &&p_fv[i]!=vb){vc=p_fv[i];break;}

    minusVec(v_begin(va),v_begin(vc),w);
    minusVec(v_begin(vb),v_begin(vc),u);


    cross( w, u,v );
    double cotTheta2 = dot( w, u ) / normVec(v);
    if( cotTheta2 < thr )cotTheta2 = thr;

    return (cotTheta1+cotTheta2)/2.;



}
//#include "my_mesh.inl"
