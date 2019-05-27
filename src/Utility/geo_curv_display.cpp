#include "geo_curv.h"
#include <stdio.h>
#include<iostream>
#include<fstream>
#include<limits>
#include <functional>
#include <eigen3/Eigen/Geometry>
#include<set>

double controlBallsize = 1.;
double controlTubesize = 1.;
void ComputeAngleAxisMatrix(Eigen::AngleAxisd &m_aa,double *oriVec,double *newVec){
    Eigen::Vector3d oriVecE,newVecE,nc;
    for(int j=0;j<3;++j)oriVecE(j) = oriVec[j];
    for(int j=0;j<3;++j)newVecE(j) = newVec[j];
    nc = oriVecE.cross(newVecE);
    if(len(nc.data())<1e-6){nc(0) = 1.0;nc(1) = 0.0;nc(2) = 0.0;}
    nc.normalize();
    double angle = acos(cosine(oriVec,newVec));
    m_aa = Eigen::AngleAxisd(angle,nc);
}

void RotateVertices(vector<double>&inv,double *oriVec,double *newVec,vector<double>&outv){

    int nvs = inv.size()/3;
    Eigen::AngleAxisd t;
    vector<Eigen::Vector3d>ori(nvs);
    vector<Eigen::Vector3d>tranformed(nvs);

    ComputeAngleAxisMatrix(t,oriVec,newVec);
    for(int j = 0;j<nvs;++j){
        auto p_spv = inv.data()+j*3;
        for(int k=0;k<3;++k)ori[j](k) = p_spv[k];
    }
    for(int j = 0;j<nvs;++j)tranformed[j] = t*ori[j];

    outv.resize(inv.size());
    for(int j = 0;j<nvs;++j){
        auto p_spv = outv.data()+j*3;
        for(int k=0;k<3;++k)p_spv[k] = tranformed[j](k);
    }

}

void ThresholdColoring(double degree,uchar *p_rgb);
namespace n_rf {




void Curve::BuildDisplay(bool isrescale, bool isuseball, bool isshowtangent, bool isspecifiedcolor){
    setparameters();
    //scale = 3;
    //cout<<"Curve display: "<<n_vertices<<' '<<n_edges<<endl;
    double sphere_scale = 0.20;
    double tube_scale = 0.20;
    //cout<<"isrescale "<<isrescale<<endl;
    //isrescale = false;
    //bool isuseball = false;
    int nonsize = 12;


    if(isrescale){
        double xmin = 999, ymin=999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
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

        //for(auto& v:vertices)v/=largestdis;
        for (int j = 0; j<n_vertices*3; j+=3)
            for(int i = 0; i<3; i++)
                vertices[j+i]=(vertices[j+i]-centers[i])/largestdis;


    }

    avedis = 0.0;
    for(int i =0; i< n_edges; i++){
        auto p_v = ev_begin(i);
        avedis+=VerticesDistance(*p_v,*(p_v+1));
    }
    avedis/=n_edges;


    scale=avedis*disscale * controlTubesize;
    avedis/=0.5;
    sphere_scale*=avedis*disscale * controlBallsize;
    tube_scale*=avedis*disscale * controlTubesize;
    if(isuseball)if(sphere_scale==sphere_scale)sphere.ReScale_uniform(sphere_scale);
    //cout<<"isuseball: "<<avedis<<' '<<sphere_scale<<endl;




    display_vertices = vertices;
    //display_edges = edges;
    display_faces.clear();
    //rvector = vertices;

    display_vnormals=display_vertices;



    double *p_spv = NULL, *p_v = NULL,*p_vn = NULL;
    uint *p_sf = NULL;
    auto n_fs = sphere.n_faces*3;
    auto n_vs = sphere.n_vertices*3;
    //cout<<"ppppp: "<<display_vertices.size()<<endl;

    if(isuseball)for(int iv = 0;iv<n_vertices;++iv){

        int offset = display_vertices.size()/3;
        p_v = v_begin(iv);
        if(vv_num(iv)>nonsize){for(int j = 0;j<sphere.n_vertices;++j){
                p_spv = sphere.v_begin(j);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]*3);
            }
        }else for(int j = 0;j<sphere.n_vertices;++j){
            p_spv = sphere.v_begin(j);
            for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]);
        }
        p_sf = sphere.fv_begin(0);
        for(int j = 0;j<n_fs;++j){
            display_faces.push_back(offset+p_sf[j]);
        }
        p_vn = sphere.vnor_begin(0);
        for(int j = 0;j<n_vs;++j){
            display_vnormals.push_back(p_vn[j]);
        }
    }
    double vec[3],nz[3] = {0,0,1};
    Eigen::AngleAxisd t;
    vector<Eigen::Vector3d>ori(cylinder.n_vertices);
    vector<Eigen::Vector3d>tranformed(cylinder.n_vertices);
    for(int j = 0;j<cylinder.n_vertices;++j){
        p_spv = cylinder.v_begin(j);
        for(int k=0;k<3;++k)ori[j](k) = p_spv[k];
    }
    vector<Eigen::Vector3d>orinor(cylinder.n_vertices);
    for(int j = 0;j<cylinder.n_vertices;++j){
        p_spv = cylinder.vnor_begin(j);
        for(int k=0;k<3;++k)orinor[j](k) = p_spv[k];
    }

    if(isuseball)for(int iv = 0;iv<n_edges;++iv){
        int offset = display_vertices.size()/3;
        auto p_ev = ev_begin(iv);
        minusVec(v_begin(p_ev[0]),v_begin(p_ev[1]),vec);
        normalize(vec);
        ComputeAngleAxisMatrix(t,nz,vec);
        Eigen::Transform<double,3,Eigen::Affine> t_m(t);
        _VerticesMidpoint(v_begin(p_ev[0]),v_begin(p_ev[1]),vec);
        //t_m.translate(Eigen::Vector3d(vec[0],vec[1],10*vec[2]));
        double d = _VerticesDistance(v_begin(p_ev[0]),v_begin(p_ev[1]));
        t_m.scale(Eigen::Vector3d(tube_scale,tube_scale,d));
        for(int j = 0;j<cylinder.n_vertices;++j){
            tranformed[j] = t_m*ori[j];
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+vec[k]);
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            tranformed[j] = t_m*orinor[j];
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            for(int k=0;k<3;++k)display_vnormals.push_back(tranformed[j](k));
        }

        p_sf = cylinder.fv_begin(0);
        for(int j = 0;j<cylinder.n_faces*3;++j){
            display_faces.push_back(offset+p_sf[j]);
        }


    }

    n_ver_f1 = display_vertices.size();
    //cout<<"fff: "<<display_faces.size()<<endl;

    n_face_f1 = display_faces.size();

    n_vnor_f1 = display_vnormals.size();






    unsigned char red[4] = {225,0,0,255};
    unsigned char blue[4] = {0,0,225,255};
    unsigned char green[4] = {0,255,0,255};
    unsigned char gray[4] = {215,215,215,255};
    unsigned char ccccoooollll[4] = {225,225,225,255};
    unsigned char* pcolor = gray;
    if(isspecifiedcolor){
        ThresholdColoring(colordegree,ccccoooollll);
        pcolor = ccccoooollll;
    }
    v_color.clear();
    for (int i = 0; i < n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(red[j]);

    if(isuseball){
        for(int iv = 0;iv<n_vertices;++iv){
            if(vv_num(iv)>nonsize)for(int k =0;k<sphere.n_vertices;++k)for(int j =0;j<4;++j)v_color.push_back(red[j]);
            else for(int k =0;k<sphere.n_vertices;++k)for(int j =0;j<4;++j)v_color.push_back(pcolor[j]);
        }
        //for (int i = 0; i <  sphere.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
        //for (int i = sphere.n_vertices; i < n_vertices * sphere.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
        for (int i = 0; i < n_edges * cylinder.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(pcolor[j]);
    }


    display_edges.clear();
    if(isshowtangent &&tangent.size()==vertices.size()){
        int offset = display_vertices.size()/3;
        double lenscale = 0.2;
        double newvv[3];
        for(int i=0;i<n_vertices;++i){
            weightedAddVec(1.0,lenscale,v_begin(i),t_begin(i),newvv);
            for(int j=0;j<3;++j)display_vertices.push_back(newvv[j]);
            display_edges.push_back(i);
            display_edges.push_back(i+offset);
            for(int j=0;j<3;++j)display_vnormals.push_back(newvv[3]);
            for(int j =0;j<4;++j)v_color.push_back(red[j]);
        }




    }



    n_col_f1 = v_color.size();
    isbuild = true;


}

void Curve::BuildDisplay(vector<double>&vcolordegree,vector<double>&ecolordegree,int method){




    setparameters();
    //scale = 3;
    //cout<<"Curve display: "<<n_vertices<<' '<<n_edges<<endl;
    double sphere_scale = 0.15;
    //cout<<"isrescale "<<isrescale<<endl;
    //isrescale = false;
    //bool isuseball = false;
    int nonsize = 12;

    avedis = 0.0;
    for(int i =0; i< n_edges; i++){
        auto p_v = ev_begin(i);
        avedis+=VerticesDistance(*p_v,*(p_v+1));
    }
    avedis/=n_edges;


    if(method!=0)avedis = 1;
    scale=avedis*disscale * controlTubesize;
    avedis/=0.5;
    sphere_scale*=avedis*disscale * controlBallsize;
    if(sphere_scale==sphere_scale)sphere.ReScale_uniform(sphere_scale);
    //cout<<"isuseball: "<<avedis<<' '<<sphere_scale<<endl;




    display_vertices = vertices;
    //display_edges = edges;
    display_faces.clear();
    //rvector = vertices;

    display_vnormals=display_vertices;



    double *p_spv = NULL, *p_v = NULL,*p_vn = NULL;
    uint *p_sf = NULL;
    auto n_fs = sphere.n_faces*3;
    auto n_vs = sphere.n_vertices*3;
    //cout<<"ppppp: "<<display_vertices.size()<<endl;

    for(int iv = 0;iv<n_vertices;++iv){

        int offset = display_vertices.size()/3;
        p_v = v_begin(iv);
        if(vv_num(iv)>nonsize){for(int j = 0;j<sphere.n_vertices;++j){
                p_spv = sphere.v_begin(j);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]*2);
            }
        }else for(int j = 0;j<sphere.n_vertices;++j){
            p_spv = sphere.v_begin(j);
            for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]);
        }
        p_sf = sphere.fv_begin(0);
        for(int j = 0;j<n_fs;++j){
            display_faces.push_back(offset+p_sf[j]);
        }
        p_vn = sphere.vnor_begin(0);
        for(int j = 0;j<n_vs;++j){
            display_vnormals.push_back(p_vn[j]);
        }
    }
    double vec[3],nz[3] = {0,0,1};
    Eigen::AngleAxisd t;
    vector<Eigen::Vector3d>ori(cylinder.n_vertices);
    vector<Eigen::Vector3d>tranformed(cylinder.n_vertices);
    for(int j = 0;j<cylinder.n_vertices;++j){
        p_spv = cylinder.v_begin(j);
        for(int k=0;k<3;++k)ori[j](k) = p_spv[k];
    }
    vector<Eigen::Vector3d>orinor(cylinder.n_vertices);
    for(int j = 0;j<cylinder.n_vertices;++j){
        p_spv = cylinder.vnor_begin(j);
        for(int k=0;k<3;++k)orinor[j](k) = p_spv[k];
    }

    double cyclinder_scale = avedis*disscale * controlTubesize*0.2;
    for(int iv = 0;iv<n_edges;++iv){
        int offset = display_vertices.size()/3;
        auto p_ev = ev_begin(iv);
        minusVec(v_begin(p_ev[0]),v_begin(p_ev[1]),vec);
        normalize(vec);
        ComputeAngleAxisMatrix(t,nz,vec);
        Eigen::Transform<double,3,Eigen::Affine> t_m(t);
        _VerticesMidpoint(v_begin(p_ev[0]),v_begin(p_ev[1]),vec);
        //t_m.translate(Eigen::Vector3d(vec[0],vec[1],10*vec[2]));
        double d = _VerticesDistance(v_begin(p_ev[0]),v_begin(p_ev[1]));
        t_m.scale(Eigen::Vector3d(cyclinder_scale,cyclinder_scale,d));
        for(int j = 0;j<cylinder.n_vertices;++j){
            tranformed[j] = t_m*ori[j];
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+vec[k]);
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            tranformed[j] = t_m*orinor[j];
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            for(int k=0;k<3;++k)display_vnormals.push_back(tranformed[j](k));
        }

        p_sf = cylinder.fv_begin(0);
        for(int j = 0;j<cylinder.n_faces*3;++j){
            display_faces.push_back(offset+p_sf[j]);
        }


    }

    n_ver_f1 = display_vertices.size();
    //cout<<"fff: "<<display_faces.size()<<endl;

    n_face_f1 = display_faces.size();

    n_vnor_f1 = display_vnormals.size();


    unsigned char red[4] = {225,0,0,255};
    unsigned char blue[4] = {0,0,225,255};
    unsigned char green[4] = {0,255,0,255};
    unsigned char gray[4] = {225,225,225,255};
    unsigned char ccccoooollll[4] = {225,225,225,255};
    v_color.clear();
    for (int i = 0; i < n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(red[j]);

    {
        for(int iv = 0;iv<n_vertices;++iv){
            ThresholdColoring(vcolordegree[iv],ccccoooollll);
            for(int k =0;k<sphere.n_vertices;++k)for(int j =0;j<4;++j)v_color.push_back(ccccoooollll[j]);
        }
        //for (int i = 0; i <  sphere.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
        //for (int i = sphere.n_vertices; i < n_vertices * sphere.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
        for (int ie = 0; ie < n_edges; ie++){
            ThresholdColoring(ecolordegree[ie],ccccoooollll);
            for(int k =0;k<cylinder.n_vertices;++k)for(int j =0;j<4;++j)v_color.push_back(ccccoooollll[j]);
        }
    }


    display_edges.clear();

    n_col_f1 = v_color.size();
    isbuild = true;




}


/*******************************************************/

void PartialCurve::BuildDisplay(bool isrescale, bool isuseball, bool isShowExist, bool isShowMerge){

    DisplayCurve.BuildDisplay(isrescale,isuseball);
    MergeCurve.BuildDisplay(isrescale,isuseball);


    //    vector<double>* getDisplayVertices(){return DisplayCurve.getDisplayVertices();}
    //    vector<double>* getDisplayVerticesNormal(){return DisplayCurve.getDisplayVerticesNormal();}
    //    vector<uint>* getDisplayEdges(){return DisplayCurve.getDisplayEdges();}
    //    vector<uint>* getDisplayFaces(){return DisplayCurve.getDisplayFaces();}
    //    vector<unsigned char>* getDisplayColor(){return DisplayCurve.getDisplayColor();}

    //        vector<double>* getDisplayVertices(){return &display_vertices;}
    //        vector<double>* getDisplayVerticesNormal(){return &display_vnormals;}
    //        vector<uint>* getDisplayEdges(){return &display_edges;}
    //        vector<uint>* getDisplayFaces(){return &display_faces;}
    //        vector<unsigned char>* getDisplayColor(){return &v_color;}

    display_vertices.clear();
    display_vnormals.clear();
    display_edges.clear();
    display_faces.clear();
    display_colors.clear();
    if(isShowExist){


        display_vertices = *(DisplayCurve.getDisplayVertices());
        display_vnormals = *(DisplayCurve.getDisplayVerticesNormal());
        display_edges = *(DisplayCurve.getDisplayEdges());
        display_faces = *(DisplayCurve.getDisplayFaces());
        display_colors = *(DisplayCurve.getDisplayColor());

    }

    if(isShowMerge){

        uint offset = display_vertices.size()/3;
        display_vertices.insert(display_vertices.end(), MergeCurve.getDisplayVertices()->begin(), MergeCurve.getDisplayVertices()->end());
        display_vnormals.insert(display_vnormals.end(), MergeCurve.getDisplayVerticesNormal()->begin(), MergeCurve.getDisplayVerticesNormal()->end());


        auto &mfaces  = *(MergeCurve.getDisplayFaces());
        for(auto &a:mfaces)display_faces.push_back(a+offset);
        auto &medges = *(MergeCurve.getDisplayEdges());
        for(auto &a:medges)display_edges.push_back(a+offset);

        uchar mcolor[4] = {255,0,0,255};
        int nc = MergeCurve.getDisplayColor()->size()/4;
        for(int i=0;i<nc;++i)for(int j=0;j<4;++j)display_colors.push_back(mcolor[j]);


    }

    if(!isEnableFrames)return;

    {
        FrameCurve.BuildDisplay(isrescale,isuseball);

        uint offset = display_vertices.size()/3;
        display_vertices.insert(display_vertices.end(), FrameCurve.getDisplayVertices()->begin(), FrameCurve.getDisplayVertices()->end());
        display_vnormals.insert(display_vnormals.end(), FrameCurve.getDisplayVerticesNormal()->begin(), FrameCurve.getDisplayVerticesNormal()->end());


        auto &mfaces  = *(FrameCurve.getDisplayFaces());
        for(auto &a:mfaces)display_faces.push_back(a+offset);
        auto &medges = *(FrameCurve.getDisplayEdges());
        for(auto &a:medges)display_edges.push_back(a+offset);

        uchar mcolor[4] = {0,255,0,255};
        int nc = FrameCurve.getDisplayColor()->size()/4;
        for(int i=0;i<nc;++i)for(int j=0;j<4;++j)display_colors.push_back(mcolor[j]);

        cout<<"build display"<<endl;

    }

    //cout<<"display_colors.n_vertices "<<display_colors.size()<<endl;
    //cout<<"MergeCurve.getDisplayColor().n_vertices "<<MergeCurve.getDisplayColor()->size()<<endl;
}




}//n_rf
