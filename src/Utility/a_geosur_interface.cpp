
#include "geo_sur.h"
#include <stdio.h>
#include<iostream>
#include<fstream>
#include<limits>
#include <functional>
//#include <eigen3/Eigen/Geometry>
#include<set>
#include<sys/time.h>
namespace n_rf {


int Surface::Initialize(string filename){
    clearup();
    bool rere = false;
    if(0){

        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/hand_in_tri.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/botijo_in_tri.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/elk_in_tri.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/fertility_tri.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/rockerarm_tri.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/cube2.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/torus.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/bunny.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/bunnyrefine.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/car.off");


        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/bunny.eobj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/doubletorus.eobj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/Sphere.eobj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/bbbnnnnnn.eobj");


        //ReadFile("/Users/Research/Geometry/RMFpro/Data/inputmodels/eight.obj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/chosen/bunny_remesh.obj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/chosen/sphere.obj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/chosen/221.off");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/chosenOpt/eight.objx");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/chosnOptPrin/venus.objx");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/inputmodels/armchair.obj");

        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/doubletorus.objx");
        //ReadFile("/Users/Research/Geometry/stripe/output.objx");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/fertility_tri3.objx");

        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/ring3.obj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/Mesh/bunny.eobj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/TestOptPrin/hand.objx");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/TestObj/bunny.obj");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/smipic/Eigen/armchair.objx");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/smipic/GlobalOpt/torus.objx");
        //ReadFile("/Users/Research/Geometry/RMFpro/Data/smipic/rockerarm/r40.objx");
        //GenerateToy();

        ReadFile("/Users/Research/Geometry/MM/result/toy_3col/OutputSuf.suf4");

        rere = true;
    }
    else{
        rere = ReadFile(filename);
    }
    if(!rere)return -1;
    setparameters();
    ReScale_uniform(1.0);
    isbuilddisp = false;
    cout<<"read finish!"<<endl;


    BuildNeighborTable();
    ReOrientFaces();
    ComputeEdgefromFace();

    ComputeEdgeLength();
    ComputeDefectAngle();
    ComputeFaceNormal();

    BuildFacesCenter();
    ComputeArea();


    BuildPickID();


    reinitflag = true;
    isbuilddisp = false;





//    for(int i =0;i<n_vertices;++i){
//        auto p_v = v_begin(i);
//        for(int j=0;j<3;++j)cout<<p_v[j]<<' ';cout<<endl;
//    }
//    for(int i =0;i<n_faces;++i){
//        auto p_v = fv_begin(i);
//        for(int j=0;j<3;++j)cout<<p_v[j]<<' ';cout<<endl;
//    }
    cout<<"Mesh Initialized!"<<endl;

    return 0;




}


void Surface::BuildDisplay(infoSurfDisp info){


    BuildDisplay(colormethod,colordegree, info.isField,info.isNormal,info.isSurface, info.isWire,info.isSingularity,info.ismark,info.length,info.width,info.upnormal);


}



bool Surface::ReadFile(string filename){


    string filepath;
    string modname;
    string extname;

    SplitFileName(filename,filepath,modname,extname);
    cout<<modname<<' '<<extname<<' '<<filepath<<endl;

    bool issuc = false;

    if(extname==".suf"){issuc = readSufFile(filename);}
    else issuc = readfile(filename);

    if(issuc){
        modelname = modname;
        prepath = filepath;
        ext=extname;
    }

    return issuc;
}




bool Surface::SaveInterface(string filename){
    //filename="/Users/Research/Geometry/stripe/outputf.objx";
    //filename="/Users/Research/Geometry/RMFpro/Data/smipic/Eigen/outputf.objx";
    //saveObjx(filename);
    filename="/Users/Research/Geometry/MM/QTProgram/222.obj";
    saveObjFile(filename);

    //filename="/Users/Research/Geometry/CrestCODE/PLY/cube.ply2";
    //savePLY2File(filename);

    //filename="/Users/Research/Geometry/RMFpro/Data/Test/cube.objx";
    //saveObjx(filename);
    //saveObjFile(filename);

}

void Surface::testSlicer(int thres){
    //JacobianEnergy();
    //double n_thres =(*(max_element(v_jacobian_energy.begin(),v_jacobian_energy.end())))/((double)thres/2.);
    //    cout<<*(max_element(v_jacobian_energy.begin(),v_jacobian_energy.end()))<<endl;
    //WeightColor(n_thres,v_jacobian_energy,weighted_color);
    //double n_thres =((double)thres/500.)*(*(max_element(vertices_curvature.begin(),vertices_curvature.end())));
    //WeightColor(thres/500.0,vertices_curvature,weighted_color);
    cout<<"thres: "<<thres<<endl;

}





}//n_rf
