#include "geo_curv.h"
#include "readers.h"
#include <stdio.h>
#include<iostream>
#include<fstream>
#include<limits>
#include <functional>
#include <eigen3/Eigen/Geometry>
#include<set>
namespace n_rf {

bool Curve::isload = false;
Mesh Curve::sphere;
Mesh Curve::cylinder;
Mesh Curve::cone;
Curve::Curve():n_vertices(0),n_edges(0),isbuild(false),isSettangent(false),isBuildFrenet(false)
{

    load();

}
void Curve::reset(){
    vertices.clear(); edges.clear(); tangent.clear();

    display_edges.clear();display_vertices.clear();v_color.clear();display_faces.clear();

    n_vertices = n_edges = 0;


    isbuild = false;
    isSettangent = false;
    isBuildFrenet = false;

}
void Curve::setparameters(){
    n_vertices = vertices.size()/3;
    n_edges = edges.size()/2;
}


bool Curve::ImportCurve(vector<double>& vertices_in, vector<uint> &edges_in){

    reset();

    //cout<<"importCurve"<<endl;
    vertices = vertices_in;
    edges = edges_in;
    //disscale = 1.0;

    setparameters();
    ValidationCheck();

    BuildEdges(false);

}

bool Curve::ImportCurve(vector<double>& vertices_in, vector<uint>& edges_in,double lscale, double *pcenter){
    ImportCurve(vertices_in, edges_in);
    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-pcenter[i])*lscale;

}

bool Curve::AddCurve(Curve& a){
    setparameters();
    int offset = n_vertices;
    for(auto v:a.vertices)vertices.push_back(v);

    for(auto e:a.edges)edges.push_back(e+offset);
    setparameters();

}
bool Curve::ImportCurve(string filename){

    vector<double>vertices_in;
    vector<uint>edges_in;
    ifstream fin(filename.data());
    if(fin.fail())cout<<"error"<<endl;

    int numv,nume;
    fin>>numv;fin>>nume;
    vertices_in.resize(numv*3);
    edges_in.resize(nume*2);
    for(auto &a:vertices_in)fin>>a;
    for(auto &a:edges_in)fin>>a;
    ImportCurve(vertices_in,edges_in);

}

bool Curve::ReadCurve(string filename){
    cout<<"Reading: "<<filename<<endl;
    ifstream fin(filename.data());
    if(fin.fail())cout<<"error"<<endl;
    reset();
    fin>>n_vertices;fin>>n_edges;
    if(n_vertices==n_edges)isloop = true;
    else if(n_vertices==n_edges+1)isloop = false;
    else {cout<<"Not a curve!"<<endl;return false;}

    cout<<"read: "<<n_vertices<<' '<<n_edges<<endl;
    vertices.resize(n_vertices*3);
    for(int i =0;i<n_vertices*3;++i)fin>>vertices[i];

    edges.resize(n_edges*2);
    bool hasedges = true;
    for(int i =0;i<n_edges*2;++i){
        fin>>edges[i];
        if(fin.eof()){hasedges = false;break;}
    }
    fin.close();
    if(hasedges)SortLoopEdges();
    BuildEdges();
    tangent.clear();

    setparameters();

    EstimateFrenetFrame();
    //BuildEdges(false);

    //BuildDisplay(false);

    return true;
}
bool Curve::ReadCurve(ifstream &fin, int p0, int p1, int p2){
    n_vertices = p0;
    if(p1==0)isloop = false;
    else if(p1==1)isloop = true;
    else {cout<<"Not a curve!"<<endl;return false;}

    vertices.resize(n_vertices*3);
    for(int i =0;i<n_vertices*3;++i)fin>>vertices[i];


    edges.clear();
    for(int i =0;i<n_vertices-1;i++){
        edges.push_back(i);
        edges.push_back(i+1);
    }
    if(isloop){edges.push_back(n_vertices-1);edges.push_back(0);}


    tangent.clear();
    //BuildDisplay();
    setparameters();
    EstimateFrenetFrame();
    BuildEdges(false);



    return true;


}
bool Curve::SaveCurve(string filename){
    ofstream fout(filename.data());
    if(fout.fail())cout<<"error"<<endl;
    fout<<n_vertices<<' '<<n_edges<<endl;
    for(int i =0; i<n_vertices;++i){
        auto p_v = v_begin(i);
        fout<<p_v[0]<<' '<<p_v[1]<<' '<<p_v[2]<<endl;

    }
    fout.close();
    return true;
}

bool Curve::load(){
    if(!isload){
//        cout<<"a"<<endl;
//        sphere.createToy(3);
//        cout<<"2"<<endl;
//        cylinder.createToy(2);

//        cone.createToy(1);
//        cone.ReScale(2.5,2.5,0.5);
//        cylinder.ReScale(0.5,0.5,1.0);

//        sphere.BuildNeighborTable();
//        cylinder.BuildNeighborTable();
//        cone.BuildNeighborTable();

//        sphere.ComputeFaceNormal(true);
//        cylinder.ComputeFaceNormal(true);
//        cone.ComputeFaceNormal(true);


        isload = true;
    }
    //sphere.readfile(string("/Users/Research/Geometry/RMFpro/cube.off"));
}

void Curve::BuildEdges(bool isseq){

    if(isseq){
        edges.clear();
        for(int i =0;i<n_vertices-1;i++){
            edges.push_back(i);
            edges.push_back(i+1);
        }
        if(isloop){edges.push_back(n_vertices-1);edges.push_back(0);n_edges = n_vertices;}
        else n_edges = n_vertices-1;
    }

    edge_len.resize(n_edges);
    inverse_edge_len.resize(n_edges);
    double evec[3];
    for(uint i = 0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        vectorize(p_ev[0],p_ev[1],evec);
        edge_len[i]=normVec(evec);
        inverse_edge_len[i] = 1./edge_len[i];

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


    //for(a:vertices2edges_accumulate)cout<<a<<' ';cout<<endl;
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

    vertices2edgesPN.resize(vertices2edges.size(),false);
    for(uint i=0;i<n_vertices;++i){
        auto p_vepn = vepn_begin(i);
        for(auto p_ve = ve_begin(i);p_ve!=ve_end(i);++p_ve,++p_vepn){
            if((*ev_begin(*p_ve))==i)(*p_vepn)=true;
            else (*p_vepn)=false;

        }
    }


    vertices2vertices.resize(vertices2edges.size());
    for(int i=0;i<n_vertices;++i){
        auto p_vepn = vepn_begin(i);
        auto p_vv = vv_begin(i);
        for(auto p_ve = ve_begin(i);p_ve!=ve_end(i);++p_ve,++p_vepn,++p_vv)
            if(*p_vepn)(*p_vv) = ev_begin(*p_ve)[1];
            else (*p_vv) = ev_begin(*p_ve)[0];
    }

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

    vertices_len.resize(n_vertices);
    for(int i=0;i<n_vertices;++i){
        auto p_ve = ve_begin(i);
        auto venum = ve_num(i);
        if(venum==2)vertices_len[i] = (elen_begin(p_ve[0])+elen_begin(p_ve[1]))/2.;
        else if(venum==1)vertices_len[i] = elen_begin(p_ve[0]);
    }




}

bool Curve::ValidationCheck(){
    for(auto a:edges)assert(a<n_vertices);
    if(0){
        double baa=0;
        for(int i=0;i<n_edges;++i){
            auto p_ev = ev_begin(i);
            if((baa = VerticesDistance(p_ev[0],p_ev[1]))<1e-6){
                cout<<"curve degenerate edge: "<<baa<<' '<<p_ev[0]<<' '<<p_ev[1]<<endl;
            }
        }
        cout<<"Curve validation check"<<endl;
    }

}


void Curve::EliminateDuplication(double thres){
    setparameters();
    vector<double>newvertices;
    vector<bool>iskeep(n_vertices,true);


    for(int i=0;i<n_vertices-1;++i){
        if(VerticesDistance(i,i+1)<thres)iskeep[i+1]=false;
    }
    if(isloop)if(VerticesDistance(0,n_vertices-1)<thres)iskeep[n_vertices-1]=false;
    int count = 0;
    for(int i=0;i<n_vertices-1;++i){
        if(iskeep[i]){
            count++;
            auto p_o = v_begin(i);
            newvertices.push_back(p_o[0]);
            newvertices.push_back(p_o[1]);
            newvertices.push_back(p_o[2]);
        }
    }
    vertices=newvertices;
    cout<<"Duplication "<<n_vertices<<' '<<count<<endl;
    setparameters();
    BuildEdges();



}


bool Curve::isBuildDisplay(){return isbuild;}





bool Curve::SortLoopEdges(){

    auto old_edges = edges;
    auto old_vertices = vertices;
    vector< vector<int> >neighborhood(n_vertices);
    for (int i = 0; i < edges.size(); i += 2){
        neighborhood[edges[i]].push_back(edges[i + 1]);
        neighborhood[edges[i + 1]].push_back(edges[i]);
    }
    if(isloop)for (auto &v_v : neighborhood)if (v_v.size() != 2){ cout << "Not a closed polygon!" << endl; return false; }
    //for (auto v_v : neighborhood)if (v_v[0] == v_v[1]){ cout << "Not a closed polygon!" << endl; }
    //for (auto v_v : neighborhood){ cout << v_v[0]<<"  "<< v_v[1] << endl; }
    vector<int>orderlist(n_vertices, -1);

    //map<pair<uint, uint>, uint>edges_index;
    vector< vector<int> >edges_index2(n_vertices, vector<int>(n_vertices));

    int seed_index = 0, pre_seed_index = -1, order = 0;
    if(!isloop){
        for(int i=0;i<neighborhood.size();++i)if (neighborhood[i].size() == 1){seed_index=i;break;}
    }
    for (int i = 0; i < n_vertices; ++i){
        orderlist[seed_index] = order;
        auto& v_v = neighborhood[seed_index];
        for (auto next_seed : v_v)if (next_seed != pre_seed_index){
            //edges_index[make_pair(seed_index, next_seed)] = order;
            //edges_index[make_pair(next_seed, seed_index)] = order;
            edges_index2[next_seed][seed_index] = order;
            edges_index2[seed_index][next_seed] = order;
            pre_seed_index = seed_index;
            seed_index = next_seed;
            break;
        }
        ++order;
    }
    //cout << order << endl;
    //for (auto new_order : orderlist)cout << new_order << endl;
    for (auto new_order : orderlist)if (new_order == -1){ cout << "Not a closed polygon!" << endl; return false; }
    for (int i = 0; i < n_vertices; i++){
        for (int j = 0; j < 3; ++j)
            vertices[3 * orderlist[i] + j] = old_vertices[3 * i + j];
    }
    for (int i = 0; i < edges.size(); i += 2){
        //auto e_key = make_pair(old_edges[i], old_edges[i + 1]);
        //auto new_index = 2 * edges_index[e_key];
        auto new_index = 2 * edges_index2[old_edges[i]][old_edges[i + 1]];
        auto p1 = orderlist[old_edges[i]];
        auto p2 = orderlist[old_edges[i + 1]];
        if (p1 < p2){ edges[new_index] = p1; edges[new_index + 1] = p2; }
        else { edges[new_index] = p2; edges[new_index + 1] = p1; }
        if (edges[new_index] == 0 && edges[new_index + 1] == n_vertices - 1)swap(edges[new_index], edges[new_index + 1]);
    }

    BuildEdges(false);

}

bool Curve::EstimateFirstDerivativeO5(vector<double>& invect, vector<double>& derivative, int interval,bool isloop){
    if(invect.size()==0 || invect.size()%interval!=0)return false;

    derivative.resize(invect.size());

    auto singleinterval = [this](vector<double>&invec, vector<double>&outvec,bool isloop){
        int vecN = invec.size();
        outvec.resize(vecN);
        if(isloop){
            outvec[0] = 1/12.0*invec[vecN-2] - 8/12.0*invec[vecN-1] + 8/12.0*invec[1] - 1/12.0*invec[2];
            outvec[1] = 1/12.0*invec[vecN-1] - 8/12.0*invec[0] + 8/12.0*invec[2] - 1/12.0*invec[3];
            outvec[vecN-1] = 1/12.0*invec[vecN-3] - 8/12.0*invec[vecN-2] + 8/12.0*invec[0] - 1/12.0*invec[1];
            outvec[vecN-2] = 1/12.0*invec[vecN-4] - 8/12.0*invec[vecN-3] + 8/12.0*invec[vecN-1] - 1/12.0*invec[0];
        }else{
            outvec[0] = -25/12.0*invec[0] +48/12.0*invec[1]-36/12.0*invec[2]+16/12.0*invec[3]-3/12.0*invec[4];
            outvec[1] = (-3/12.0*invec[0] -10/12.0*invec[1]+18/12.0*invec[2]-6/12.0*invec[3]+1/12.0*invec[4]);
            outvec[vecN-1] = -(-25/12.0*invec[vecN-1] +48/12.0*invec[vecN-2]-36/12.0*invec[vecN-3]+16/12.0*invec[vecN-4]-3/12.0*invec[vecN-5]);
            outvec[vecN-2] = -(-3/12.0*invec[vecN-1] -10/12.0*invec[vecN-2]+18/12.0*invec[vecN-3]-6/12.0*invec[vecN-4]+1/12.0*invec[vecN-5]);
        }

        for(int i=2;i<vecN-2;i++){
            outvec[i] = 1/12.0*invec[i-2] - 8/12.0*invec[i-1] + 8/12.0*invec[i+1] - 1/12.0*invec[i+2];
        }

    };
    int vecNN = invect.size()/interval;
    vector<double>invec(vecNN);
    vector<double>outvec(vecNN);
    for(int i=0;i<interval;++i){
        for(int j=0;j<vecNN;++j)invec[j] = invect[i+j*interval];
        singleinterval(invec,outvec,isloop);
        for(int j=0;j<vecNN;++j)derivative[i+j*interval] = outvec[j];
    }

    return true;
}
void Curve::BuildTangent(){

    if(n_vertices==0||n_edges==0)return;
    tangent.resize(n_vertices*3);
    EstimateFirstDerivativeO5(vertices,tangent,3,isloop);
    for(int i=0;i<n_vertices;++i)normalize(t_begin(i));
    for(int i=0;i<tangent.size();++i)if(tangent[i]!=tangent[i])exit(121);
    isSettangent = true;
}
void Curve::EstimateFrenetFrame(){

    if(n_vertices==0||n_edges==0)return;

    //BuildTangent();

    tangent.resize(n_vertices*3);
    EstimateFirstDerivativeO5(vertices,tangent,3,isloop);
    linespeed.resize(n_vertices);
    for(int i=0;i<n_vertices;++i)linespeed[i] = len(t_begin(i));
    //for(int i=0;i<n_vertices;++i)cout<<linespeed[i]<<' ';
    //cout<<endl;



    for(int i=0;i<n_vertices;++i)normalize(t_begin(i));
    for(int i=0;i<tangent.size();++i)if(tangent[i]!=tangent[i])exit(121);
    isSettangent = true;

    Frenetnormal.resize(n_vertices*3);
    EstimateFirstDerivativeO5(tangent,Frenetnormal,3,isloop);

    curvature.resize(n_vertices);
    for(int i=0;i<n_vertices;++i)curvature[i] = len(nor_begin(i))/linespeed[i];

    normCurvature.resize(n_vertices);
    double cmaxC = 3.5/(*max_element(curvature.begin(),curvature.end()));
    for(int i=0;i<n_vertices;++i)normCurvature[i] =min(1., curvature[i] *cmaxC);
    CurveScalarFunctionSmoothing(normCurvature,20);

    //for(int i=0;i<n_vertices;++i)cout<<curvature[i]<<' ';
    //cout<<endl;

    for(int i=0;i<n_vertices;++i)normalize(nor_begin(i));
    for(int i=0;i<tangent.size();++i)if(Frenetnormal[i]!=Frenetnormal[i])exit(121);


    Frenetbinormal.resize(n_vertices*3);
    for(int i=0;i<n_vertices;++i)cross(t_begin(i),nor_begin(i),binor_begin(i));
    for(int i=0;i<n_vertices;++i)normalize(binor_begin(i));


    vector<double>dbinormal(n_vertices*3);
    torsion.resize(n_vertices);
    EstimateFirstDerivativeO5(Frenetbinormal,dbinormal,3,isloop);

    for(int i=0;i<n_vertices;++i)torsion[i] = -dot(nor_begin(i),&(dbinormal[i*3]));
    for(int i=0;i<n_vertices;++i)torsion[i] /=linespeed[i];

    if(1)for(int i=0;i<n_vertices;++i)inversevec(binor_begin(i),binor_begin(i));
    //for(int i=0;i<n_vertices;++i)cout<<torsion[i]<<' ';
    //cout<<endl;





    isBuildFrenet = true;
    isSettangent = true;
    //cout<<"EstimateFrenetFrame"<<endl;
    //for(auto j:Frenetbinormal)cout<<j<<' ';cout<<endl;


}
void Curve::RescaleUniform(){
    double xmin = 999e8, ymin=999e8,zmin=999e8,xmax=-999e8,ymax=-999e8,zmax=-999e8;
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
    //double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    vector<double>dis(3);
    dis[0] = (xmax-xmin)/2;dis[1] = (ymax-ymin)/2;dis[2] = (zmax-zmin)/2;
    //cout<<centers[0]<<"  "<<centers[1]<<"  "<<centers[2]<<"   "<<largestdis<<endl;

    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-centers[i])/dis[i];

}

void Curve::Rescale(double &lscale, double *pcenter){

    double xmin = 999e8, ymin=999e8,zmin=999e8,xmax=-999e8,ymax=-999e8,zmax=-999e8;
    for(int i = 0; i<3; i++)pcenter[i] = 0;
    for (int i = 0; i<n_vertices; i++){
        auto point = v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);
        pcenter[0]+=point[0];pcenter[1]+=point[1];pcenter[2]+=point[2];
        //if(centers[0]!=centers[0]){cout<<centers[0]<<"  "<<point[0]<<" "<<i/dim<<endl;}
    }
    pcenter[0]/=n_vertices;pcenter[1]/=n_vertices;pcenter[2]/=n_vertices;
    lscale = 1/(max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2);


    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-pcenter[i])*lscale;




}
void Curve::CurveScalarFunctionSmoothing(vector<double>&scalarF,int maxiter){


    vector<double>buffer1(scalarF.size());
    vector<double>buffer2(scalarF.size());
    vector<double>*pre = &buffer1,*cur = &buffer2;
    *pre = scalarF;

    double kaipa = 0.1,lamda = 0.5;
    auto gama = 1 / (kaipa - 1 / lamda);

    auto oneiter = [this](double coef,vector<double>*pre,vector<double>*cur){
        double ocoef = 1- coef;
        for(int i =0;i<n_vertices;++i){
            double aaa = 0;
            for(auto p_vv = vv_begin(i);p_vv!=vv_end(i);++p_vv){
                aaa+=pre->at(*p_vv);
            }
            aaa/=(vv_num(i));
            cur->at(i) = ocoef * pre->at(i) + coef * aaa;
        }
    };

    for(int iter=0;iter<maxiter;++iter){

        oneiter(lamda,pre,cur);
        swap(pre,cur);
        oneiter(gama,pre,cur);
        swap(pre,cur);

    }
    scalarF=(*pre);




}

class Triple{
public:
    float cost;
    int kid, father;
    Triple(float c, int k, int f) :cost(c), kid(k), father(f){}
    Triple() :cost(-1), kid(-1), father(-1){}
    bool operator < (const Triple& T)const{ return this->cost < T.cost; }
    bool operator >(const Triple& T)const{ return this->cost > T.cost; }
};

void Curve::CurveNetAnalayze(vector<uint>&edgesMat, vector<int> &edges2Cells,int _nCell){


    vector<int>v2vnum(n_vertices);

    for(int i=0;i<n_vertices;++i)v2vnum[i] = vv_num(i);

    vector<bool>keypoints(n_vertices,false);
    vector< int >keypoint_inverse_ids(n_vertices,-1);
    keyVer.clear();
    int kind = 0;
    for(int i=0;i<n_vertices;++i)if(v2vnum[i]>2){
        keyVer.push_back(i);keypoints[i] = true;
        keypoint_inverse_ids[i] = kind++;
    }

    int n_c_edges = 0;
    vector<bool>visited_edges(n_edges, false);
    //vector< int >superedges;
    vector<bool>visited_vertices(n_vertices, false);

    while (true){
        //bool a_end = true, b_end = true;
        set<uint>eL;
        set<int>eC;
        int seed_edge = -1;
        for (uint i = 0; i < n_edges; ++i)if (!visited_edges[i]){
            seed_edge = i;
            eL.insert(edgesMat[seed_edge*2]);eL.insert(edgesMat[seed_edge*2+1]);
            eC.insert(edges2Cells[seed_edge*2]);eC.insert(edges2Cells[seed_edge*2+1]);
            break;
        }
        if (seed_edge == -1)break;
        vector<uint>c_edge_unit;
        c_edge_unit.push_back(seed_edge);
        visited_edges[seed_edge] = true;
        auto p_e = &(edges[seed_edge * 2]);
        if(keypoints[p_e[0]]&&keypoints[p_e[1]]){
            cout<<seed_edge<<endl;
        }
        for (int j = 0; j < 2; ++j){
            auto v_a = p_e[j];
            auto e_a = seed_edge;
            while (!keypoints[v_a] && !visited_vertices[v_a]){
                visited_vertices[v_a] = true;
                for (auto p_v2e = ve_begin(v_a); p_v2e != ve_end(v_a); ++p_v2e){
                    if (*p_v2e == e_a)continue;
                    if(visited_edges[*p_v2e])continue;
                    e_a = *p_v2e;
                    for (int i = 0; i < 2; ++i)if (edges[e_a * 2 + i] != v_a){
                        v_a = edges[e_a * 2 + i]; break;
                    }
                    c_edge_unit.push_back(e_a);
                    visited_edges[e_a] = true;
                    eL.insert(edgesMat[e_a*2]);eL.insert(edgesMat[e_a*2+1]);
                    eC.insert(edges2Cells[e_a*2]);eC.insert(edges2Cells[e_a*2+1]);
                    break;
                }
            }
            superedges.push_back(keypoint_inverse_ids[v_a]);
        }
        seg2edges.push_back(c_edge_unit);
        assert(eL.size()==2);
        assert(eC.size()==2);
        for(auto a:eL)seg2Mat.push_back(a);
        for(auto a:eC)seg2Cells.push_back(a);

    }

    nSeg = seg2edges.size();

    SegNSeg.resize(nSeg);
    vector<set<int>>SegNSegSet(nSeg);
    vector<vector<int>>keyP2Seg(keyVer.size());
    for(int i=0;i<nSeg;++i){
        if(superedges[i*2]==-1){
            assert(superedges[i*2]==superedges[i*2+1]);continue;
        }
        keyP2Seg[superedges[i*2]].push_back(i);
        keyP2Seg[superedges[i*2+1]].push_back(i);
    }

    for(auto &a:keyP2Seg)for(auto b:a)for(auto c:a)if(c!=b)SegNSegSet[b].insert(c);
    for(int i=0;i<nSeg;++i)for(auto a:SegNSegSet[i])SegNSeg[i].push_back(a);








    nCell = _nCell;
    //nCell = *max_element(edges2Cells.begin(),edges2Cells.end())+1;
    cell2Seg.clear();
    cell2Seg.resize(nCell);
    for(int i=0;i<seg2edges.size();++i){
        cell2Seg[seg2Cells[i*2]].push_back(i);
        cell2Seg[seg2Cells[i*2+1]].push_back(i);
    }
    cout << "c_edges: " << seg2edges.size() <<' '<<nCell<<endl;

    for(int i=0;i<nCell;++i){
        cout<<i<<": "<<endl;
        for(auto a:cell2Seg[i])cout<<a<<' ';cout<<endl;
    }

    //for(auto a:seg2Cells)assert(a==2);

    //return;
    seg2ver.clear();
    seg2ver.resize(nSeg);
    for(int i=0;i<nSeg;++i){
        int beVind = edges[2*seg2edges[i][0]];
        if(superedges[i*2]!=-1)beVind = keyVer[superedges[i*2]];

        vector<uint>newOrderEdge;
        vector<uint>newOrderEdgeMat;
        vector<bool>pickV(n_vertices,false);
        vector<bool>lvisited_edges(n_edges,false);
        vector<bool>validE(n_edges,false);

        for(auto a:seg2edges[i])validE[a] = true;


        while(!pickV[beVind]){
            pickV[beVind] =true;
            seg2ver[i].push_back(beVind);
            for (auto p_v2e = ve_begin(beVind); p_v2e != ve_end(beVind); ++p_v2e){
                auto e_a = *p_v2e;
                if(!validE[e_a])continue;
                if(lvisited_edges[e_a])continue;
                lvisited_edges[e_a] = true;
                int k = 0;
                for (; k < 2; ++k)if (edges[e_a * 2 + k] != beVind){
                    beVind = edges[e_a * 2 + k]; break;
                }
                if(k==0){
                    swap(edges[e_a * 2],edges[e_a * 2+1]);
                    newOrderEdge.push_back(edges[e_a * 2+1]);newOrderEdge.push_back(edges[e_a * 2]);
                    newOrderEdgeMat.push_back(edgesMat[e_a * 2+1]);newOrderEdgeMat.push_back(edgesMat[e_a * 2]);
                }else{
                    newOrderEdge.push_back(edges[e_a * 2]);newOrderEdge.push_back(edges[e_a * 2+1]);
                    newOrderEdgeMat.push_back(edgesMat[e_a * 2]);newOrderEdgeMat.push_back(edgesMat[e_a * 2+1]);
                }
                break;
            }

        }
        seg2Mat[2*i] = newOrderEdgeMat[0];
        seg2Mat[2*i+1] = newOrderEdgeMat[1];
        for(int j=0;j<seg2edges[i].size();++j)assert(newOrderEdgeMat[0]==newOrderEdgeMat[2*j]);
        for(int j=0;j<seg2edges[i].size();++j)assert(newOrderEdgeMat[1]==newOrderEdgeMat[2*j+1]);






    }







}
void Curve::CurveNetSimplification(int s,vector<double>&out_ver,vector<uint>&out_edges,vector<uint>&outMat,int deleteIso){

    out_ver.clear();
    out_edges.clear();
    outMat.clear();

    int vind = 0;
    vector<int>vInd(n_vertices,-1);
    for(int i=0;i<nSeg;++i){
        auto p_seg2M = seg2Mat.data()+2*i;
        auto &seg2V = seg2ver[i];
        int stv = seg2V[0];
        if(seg2V.size()<deleteIso*s && superedges[i*2]==-1)continue;
        if(seg2V.size()<s){
            auto& seg2E = seg2edges[i];
            for(auto a:seg2E){
                out_edges.push_back(edges[a*2]);out_edges.push_back(edges[a*2+1]);
                outMat.push_back(p_seg2M[0]);outMat.push_back(p_seg2M[1]);
            }
            continue;
        }
        for(int j=s;j<seg2V.size();j+=s){
            out_edges.push_back(stv);out_edges.push_back(seg2V[j]);
            stv = seg2V[j];
            outMat.push_back(p_seg2M[0]);outMat.push_back(p_seg2M[1]);
        }
        if(superedges[i*2]==-1){
            out_edges.push_back(stv);out_edges.push_back(seg2V[0]);
            outMat.push_back(p_seg2M[0]);outMat.push_back(p_seg2M[1]);

        }else{
            auto endKeyV = seg2V[seg2V.size()-1];
            if(stv!=endKeyV){
                out_edges.push_back(stv);out_edges.push_back(endKeyV);
                outMat.push_back(p_seg2M[0]);outMat.push_back(p_seg2M[1]);
            }

        }
    }
    for(auto &a : out_edges)vInd[a] = 0;
    for(int j=0;j<n_vertices;++j)if(vInd[j]==0){
        auto p_v = v_begin(j);
        for(int i=0;i<3;++i)out_ver.push_back(p_v[i]);
        vInd[j]=vind++;
    }
    for(auto &a : out_edges)a = vInd[a];


}

void BuildRotationMinimizingFramesCore(int numV, double *v, double *tant, double *vrmf, double *initialRV, double *endRV){

    double vec[3];

    auto t_begin = [&tant](int ind){
        return tant+ind*3;
    };

    auto rr_begin = [&vrmf](int ind){
        return vrmf+ind*3;
    };
    auto v_begin = [&v](int ind){
        return v+ind*3;
    };

    if(initialRV==NULL&&endRV==NULL)return;




    if(initialRV!=NULL){
        product(dot(initialRV,t_begin(0)),t_begin(0),vec);
        minusVec(initialRV,vec,vec);
        normalize(vec);
    }else{
        double initialRVFake[3]={1,0,0};
        if(dot(t_begin(0),initialRVFake)>1-5e-2)initialRVFake[1]=1;
        product(dot(initialRV,t_begin(0)),t_begin(0),vec);
        minusVec(initialRV,vec,vec);
        normalize(vec);
    }

    //rvectorRMF.resize(3*n_vertices);
    for(int i =0;i<3;++i)vrmf[i] = vec[i];
    normalize(rr_begin(0));
    double v1[3],c1,rl[3],tmp[3],tl[3],v2[3],c2;
    for(int i = 1;i<numV;++i){
        minusVec(v_begin(i-1),v_begin(i),v1);
        c1 = 2/dot(v1,v1);
        product(c1*dot(v1,rr_begin(i-1)),v1,tmp);
        minusVec(rr_begin(i-1),tmp,rl);
        product(c1*dot(v1,t_begin(i-1)),v1,tmp);
        minusVec(t_begin(i-1),tmp,tl);
        minusVec(t_begin(i),tl,v2);
        c2 = 2/dot(v2,v2);
        product(c2*dot(v2,rl),v2,tmp);
        minusVec(rl,tmp,rr_begin(i));
    }

    for(int i=0;i<numV;++i)normalize(rr_begin(i));

    if(endRV==NULL)return;

    //return;
    double angle =angleNor(rr_begin(numV-1),endRV);
    vector<Eigen::Vector3d>ori_vec(numV);
    vector<Eigen::Vector3d>rot_vec(numV);
    cout<<"angle: "<<angle<<endl;
    //double vec[3];
    cross(rr_begin(numV-1),endRV,vec);
    if(dot(vec,t_begin(numV-1))<0)angle = -angle;
    angle/=(numV-1);

    vector<Eigen::Vector3d>t_vec(numV);
    for(int i = 0;i<numV;++i)
    {
        auto p = rr_begin(i);
        for(int j =0;j<3;++j)ori_vec[i](j) = p[j];
    }
    for(int i = 0;i<numV;++i)
    {
        auto p = t_begin(i);
        for(int j =0;j<3;++j)t_vec[i](j) = p[j];
    }

    if(initialRV!=NULL){
        for(int i = 0;i<numV;++i){
            Eigen::AngleAxisd t(angle*i,t_vec[i]);
            rot_vec[i] = t*ori_vec[i];
        }
    }else{
        for(int i = 0;i<numV;++i){
            Eigen::AngleAxisd t(angle*(numV-i-1),t_vec[i]);
            rot_vec[i] = t*ori_vec[i];
        }
    }

    for(int i = 0;i<numV;++i)
    {
        auto p = rr_begin(i);
        for(int j =0;j<3;++j)p[j] = rot_vec[i](j);
    }

}
void Curve::CurveNormalParalellTransport(vector<double>&oriNor,vector<double>&rmf){


    rmf.resize(oriNor.size());
    auto p_outrmf = rmf.data();
    for(int i=0;i<nSeg;++i){
        vector<double>segverices;
        vector<double>segtangent;
        vector<double>segnor;

        auto &seg2V = seg2ver[i];

        for(int j=0;j<seg2V.size();j++){
            auto p_v = v_begin(seg2V[j]);
            for(int k=0;k<3;++k)segverices.push_back(p_v[k]);
        }
        bool isloopseg = superedges[i*2]==-1;

        if(0)EstimateFirstDerivativeO5(segverices,segtangent,3,isloopseg);
        else{
            segtangent.resize(segverices.size());
            auto p_tang = segtangent.data();
            for(int j=0;j<seg2V.size()-1;j++){
                minusVec(v_begin(seg2V[j+1]),v_begin(seg2V[j]),p_tang+j*3);
            }
            copyVec(p_tang+(seg2V.size()-2)*3,p_tang+(seg2V.size()-1)*3);
            for(int j=0;j<seg2V.size();j++)normalize(p_tang+j*3);
        }

        segnor.resize(segverices.size());
        BuildRotationMinimizingFramesCore(seg2V.size(), segverices.data(), segtangent.data(), segnor.data(), oriNor.data()+seg2V[0]*3, oriNor.data()+seg2V[seg2V.size()-1]*3);



        auto p_rmf = segnor.data();
        auto p_tang = segtangent.data();
        for(int j=0;j<seg2V.size();j++){
            if(1)copyVec(p_rmf+j*3,p_outrmf+seg2V[j]*3);
            else copyVec(p_tang+j*3,p_outrmf+seg2V[j]*3);
        }

    }




}

int PartialCurve::ImportCurve(vector<double>& vertices_in, vector<uint>& edges_in,vector<uint> &edgesMat,vector<int>&edges2Cells,int nCells){


    mixCurve.ImportCurve(vertices_in,edges_in);
    cout<<"22"<<endl;
    mixCurve.CurveNetAnalayze(edgesMat,edges2Cells,nCells);


    return Processing();



}
void PartialCurve::Relocation(int visitStep){

    if(visitStep<0||visitStep>=mixCurve.nCell)return;
    if(mixCurve.n_vertices==0||mixCurve.n_edges==0)return;
    DisplayCurve.reset();

    auto &showSeq = disSeg[visitStep];

    vector<int>newVind(mixCurve.n_vertices,-1);
    vector<uint>newedges;
    vector<double>newVpos;

    for(auto a:showSeq)for(auto b:mixCurve.seg2edges[a]){
        auto p_ev = mixCurve.ev_begin(b);
        for(int j=0;j<2;++j){
            newedges.push_back(p_ev[j]);
            newVind[p_ev[j]] = 0;
        }
    }
    int vind = 0;
    for(int i=0;i<mixCurve.n_vertices;++i)if(newVind[i]==0){
        auto p_v = mixCurve.v_begin(i);
        newVind[i] = vind;++vind;
        for(int j=0;j<3;++j)newVpos.push_back(p_v[j]);
    }
    for(auto &a:newedges)a = newVind[a];

    DisplayCurve.ImportCurve(newVpos,newedges);

    /*-----------------------------------------*/
    /*-----------------------------------------*/

    newedges.clear();newVpos.clear();
    newVind.clear();newVind.resize(mixCurve.n_vertices,-1);
    auto &mergeSeq = dismergeSeg[visitStep];
    for(auto a:mergeSeq)for(auto b:mixCurve.seg2edges[a]){
        auto p_ev = mixCurve.ev_begin(b);
        for(int j=0;j<2;++j){
            newedges.push_back(p_ev[j]);
            newVind[p_ev[j]] = 0;
        }
    }
    vind = 0;
    for(int i=0;i<mixCurve.n_vertices;++i)if(newVind[i]==0){
        auto p_v = mixCurve.v_begin(i);
        newVind[i] = vind;++vind;
        for(int j=0;j<3;++j)newVpos.push_back(p_v[j]);
    }
    for(auto &a:newedges)a = newVind[a];

    MergeCurve.ImportCurve(newVpos,newedges);
    //cout<<"MergeCurve.n_vertices "<<MergeCurve.n_vertices<<endl;
    //cout<<"DisplayCurve.n_vertices "<<DisplayCurve.n_vertices<<endl;


    /*-----------------------------------------*/
    /*-----------------------------------------*/

    if(!isEnableFrames)return;
    int nvf = frameVs.size()/3;
    newedges.clear();newVpos.clear();
    newVind.clear();newVind.resize(nvf,-1);
    vector<bool>enableF(frameF2Es.size(),false);
    for(int i=0;i<visitStep+1;++i){
        int vic = mergingorder[i];
        for(auto a:frameCell2Fs[vic])for(auto b:frameF2Es[a])newVind[b] = 0;
        for(auto a:frameCell2Fs[vic])enableF[a]=true;
    }
    vind = 0;
    for(int i=0;i<nvf;++i)if(newVind[i]==0){
        auto p_v = &(frameVs[3*i]);
        newVind[i] = vind;++vind;
        for(int j=0;j<3;++j)newVpos.push_back(p_v[j]);
    }

    for(int i=0;i<frameF2Es.size();++i)if(enableF[i]){
        for(auto b:frameF2Es[i])newedges.push_back(newVind[b]);
    }

    FrameCurve.ImportCurve(newVpos,newedges);




}
int PartialCurve::Processing(){

    //random_shuffle();

    //    vector< vector<uint> >seg2edges;
    //    vector< uint >superedges;
    //    vector< uint >superedgesMat;

    //    vector< vector<uint> >cell2SuperEdges;
    //    vector< uint >seg2Cells;

    cout<<"Processing"<<endl;
    int nCell = mixCurve.nCell;
    int nSeg = mixCurve.nSeg;
    auto &cell2Seg = mixCurve.cell2Seg;
    auto &seg2Cells = mixCurve.seg2Cells;

    mergingorder.resize(nCell);
    for(int i=0;i<nCell;++i)mergingorder[i] = i;
    //random_shuffle(mergingorder.begin(),mergingorder.end());

    mergingorder[nCell-1] = 0;
    for(int i=0;i<nCell-1;++i)mergingorder[i] = i+1;

    //    vector<vector<bool>>existSeg;
    //    vector<vector<uint>>disSeg;

    existSeg.resize(nCell,vector<bool>(nSeg,false));
    disSeg.resize(nCell);

    mergeSeg.resize(nCell,vector<bool>(nSeg,false));
    dismergeSeg.resize(nCell);


    for(auto a:cell2Seg[mergingorder[0]])existSeg[0][a] = true;
    disSeg[0] = cell2Seg[mergingorder[0]];
    for(int i=1;i<nCell;++i){
        vector<bool>newInSeg(nSeg,false);
        for(auto a:cell2Seg[mergingorder[i]])newInSeg[a] = true;
        for(int j=0;j<nSeg;++j)existSeg[i][j] = existSeg[i-1][j]^newInSeg[j];
        for(int j=0;j<nSeg;++j)if(existSeg[i][j]) disSeg[i].push_back(j);

        for(int j=0;j<nSeg;++j)mergeSeg[i][j] = existSeg[i-1][j] & newInSeg[j];
        for(int j=0;j<nSeg;++j)if(mergeSeg[i][j]) dismergeSeg[i].push_back(j);
    }

    for(int i=0;i<nCell;++i){

        cout<<mergingorder[i]<<": ";
        for(auto a:disSeg[i])cout<<a<<' '; cout<<endl;

        for(auto a:dismergeSeg[i])cout<<a<<' '; cout<<endl;

    }



    return nCell;




}



void PartialCurve::ReadFrameCurves(string prefix, string fvendfix,string ffendfix, string c2fendfix){
    //    vector<float>frameVs;
    //    vector<vector<int>>frameFs;
    //    vector<vector<int>>frameCell2Fs;
    readVecFile(prefix+fvendfix,frameVs);
    readVVecFile(prefix+ffendfix,frameF2Es);
    readVVecFile(prefix+c2fendfix,frameCell2Fs);

    vector<uint>tmpedges({0,1});
    vector<double>tmpV;
    for(auto a:frameVs)tmpV.push_back(a);
    for(auto a:mixCurve.vertices)tmpV.push_back(a);
    FrameCurve.ImportCurve(tmpV,tmpedges);
    FrameCurve.RescaleUniform();

    memcpy(frameVs.data(),FrameCurve.vertices.data(),frameVs.size()*sizeof(double));
    memcpy(mixCurve.vertices.data(),FrameCurve.vertices.data()+frameVs.size(),mixCurve.vertices.size()*sizeof(double));

    FrameCurve.disscale = 0.2;
    isEnableFrames = true;


}


void CrossSections::RescaleToUniform(){

    double xmin = 10e20, ymin=10e20,zmin=10e20,xmax=-10e20,ymax=-10e20,zmax=-10e20;
    double centers[3] = {0,0,0};
    int n_vertices = 0;
    for(Contour &a:this->CSs){
        n_vertices+=a.vertices.size()/3;
        for (int i = 0,en = a.vertices.size()/3; i<en; i++){
            auto point = a.vertices.data()+i*3;
            xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
            ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
            zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);
            centers[0]+=point[0];centers[1]+=point[1];centers[2]+=point[2];
            //if(centers[0]!=centers[0]){cout<<centers[0]<<"  "<<point[0]<<" "<<i/dim<<endl;}
        }
    }
    centers[0]/=n_vertices;centers[1]/=n_vertices;centers[2]/=n_vertices;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    //cout<<centers[0]<<"  "<<centers[1]<<"  "<<centers[2]<<"   "<<largestdis<<endl;

    double lscale = 1/largestdis;
    double pcenter[3];
    for(int i=0;i<3;++i)pcenter[i] = centers[i];

    for(Contour &a:this->CSs){
        a.plane_para[3] = (a.plane_para[3]+dot(a.plane_para,pcenter))*lscale;
        int n_vertices = a.vertices.size()/3;
        for (int j = 0; j<n_vertices*3; j+=3)
            for(int i = 0; i<3; i++)
                a.vertices[j+i]=(a.vertices[j+i]-pcenter[i])*lscale;

    }



}



}//n_rf

