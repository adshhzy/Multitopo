
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
#include<set>
extern double controlBallsize;

void ThresholdColoring(double degree,uchar *p_rgb);
void ThresholdColoring(double degree,double s,double v,uchar *p_rgb);
namespace n_rf {

void Surface::SetDisplayTransparency(uchar alpha){

    int vnum = display_vcolor.size();
    for(int i=3;i<vnum;i+=4)display_vcolor[i] = alpha;



}

void Surface::BuildDisplay(int colormethod, double colordegree,bool isfield, bool isnormal,bool issurf, bool iswire, bool issingularity, bool ismark, int length, int width,int upnormal){



    if(n_vertices==0)return;
    if(n_faces==0)return;



    double scale = 0.0015 * length;
    double scalevn = 0.0005;
    double scaleFieldupn = 0.0015/30*upnormal;
    double scalew = 0.0045/30*width;

    unsigned char dark[4] = {0,0,0,255};
    unsigned char red[4] = {250,0,0,255};
    unsigned char blue[4] = {0,0,255,255};
    unsigned char gray[4] = {225,225,225,255};
    unsigned char colllour[4];colllour[3] = 255;

    bool isUseDoubleColorSingleSide = colormethod == 0;
    bool isUseDoubleColorDoubleSide = colormethod == 2;
    bool isUseBackSideColorDoubleSide = colormethod == 3;
    ThresholdColoring(colordegree,colllour);
    //cout<<colordegree<<" colordegree: "<<int(colllour[0])<<' '<<int(colllour[1])<<endl;


    //cout<<"weighted_fcolor.size(): "<<weighted_fcolor.size()<<endl;
    if(isFaceRenderMode != isfield){isbuilddisp=false;isFaceRenderMode = isfield;}
    if(ismarkC != ismark){isresetColor=true;ismarkC = ismark;}
    if(!isbuilddisp || reinitflag){

        if(isFaceRenderMode){
            //cout<<"mmmmmmmmmmmmmm"<<endl;
            display_vertices.resize(n_faces*9*2);
            display_normal.resize(n_faces*9*2);
            auto p_dnor = display_normal.data();
            auto p_fnd = faces_normal.data();
            if(displayInvertfaceNormal.size()==faces_normal.size())p_fnd = displayInvertfaceNormal.data();
            for(int i=0;i<n_faces;++i){
                auto p_fn= p_fnd+i*3;
                auto p_cdnor = p_dnor+i*9;
                for(int j=0;j<3;++j)copyVec(p_fn,p_cdnor+j*3);
            }
            p_dnor = display_normal.data() + n_faces*9;
            for(int i=0;i<n_faces;++i){
                auto p_fn= p_fnd+i*3;
                auto p_cdnor = p_dnor+i*9;
                for(int j=0;j<3;++j)inversevec(p_fn,p_cdnor+j*3);
            }
            auto p_dv = display_vertices.data();
            for(int i=0;i<n_faces;++i){
                auto p_fv= fv_begin(i);
                auto p_cdv = p_dv+i*9;
                auto p_fn= p_fnd+i*3;
                //for(int j=0;j<3;++j)copyVec(v_begin(p_fv[j]),p_cdv+j*3);
                for(int j=0;j<3;++j)weightedAddVec(0.00001,1.,p_fn,v_begin(p_fv[j]),p_cdv+j*3);
            }
            p_dv = display_vertices.data() + n_faces*9;
            for(int i=0;i<n_faces;++i){
                auto p_fv= fv_begin(i);

                auto p_cdv = p_dv+i*9;
                //for(int j=0;j<3;++j)copyVec(v_begin(p_fv[j]),p_cdv+j*3);
                auto p_fn= p_fnd+i*3;
                //for(int j=0;j<3;++j)copyVec(v_begin(p_fv[j]),p_cdv+j*3);
                for(int j=0;j<3;++j)weightedAddVec(-0.00001,1.,p_fn,v_begin(p_fv[j]),p_cdv+j*3);
            }
            //for(int i=0;i<n_faces*9;++i)*(p_dv+i)+=1.01;

        }else{

            display_vertices = vertices;
            display_normal = vertices_normal;
        }


        reinitflag = false;
        isresetColor = true;

        //ne = display_edges.size();
        //nc = display_vcolor.size();
        nv = display_vertices.size();
        //nf = display_faces.size();
        nn = display_normal.size();



    }

    if(!isbuilddisp||isresetColor){


        if(isFaceRenderMode){
            //cout<<"adasdsadsadasdasdasdasdasdsadasda"<<endl;
            display_vcolor.resize(n_faces*4*3*2);
            int offset = n_faces*4*3;
            //cout<<"weighted_fcolor.size()==n_faces*4: "<<(weighted_fcolor.size()==n_faces*4)<<endl;
            if(ismark && weighted_fcolor.size()==n_faces*4){
                //cout<<"adsajkjljfkdlsjglkfjd"<<endl;


                if(isUseDoubleColorSingleSide || isUseDoubleColorDoubleSide || isUseBackSideColorDoubleSide){
                    for(int i=0;i<n_faces;++i){
                        auto ind = i*12;
                        auto ind2 = i*4;
                        for(int j=0;j<3;++j)for(int k=0;k<4;++k)display_vcolor[ind+j*4+k] = weighted_fcolor[ind2+k];
                    }
                }else{
                    for(int i = 0;i<n_faces*3;++i){
                        //for(int j=0;j<4;++j)display_vcolor.push_back(green[j]);
                        for(int j=0;j<4;++j)display_vcolor[i*4+j] = colllour[j];
                    }
                }


                if(isUseDoubleColorDoubleSide){
                    for(int i=0;i<n_faces;++i){
                        auto ind = i*12 + offset;
                        auto ind2 = i*4;
                        for(int j=0;j<3;++j)for(int k=0;k<4;++k)display_vcolor[ind+j*4+k] = weighted_fcolor[ind2+k];
                    }
                }else if(isUseBackSideColorDoubleSide && weighted_fcolorBackside.size()==n_faces*4){
                    for(int i=0;i<n_faces;++i){
                        auto ind = i*12 + offset;
                        auto ind2 = i*4;
                        for(int j=0;j<3;++j)for(int k=0;k<4;++k)display_vcolor[ind+j*4+k] = weighted_fcolorBackside[ind2+k];
                    }

                }else{
                    for(int i = 0;i<n_faces*3;++i){
                        //for(int j=0;j<4;++j)display_vcolor.push_back(green[j]);
                        for(int j=0;j<4;++j)display_vcolor[offset+i*4+j] = colllour[j];
                    }
                }


            }else{
                //cout<<"adsajkjljfkdlsjglkfjd"<<endl;
                if(isUseDoubleColorDoubleSide){
                    for(int i=0;i<n_faces;++i){
                        auto ind = i*12;
                        auto ind2 = i*4;
                        for(int j=0;j<3;++j)for(int k=0;k<4;++k)display_vcolor[ind+j*4+k] = weighted_fcolor[ind2+k];
                    }
                }else if(isUseBackSideColorDoubleSide && weighted_fcolorBackside.size()==n_faces*4){
                    for(int i=0;i<n_faces;++i){
                        auto ind = i*12;
                        auto ind2 = i*4;
                        for(int j=0;j<3;++j)for(int k=0;k<4;++k)display_vcolor[ind+j*4+k] = weighted_fcolorBackside[ind2+k];
                    }


                }else{
                    for(int i = 0;i<n_faces*3;++i){
                        //for(int j=0;j<4;++j)display_vcolor.push_back(green[j]);
                        for(int j=0;j<4;++j)display_vcolor[i*4+j] = colllour[j];
                    }
                }

                //cout<<"adsajkjljfkdlsjglkfjd "<<weighted_fcolor.size()<<endl;
                if(isUseDoubleColorSingleSide || isUseDoubleColorDoubleSide){
                    for(int i=0;i<n_faces;++i){
                        auto ind = i*12 + offset;
                        auto ind2 = i*4;
                        for(int j=0;j<3;++j)for(int k=0;k<4;++k)display_vcolor[ind+j*4+k] = weighted_fcolor[ind2+k];
                    }
                }else{
                    for(int i = 0;i<n_faces*3;++i){
                        //for(int j=0;j<4;++j)display_vcolor.push_back(green[j]);
                        for(int j=0;j<4;++j)display_vcolor[offset+i*4+j] = colllour[j];
                    }
                }

                //cout<<"adsajkjljfkdlsjglkfjd"<<endl;
            }
        }else{
            if(ismark && weighted_color.size()==n_vertices*4)display_vcolor = weighted_color;
            //else if(constrainActivated && weighted_color_constrain.size()==n_vertices*4)display_vcolor = weighted_color_constrain;
            else{
                display_vcolor.resize(n_vertices*4);
                for(int i = 0;i<n_vertices;++i){
                    //for(int j=0;j<4;++j)display_vcolor.push_back(green[j]);
                    for(int j=0;j<4;++j)display_vcolor[i*4+j] = colllour[j];
                }
            }
        }
        nc = display_vcolor.size();
        isresetColor=false;

    }

    display_vertices.resize(nv);
    //display_faces.resize(nf);
    display_normal.resize(nn);
    //display_edges.resize(ne);
    display_vcolor.resize(nc);
    display_field_dot.clear();


    //cout<<"adasdsadsadasdasdasdasdasdsadasda"<<endl;

    if(issurf){
        if(isFaceRenderMode){
            display_faces.resize(n_faces*6);
            int n3f = n_faces*3;

            auto p_df2 = display_faces.data()+n3f;
            auto p_df1 = display_faces.data();
            for(int i=0;i<n3f;++i)p_df1[i]=i;
            //for(int i=0;i<n3f;++i)p_df2[i]=i;
            for(int i=0;i<n_faces;++i)rotateTriVec3D(p_df1+i*3,p_df2+i*3);
            //cout<<p_df[0]<<' '<<p_df[1]<<' '<<p_df[2]<<' '<<endl;
            for(int i=0;i<n3f;++i)p_df2[i]+=n3f;
        }else display_faces = faces2vertices;
    }else display_faces.clear();
    if (iswire) {
        scalevn = 0;
        auto bias = display_vertices.size()/3;
        display_edges = edges;
        //display_edges = MSTedges;
        for(int i = 0;i<n_vertices;++i){
            auto p_v = v_begin(i);
            auto p_n = vnor_begin(i);
            for(int j=0;j<3;++j)display_vertices.push_back(p_v[j]+p_n[j]*scalevn);

        }


        for(auto &a:display_edges)a+=bias;
        //display_vcolor.reserve(bias*4+display_vcolor.size());
        for(int i = 0;i<n_vertices;++i){
            //for(int j=0;j<4;++j)display_vcolor.push_back(green[j]);
            for(int j=0;j<4;++j)display_vcolor.push_back(dark[j]);
        }
        display_normal.insert(display_normal.end(), vertices_normal.begin(), vertices_normal.end());
    }else display_edges.clear();

    //cout<<"disp: "<<display_vertices.size()<<' '<<display_edges.size()<<' '<<vertices.size()<<' '<<nv<<endl;







    auto drawsingularity=[this](vector<uint>* p_sing,double *p_pos,uchar *colour,double scale){
        int unmofSin = p_sing->size();
        double pos[3];
        for(int i=0;i<unmofSin;++i){
            int offset = display_vertices.size()/3;
            auto p_v = p_pos+3*(p_sing->at(i));
            for(int j = 0;j<sphere.n_vertices;++j){
                auto p_spv = sphere.v_begin(j);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]*scale);
            }
            for(int j = 0;j<sphere.n_vertices;++j){
                auto p_vn = sphere.vnor_begin(j);
                for(int k=0;k<3;++k)display_normal.push_back(p_vn[k]);
            }
            auto p_sf = sphere.fv_begin(0);
            for(int j = 0;j<sphere.n_faces*3;++j){
                display_faces.push_back(offset+p_sf[j]);
            }
            for(int j = 0;j<sphere.n_vertices;++j){
                for(int k=0;k<4;++k)display_vcolor.push_back(colour[k]);
            }
        }
    };

    auto drawsingularityVer=[this](int n_ver,double *p_pos,uchar *colour,double scale){
        int unmofSin = n_ver;
        double pos[3];
        for(int i=0;i<unmofSin;++i){
            int offset = display_vertices.size()/3;
            auto p_v = p_pos+3*i;
            for(int j = 0;j<sphere.n_vertices;++j){
                auto p_spv = sphere.v_begin(j);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]*scale);
            }
            for(int j = 0;j<sphere.n_vertices;++j){
                auto p_vn = sphere.vnor_begin(j);
                for(int k=0;k<3;++k)display_normal.push_back(p_vn[k]);
            }
            auto p_sf = sphere.fv_begin(0);
            for(int j = 0;j<sphere.n_faces*3;++j){
                display_faces.push_back(offset+p_sf[j]);
            }
            for(int j = 0;j<sphere.n_vertices;++j){
                for(int k=0;k<4;++k)display_vcolor.push_back(colour[i*4+k]);
            }
        }
    };

    auto drawsingularityVer2=[this](int n_ver,double *p_pos,uchar *colour,double scale){
        int unmofSin = n_ver;
        double pos[3];
        for(int i=0;i<unmofSin;++i){
            int offset = display_vertices.size()/3;
            auto p_v = p_pos+3*i;
            for(int j = 0;j<sphere.n_vertices;++j){
                auto p_spv = sphere.v_begin(j);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]*scale);
            }
            for(int j = 0;j<sphere.n_vertices;++j){
                auto p_vn = sphere.vnor_begin(j);
                for(int k=0;k<3;++k)display_normal.push_back(p_vn[k]);
            }
            auto p_sf = sphere.fv_begin(0);
            for(int j = 0;j<sphere.n_faces*3;++j){
                display_faces.push_back(offset+p_sf[j]);
            }
            for(int j = 0;j<sphere.n_vertices;++j){
                for(int k=0;k<4;++k)display_vcolor.push_back(colour[k]);
            }
        }
    };


    auto displayfield2=[this](double *pcent,double *pfield, double *pnor, uint num,double scale,double scaleFieldupn,uchar *colour){
        uint nnv = display_vertices.size()/3;
        for(uint i = 0;i<num;++i){
            auto p_fcnet = pcent+i*3;
            auto p_fn = pnor+i*3;
            auto p_fvec = pfield+i*3;

            for(int j =0;j<3;++j)display_vertices.push_back(p_fcnet[j] + p_fn[j]*scaleFieldupn);
            for(int j =0;j<3;++j)display_vertices.push_back(p_fcnet[j] + p_fn[j]*scaleFieldupn + p_fvec[j]*scale);
            for(int j =0;j<3;++j)display_normal.push_back(p_fn[j]);
            for(int j =0;j<3;++j)display_normal.push_back(p_fn[j]);

            for(int j =0;j<4;++j)display_vcolor.push_back(colour[j]);
            for(int j =0;j<4;++j)display_vcolor.push_back(colour[j]);

        }

        for(uint i = 0;i<num*2;++i)display_edges.push_back(nnv+i);
    };


    auto displayfield=[this](double *pcent,double *pfield, double *pofield, double *pnor, uint num,double scale,double scalew,double scaleFieldupn,uchar *colour){
        //uint nnv = display_vertices.size()/3;
        double vec[3],vecw[3];
        double find[6] = {0,2,1,1,2,3};
        display_field_dot.clear();
        for(uint i = 0;i<num;++i){
            if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;

            auto ind = display_vertices.size()/3;
            auto p_fcnet = pcent+i*3;
            auto p_fn = pnor+i*3;
            auto p_fvec = pfield+i*3;
            auto p_ofvec = pofield+i*3;

            for(int j =0;j<3;++j)vec[j] = p_fcnet[j] + p_fn[j]*scaleFieldupn;
            for(int j =0;j<3;++j)vecw[j] = p_ofvec[j]*scalew;

            for(int j =0;j<3;++j)display_field_dot.push_back(vec[j]);

            for(int j =0;j<3;++j)display_vertices.push_back(vec[j] + vecw[j]);
            for(int j =0;j<3;++j)display_vertices.push_back(vec[j] - vecw[j]);

            for(int j =0;j<3;++j)vec[j] +=  p_fvec[j]*scale;
            for(int j =0;j<3;++j)display_vertices.push_back(vec[j] + vecw[j]);
            for(int j =0;j<3;++j)display_vertices.push_back(vec[j] - vecw[j]);

            for(int k=0;k<4;++k)for(int j =0;j<3;++j)display_normal.push_back(p_fn[j]);
            for(int k=0;k<4;++k)for(int j =0;j<4;++j)display_vcolor.push_back(colour[j]);

            for(int j=0;j<6;++j)display_faces.push_back(ind+find[j]);



        }



    };


    if(isnormal){
        if(vertices_normal.size()!=0 )displayfield2(v_begin(0),vnor_begin(0),vnor_begin(0),n_vertices, scale,scaleFieldupn,blue);
    }

    if(issingularity){
        if(weighted_color.size()==n_vertices*4)drawsingularityVer(n_vertices,vertices.data(),weighted_color.data(),scale*0.2);

        if(hidden_ctrV.size()!=0){
            int nctrV = hidden_ctrV.size()/3;
            drawsingularityVer2(nctrV,hidden_ctrV.data(),gray,scale*0.13*controlBallsize);
            int offset = display_vertices.size()/3;
            for(auto a:hidden_ctrV)display_vertices.push_back(a);
            for(auto a:hidden_ctrVnor)display_normal.push_back(a);
            for(int i=0;i<nctrV;++i)for(int j=0;j<4;++j)display_vcolor.push_back(red[j]);
            for(auto a:hidden_ctrE)display_edges.push_back(offset+a);
        }

    }


    //if(ismark)if(markbufferv.size()!=0)drawsingularity(&markbufferv, v_begin(0),blue,0.01);
    isbuilddisp = true;

    //cout<<"display_faces: "<<display_faces.size()<<endl;

}



void Surface::sparseSampling(int a){
    if(n_vertices==0)return;

    sparsecoff = a;

    if(sparsecoff==0){
        sparseShowField.clear();
        sparseShowField.resize(n_vertices,true);
        return;
    }

    list<uint>ll1,ll2;

    auto issample = [this,a,&ll1,&ll2](uint ind){

        if(vertices_flinge[ind])return false;

        auto p_pre = &ll1;auto p_cur = &ll2;

        uint *pp;
        int num;
        p_pre->clear();p_cur->clear();
        p_pre->push_back(ind);

        for(int i=0;i<a;++i){
            while(!p_pre->empty()){
                ind = p_pre->front();
                p_pre->pop_front();
                num = vv_num(ind);
                pp=vv_begin(ind);
                for(int j=0;j<num;++j){
                    if(sparseShowField[pp[j]])return false;
                    else p_cur->push_back(pp[j]);
                }

            }
            swap(p_pre,p_cur);
        }

        return true;


    };

    sparseShowField.clear();
    sparseShowField.resize(n_vertices,false);
    vector<bool>visited(n_vertices,false);
    uint seed = 100;
    list<uint>ll3,*p_ppre,*p_pcur;
    list<uint>ll4;



    p_ppre = &ll3;p_pcur = &ll4;

    p_ppre->clear();p_pcur->clear();
    p_ppre->push_back(seed);
    visited[seed] = true;
    int iter = 0,knn = 0;
    while(true){
        ++iter;

        while(!p_ppre->empty()){
            auto ind = p_ppre->front();
            p_ppre->pop_front();

            auto vvum = vv_num(ind);
            auto p_vv = vv_begin(ind);
            if(issample(ind)){sparseShowField[ind] = true;++knn;}
            for(int i=0; i<vvum; ++i){
                if(!visited[p_vv[i]]){
                    visited[p_vv[i]] = true;
                    p_pcur->push_back(p_vv[i]);
                }
            }
        }
        //cout<<p_pcur->size()<<' ';cout<<endl;
        if(p_pcur->empty())break;
        swap(p_ppre,p_pcur);


    }

    //cout<<"sparse: "<<iter<<' '<<knn<<endl;
    //for(auto a:visited)cout<<a<<' ';cout<<endl;
}

void Surface::GetPerFaceColorDegree(vector<double>&facedegree){
    invert_faces_colordegree = facedegree;
    //    invert_vertices_colordegree.resize(n_vertices);
    //    for(int i=0;i<n_faces;++i){
    //        auto p_fv = fv_begin(i);
    //        for(int j=0;j<3;++j){invert_vertices_colordegree[p_fv[j]] = invert_faces_colordegree[i];}
    //    }

    //    weighted_color.resize(n_vertices*4);
    //    auto p_weighted = weighted_color.data();
    //    for(int i=0;i<n_vertices;++i){
    //        ThresholdColoring(invert_vertices_colordegree[i],p_weighted+4*i);
    //        weighted_color[4*i+3]=255;
    //    }

    weighted_fcolor.resize(n_faces*4);
    auto p_fweighted = weighted_fcolor.data();
    for(int i=0;i<n_faces;++i){
        weighted_fcolor[4*i+3]=255;
        ThresholdColoring(invert_faces_colordegree[i],p_fweighted+4*i);

    }

}

void Surface::GetPerFaceColorDegreeBackside(vector<double>&facedegree){
    invert_faces_colordegree_backside = facedegree;
    //    invert_vertices_colordegree.resize(n_vertices);
    //    for(int i=0;i<n_faces;++i){
    //        auto p_fv = fv_begin(i);
    //        for(int j=0;j<3;++j){invert_vertices_colordegree[p_fv[j]] = invert_faces_colordegree[i];}
    //    }

    //    weighted_color.resize(n_vertices*4);
    //    auto p_weighted = weighted_color.data();
    //    for(int i=0;i<n_vertices;++i){
    //        ThresholdColoring(invert_vertices_colordegree[i],p_weighted+4*i);
    //        weighted_color[4*i+3]=255;
    //    }

    weighted_fcolorBackside.resize(n_faces*4);
    auto p_fweighted = weighted_fcolorBackside.data();
    for(int i=0;i<n_faces;++i){
        weighted_fcolorBackside[4*i+3]=255;
        ThresholdColoring(invert_faces_colordegree_backside[i],p_fweighted+4*i);

    }

}

void Surface::GetPerVertexColorDegree(vector<double>&verticesdegree){

    invert_vertices_colordegree = verticesdegree;

    weighted_color.resize(n_vertices*4);
    auto p_weighted = weighted_color.data();
    for(int i=0;i<n_vertices;++i){
        ThresholdColoring(invert_vertices_colordegree[i],p_weighted+4*i);
        weighted_color[4*i+3]=255;
    }


}

void Surface::SetDisplayNormal(vector<bool>&isinvertFnormal){

    displayInvertfaceNormal = faces_normal;
    for(int i=0;i<n_faces;++i)if(isinvertFnormal[i]){
        auto p_fn = displayInvertfaceNormal.data()+i*3;
        inversevec(p_fn,p_fn);

    }

}

void Surface::WeightColor(const double thres,const vector<double>&weights,vector<uchar>&out_color){

    int n_para = weights.size();
    out_color.resize(n_para*4);
    double hdegree = 220;
    bool useHSV = true;
    vector<double>huehue(n_para);
    if(useHSV){
        for(int i =0;i<n_para;++i)huehue[i]=hdegree*(1-min(1.,fabs(weights[i])/(double)thres));
        meshScalarfunctionSmoothing(huehue,20);
    }
    for(int i =0;i<n_para;++i){
        auto a = fmin(fabs(weights[i]),(double)thres);
        if(!useHSV){
            a=min(255.,a/(double)thres*255);
            weighted_color[i*4] = uchar(a);
            weighted_color[i*4+1] = 0;
            weighted_color[i*4+2] = 255-uchar(a);
            weighted_color[i*4+3] = 255;
        }else{

            ThresholdColoring(huehue[i],0.8,0.8,&(weighted_color[i*4]));
            weighted_color[i*4+3] = 255;
        }
    }
    isresetColor=true;
}


void Surface::meshScalarfunctionSmoothing(vector<double>& scalarF,int maxiter){

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


        //        for(int i =0;i<n_vertices;++i){
        //            double aaa = 0;
        //            for(auto p_vv = vv_begin(i);p_vv!=vv_end(i);++p_vv){
        //                aaa+=pre->at(*p_vv);
        //            }
        //            aaa/=(vv_num(i));
        //            cur->at(i) = coef * pre->at(i) + ocoef * aaa;
        //        }
        //        swap(pre,cur);


    }
    scalarF=(*pre);





}


}//n_rf


struct rgb{
    double r;       // percent
    double g;       // percent
    double b;       // percent
    rgb(double r,double g,double b):r(r),g(g),b(b){}
    rgb(){}
};

struct hsv{
    double h;       // angle in degrees
    double s;       // percent
    double v;       // percent
    hsv(double h,double s,double v):h(h),s(s),v(v){}
    hsv(){}
};

hsv rgb2hsv(rgb &in)
{
    hsv         out;
    double      min, max, delta;

    min = in.r < in.g ? in.r : in.g;
    min = min  < in.b ? min  : in.b;

    max = in.r > in.g ? in.r : in.g;
    max = max  > in.b ? max  : in.b;

    out.v = max;                                // v
    delta = max - min;
    if (delta < 0.00001)
    {
        out.s = 0;
        out.h = 0; // undefined, maybe nan?
        return out;
    }
    if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
        out.s = (delta / max);                  // s
    } else {
        // if max is 0, then r = g = b = 0
        // s = 0, v is undefined
        out.s = 0.0;
        out.h = NAN;                            // its now undefined
        return out;
    }
    if( in.r >= max )                           // > is bogus, just keeps compilor happy
        out.h = ( in.g - in.b ) / delta;        // between yellow & magenta
    else
        if( in.g >= max )
            out.h = 2.0 + ( in.b - in.r ) / delta;  // between cyan & yellow
        else
            out.h = 4.0 + ( in.r - in.g ) / delta;  // between magenta & cyan

    out.h *= 60.0;                              // degrees

    if( out.h < 0.0 )
        out.h += 360.0;

    return out;
}


rgb hsv2rgb(hsv &in)
{
    double      hh, p, q, t, ff;
    long        i;
    rgb         out;

    if(in.s <= 0.0) {       // < is bogus, just shuts up warnings
        out.r = in.v;
        out.g = in.v;
        out.b = in.v;
        return out;
    }
    hh = in.h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.v * (1.0 - in.s);
    q = in.v * (1.0 - (in.s * ff));
    t = in.v * (1.0 - (in.s * (1.0 - ff)));

    switch(i) {
    case 0:
        out.r = in.v;
        out.g = t;
        out.b = p;
        break;
    case 1:
        out.r = q;
        out.g = in.v;
        out.b = p;
        break;
    case 2:
        out.r = p;
        out.g = in.v;
        out.b = t;
        break;

    case 3:
        out.r = p;
        out.g = q;
        out.b = in.v;
        break;
    case 4:
        out.r = t;
        out.g = p;
        out.b = in.v;
        break;
    case 5:
    default:
        out.r = in.v;
        out.g = p;
        out.b = q;
        break;
    }
    return out;
}


void ThresholdColoring(double degree,double s,double v,uchar *p_rgb){

    hsv hhh(degree,s,v);
    rgb rrr = hsv2rgb(hhh);
    p_rgb[0] = uchar(rrr.r*255.);
    p_rgb[1] = uchar(rrr.g*255.);
    p_rgb[2] = uchar(rrr.b*255.);
    //weighted_color[i*4+3] = 255;

}
vector<vector<uchar>>HHHH_tencolors_HHHH({
                                             {150,150,150},
                                             {255,255,0},
                                             {255,0,0},
                                             {0,255,0},
                                             {0,0,255},
                                             {0,255,255},
                                             {255,0,255},
                                             {255,128,0},
                                             {255,128,128},
                                             {128,0,128},
                                             {255,217,217},
                               });

/*
                                   {150,150,150},
                                   {255,255,0},
                                   {255,0,0},
                                   {0,255,0},
                                   {0,0,255},
                                   {0,255,255},
                                   {255,0,255},
                                   {255,128,0},
                                   {255,128,128},
                                   {128,0,128},
                                   {255,217,217},


                                   {150,150,150},
                                   {255,0,0},
                                   {0,0,255},
                                   {0,255,0},
                                   {255,255,0},
                                   {0,255,255},
                                   {255,0,255},
                                   {255,128,0},
                                   {255,128,128},
                                   {128,0,128},
                                   {255,217,217},


*/
void ThresholdColoring(double degree,uchar *p_rgb){
    if(degree == -10000.){//10000
        p_rgb[0] = 140;p_rgb[1] = 140;p_rgb[2] = 140;
        return;
    }else if(degree<-100){
        //p_rgb[0] = 40;p_rgb[1] = 40;p_rgb[2] = 40;
        p_rgb[0] = 0;p_rgb[1] = 0;p_rgb[2] = 0;
    }else if(degree<0){
        p_rgb[0] = 150;p_rgb[1] = 150;p_rgb[2] = 150;
        //p_rgb[3] = 100;
    }else {
        if(0)ThresholdColoring(degree,0.8,0.8,p_rgb);
        else{
            int ind = int(degree/340.*7);
            if(ind >10)ind=10;
            for(int i=0;i<3;++i)p_rgb[i] = HHHH_tencolors_HHHH[ind][i];
            rgb rgb1(p_rgb[0],p_rgb[1],p_rgb[2]);
            hsv hsv1 = rgb2hsv(rgb1);
            ThresholdColoring(hsv1.h,0.8,0.8,p_rgb);

        }
    }
}

