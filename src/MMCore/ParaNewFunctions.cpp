#include "ParallelArrangement.h"
#include <Eigen/Dense>
#include "Util.h"
#include "tetgen.h"
#include <set>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iostream>
#include <string>

#include "../Utility/readers.h"
#include "../Utility/geo_sur.h"
#include "../Utility/a_multisur.h"
// using Eigen::Vector3f;
// using namespace cimg_library;




void ParallelArrangement::loadInfoFromSpaceDivision(Arrangement &ar, int MMprop_f2vmethod, string outfoldername){


    this->outfoldername = outfoldername;
    _framef2Cs = ar._framef2Cs;
    _frameVstack = ar._frameVstack;
    _frameVs = ar._frameVs;
    _frameFs = ar._frameFs;
    _frameCells = ar._frameCells;
    _frameF2PlaneParas = ar._frameF2PlaneParas;
    _Cs2PlaneParas = ar._Cs2PlaneParas;
    _nCell = ar._nCell;
    _nCrossSection = ar._nCrossSection;
    _nFac = _framef2Cs.size();


    _globalEdges = ar._globalEdges;
    _globalEdgesMat = ar._globalEdgesMat;
    if(0){
        int mappp[3] = {1,2,0};
        vector<int>maxMmapping(100,-1);
        for(auto &a:_globalEdgesMat)for(auto &b:a)maxMmapping[b] = 0;
        int tttind = 0;
        for(auto &a:maxMmapping)if(a==0)a=tttind++;
        for(auto &a:_globalEdgesMat)for(auto &b:a)b = mappp[maxMmapping[b]];
        //for(auto &a:_globalEdgesMat)for(auto &b:a)b = maxMmapping[b];
    }


    _framefisCs.resize(_framef2Cs.size());
    for(int i=0;i<_nFac;++i)if(_framef2Cs[i]<0)_framefisCs[i] = false;else _framefisCs[i] = true;

    _frameF2Cell.clear();
    _frameF2Cell.resize(_nFac);
    for(int i=0;i<_frameCells.size();++i)for(auto &b:_frameCells[i])_frameF2Cell[b].push_back(i);
    for(auto &a:_frameF2Cell)assert(a.size()<=2);

    tetgenio &out = ar.out;

    cout<<"number of Vertices: "<<out.numberofpoints<<" number of Tets: "<<out.numberoftetrahedra<<endl;

    /***************************************/

    vector< vector<int> >FrameEdges_cluster(_frameFs.size());
    for(int i=0;i<_frameFs.size();++i){
        for(int j=0;j<_frameFs[i].size()-1;++j){
            FrameEdges_cluster[i].push_back(_frameFs[i][j]);
            FrameEdges_cluster[i].push_back(_frameFs[i][j+1]);
        }
        FrameEdges_cluster[i].push_back(_frameFs[i][_frameFs[i].size()-1]);
        FrameEdges_cluster[i].push_back(_frameFs[i][0]);
    }

    string vvfname = outfoldername+string("cellsframe");
    writeVVecFile(vvfname+string("_fe"),FrameEdges_cluster);
    writeVecFile(vvfname+string("_fv"),_frameVstack);
    writeVVecFile(vvfname+string("_cf"),_frameCells);



    /***************************************/




    auto pointMarkerT = [](int f){
        if(f==-1)return -2;//framework points
        if(f==-10000)return -1;//contour points shared by multiple cross-section
        if(f==0)return -3;//new inside points
        if(f<-99)return -f-100;//normal contour points on which cross-section
        if(f>99)return f-100+10000;//new points on which cross-section

        return -100000;//error

    };

    // Extract the output tets from tetgen data structure to our own data
    // structure.
    cout<<"Extract the output tets from tetgen data structure to our own data structure"<<endl;
    _tetVs.resize(out.numberofpoints, vector<float>(3));
    _tetVMarkers.resize(out.numberofpoints);
    _tetVMarkers_overlap.resize(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; ++i) {
        for (int j = 0; j < 3; ++j) {
            _tetVs[i][j] = out.pointlist[i * 3 + j];
        }
    }
    for (int i = 0; i < out.numberofpoints; ++i) {
        _tetVMarkers[i] = pointMarkerT(out.pointmarkerlist[i]);
    }
    _frameF2tetFs.resize(_frameFs.size());
    cout<<"_frameFs.size(): "<<_frameFs.size()<<endl;
    vector<int> oneTetF(3);
    cout<<"step 1"<<endl;
    for (int i = 0; i < out.numberoftrifaces; ++i) {
        int marker = out.trifacemarkerlist[i]-100;
        //cout<<marker<<' ';
        for (int j = 0; j < 3; ++j) {
            oneTetF[j] = out.trifacelist[i * 3 + j];
            _tetVMarkers_overlap[oneTetF[j]].insert(marker);
        }
        sort(oneTetF.begin(), oneTetF.end());
        _frameF2tetFs[marker].push_back(oneTetF);
    }

    //for (int i = 0; i < _frameF2tetFs.size(); ++i)cout<<_frameF2tetFs[i].size()<<endl;
    //cout<<endl;
    cout<<"step 2"<<endl;
    _cell2tetTs.resize(_frameCells.size());
    vector<unordered_set<int>> cellVsSet(_frameCells.size());
    assert(out.numberoftetrahedronattributes == 1);
    vector<int> oneTetT(4);
    for (int i = 0; i < out.numberoftetrahedra; ++i) {
        int marker = -out.tetrahedronattributelist[i] - 1;//cell index
        for (int j = 0; j < 4; ++j) {
            oneTetT[j] = out.tetrahedronlist[i * 4 + j];
            cellVsSet[marker].insert(oneTetT[j]);
        }
        sort(oneTetT.begin(), oneTetT.end());
        _cell2tetTs[marker].push_back(oneTetT);
    }
    cout<<"3"<<endl;
    _cell2tetVIds.resize(_frameCells.size());
    for (int i = 0; i < _frameCells.size(); ++i) {
        _cell2tetVIds[i].resize(cellVsSet[i].size());
        copy(cellVsSet[i].begin(), cellVsSet[i].end(), _cell2tetVIds[i].begin());
        sort(_cell2tetVIds[i].begin(),_cell2tetVIds[i].end());
    }
    _cell2MapTetV2CellV.resize(_nCell);
    for (int i = 0; i < _nCell; ++i) {
        auto& vMap = _cell2MapTetV2CellV[i];
        for (int j = 0; j < _cell2tetVIds[i].size(); ++j) {
            int tetVi = _cell2tetVIds[i][j];
            vMap[tetVi] = j;
        }
    }

    cout<<"out.numberofedges: "<<out.numberofedges<<endl;
    _tetEdges.resize(out.numberofedges*2);//_tetEdgeMarkers.resize(out.numberofedges);
    for(int i=0;i<out.numberofedges*2;++i)_tetEdges[i] = out.edgelist[i];
    //for(int i=0;i<out.numberofedges;++i)_tetEdgeMarkers[i] = out.edgemarkerlist[i];

    cout<<"----------------------- Label constrainted tetVs -----------------------"<<endl;
    // ----------------------- Label constrainted tetVs -----------------------
    // By default, everything is labeled as background: label 0.

    CtrEdgesRetracing();

    _constTetV2Label.resize(_tetVs.size(), 0);


    vector<int>OutfacesMat,ReOrientatedFaces;
    if(MMprop_f2vmethod==2 || MMprop_f2vmethod==4){
        HandleFaceletLabels(MMprop_f2vmethod, OutfacesMat,ReOrientatedFaces);
        CtrNetWorkReconstruction(ReOrientatedFaces,OutfacesMat,_constTetV2Label);
    }else{
        HandleFaceletLabels2(MMprop_f2vmethod, OutfacesMat,ReOrientatedFaces);
        CtrNetWorkReconstruction2(ReOrientatedFaces,OutfacesMat,_constTetV2Label);
        OutputCellbox2(ReOrientatedFaces,OutfacesMat);
    }

    //return;
    //MingUtility::printVector(_constTetV2Label);

    cout<<"----------------------- Get per cell labeling -----------------------"<<endl;
    // --------------------- Get per cell labeling ---------------------
    // Each cell might not use all the materials. To make the code easier, make
    // the labels consecutive in each cell.
    _cell2labeling.resize(_nCell);
    _cell2labels.resize(_nCell);
    // _tetLabel2cellLabel.resize(_nMat);
    _cell2MapTetL2CellL.resize(_nCell);
    for (int i = 0; i < _nCell; ++i) {
        // Find all the unique labels in cell i.
        unordered_set<int> cellLsSet;
        auto& oneCellLs = _cell2labels[i];
        for (auto vi : _cell2tetVIds[i])if(_tetVMarkers_overlap[vi].size()>0) {
            cellLsSet.insert(_constTetV2Label[vi]);
        }
        oneCellLs.resize(cellLsSet.size());
        copy(cellLsSet.begin(), cellLsSet.end(), oneCellLs.begin());
        sort(oneCellLs.begin(), oneCellLs.end());
        // Map the unique labels in cell i to consecutive numbers.
        auto& labelMap = _cell2MapTetL2CellL[i];
        for (int j = 0; j < oneCellLs.size(); ++j) {
            labelMap[oneCellLs[j]] = j;
        }
        // Use the new label to extract the initial labeling (mainly the boundary
        // constraint) for cell i.
        auto& oneCellLabeling = _cell2labeling[i];
        oneCellLabeling.resize(_cell2tetVIds[i].size());
        for (int j = 0; j < _cell2tetVIds[i].size(); ++j) {
            int vi = _cell2tetVIds[i][j];
            oneCellLabeling[j] = labelMap[_constTetV2Label[vi]];
        }
    }






    // ----------------------- Cross Section -----------------------
    // TODO this is for extracting cross section curves. Ideally we want to
    // extract the curves in ParallelArrangement class. However, I didn't have
    // enough time to move the feature from class Mesh, which I implemented along
    // ago, to this newly added class. MaybeI will change it in the future, or
    // not...
    int nFrameF = _frameFs.size();
    _frameF2ActiveFs.resize(nFrameF);
    _frameF2CrvVs.resize(nFrameF);
    _frameF2CrvEs.resize(nFrameF);

    for(int i=0;i<_nCell;++i){
        vector<int>tmpvec;
        for(auto &a:_frameCells[i])if(_framef2Cs[a]>=0)tmpvec.push_back(a);
        _frameCell2csFrameIDs.push_back(tmpvec);
        cout<<i<<": ";for(auto a:tmpvec)cout<<a<<' ';cout<<endl;
    }



    _isFrameFReported.resize(nFrameF, true);


    vector<vector<int>> ePair;
    MingUtility::combination(3, 2, ePair);
    for (int fi = 0; fi < nFrameF; ++fi) {
        if(_framef2Cs[fi]<0)continue;
        _isFrameFReported[fi] = false;
        map<vector<int>, vector<int>> edges;
        for (const auto& oneF : _frameF2tetFs[fi]) {
            // Check activeness of oneF.
            if (!(_constTetV2Label[oneF[0]] == _constTetV2Label[oneF[1]] &&
                  _constTetV2Label[oneF[0]] == _constTetV2Label[oneF[2]])) {
                // Check activeness of oneE.
                for (const auto& oneP : ePair) {
                    vector<int> oneE = {oneF[oneP[0]], oneF[oneP[1]]};
                    sort(oneE.begin(), oneE.end());
                    if (_constTetV2Label[oneE[0]] != _constTetV2Label[oneE[1]]) {
                        edges[oneE].push_back(_frameF2ActiveFs[fi].size());
                    }
                }
                _frameF2ActiveFs[fi].push_back(oneF);
            }
        }
        for (const auto& oneE : edges) {
            //assert(oneE.second.size() == 2);
            _frameF2CrvEs[fi].push_back(oneE.second);
        }
    }

}

void ParallelArrangement::HandleFaceletLabels(int f2vmethod,vector<int>&OutfacesMat,vector<int>&ReOrientatedFaces){

    _constTetV2Label.resize(_tetVs.size(), 0);
    for(int i=0;i<_frameF2tetFs.size();++i){
        vector<double>vpos;
        vector<int>newVind(_tetVs.size(),-1);
        vector<int>inverseFv2Tv;
        vector<unsigned int>newFaces;
        vector<unsigned int>newEdges;
        int accInd = 0;
        for(auto &a:_frameF2tetFs[i])for(auto b:a)newVind[b] = 0;
        //if(i==0)for(auto &a:_frameF2tetFs[i])for(auto b:a)cout<<b<<' ';cout<<endl;
        for(int j=0;j<_tetVs.size();++j)if(newVind[j]==0){
            newVind[j]=accInd;++accInd;
            for(int k=0;k<3;++k)vpos.push_back(_tetVs[j][k]);
            inverseFv2Tv.push_back(j);
        }

        for(int j=0;j<_frameF2tetFs[i].size();++j)for(int k=0;k<3;++k){
            newFaces.push_back((uint)newVind[_frameF2tetFs[i][j][k]]);
            //newFaces.push_back((uint)_frameF2tetFs[i][j][k]);
        }
        for(auto a:_globalEdges[i])newEdges.push_back((uint)newVind[a]);

        string filename = outfoldername+to_string(i);

        n_rf::Surface planesurface;
        planesurface.importSurface(vpos,newFaces,true);

        vector<int>facesMat(planesurface.n_faces,0),verticesMat(planesurface.n_vertices,0);
        if(_framef2Cs[i]>=0 && newEdges.size()>0){
            double plnormal[3];
            for(int j=0;j<3;++j)plnormal[j] = _frameF2PlaneParas[i][j];
            if(MyUtility::dot(plnormal,planesurface.vnor_begin(0))<0)planesurface.InverseFaces();
            //if(f2vmethod==4)f2vmethod = 40000+_framef2Cs[i];
            planesurface.MMpropagation(f2vmethod,newEdges,_globalEdgesMat[i],facesMat,verticesMat);
            for(auto fv:planesurface.faces2vertices)ReOrientatedFaces.push_back(inverseFv2Tv[fv]);
            for(auto f:facesMat)OutfacesMat.push_back(f);
        }

        for(int j=0;j<accInd;++j){
            if(_constTetV2Label[inverseFv2Tv[j]]!=0)assert( _constTetV2Label[inverseFv2Tv[j]] == verticesMat[j]);
            _constTetV2Label[inverseFv2Tv[j]] = verticesMat[j];
        }


        writeObjFile(filename,planesurface.vertices,planesurface.faces2vertices);
        writeContourEdgeTxtFile(filename,newEdges);
        writeVecFile(filename+"_f",facesMat);
        writeVecFile(filename+"_v",verticesMat);


    }

    _nMat = *max_element(_constTetV2Label.begin(),_constTetV2Label.end())+1;

    //for(int i=0;i<_frameF2tetFs.size();++i)cout<<"_framef2Cs: "<<_framef2Cs[i]<<endl;

}
void ParallelArrangement::HandleFaceletLabels2(int f2vmethod,vector<int>&OutfacesMat,vector<int>&ReOrientatedFaces){

    _constTetV2Label.resize(_tetVs.size(), 0);
    vector<vector<int>>Cs2frame(_nCrossSection);
    for(int i=0;i<_frameF2tetFs.size();++i)if(_framef2Cs[i]>=0)Cs2frame[_framef2Cs[i]].push_back(i);
    for(int i=0;i<_nCrossSection;++i){
        vector<double>vpos;
        vector<int>newVind(_tetVs.size(),-1);
        vector<int>inverseFv2Tv;
        vector<unsigned int>newFaces;
        vector<unsigned int>newEdges;
        vector<int>newEdgesMat;
        int accInd = 0;
        for(auto b:Cs2frame[i])for(auto &a:_frameF2tetFs[b])for(auto b:a)newVind[b] = 0;
        //if(i==0)for(auto &a:_frameF2tetFs[i])for(auto b:a)cout<<b<<' ';cout<<endl;
        for(int j=0;j<_tetVs.size();++j)if(newVind[j]==0){
            newVind[j]=accInd;++accInd;
            for(int k=0;k<3;++k)vpos.push_back(_tetVs[j][k]);
            inverseFv2Tv.push_back(j);
        }

        for(auto b:Cs2frame[i])for(int j=0;j<_frameF2tetFs[b].size();++j)for(int k=0;k<3;++k){
            newFaces.push_back((uint)newVind[_frameF2tetFs[b][j][k]]);
            //newFaces.push_back((uint)_frameF2tetFs[i][j][k]);
        }
        for(auto b:Cs2frame[i])for(auto a:_globalEdges[b])newEdges.push_back((uint)newVind[a]);
        for(auto b:Cs2frame[i])for(auto a:_globalEdgesMat[b])newEdgesMat.push_back(a);

        string filename = outfoldername+to_string(i);

        n_rf::Surface planesurface;
        planesurface.importSurface(vpos,newFaces,true);

        vector<int>facesMat(planesurface.n_faces,0),verticesMat(planesurface.n_vertices,0);
        {
            double plnormal[3];
            for(int j=0;j<3;++j)plnormal[j] = _Cs2PlaneParas[i][j];
            if(MyUtility::dot(plnormal,planesurface.vnor_begin(0))<0)planesurface.InverseFaces();
            planesurface.MMpropagation(f2vmethod,newEdges,newEdgesMat,facesMat,verticesMat);
            for(auto fv:planesurface.faces2vertices)ReOrientatedFaces.push_back(inverseFv2Tv[fv]);
            for(auto f:facesMat)OutfacesMat.push_back(f);
        }

        for(int j=0;j<accInd;++j){
            if(_constTetV2Label[inverseFv2Tv[j]]!=0)assert( _constTetV2Label[inverseFv2Tv[j]] == verticesMat[j]);
            _constTetV2Label[inverseFv2Tv[j]] = verticesMat[j];
        }


        writeObjFile(filename,planesurface.vertices,planesurface.faces2vertices);
        writeContourEdgeTxtFile(filename,newEdges);
        writeVecFile(filename+"_f",facesMat);
        writeVecFile(filename+"_v",verticesMat);


    }

    _nMat = *max_element(_constTetV2Label.begin(),_constTetV2Label.end())+1;

    //for(int i=0;i<_frameF2tetFs.size();++i)cout<<"_framef2Cs: "<<_framef2Cs[i]<<endl;

}

void BFS(vector<int>&edges,int nv,int st,int ed,vector<int>&seq){

    vector< vector<int> >nne(nv);
    for(int i=0;i<edges.size()/2;++i){
        int v1 = edges[i*2], v2 = edges[i*2+1];
        nne[v1].push_back(v2);
        nne[v2].push_back(v1);
    }

    list<pair<int,int> >flinge;
    vector<bool>visited(nv,false);
    flinge.push_back(make_pair(st,-1));
    vector<int>actionlist;
    vector<int>ansestor(nv,-2);
    ansestor[st] = 0;


    bool done = false;

    while(!flinge.empty()){
        pair<int,int> curI = flinge.front();
        flinge.pop_front();
        int cur = curI.first;
        if(visited[cur])continue;
        visited[cur] = true;
        ansestor[cur] = curI.second;
        if(cur==ed){done=true;break;}

        for(auto vv: nne[cur]){
            if(visited[vv])continue;
            flinge.push_back(make_pair(vv,cur));
        }
    }
    assert(done==true);
    vector<int>iseq;

    while(ansestor[ed]!=-1){
        iseq.push_back(ed);
        ed = ansestor[ed];
    }

    iseq.push_back(st);
    for(int i=iseq.size()-1;i>=0;--i)seq.push_back(iseq[i]);


}

void ParallelArrangement::CtrEdgesRetracing(){


    //    vector<vector<int>>_globalEdges;
    //    vector<vector<int>>_globalEdgesMat;
    //    vector<int>_tetEdges;
    //    vector<int>_tetEdgeMarkers;
    //    vector<int>_tetVMarkers;
    vector< vector<int> >globalEdgesPool(_frameFs.size());
    vector< vector<int> >globalEdgesMarkerPool(_frameFs.size());


    cout<<"_tetEdges.size(): "<<_tetEdges.size()<<endl;

    for(int i = 0;i<_tetEdges.size()/2;++i){
        int e1 = _tetEdges[i*2], e2 = _tetEdges[i*2+1];
        int m1 = _tetVMarkers[e1],m2 = _tetVMarkers[e2];
        if(m1>10000-1)m1-=10000;if(m2>10000-1)m2-=10000;

        if((m1>-1 || m2>-1) && m1>-2 && m2 >-2){
            if(m1>-1 && m2>-1)if(m1!=m2)continue;
            int ff=-10;
            if(m1>-1)ff=m1;else ff = m2;
            globalEdgesPool[ff].push_back(e1);globalEdgesPool[ff].push_back(e2);
            globalEdgesMarkerPool[ff].push_back( _tetVMarkers[e1]);globalEdgesMarkerPool[ff].push_back( _tetVMarkers[e2]);
        }
    }
    for(int i=0;i<globalEdgesPool.size();++i){
        bool isspliteCtr = false;
        for(int j=0;j<globalEdgesMarkerPool[i].size()/2;++j){
            int m1 = globalEdgesMarkerPool[i][j*2],m2 = globalEdgesMarkerPool[i][j*2+1];
            if((m1>10000-1 && m2<10000 && m2>-1) || (m2>10000-1 && m1<10000 && m1>-1)){isspliteCtr = true;break;}
        }
        if(!isspliteCtr){globalEdgesPool[i].clear();globalEdgesMarkerPool[i].clear();}
    }



    for(int i=0;i<_globalEdges.size();++i){
        if(globalEdgesPool[i].size()==0)continue;

        vector<bool>isexist(_globalEdges[i].size()/2,false);

        vector<int>missingEdge;
        vector<vector<int> >completionE;
        for(int j = 0;j<isexist.size();++j){
            int v1 = _globalEdges[i][j*2], v2 = _globalEdges[i][j*2+1];
            for(int k=0;k<globalEdgesPool.size()/2;++k){
                int vv1 = globalEdgesPool[i][k*2], vv2 = globalEdgesPool[i][k*2+1];
                if((v1==vv1 && v2 ==vv2)||(v1==vv2 && v2 ==vv1)){
                    isexist[j] = true;
                    globalEdgesPool[i].erase(globalEdgesPool[i].begin()+k*2);
                    globalEdgesPool[i].erase(globalEdgesPool[i].begin()+k*2);
                    globalEdgesMarkerPool[i].erase(globalEdgesMarkerPool[i].begin()+k*2);
                    globalEdgesMarkerPool[i].erase(globalEdgesMarkerPool[i].begin()+k*2);
                    break;
                }
            }
            if(!isexist[j])missingEdge.push_back(j);
        }
        if(missingEdge.size()==0)continue;
        vector<int>tmpind(_tetVs.size(),-1);
        vector<int>invtmpind;
        for(auto a:globalEdgesPool[i]){tmpind[a]=0;}
        int inddd = 0;
        for(int j=0;j<tmpind.size();++j)if(tmpind[j]==0){tmpind[j]=inddd++;invtmpind.push_back(j);}
        vector<int>tmpedge = globalEdgesPool[i];
        for(auto&a:tmpedge)a=tmpind[a];
        vector<vector<int>>seqs(missingEdge.size());
        for(int j=0;j< missingEdge.size();++j){
            auto a = missingEdge[j];
            int v1 = tmpind[_globalEdges[i][a*2]], v2 = tmpind[_globalEdges[i][a*2+1]];
            BFS(tmpedge,inddd,v1,v2,seqs[j]);
            for(auto &b:seqs[j])b = invtmpind[b];
        }

        vector<int>newEBuf;
        vector<int>newEMBuf;int mi = 0;
        for(int j = 0;j<isexist.size();++j){
            if(isexist[j]){
                newEBuf.push_back(_globalEdges[i][j*2]);
                newEBuf.push_back(_globalEdges[i][j*2+1]);
                newEMBuf.push_back(_globalEdgesMat[i][j*2]);
                newEMBuf.push_back(_globalEdgesMat[i][j*2+1]);
            }else{
                for(int k=0;k<seqs[mi].size()-1;++k){
                    newEBuf.push_back(seqs[mi][k]);
                    newEBuf.push_back(seqs[mi][k+1]);
                    newEMBuf.push_back(_globalEdgesMat[i][j*2]);
                    newEMBuf.push_back(_globalEdgesMat[i][j*2+1]);
                }
                ++mi;
            }
        }
        _globalEdges[i] = newEBuf;
        _globalEdgesMat[i] = newEMBuf;
        //        cout<<i<<": ";
        //        for(auto &a:seqs){
        //            for(auto b:a)cout<<b<<' ';
        //            cout<<"||";
        //        }
        //        cout<<endl;




    }
    cout<<"Edges Retracing end!"<<endl;
    return;
    for(int i=0;i<globalEdgesPool.size();++i){
        if(globalEdgesPool[i].size()==0)continue;
        cout<<i<<" : ";
        for(int j=0;j<globalEdgesMarkerPool[i].size()/2;++j)
            cout<<globalEdgesPool[i][j*2]<<' '<<globalEdgesPool[i][j*2+1]<<" || ";cout<<endl;
        for(int j=0;j<globalEdgesMarkerPool[i].size()/2;++j)
            cout<<globalEdgesMarkerPool[i][j*2]<<' '<<globalEdgesMarkerPool[i][j*2+1]<<" || ";cout<<endl;
        for(auto a: _globalEdges[i])cout<<a<<' ';cout<<endl;
    }



}

void ParallelArrangement::OutputCellbox2(vector<int> &ReOrientatedFaces, vector<int> &facesMat){


    vector<double>vpos;
    vector<int>newVind(_tetVs.size(),-1);
    vector<int>inverseFv2Tv;
    vector<unsigned int>newFaces;
    vector<unsigned int>newEdges;
    vector<int>frameVMat;
    vector<int>faces2Cs;
    vector<int>faces2Cell;

    for(auto a:ReOrientatedFaces)newFaces.push_back(uint(a));





    vector<vector<int>>Cs2frame(_nCrossSection);
    for(int i=0;i<_frameF2tetFs.size();++i)if(_framef2Cs[i]>=0)Cs2frame[_framef2Cs[i]].push_back(i);
    for(int b=0;b<_nCrossSection;++b){
        for(auto i:Cs2frame[b]){
            for(int j=0;j<_frameF2tetFs[i].size();++j)faces2Cs.push_back(_framef2Cs[i]);
            assert(_frameF2Cell[i].size()==2);
            for(int j=0;j<_frameF2tetFs[i].size();++j){
                faces2Cell.push_back(_frameF2Cell[i][0]); faces2Cell.push_back(_frameF2Cell[i][1]);
            }
        }
    }



    for(int i=0;i<_frameF2tetFs.size();++i)if(_framef2Cs[i]<0){
        for(int j=0;j<_frameF2tetFs[i].size();++j)faces2Cs.push_back(_framef2Cs[i]);
        for(int j=0;j<_frameF2tetFs[i].size();++j)facesMat.push_back(0);
        for(auto &a:_frameF2tetFs[i])for(auto b:a)newFaces.push_back(b);
        if(_frameF2Cell[i].size()==2){
            for(int j=0;j<_frameF2tetFs[i].size();++j){
                faces2Cell.push_back(_frameF2Cell[i][0]); faces2Cell.push_back(_frameF2Cell[i][1]);
            }
        }else{
            for(int j=0;j<_frameF2tetFs[i].size();++j){
                faces2Cell.push_back(_frameF2Cell[i][0]); faces2Cell.push_back(-1);
            }
        }
    }

    int accInd = 0;
    for(auto a:newFaces)newVind[a] = 0;

    for(int j=0;j<_tetVs.size();++j)if(newVind[j]==0){
        newVind[j]=accInd;++accInd;
        for(int k=0;k<3;++k)vpos.push_back(_tetVs[j][k]);
    }

    for(auto &a:newFaces)a = newVind[a];
    for(int i=0;i<_frameF2tetFs.size();++i)for(auto a:_globalEdges[i])newEdges.push_back((uint)newVind[a]);


    string filename = outfoldername+to_string(200000);
    writeObjFile(filename,vpos,newFaces);
    writeContourEdgeTxtFile(filename,newEdges);
    writeVecFile(filename+"_f",facesMat);
    writeVecFile(filename+"_v",frameVMat);
    writeVecFile(filename+"_fc",faces2Cs);
    writeVecFile(filename+"_fce",faces2Cell);


}

void ParallelArrangement::CtrNetWorkReconstruction2(vector<int> &ReOrientatedFaces, vector<int>&facesMat, vector<int>&verticesMat){


    vector<double>vpos;
    vector<int>newVind(_tetVs.size(),-1);
    vector<int>inverseFv2Tv;
    vector<unsigned int>newFaces;
    vector<unsigned int>newEdges;
    vector<int>frameVMat;
    vector<int>faces2Cs;
    vector<int>faces2Cell;

    int accInd = 0;
    for(auto a:ReOrientatedFaces)newVind[a] = 0;

    vector<vector<int>>Cs2frame(_nCrossSection);
    for(int i=0;i<_frameF2tetFs.size();++i)if(_framef2Cs[i]>=0)Cs2frame[_framef2Cs[i]].push_back(i);
    for(int b=0;b<_nCrossSection;++b){
        for(auto i:Cs2frame[b]){
            for(int j=0;j<_frameF2tetFs[i].size();++j)faces2Cs.push_back(_framef2Cs[i]);
            assert(_frameF2Cell[i].size()==2);
            for(int j=0;j<_frameF2tetFs[i].size();++j){
                faces2Cell.push_back(_frameF2Cell[i][0]); faces2Cell.push_back(_frameF2Cell[i][1]);
            }

        }
    }

    //if(i==0)for(auto &a:_frameF2tetFs[i])for(auto b:a)cout<<b<<' ';cout<<endl;

    mergeCsV2Tv.clear();
    for(int j=0;j<_tetVs.size();++j)if(newVind[j]==0){
        newVind[j]=accInd;++accInd;mergeCsV2Tv.push_back(j);
        for(int k=0;k<3;++k)vpos.push_back(_tetVs[j][k]);
        inverseFv2Tv.push_back(j);
        frameVMat.push_back(verticesMat[j]);
    }

    for(auto a:ReOrientatedFaces)newFaces.push_back((uint)newVind[a]);

    //    for(int i=0;i<_frameF2tetFs.size();++i){
    //        if(_framef2Cs[i]<0 || _globalEdges[i].size()<=0)continue;
    //        for(int j=0;j<_frameF2tetFs[i].size();++j)for(int k=0;k<3;++k){
    //            newFaces.push_back((uint)newVind[_frameF2tetFs[i][j][k]]);
    //            //newFaces.push_back((uint)_frameF2tetFs[i][j][k]);
    //        }
    //    }
    for(int i=0;i<_frameF2tetFs.size();++i)for(auto a:_globalEdges[i])newEdges.push_back((uint)newVind[a]);




    string filename = outfoldername+to_string(100000);


    mergeCs.importSurface(vpos,newFaces,frameVMat,faces2Cell,faces2Cs);
    writeObjFile(filename,mergeCs.vertices,mergeCs.faces2vertices);
    writeContourEdgeTxtFile(filename,newEdges);
    writeVecFile(filename+"_f",facesMat);
    writeVecFile(filename+"_v",frameVMat);
    writeVecFile(filename+"_fc",faces2Cs);
    writeVecFile(filename+"_fce",faces2Cell);

    newCtrNet.ImportCurve(mergeCs.hidden_ctrV,mergeCs.hidden_ctrE,mergeCs.hidden_ctrELable,mergeCs.hidden_ctrE2Cells,_nCell);




    auto& seg2edges = newCtrNet.mixCurve.seg2edges;
    seg_dualedges2Tv.resize(seg2edges.size());

    for(int i=0;i<seg2edges.size();++i){
        auto &dualedges2Tv = seg_dualedges2Tv[i];
        dualedges2Tv.clear();
        for(auto a:seg2edges[i]){
            auto v1 = mergeCsV2Tv[mergeCs.hidden_ctrActiveE2V[a*2]];
            auto v2 = mergeCsV2Tv[mergeCs.hidden_ctrActiveE2V[a*2+1]];
            if(v1>v2)swap(v1,v2);
            dualedges2Tv.push_back(v1);
            dualedges2Tv.push_back(v2);
        }

    }

    auto& ActiveFV = mergeCs.hidden_ctrActiveF2V;
    ActiveFV2Tv.resize(ActiveFV.size());

    for(int i=0;i<ActiveFV.size();++i){
        ActiveFV2Tv[i] = mergeCsV2Tv[ActiveFV[i]];
    }

    ActiveFCP.resize(ActiveFV.size());
    ActiveFCP_setup.clear();ActiveFCP_setup.resize(ActiveFV.size()/3,0);
    for(int i=0;i<ActiveFCP_setup.size();++i){
        sort(ActiveFV2Tv.begin()+i*3,ActiveFV2Tv.begin()+i*3+3);
    }

    //return;

    auto& ActiveNonE2V = mergeCs.hidden_ctrActiveNonE2V;
    ActiveNonEV2Tv.resize(ActiveNonE2V.size());

    for(int i=0;i<ActiveNonEV2Tv.size();++i){
        ActiveNonEV2Tv[i] = mergeCsV2Tv[ActiveNonE2V[i]];
    }
    for(int i=0;i<ActiveNonEV2Tv.size()/2;++i){
        sort(ActiveNonEV2Tv.begin()+i*2,ActiveNonEV2Tv.begin()+i*2+2);
    }



}

void ParallelArrangement::CtrNetWorkReconstruction(vector<int> &ReOrientatedFaces, vector<int>&facesMat, vector<int>&verticesMat){


    vector<double>vpos;
    vector<int>newVind(_tetVs.size(),-1);
    vector<int>inverseFv2Tv;
    vector<unsigned int>newFaces;
    vector<unsigned int>newEdges;
    vector<int>frameVMat;
    vector<int>faces2Cs;
    vector<int>faces2Cell;

    int accInd = 0;
    for(auto a:ReOrientatedFaces)newVind[a] = 0;



    for(int i=0;i<_frameF2tetFs.size();++i){

        //if(_framef2Cs[i]<0 || _globalEdges[i].size()<=0)continue;
        if(_framef2Cs[i]<0 || _globalEdges[i].size()<=0)continue;

        for(int j=0;j<_frameF2tetFs[i].size();++j)faces2Cs.push_back(_framef2Cs[i]);
        assert(_frameF2Cell[i].size()==2);
        for(int j=0;j<_frameF2tetFs[i].size();++j){
            faces2Cell.push_back(_frameF2Cell[i][0]); faces2Cell.push_back(_frameF2Cell[i][1]);
        }

    }
    //if(i==0)for(auto &a:_frameF2tetFs[i])for(auto b:a)cout<<b<<' ';cout<<endl;

    mergeCsV2Tv.clear();
    for(int j=0;j<_tetVs.size();++j)if(newVind[j]==0){
        newVind[j]=accInd;++accInd;mergeCsV2Tv.push_back(j);
        for(int k=0;k<3;++k)vpos.push_back(_tetVs[j][k]);
        inverseFv2Tv.push_back(j);
        frameVMat.push_back(verticesMat[j]);
    }

    for(auto a:ReOrientatedFaces)newFaces.push_back((uint)newVind[a]);

    //    for(int i=0;i<_frameF2tetFs.size();++i){
    //        if(_framef2Cs[i]<0 || _globalEdges[i].size()<=0)continue;
    //        for(int j=0;j<_frameF2tetFs[i].size();++j)for(int k=0;k<3;++k){
    //            newFaces.push_back((uint)newVind[_frameF2tetFs[i][j][k]]);
    //            //newFaces.push_back((uint)_frameF2tetFs[i][j][k]);
    //        }
    //    }
    for(int i=0;i<_frameF2tetFs.size();++i)for(auto a:_globalEdges[i])newEdges.push_back((uint)newVind[a]);




    string filename = outfoldername+to_string(100000);


    mergeCs.importSurface(vpos,newFaces,frameVMat,faces2Cell,faces2Cs);
    writeObjFile(filename,mergeCs.vertices,mergeCs.faces2vertices);
    writeContourEdgeTxtFile(filename,newEdges);
    writeVecFile(filename+"_f",facesMat);
    writeVecFile(filename+"_v",frameVMat);
    writeVecFile(filename+"_fc",faces2Cs);
    writeVecFile(filename+"_fce",faces2Cell);

    newCtrNet.ImportCurve(mergeCs.hidden_ctrV,mergeCs.hidden_ctrE,mergeCs.hidden_ctrELable,mergeCs.hidden_ctrE2Cells,_nCell);




    auto& seg2edges = newCtrNet.mixCurve.seg2edges;
    seg_dualedges2Tv.resize(seg2edges.size());

    for(int i=0;i<seg2edges.size();++i){
        auto &dualedges2Tv = seg_dualedges2Tv[i];
        dualedges2Tv.clear();
        for(auto a:seg2edges[i]){
            auto v1 = mergeCsV2Tv[mergeCs.hidden_ctrActiveE2V[a*2]];
            auto v2 = mergeCsV2Tv[mergeCs.hidden_ctrActiveE2V[a*2+1]];
            if(v1>v2)swap(v1,v2);
            dualedges2Tv.push_back(v1);
            dualedges2Tv.push_back(v2);
        }

    }

    auto& ActiveFV = mergeCs.hidden_ctrActiveF2V;
    ActiveFV2Tv.resize(ActiveFV.size());

    for(int i=0;i<ActiveFV.size();++i){
        ActiveFV2Tv[i] = mergeCsV2Tv[ActiveFV[i]];
    }

    ActiveFCP.resize(ActiveFV.size());
    ActiveFCP_setup.clear();ActiveFCP_setup.resize(ActiveFV.size()/3,0);
    for(int i=0;i<ActiveFCP_setup.size();++i){
        sort(ActiveFV2Tv.begin()+i*3,ActiveFV2Tv.begin()+i*3+3);
    }


    auto& ActiveNonE2V = mergeCs.hidden_ctrActiveNonE2V;
    ActiveNonEV2Tv.resize(ActiveNonE2V.size());

    for(int i=0;i<ActiveNonEV2Tv.size();++i){
        ActiveNonEV2Tv[i] = mergeCsV2Tv[ActiveNonE2V[i]];
    }
    for(int i=0;i<ActiveNonEV2Tv.size()/2;++i){
        sort(ActiveNonEV2Tv.begin()+i*2,ActiveNonEV2Tv.begin()+i*2+2);
    }



}

void ParallelArrangement::OutputLoops(int loopi, vector<double>&loopVs, vector<uint>&loopEs, int &label){
    loopVs.clear();
    loopEs.clear();
    if(loopi<0 || loopi>=_nLoops)return;
    int nv = 0;
    //newCtrNet.mixCurve.cell2Seg[celli];
    auto &pickSeg = compactloop2Segs[loopi];
    auto &seg2edges = newCtrNet.mixCurve.seg2edges;
    vector<bool>pickedges(newCtrNet.mixCurve.n_edges,false);
    vector<int>pickVs(newCtrNet.mixCurve.n_vertices,-1);

    for(auto a:pickSeg)for(auto b:seg2edges[a])pickedges[b] = true;
    for(int i=0;i<pickedges.size();++i)if(pickedges[i]){
        auto p_ev = newCtrNet.mixCurve.ev_begin(i);
        pickVs[p_ev[0]] = 0;
        pickVs[p_ev[1]] = 0;
        loopEs.push_back(p_ev[0]);loopEs.push_back(p_ev[1]);
    }
    for(int i=0;i<pickVs.size();++i)if(pickVs[i]==0){
        pickVs[i] = nv++;
        auto p_v = newCtrNet.mixCurve.v_begin(i);
        for(int j=0;j<3;++j)loopVs.push_back(p_v[j]);
    }
    for(auto &a:loopEs)a = pickVs[a];

    label = compactloop2label[loopi];

    //    cout<<"_nLoops: "<<_nLoops<<endl;
    //    for(auto a:pickSeg)cout<<a<<' ';cout<<endl;
    //    for(auto a:pickSeg){
    //        for(auto b:seg2edges[a])cout<<b<<' ';cout<<endl;
    //    }
}

int ParallelArrangement::MakeupLoop2SegsLists(){

    _nLoops = 0;

    loop2compactloopN.clear();
    loop2compactloopN.resize(_nCell);
    compactloop2Segs.clear();
    compactloop2label.clear();
    vector<vector<int>>tmp_compactloop2Segs;
    vector<int>tmp_compactloop2label;
    for(int i=0;i<_nCell;++i){
        for(int j=0;j<loop2Segs[i].size();++j){
            loop2compactloopN[i].push_back(_nLoops++);
            tmp_compactloop2Segs.push_back(loop2Segs[i][j]);
            tmp_compactloop2label.push_back(_bdLoops2lable[i][j]);
        }
    }

    for(auto &a:tmp_compactloop2Segs)sort(a.begin(),a.end());
    vector<int>newMergeLoop(_nLoops,-1);
    _nLoops = 0;


    for(int i=0;i<tmp_compactloop2Segs.size();++i)if(newMergeLoop[i]==-1){
        newMergeLoop[i] = _nLoops;
        compactloop2Segs.push_back(tmp_compactloop2Segs[i]);
        compactloop2label.push_back(tmp_compactloop2label[i]);
        for(int j=i+1;j<tmp_compactloop2Segs.size();++j){
            if(tmp_compactloop2label[i]==tmp_compactloop2label[j] && tmp_compactloop2Segs[i]==tmp_compactloop2Segs[j])newMergeLoop[j] = _nLoops;
        }
        _nLoops++;
    }

    for(auto &a:loop2compactloopN)for(auto &b:a)b = newMergeLoop[b];



}
void ParallelArrangement::MapLoopToSegs(int celli, vector< vector<int> >&_bdLoopE2V){

    if(loop2Segs.size()!=_nCell)loop2Segs.resize(_nCell);
    assert(celli<_nCell);

    auto& cloop2Segs = loop2Segs[celli];
    cloop2Segs.resize(_bdLoopE2V.size());

    auto &segs = newCtrNet.mixCurve.cell2Seg[celli];
    // auto &


    vector<int>usesegs(segs.size(),0);
    auto isContain = [this](int v1, int v2, vector<int>&candidateEvs){
        int nE = candidateEvs.size()/2;
        bool re = false;
        for(int i=0;i<nE;++i){
            if(candidateEvs[i*2] ==v1 && candidateEvs[i*2+1]==v2){
                re = true;
                break;
            }
        }
        return re;

    };
    auto &cur_cell2tetVIds = _cell2tetVIds[celli];
    for(int i=0;i<_bdLoopE2V.size();++i){
        auto mapbdloopEV2TetV = _bdLoopE2V[i];
        for(auto &a:mapbdloopEV2TetV)a = cur_cell2tetVIds[a];
        for(int j = 0;j<mapbdloopEV2TetV.size()/2;++j)if(mapbdloopEV2TetV[j*2]>mapbdloopEV2TetV[j*2+1])swap(mapbdloopEV2TetV[j*2],mapbdloopEV2TetV[j*2+1]);
        cloop2Segs[i].clear();



        for(int j=0;j<segs.size();++j){
            auto &seg2Tv = seg_dualedges2Tv[segs[j]];
            int e_num = seg2Tv.size()/2;
            bool isc = false;
            for(int k=0;k<e_num;++k){
                isc = isContain(seg2Tv[k*2],seg2Tv[k*2+1],mapbdloopEV2TetV);
                if(!isc)break;
            }

            if(isc){
                cloop2Segs[i].push_back(segs[j]);
                usesegs[j]++;
            }
        }


    }

    for(auto &a:cloop2Segs)compactloop2Segs.push_back(a);
    for(auto a:usesegs)assert(a==2);


    //vector< vector<uint> >seg2edges;
    //vector< vector<uint> >cell2Seg;

    //hidden_ctrActiveE2V



}

void ParallelArrangement::MapActiveFCP(int celli, vector< int >&_ActFV, vector<vector<float>>&_ActFCP, vector<int> &_ActNonE2V,bool isForMainpipe){


    if(FbdCtrV2ActiveF_percell.size()!=_nCell)FbdCtrV2ActiveF_percell.resize(_nCell);
    if(EbdCtrV2ActiveNonE_percell.size()!=_nCell)EbdCtrV2ActiveNonE_percell.resize(_nCell);

    if(numofActiveF_percell.size()!=_nCell)numofActiveF_percell.resize(_nCell);
    if(numofActiveNonE_percell.size()!=_nCell)numofActiveNonE_percell.resize(_nCell);


    vector<int>ActFTV;
    vector<int>ActNonE2TV;
    auto &cur_cell2tetVIds = _cell2tetVIds[celli];
    for(auto a:_ActFV)ActFTV.push_back(cur_cell2tetVIds[a]);
    for(auto a:_ActNonE2V)ActNonE2TV.push_back(cur_cell2tetVIds[a]);

    int nf = ActFTV.size()/3;
    int ne =  ActNonE2TV.size()/2;
    assert(_ActFCP.size()==nf+ne);


    numofActiveF_percell[celli] = nf;
    numofActiveNonE_percell[celli] = ne;

    for(int i=0;i<nf;++i){
        sort(ActFTV.begin()+i*3,ActFTV.begin()+i*3+3);
    }
    for(int i=0;i<ne;++i){
        sort(ActNonE2TV.begin()+i*2,ActNonE2TV.begin()+i*2+2);
    }


    auto matchF = [this](int *pf, vector<int>&candidateFvs){
        int nff = candidateFvs.size()/3;
        int re = -1;
        for(int kd=0;kd<nff;++kd){
            int ind = kd*3;
            if(candidateFvs[ind] ==pf[0] && candidateFvs[ind+1]==pf[1] && candidateFvs[ind+2]==pf[2]){
                re = kd;
                break;
            }
        }
        return re;

    };
    auto matchE = [this](int *pE, vector<int>&candidateEvs){
        int nev = candidateEvs.size()/2;
        int re = -1;
        for(int kd=0;kd<nev;++kd){
            int ind = kd*2;
            if(candidateEvs[ind] ==pE[0] && candidateEvs[ind+1]==pE[1] || candidateEvs[ind] ==pE[1] && candidateEvs[ind+1]==pE[0] ){
                re = kd;
                break;
            }
        }
        return re;

    };

    auto &FbdCtrV2ActiveF = FbdCtrV2ActiveF_percell[celli];
    FbdCtrV2ActiveF.resize(nf);

    for(int i=0;i<nf;++i){
        FbdCtrV2ActiveF[i] = matchF(ActFTV.data()+i*3,ActiveFV2Tv);
    }

    auto &EbdCtrV2ActiveNonE = EbdCtrV2ActiveNonE_percell[celli];
    EbdCtrV2ActiveNonE.resize(ne);

    for(int i=0;i<ne;++i){
        EbdCtrV2ActiveNonE[i] = matchE(ActNonE2TV.data()+i*2,ActiveNonEV2Tv);
    }



    auto &mapF2Vpos = mergeCs.mapActiveF2Vpos;
    for(int i=0;i<nf;++i)if(FbdCtrV2ActiveF[i]!=-1){


        copyVec(mapF2Vpos[FbdCtrV2ActiveF[i]].data(),_ActFCP[i].data());
        //        if(matchFFF!=0){
        //            copyVec(ActiveFCP.data()+FbdCtrV2ActiveF[i]*3,_ActFCP[i].data());
        //        }else{
        //            copyVec(_ActFCP[i].data(),ActiveFCP.data()+FbdCtrV2ActiveF[i]*3);
        //        }
        if(isForMainpipe){
            auto &matchFFF = ActiveFCP_setup[FbdCtrV2ActiveF[i]];
            matchFFF++;
            assert(matchFFF<3);
        }
    }

    auto &mapNonE2Vpos = mergeCs.mapActiveNonE2Vpos;
    for(int i=0;i<ne;++i){
        copyVec(mapNonE2Vpos[EbdCtrV2ActiveNonE[i]].data(),_ActFCP[i+nf].data());
    }

    //cout<<1<<endl;




}


void ParallelArrangement::GetCellTopo(int celli, vector<CellTopology>&cellTopo, vector<vector<int> > &label2bdLoops, vector<int> &cellNexplored2Topogroup){


    if(_cellTopo.size()!=_nCell)_cellTopo.resize(_nCell);
    if(_cellTopo_ori.size()!=_nCell)_cellTopo_ori.resize(_nCell);

    if(_label2bdLoops.size()!=_nCell)_label2bdLoops.resize(_nCell);
    if(_bdLoops2lable.size()!=_nCell)_bdLoops2lable.resize(_nCell);
    if(_cellNexplored2Topogroup.size()!=_nCell)_cellNexplored2Topogroup.resize(_nCell);

    _cellTopo[celli].clear();
    _cellNexplored2Topogroup[celli].clear();
    _cellTopo[celli] = cellTopo;
    _cellTopo_ori[celli] = cellTopo;
    _cellNexplored2Topogroup[celli] = cellNexplored2Topogroup;
    //_label2bdLoops[celli] = label2bdLoops;

    _label2bdLoops[celli].clear();
    _label2bdLoops[celli].resize(_nMat);
    auto& maplabels = _cell2labels[celli];

    for(int i=0;i<label2bdLoops.size();++i)_label2bdLoops[celli][maplabels[i]] = label2bdLoops[i];
    auto &Loops2lable = _bdLoops2lable[celli];
    Loops2lable.clear();
    Loops2lable.resize(loop2Segs[celli].size(),-1);

    for(int i=0;i<_nMat;++i)for(auto a:_label2bdLoops[celli][i]){assert(Loops2lable[a]==-1);Loops2lable[a] = i;}

    for(auto a:Loops2lable)assert(a!=-1);

}
void ParallelArrangement::PrintPickInfo(int celli, int topoi){
    if(celli<0 || celli>=_cellTopo.size())return;

    if(topoi<0 || topoi>_cellTopo[celli].size())return;
    CellTopology& cT = _cellTopo[celli][topoi];

    cout<<topoi<<" group: "<<_cellNexplored2Topogroup[celli][topoi]<<" junction points: "<<cT._junctionpoints<<" score: "<<cT._score<<endl;

}
void ParallelArrangement::GetCell(int celli, vector<vector<float>>& cellVs,vector<int>& cellVMarkers,
                                  vector<vector<int>>& cellTs, vector<int>& cellInitLabel,
                                  int& cellNMat) {
    if (celli >= _nCell) return;
    cellTs = _cell2tetTs[celli];
    cellInitLabel = _cell2labeling[celli];
    cellNMat = _cell2labels[celli].size();
    int nV = _cell2tetVIds[celli].size();
    cellVs.resize(nV);
    cellVMarkers.resize(nV);
    for (int i = 0; i < nV; ++i) {
        cellVs[i] = _tetVs[_cell2tetVIds[celli][i]];
    }
    for (int i = 0; i < nV; ++i) {
        cellVMarkers[i] = 0;
        for(auto a:_tetVMarkers_overlap[_cell2tetVIds[celli][i]])if(_framef2Cs[a]>-1)++cellVMarkers[i];
    }
    for (int i = 0; i < nV; ++i) {
        if(cellVMarkers[i] >1) cellVMarkers[i]=-1;
        else cellVMarkers[i] = 0;
    }
    for (int ti = 0; ti < cellTs.size(); ++ti) {
        for (int i = 0; i < 4; ++i) {
            cellTs[ti][i] = _cell2MapTetV2CellV[celli][cellTs[ti][i]];
        }
    }
}




typedef vector<vector<rgb>> imagecols;
void ParallelArrangement::MakeContourfromImage(const vector<string> &filenames, string outfilename, double zdis){

    _nCrossSection = filenames.size();
    int width = 0;
    int height = 0;

    // --------------------- Load cross sections ----------------------
    // Load all cross section images and save the rgb value of each pixel.
    vector<imagecols> image2cols(_nCrossSection);
    set<rgb> colsSet;
//    for (int i = 0; i < _nCrossSection; ++i) {
//        cimg_library::CImg<float> cimg(filenames[i].c_str());
//        auto& oneimage = image2cols[i];
//        width = cimg.width();
//        height = cimg.height();
//        oneimage.resize(width, vector<rgb>(height, vector<int>(3)));
//        for (int x = 0; x < width; ++x) {
//            oneimage[x].resize(height);
//            for (int y = 0; y < height; ++y) {
//                for (int c = 0; c < 3; ++c) {
//                    oneimage[x][y][c] = cimg(x, y, 0, c);
//                }
//                // MingUtility::printVector(oneimage[x][y]);
//                colsSet.insert(oneimage[x][y]);
//            }
//        }
//    }
    // cout << "width=" << width << endl;
    // cout << "height=" << height << endl;
    // Sort all the colors used in this data set, the index of each color is the
    // label of corresponding material in follow up algorithm.
    // TODO here I assume the background color is {0,0,0}, so it will be the 0
    // label. Let the background color be the first label is important because the
    // cell facets with out image constraint is labeled as 0 in my code. If other
    // colors are used as the background, change the sort to a mammual indexing
    // and make the background color to be the first one in _cols. Or change the
    // rgb2label map.
    _cols.resize(colsSet.size());
    copy(colsSet.begin(), colsSet.end(), _cols.begin());
    sort(_cols.begin(), _cols.end());
    map<rgb, int> rgb2label;
    for (int i = 0; i < _cols.size(); ++i) {
        rgb2label[_cols[i]] = i;
    }
    //rgb2label[_cols[_cols.size()-1]] = 2;
    _nMat = _cols.size();

    vector<vector<vector<int>>>fbuffer(_nCrossSection);

    for (int i = 0; i < _nCrossSection; ++i) {

        auto& oneimage = image2cols[i];
        auto& onebuffer = fbuffer[i];


        onebuffer.resize(width);
        for (int x = 0; x < width; ++x) {
            onebuffer[x].resize(height,-1);
            for (int y = 0; y < height; ++y) {
                onebuffer[x][y] =  rgb2label[oneimage[x][y]];
            }
        }
    }



    n_rf::CrossSections CS;
    CS.ReadCrossSectionsfromImages(fbuffer,zdis);

    CS.ManipulateCrossSections();


    CS.WriteCrossSections(outfilename);


    //CS.RescaleToUniform();
    //string fpath("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/Draw/");
    //CS.WriteCrossSectionsToFCM(fpath+string("tctrloop.txt"),3);

    //CS.WriteCrossSectionsToFCMtCtr(fpath+string("tctr.txt"));











}

void ParallelArrangement::CreateTotalSurf(vector<vector<double>>&TopoVs, vector<vector<uint>>&TopoFs, vector<vector<int>>&TopoFMs, vector<vector<uint> > &TopoCtrs,
                                          bool isOutputFiles,n_rf::SufStructure* pSuf){

    //    if(FbdCtrV2ActiveF_percell.size()!=_nCell)FbdCtrV2ActiveF_percell.resize(_nCell);
    //    if(EbdCtrV2ActiveFpair_percell.size()!=_nCell)EbdCtrV2ActiveFpair_percell.resize(_nCell);


    //return;
    int totalv = 0;
    vector<double>Vs;
    vector<uint>Fs;
    vector<int>FMs;
    vector<int>Ctrs;




    int VInd = 0;
    vector<int>pickActiveF(ActiveFV2Tv.size()/3,-1);
    //for(auto &b:FbdCtrV2ActiveF_percell)for(auto a:b)pickActiveF[a] = 0;
    for(int i=0;i<pickActiveF.size();++i){
        pickActiveF[i] = VInd++;
        //auto p_v = ActiveFCP.data()+i*3;
        //for(int j=0;j<3;++j)ctrVs.push_back(p_v[j]);
    }

    cout<<VInd<<endl;


    vector<int>mapActiveNonE2PInd(mergeCs.mapActiveNonE2Vpos.size());
    for(auto &a:mapActiveNonE2PInd)a = VInd++;
    //    for(auto &a:mergeCs.mapActiveNonE2Vpos){
    //        mapActiveNonE2PInd[a.first] = VInd++;
    //    }
    //for(auto &b:EbdCtrV2ActiveFpair_percell)for(auto a:b)pickActiveF[a] = 0;


    int ctrVnum = VInd;
    vector<double>ctrVs(VInd*3);
    vector<bool>ispickctrVs(VInd,false);

    cout<<VInd<<endl;

    vector<vector<int>>mappingsToGlobal(_nCell);


    for(int i=0;i<_nCell;++i){
        auto &FbdCtrV2ActiveF  = FbdCtrV2ActiveF_percell[i];
        auto &EbdCtrV2ActiveNonE  = EbdCtrV2ActiveNonE_percell[i];


        auto &Vpos = TopoVs[i];
        auto &F2V = TopoFs[i];
        auto &FM = TopoFMs[i];
        //auto &Ctr = TopoCtrs[i];
        int nv = Vpos.size()/3;
        vector<int>&vindtable = mappingsToGlobal[i];
        vindtable.clear();
        vindtable.resize(nv,-1);

        int nf = numofActiveF_percell[i];
        int ne = numofActiveNonE_percell[i];
        for(int j=0;j<nf;++j){

            vindtable[j] = pickActiveF[FbdCtrV2ActiveF[j]];
            if(!ispickctrVs[vindtable[j]]){
                ispickctrVs[vindtable[j]] = true;
                copyVec(Vpos.data()+j*3,ctrVs.data()+vindtable[j]*3);
            }
        }
        for(int j=0;j<ne;++j){
            int indj = j+nf;
            vindtable[indj] = mapActiveNonE2PInd[EbdCtrV2ActiveNonE[j]];
            if(!ispickctrVs[vindtable[indj]]){
                ispickctrVs[vindtable[indj]] = true;
                copyVec(Vpos.data()+indj*3,ctrVs.data()+vindtable[indj]*3);
            }
        }
        for(int j=nf+ne;j<nv;++j){
            auto p_v = Vpos.data()+j*3;
            for(int k=0;k<3;++k){
                Vs.push_back(p_v[k]);
            }
            vindtable[j] = VInd++;
        }

        for(auto a:F2V)Fs.push_back(vindtable[a] );
        //for(auto a:Ctr)Ctrs.push_back(vindtable[a] );
        FMs.insert(FMs.end(),FM.begin(),FM.end());
    }

    //writeSufFile(string filename, const vector<double>&vertices, const vector<unsigned int>&faces2vertices, const vector<int>&facesMat, const vector<int> &CtrEdges);

    ctrVs.insert(ctrVs.end(),Vs.begin(),Vs.end());

    if(isOutputFiles)writeSufFile(string("../output"),ctrVs,Fs,FMs,mergeCs.hidden_writeCtrE);
    //writeOffFile(string("../output"),ctrVs,Fs);
    //writeObjFile(string("../output"),ctrVs,Fs);
    if(pSuf!=NULL)pSuf->Construct(ctrVs,Fs,FMs,mergeCs.hidden_writeCtrE,mappingsToGlobal);


    vector<vector<int>>CtrEperCs(_nCrossSection);
    vector<vector<int>>CtrEMperCs(_nCrossSection);
    vector<vector<uint>>uCtrEperCs(_nCrossSection);

    auto &hidden_writeCtrE = mergeCs.hidden_writeCtrE;
    auto &hidden_ctrELable = mergeCs.hidden_ctrELable;
    for(int i=0;i<mergeCs.hidden_CtrE2CSs.size();++i){
        int c = mergeCs.hidden_CtrE2CSs[i];
        CtrEperCs[c].push_back(hidden_writeCtrE[i*2]);
        CtrEperCs[c].push_back(hidden_writeCtrE[i*2+1]);
        CtrEMperCs[c].push_back(hidden_ctrELable[i*2]);
        CtrEMperCs[c].push_back(hidden_ctrELable[i*2+1]);
        uCtrEperCs[c].push_back(hidden_writeCtrE[i*2]);
        uCtrEperCs[c].push_back(hidden_writeCtrE[i*2+1]);

    }
    vector<float>ctrctrVs;
    for(int i=0;i<ctrVnum*3;++i)ctrctrVs.push_back(ctrVs[i]);
    vector<vector<double>>Cs2PlaneParas(_Cs2PlaneParas.size());
    for(int i=0;i<_Cs2PlaneParas.size();++i)for(auto b:_Cs2PlaneParas[i])Cs2PlaneParas[i].push_back(b);

    if(isOutputFiles){
        n_rf::CrossSections CS;
        CS.ReadCrossSections(ctrVs,Cs2PlaneParas,uCtrEperCs,CtrEMperCs);
        CS.WriteCrossSections(outfoldername_meta + string("ctr"));
        writeCtrGraphFile(outfoldername_meta + string("ctr"),ctrctrVs,CtrEperCs,CtrEMperCs,_Cs2PlaneParas);
        OutputNewCtrs(outfoldername_meta + string("ctr"),mergeCs.hidden_ctrV,mergeCs.hidden_ctrE);
        OutputMappingInfo(outfoldername_meta);
    }

//    if(isOutputFiles){
//        n_rf::CrossSections CS;
//        CS.ReadCrossSections(ctrVs,Cs2PlaneParas,uCtrEperCs,CtrEMperCs);
//        CS.WriteCrossSections(string("../mappingInfo/ctr"));
//        writeCtrGraphFile(string("../mappingInfo/ctr"),ctrctrVs,CtrEperCs,CtrEMperCs,_Cs2PlaneParas);
//        OutputNewCtrs(string("../mappingInfo/ctr"),mergeCs.hidden_ctrV,mergeCs.hidden_ctrE);
//        OutputMappingInfo(string("../mappingInfo/"));
//    }



    //    auto &hidden_CtrE = mergeCs.hidden_ctrE;
    //    for(int i=0;i<hidden_CtrE.size();++i){
    //        if(hidden_CtrE[i]!=hidden_writeCtrE[i])cout<<"adsada"<<endl;
    //    }




}

void ParallelArrangement::SetMetaFolder(string meta_folder){
    outfoldername_meta = meta_folder;
}

void ParallelArrangement::WriteToMappingInfo(){





     OutputNewCtrs(outfoldername_meta+string("ctr"),mergeCs.hidden_ctrV,mergeCs.hidden_ctrE);
     OutputMappingInfo(outfoldername_meta);




}

void ParallelArrangement::OutputNewCtrs(string filename,const vector<double>&vertices,const vector<uint>&edge2vertices){

    filename = filename + ".obj";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output Obj file " << filename << endl;
        //return false;
    }



    int n_vertices = vertices.size()/3;

    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }


    {
        int n_edges = edge2vertices.size()/2;
        for(int i=0;i<n_edges;++i){
            auto p_ev = edge2vertices.data()+i*2;
            outer << "l " << p_ev[0]+1<< " "<< p_ev[1]+1 << endl;
        }
    }


    outer.close();
    cout<<"saving finish: "<<filename<<endl;
    //return true;
}

void makemethodfolder(string foldername);
void ParallelArrangement::OutputMappingInfo(string filepath){

    makemethodfolder(filepath + "loop2Segs/");
    makemethodfolder(filepath + "topoInfo/");
    makemethodfolder(filepath + "topoInfoVVec/");

    auto& seg2edges = newCtrNet.mixCurve.seg2edges;

    string filename = filepath + "seg2edges";
    writeVVecFile(filename,seg2edges);



    for(int celli=0;celli<_nCell;++celli){
        auto& cloop2Segs = loop2Segs[celli];
        string filename = filepath + "loop2Segs/cell"+to_string(celli);
        writeVVecFile(filename,cloop2Segs);
    }
    //cout<<"0"<<endl;
    filename = filepath + "loop2label";
    writeVVecFile(filename,_bdLoops2lable);

    //cout<<"1"<<endl;
    string prepath = filepath + "topoInfo/info_";
    _cellTopo_ori;

    vector<vector<double>>score(_cellTopo_ori.size());
    for(int celli=0;celli<_cellTopo_ori.size();++celli){
        string prepath2 = prepath+to_string(celli)+string("_");
        cout<<prepath2<<endl;
        score[celli].resize(_cellTopo_ori[celli].size());
        for(int topoi = 0;topoi<_cellTopo_ori[celli].size();++topoi){
            _cellTopo_ori[celli][topoi].writetxt(prepath2+to_string(topoi));
            score[celli][topoi] = _cellTopo_ori[celli][topoi]._score;
        }

    }

    prepath = filepath + "topoInfoVVec/info_";
    writeVVecFile(prepath+"score",score);

    for(int celli=0;celli<_cellTopo_ori.size();++celli){
        string prepath2 = prepath+to_string(celli)+string("_");
        //cout<<prepath2<<endl;
        score[celli].resize(_cellTopo_ori[celli].size());
        for(int topoi = 0;topoi<_cellTopo_ori[celli].size();++topoi){
            writeVVecFile(prepath2+to_string(topoi),_cellTopo_ori[celli][topoi]._comps);
        }

    }








}
