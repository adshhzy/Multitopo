
#include "Arrangement.h"
#include "../../Utility/utility.h"
void Arrangement::constructArrangement(
        int ssvernum, int ssedgenum, int ssfacenum, int ssspacenum,
        float *ssver, int *ssedge, int **ssface, int **ssspace,
        int *ssfaceedgenum, int *ssspacefacenum, int **ssspace_planeside, int *ssface_planeindex,

        int planenum, float *pparam, float cellFrame [],

        int *ctrfvernum, int *ctrfedgenum,
        float **ctrfverpos, int **ctrfvertype, int **ctrfverval,
        int **ctrfedge, int **ctrfedgetype, int **ctrfedgeval
        ){
    /* -------------------------------------------------------------------- */
    /* ------------------- MING: construct Arrangement -------------------- */
    /* Section 1: Construct Vs, seg2V, F2seg */
    vector<vector<float> > verOnCellE;
    vector<vector<int> > verIdOnCellE;
    vector<float> oneV(3, 0);
    vector<int> oneE(2, 0);
    vector<int> oneSeg;
    vector<int> segEnds;
    vector<int> ind2ind;
    vector<vector<int> > nbVs; // will be empty after tracing the segs, do not plan to use it afterwards
    vector<vector<int> > nbVsLeftMaterial; // will be empty after tracing the segs, do not plan to use it afterwards
    set<int> edgeMats;
    int matchVi, verOnCellEi, v0, v1;
    verOnCellE.resize(ssedgenum);
    verIdOnCellE.resize(ssedgenum);
    F2seg.resize(ssfacenum);
    nF = ssfacenum;
    nSegV = 0;

    vector<int> toDelSegVs;

    for(int i=0; i<3; ++i){
        cellFrameLowerCorner[i] = cellFrame[i];
        cellFrameScale[i] = cellFrame[3+i]-cellFrame[i];
    }
    for (int fi = 0; fi < ssfacenum; ++fi){
        if (ctrfvernum[fi] == 0){
            continue;
        }
        segEnds.clear();
        ind2ind.clear();
        ind2ind.resize(ctrfvernum[fi]);
        for (int vi = 0; vi < ctrfvernum[fi]; ++vi){

            /*Step 1: store all V in current face into arrangement Vs*/
            oneV[0] = ctrfverpos[fi][3 * vi + 0];
            oneV[1] = ctrfverpos[fi][3 * vi + 1];
            oneV[2] = ctrfverpos[fi][3 * vi + 2];
            // ver is inside of the face completely
            if (ctrfvertype[fi][vi] == /*CTRTYPE_SUBFACE*/3){
                ind2ind[vi] = (int)Vs.size();
                Vs.push_back(oneV);
                V2seggV.push_back(-1);
            }
            // ver is on some cell edge
            else {
                segEnds.push_back(vi);
                verOnCellEi = ctrfverval[fi][vi];
                matchVi = findMatchVer(verOnCellE[verOnCellEi], oneV);
                // a new intersection ver is found
                if (matchVi == -1){
                    ind2ind[vi] = (int)Vs.size();
                    Vs.push_back(oneV);
                    verIdOnCellE[verOnCellEi].push_back(ind2ind[vi]);
                    verOnCellE[verOnCellEi].push_back(oneV[0]);
                    verOnCellE[verOnCellEi].push_back(oneV[1]);
                    verOnCellE[verOnCellEi].push_back(oneV[2]);
                    V2seggV.push_back(nSegV++);
                }
                // the intersection ver is processed before
                else{
                    ind2ind[vi] = verIdOnCellE[verOnCellEi][matchVi];
                }
            }
        }
        /*Step 2: trace segment out from stored edges*/

        int leftMaterial;
        nbVs.clear();
        nbVsLeftMaterial.clear();
        nbVs.resize(ctrfvernum[fi]);
        nbVsLeftMaterial.resize(ctrfvernum[fi]);
        // go through all edges, store the neighborhood
        for (int ei = 0; ei < ctrfedgenum[fi]; ++ei){
            v0 = ctrfedge[fi][4 * ei + 0];
            v1 = ctrfedge[fi][4 * ei + 1];
            nbVs[v0].push_back(v1);
            nbVs[v1].push_back(v0);
            nbVsLeftMaterial[v0].push_back(ctrfedge[fi][4 * ei + 2]);
            nbVsLeftMaterial[v1].push_back(ctrfedge[fi][4 * ei + 3]);
            edgeMats.insert(ctrfedge[fi][4 * ei + 2]);
            edgeMats.insert(ctrfedge[fi][4 * ei + 3]);
        }
        // for now, we only consider 2 materials: in/out
        matMarks = vector<int>(edgeMats.begin(), edgeMats.end());

        for(auto a:edgeMats)cout<<a<<' ';cout<<endl;
        // assumption here: segments in each face are disjoint completely
        // so 2 nbs for interior ver; 1 for end ver;
        // handle open segment first
        vector<bool> unvisited(ctrfvernum[fi], true);
        for (int end : segEnds){
            if (!nbVs[end].empty()){
                // a new open segment is found
                traceOneSeg(nbVs, nbVsLeftMaterial, ind2ind, unvisited, oneSeg, end, leftMaterial);

                // handle the bug where the newly added intersection points are too close to its nb, so tetgen will crash later
                // remove the close nbV(s)
                cleanTwoEndsOfOneSeg(oneSeg, toDelSegVs);

                F2seg[fi].push_back((int)seg2V.size());
                seg2F.push_back(fi);

                seg2V.push_back(oneSeg);
                seg2lmat.push_back(leftMaterial);
            }
        }
        cout<<"handle closed segment"<<endl;
        // handle closed segment
        int end;
        for (int vi = 0; vi < ctrfvernum[fi]; ++vi){
            // a new closed segment is found
            if (unvisited[vi]){
                end = vi;
                traceOneSeg(nbVs, nbVsLeftMaterial, ind2ind, unvisited, oneSeg, end, leftMaterial);
                F2seg[fi].push_back((int)seg2V.size());
                seg2F.push_back(fi);
                seg2V.push_back(oneSeg);
                seg2lmat.push_back(leftMaterial);
            }
        }
        cout<<"facet: "<<fi<<endl;

    }
    nSeg = (int)seg2V.size();
    for(size_t i=0; i<nSeg; ++i){
        dealDegeneratedOneSegWithTwoV(seg2V[i]);
    }

    if(!toDelSegVs.empty()){
        sort(toDelSegVs.begin(), toDelSegVs.end());
        // seg2V
        for(size_t i=0; i<seg2V.size(); ++i){
            for(size_t j=0; j<seg2V[i].size(); ++j){
                seg2V[i][j] = seg2V[i][j] - Utility::findNEleSmallerInSortedList(toDelSegVs, seg2V[i][j]);
            }
        }
        // verIdOnCellE
        for(size_t i=0; i<verIdOnCellE.size(); ++i){
            for(size_t j=0; j<verIdOnCellE[i].size(); ++j){
                verIdOnCellE[i][j] = verIdOnCellE[i][j] - Utility::findNEleSmallerInSortedList(toDelSegVs, verIdOnCellE[i][j]);
            }
        }
        // V2seggv
        vector<int> new_V2seggV(Vs.size()-toDelSegVs.size(),-1);
        for(int vi=0; vi<V2seggV.size(); ++vi){
            if(V2seggV[vi]!=-1){
                new_V2seggV[vi - Utility::findNEleSmallerInSortedList(toDelSegVs, vi)] = V2seggV[vi];
            }
        }
        V2seggV = new_V2seggV;
        // Vs
        for(int i=toDelSegVs.size()-1; i>=0; i--){
            Vs.erase(Vs.begin()+toDelSegVs[i]);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------

    cout<<"Section 2: Construct segGraph"<<endl;
    /* Section 2: Construct segGraph */
    segGraph.resize(nSegV);
    segGraph_b.resize(nSegV, BigInt(nSeg));

    for (int si = 0; si < nSeg; ++si){
        v0 = V2seggV[seg2V[si].front()];
        v1 = V2seggV[seg2V[si].back()];
        if (v0 != v1){
            segGraph[v0].push_back(si);
            segGraph[v1].push_back(si);
            segGraph_b[v0].set(si);
            segGraph_b[v1].set(si);
        }
    }

    /* Section 3: Construct cell2F, cellEmpty, cellGraph*/
    nCell = ssspacenum;
    nFront = nCell;
    nSeg = (int)seg2V.size();
    cell2F.resize(nCell);
    cell2Fside.resize(nCell);
    cellEmpty.resize(nCell);
    cellGraph.resize(nCell);
    cell2seg.resize(nCell);
    cell2seg_b.resize(ssspacenum, BigInt(nSeg));
    vector<vector<int> > cellsOfF;
    cellsOfF.resize(ssfacenum);
    bool isempty;
    for (int ci = 0; ci < nCell; ++ci){
        isempty = true;
        for (int fi = 0; fi < ssspacefacenum[ci]; ++fi){
            cell2F[ci].push_back(ssspace[ci][fi]);
            cell2Fside[ci].push_back(ssspace_planeside[ci][fi]);
            cellsOfF[ssspace[ci][fi]].push_back(ci);
            if (!F2seg[ssspace[ci][fi]].empty()){
                isempty = false;
            }
        }
        cellEmpty[ci] = isempty;
    }
    int c1, c2;
    for (int fi = 0; fi < ssfacenum; ++fi){
        F2PlaneID.push_back(ssface_planeindex[fi]);
        if (cellsOfF[fi].size() != 2) continue;
        c1 = cellsOfF[fi][0];
        c2 = cellsOfF[fi][1];
        cellGraph[c1].push_back(c2);
        cellGraph[c2].push_back(c1);

        for (int si : F2seg[fi]){
            cell2seg_b[c1].set(si);
            cell2seg_b[c2].set(si);
            cell2seg[c1].push_back(si);
            cell2seg[c2].push_back(si);
        }
    }

    vector<float> para;
    para.resize(4);
    for (int pi=0; pi<planenum; ++pi){
        para[0] = pparam[pi * 4];
        para[1] = pparam[pi * 4 + 1];
        para[2] = pparam[pi * 4 + 2];
        para[3] = pparam[pi * 4 + 3];
        planeParas.push_back(para);
    }

    /* Section 4: Construct cellVs and cellEs, for debug*/
    for (int vi = 0; vi < ssvernum; ++vi){
        oneV[0] = ssver[3 * vi + 0];
        oneV[1] = ssver[3 * vi + 1];
        oneV[2] = ssver[3 * vi + 2];
        cellVs.push_back(oneV);
    }
    list<pair<float, int> > sorted_VnDistOnCellE;
    sorted_verIdOnCellE.resize(ssedgenum);
    for (int ei = 0; ei < ssedgenum; ++ei){
        oneE[0] = ssedge[2 * ei + 0];
        oneE[1] = ssedge[2 * ei + 1];
        cellEs.push_back(oneE);

        sorted_VnDistOnCellE.clear();
        for (int vi : verIdOnCellE[ei]){
            sorted_VnDistOnCellE.push_back(make_pair(EuclideanDistSquare1Cell2Seg(oneE[0], vi), vi));
        }
        sorted_VnDistOnCellE.sort();
        sorted_verIdOnCellE[ei].reserve(sorted_VnDistOnCellE.size());
        for (auto it = sorted_VnDistOnCellE.begin(); it != sorted_VnDistOnCellE.end(); ++it){
            sorted_verIdOnCellE[ei].push_back(it->second);
        }
    }

    cellFs.resize(ssfacenum);
    for (int fi = 0; fi < ssfacenum; ++fi){
        for (int ei = 0; ei < ssfaceedgenum[fi]; ++ei){
            cellFs[fi].push_back(ssface[fi][ei]);
        }
    }
    reOrderCellEinCellF();

    findCloseCyclesInCells();

    /* -------------------------------------------------------------------- */
#if _SAVEDATA_MING
    writeArrangement((t_savedir+"arrangement_beforeCell.txt").c_str());
#endif
    moveSegPtsBackToPlane();
    moveCellPtsBackToPlane();
#if _SAVEDATA_MING
    writeArrangement((t_savedir+"arrangement_beforeCell2.txt").c_str());
#endif

#if _SAVERES_MING
    writeLuCtrGraph((t_savedir+"forLu.ctrGraph").c_str());
#endif

    cout<<"Section 5: explor topology within each cell, tetgen + random walk + critical point"<<endl;

    /*  Section 5: explor topology within each cell, tetgen + random walk + critical point*/

    tetV.resize(nCell);
    tetVh.resize(nCell);
    tetVmarker.resize(nCell);
    tetnV.resize(nCell);
    tetQ.resize(nCell);
    allGroupings.resize(nCell);
    allGroupingScore.resize(nCell);
    allGroupingIso.resize(nCell);
    cell2Loc_setIso.resize(nCell);
    cell2Loc_Global_setMap.resize(nCell);
    cell2Loc_setActEdge.resize(nCell);
    bitIntHash allSets_b;
    maxNCyc = -1;
    for (int ci = 0; ci < nCell; ++ci){
        if (cell2ncyc[ci] > maxNCyc){
            maxNCyc = cell2cyc[ci].size();
        }
    }
    nSets = 0;
    allCellTets.init(nCell);
#if _DO_PARALLEL_
#pragma omp parallel for default(shared)// num_threads(8)
#endif

#if _CG_LOG
    clock_t start, end;
    start = clock();
#endif
    int n_record = 0;
    for (int ci = 0; ci < nCell; ++ci){
        if (cellEmpty[ci] ) continue;
#if _CG_LOG

        cout << "   cell ["<< ci << "] : BEGIN ..." << endl;
#endif
        explorTopologyInCell(ci, allSets_b);
#if _CG_LOG
        cout << "   cell ["<< ci << "] : DONE!" << endl << "   -----------------------------------" << endl;
        cout << "   # of candidates: " << allGroupings[ci].size()<< endl;
        end = clock();
        cout << "      Time Used Until Now: "<< (end - start) / CLOCKS_PER_SEC << "s " << endl << "   -----------------------------------" << endl;
#endif
    }
    nSets = allSets_b.size();
    allSets.resize(nSets);
    for (auto set = allSets_b.begin(); set!=allSets_b.end(); ++set){
        allSets[set->second]=(set->first.getOnesInd());
    }
    int nTets = 0;
    for(int i=0; i<nCell; ++i){
        nTets += allCellTets.cells[i].cellTets.size();
    }
    cout << "nTets = " << nTets << endl;
#if _SAVEDATA_MING
    Utility::writeVector(allSets, (t_savedir+"allSets.txt").c_str());
#endif
}


void Arrangement::moveIntersectionPtsBackToPlane(){
    set<int> intsVs;
    unordered_map<int, set<int> > intsV2PlaneIDs;
    vector<int> vFrontBack(2);
    set<int> tmpPlaneIDs;
    for(size_t i=0; i<seg2V.size(); ++i){
        vFrontBack[0] = seg2V[i].front();
        vFrontBack[1] = seg2V[i].back();
        if(vFrontBack[0]!=vFrontBack[1]){
            for(int j=0; j<2; ++j){
                // a new intersection V is found
                if(intsVs.find(vFrontBack[j])==intsVs.end()){
                    intsVs.insert(vFrontBack[j]);
                }
                tmpPlaneIDs = intsV2PlaneIDs[vFrontBack[j]];
                tmpPlaneIDs.insert(F2PlaneID[seg2F[i]]);
                intsV2PlaneIDs[vFrontBack[j]] = tmpPlaneIDs;
            }
        }
    }
    for(auto vi : intsVs){
        Eigen::Vector3f curV(Vs[vi][0],Vs[vi][1],Vs[vi][2]);
        tmpPlaneIDs = intsV2PlaneIDs[vi];
        for(auto pi : tmpPlaneIDs){
            Eigen::Vector3f n(planeParas[pi][0],planeParas[pi][1],planeParas[pi][2]);
            float d = -planeParas[pi][3] / n.norm();
            n.normalize();
            curV = curV - (curV.dot(n)+d)*n;
        }
        Vs[vi][0] = curV[0];
        Vs[vi][1] = curV[1];
        Vs[vi][2] = curV[2];
    }
}


void Arrangement::writeLuCtrGraph(const char* filename){
    vector<vector<int> > plane2segs(planeParas.size());
    for(size_t segi=0; segi<seg2V.size(); ++segi){
        int pi = F2PlaneID[seg2F[segi]];
        if (pi < 0) continue;
        plane2segs[pi].push_back(segi);
    }
    ofstream ofs(filename, ofstream::out);
    if (!ofs.good()) {
        cout << "Can not write ctrGraph file " << filename << endl;
        return;
    }
    // Vs
    vector<vector<float> > scaledVs = Vs;
    Utility::scaleDnVers(scaledVs,unifyScaleRatio, unifyTransVec);
    ofs << "n " << scaledVs.size() << "\n";
    for( auto & oneV : scaledVs){
        ofs << "v " << oneV[0] << " " << oneV[1] << " " << oneV[2] << endl;
    }
    // Planes
    ofs << "n " << plane2segs.size() << "\n";
    for(size_t pi=0; pi<plane2segs.size(); ++pi){
        int nSegEdge = 0;
        for(auto segi : plane2segs[pi]){
            nSegEdge += (seg2V[segi].size()-1);
        }
        ofs << "p " << planeParas[pi][0] << " " << planeParas[pi][1] << " " << planeParas[pi][2] << " " << planeParas[pi][3] << " " << nSegEdge << endl;
        for(auto segi : plane2segs[pi]){
            int leftMat = seg2lmat[segi];
            int rightMat = TagMapping::getOtherMat(leftMat, matMarks);
            for(size_t vi=0; vi<seg2V[segi].size()-1; ++vi){
                ofs << "e " << seg2V[segi][vi] << " " << seg2V[segi][vi+1] << " " << leftMat << " " << rightMat << endl;
            }
        }
    }
    ofs.close();
}

void Arrangement::moveSegPtsBackToPlane(){
    for(size_t segi=0; segi<seg2V.size(); ++segi){
        for(auto vi : seg2V[segi]){
            Eigen::Vector3f curV(Vs[vi][0],Vs[vi][1],Vs[vi][2]);
            int pi = F2PlaneID[seg2F[segi]];
            Eigen::Vector3f n(planeParas[pi][0],planeParas[pi][1],planeParas[pi][2]);
            float d = -planeParas[pi][3] / n.norm();
            n.normalize();
            curV = curV - (curV.dot(n)+d)*n;
            Vs[vi][0] = curV[0];
            Vs[vi][1] = curV[1];
            Vs[vi][2] = curV[2];
        }
    }
}


void Arrangement::moveCellPtsBackToPlane(){
    vector<set<int> > CellV2PlaneIDs(cellVs.size());
    int pInd;
    for(int ci=0; ci<cell2F.size(); ++ci){
        for(auto fInd : cell2F[ci]){
            if(fInd < 0) continue;

            pInd = F2PlaneID[fInd];
            for(auto eInd : cellFs[fInd]){
                CellV2PlaneIDs[cellEs[eInd][0]].insert(pInd);
                CellV2PlaneIDs[cellEs[eInd][1]].insert(pInd);
            }
        }
    }
    for(int vi=0; vi<cellVs.size(); ++vi){
        Eigen::Vector3f curV(cellVs[vi][0],cellVs[vi][1],cellVs[vi][2]);
        for(auto pi : CellV2PlaneIDs[vi]){
            if(pi < 0 ) continue;
            Eigen::Vector3f n(planeParas[pi][0],planeParas[pi][1],planeParas[pi][2]);
            float d = -planeParas[pi][3] / n.norm();
            n.normalize();
            curV = curV - (curV.dot(n)+d)*n;
            cellVs[vi][0] = curV[0];
            cellVs[vi][1] = curV[1];
            cellVs[vi][2] = curV[2];
        }
    }
}




void outputCurve(string filename,vector<float>&vv,vector<int>&ee);
void projectAllOnEdgeVertices(vector<float> &vpos, vector<vector<int> > &v2ps, float *pparam);
void Arrangement::constructArrangementMM(
        int ssvernum, int ssedgenum, int ssfacenum, int ssspacenum,
        float *ssver, int *ssedge, int **ssface, int **ssspace,
        int *ssfaceedgenum, int *ssspacefacenum, int **ssspace_planeside, int *ssface_planeindex,

        int planenum, float *pparam, float pbbox [],

        int *ctrfvernum, int *ctrfedgenum,
        float **ctrfverpos, int **ctrfvertype, int **ctrfverval,
        int **ctrfedge, int **ctrfedgetype, int **ctrfedgeval
        ){



    /* ------------------- phase 3: partiting the space -------------------- */
    /*tss, ss: (tmp) subspace*/
    /*divide the 3D bounding box into subspaces, each ss has several oriented faces, each face is part of some plane*/
    /*here V and E are the ones from the bbox, but faces are later shared with curves*/
    //int ssvernum;				//# of all V
    //float* ssver;				//coord of all V
    //int ssedgenum;			//# of all E
    //int* ssedge;				//two ends of all E
    //int ssfacenum;			//# of all F
    //int* ssfaceedgenum;		//# of Es in each F
    //int** ssface;				//E ids of each face
    //int* ssface_planeindex;	//plane ID cooresponding each face
    //int ssspacenum;			//# of spaces
    //int* ssspacefacenum;		//# of faces in each space
    //int** ssspace;			//F ids of each space
    //int** ssspace_planeside;	//which side of the oriented plane is each face in each space

    streambuf* orig_buf = cout.rdbuf();
    cout.rdbuf(NULL);

    bool isOnEdge = true;


    // ----------------------- Generate tets ------------------------
    // Find the centroid of each cell for classifying tets in to cells.
    vector<vector<float> > cellCenter(ssspacenum, vector<float>(3,0));
    vector<vector<int> > cellVind(ssspacenum, vector<int>(0));
    for(int i=0;i<ssspacenum;++i){
        vector<bool>vpick(ssvernum,false);
        for(int j=0;j<ssspacefacenum[i];++j){
            int fff = ssspace[i][j];
            for(int k=0;k<ssfaceedgenum[fff];++k){
                int eee = ssface[fff][k];
                vpick[ssedge[2*eee]]=true;
                vpick[ssedge[2*eee+1]]=true;
            }
        }
        for(int j=0;j<ssvernum;++j)if(vpick[j])cellVind[i].push_back(j);
        for(auto a:cellVind[i]){
            for(int j=0;j<3;++j)cellCenter[i][j]+=ssver[3*a+j];
        }

        for(int j=0;j<3;++j)cellCenter[i][j]/=cellVind[i].size();
        //for(int j=0;j<3;++j)cout<<cellCenter[i][j]<<' ';cout<<endl;
        //for(auto a:cellVind[i])cout<<a<<' ';cout<<endl;
    }
    //cout<<endl;

    cout<<"ssspacenum "<<ssspacenum<<endl;
    vector<vector<int> > faceVind(ssfacenum, vector<int>(0));
    for(int i=0;i<ssfacenum;++i){
        set<int>fvinds;
        vector<int>fvindv;
        for(int j=0;j<ssfaceedgenum[i];++j){
            int eee = ssface[i][j];
            fvinds.insert(ssedge[2*eee]);
            fvinds.insert(ssedge[2*eee+1]);
            //cout<<ssedge[2*eee]<<' '<<ssedge[2*eee+1]<<' ';
        }
        //cout<<endl;
        for(auto iter=fvinds.begin();iter!=fvinds.end();++iter)fvindv.push_back(*iter);
        int curv = fvindv[0];faceVind[i].push_back(curv);
        //cout<<fvindv.size()<<endl;
        vector<bool>eflag(ssfaceedgenum[i],false);
        for(int k=0;k<fvindv.size()-1;++k){
            for(int j=0;j<ssfaceedgenum[i];++j){
                int eee = ssface[i][j];
                if(eflag[j])continue;

                if(curv==ssedge[2*eee]){
                    faceVind[i].push_back(ssedge[2*eee+1]);
                    curv = ssedge[2*eee+1];
                    eflag[j] = true;
                    break;
                }else if(curv==ssedge[2*eee+1]){
                    faceVind[i].push_back(ssedge[2*eee]);
                    curv = ssedge[2*eee];
                    eflag[j] = true;
                    break;
                }
            }
        }
        //for(auto a:faceVind[i])cout<<a<<' ';cout<<endl;


    }




    /*ctrf: contour face*/
    /*put conours into space faces constructed in last phase*/
    //int* ctrfvernum;			//# of V in each face
    //float** ctrfverpos;		//coord of V in each face
    //int** ctrfvertype;		//type of V in each face (completely in, or on the shared edge)
    //int** ctrfverval;			//?
    //int* ctrfedgenum;			//# of E in each face
    //int** ctrfedge;			//two ends of each E in each face, and left/right materials
    //int** ctrfedgetype;		//type of E in each face
    //int** ctrfedgeval;		//?
    //int** ctrfedgeancestor;	//?

    vector<int>onEVposNum(ssfacenum,0);
    vector< vector<int> >mapEV2V(ssfacenum);
    vector<float>onEVpos;
    vector<float>ctrVpos;
    vector<int>onEVfmap;
    vector<int>onEVvmap;
    vector<int>ctrfvernum_acc(ssfacenum,0);

    vector<int>ctrfverOutEnum_acc(ssfacenum,0);
    vector<int>ctrfverOutEnum(ssfacenum,0);
    vector<float>outEVpos;vector<int>outEVfaceind;


    vector<vector<int>>newInd(ssfacenum);
    vector<vector<int>>outEdges(ssfacenum);
    vector<vector<int>>outEdgesMap(ssfacenum);
    vector<int>outEdgesnum(ssfacenum);

    vector<vector<int>>onEEdges(ssfacenum);
    vector<vector<int>>onEdgesMap(ssfacenum);
    vector<int>onEdgesnum(ssfacenum);
    vector<int>onEEdge2face;

    vector<vector<int>>onEEdgesVInd(ssfacenum);
    vector<int>inverseOnEVind;

    int onEVacc = 0,outEVacc = 0;



    for(int i=0;i<ssfacenum;++i){
        newInd[i].resize(ctrfvernum[i]);
        for(int j=0;j<ctrfvernum[i];++j){
            for(int k=0;k<3;++k){
                ctrVpos.push_back(ctrfverpos[i][j*3+k]);
            }
            if(ctrfvertype[i][j]!=3){
                onEVposNum[i]++;
                newInd[i][j] = -1-onEVacc;onEVacc++;
                for(int k=0;k<3;++k){
                    onEVpos.push_back(ctrfverpos[i][j*3+k]);
                    cout<<ctrfverpos[i][j*3+k]<<' ';
                }
                onEEdge2face.push_back(i);
                mapEV2V[i].push_back(j);
                onEVfmap.push_back(i);onEVvmap.push_back(j);
                cout<<'('<<j<<')'<<" || ";
            }else{
                newInd[i][j] = outEVacc;outEVacc++;
                //newInd[i][j] = ctrfverOutEnum[i];
                ctrfverOutEnum[i]++;
                for(int k=0;k<3;++k){
                    outEVpos.push_back(ctrfverpos[i][j*3+k]);
                }
                outEVfaceind.push_back(i);


            }
        }
        cout<<i<<endl;
    }

    int totalonEv = 0,totaloutEv = 0,totalctrFv = 0;
    for(int i=0;i<ssfacenum;++i)totalonEv+=onEVposNum[i];
    for(int i=0;i<ssfacenum;++i)totalctrFv+=ctrfvernum[i];
    for(int i=1;i<ssfacenum;++i)ctrfvernum_acc[i]+=ctrfvernum_acc[i-1]+ctrfvernum[i-1];
    for(int i=1;i<ssfacenum;++i)ctrfverOutEnum_acc[i]+=ctrfverOutEnum_acc[i-1]+ctrfverOutEnum[i-1];
    for(int i=0;i<ssfacenum;++i)totaloutEv+=ctrfverOutEnum[i];

    for(int i=0;i<ssfacenum;++i){
        for(int j=0;j<ctrfedgenum[i];++j){
            int e1 = newInd[i][ctrfedge[i][4*j]],e2 = newInd[i][ctrfedge[i][4*j+1]];
            if(e1>=0 && e2>=0){
                outEdges[i].push_back(e1);
                outEdges[i].push_back(e2);
                outEdgesMap[i].push_back(j);
            }else{
                onEEdges[i].push_back(e1);
                onEEdges[i].push_back(e2);
                onEdgesMap[i].push_back(j);
            }
        }
    }

    for(int i=0;i<ssfacenum;++i)outEdgesnum[i] = outEdges[i].size()/2;
    for(int i=0;i<ssfacenum;++i)cout<<ctrfedgenum[i]-outEdgesnum[i]<<' ';cout<<endl;


    auto veclen = [](float *a,float *b){
        float r = 0;
        for(int i=0;i<3;++i)r+=(a[i]-b[i])*(a[i]-b[i]);
        return r;
    };

    vector<float>globalVpos;
    vector<float>mergeOnEVpos;


    vector<int>mergeOnEdgeVind(totalonEv,-1);
    vector<vector<int>>mergeV2faces;
    int onEmergeind = 0;
    for(int i=0; i<totalonEv;++i)if(mergeOnEdgeVind[i]==-1){
        mergeOnEdgeVind[i]=onEmergeind;++onEmergeind;
        onEEdgesVInd[onEVfmap[i]].push_back(mergeOnEdgeVind[i]);
        vector<int>tmpbatch;
        set<int>tmpset;
        tmpset.insert(ssface_planeindex[onEEdge2face[i]]);
        for(int j = i+1;j<totalonEv;++j){

            if(veclen(&(onEVpos[i*3]),&(onEVpos[j*3]))<1e-8)//1e-5
            {
                mergeOnEdgeVind[j] = mergeOnEdgeVind[i];
                tmpset.insert(ssface_planeindex[onEEdge2face[j]]);
                onEEdgesVInd[onEVfmap[j]].push_back(mergeOnEdgeVind[i]);
            }
        }
        for(int j = 0;j<3;++j)mergeOnEVpos.push_back(onEVpos[i*3+j]);
        for(auto iter = tmpset.begin();iter!=tmpset.end();++iter)tmpbatch.push_back(*iter);
        mergeV2faces.push_back(tmpbatch);
        assert(tmpbatch.size()<=2);
        for(auto a:tmpbatch)cout<<a<<' ';cout<<endl;
    }


    projectAllOnEdgeVertices(mergeOnEVpos,mergeV2faces,pparam);
    cout<<"mergeOnEdgeVind: ";for(auto a:mergeOnEdgeVind)cout<<a<<' ';cout<<endl;

    cout<<mergeOnEVpos.size()/3<<' '<<onEmergeind<<endl;

    vector<int>mergeOutEdgeVind(outEVacc,0);
    vector<vector<int> >newOnEEdges(ssfacenum);
    vector<vector<int> >newOutEEdges(ssfacenum);
    vector<vector<int> >newonEdgesMap(ssfacenum);

    for(int i=0;i<ssfacenum;++i){
        cout<<i<<' ';
        for(int j=0;j<onEEdges[i].size()/2;++j){
            int e1 = onEEdges[i][j*2],e2 = onEEdges[i][j*2+1],eOn,eOut;
            float *p_ev1,*p_ev2;
            if(e1<0) p_ev1 = &(mergeOnEVpos[3*(mergeOnEdgeVind[-e1-1])]);else p_ev1 = &(outEVpos[e1*3]);
            if(e2<0) p_ev2 = &(mergeOnEVpos[3*(mergeOnEdgeVind[-e2-1])]);else p_ev2 = &(outEVpos[e2*3]);
            cout<<veclen(p_ev1,p_ev2)<<' ';
            cout<<e1<<' '<<e2<<' ';
            if(e1<0)eOn=e1;else eOut=e1;
            if(e2<0)eOn=e2;else eOut=e2;
            if(veclen(p_ev1,p_ev2)<0){mergeOutEdgeVind[eOut] = eOn;cout<<"ZXXXXXXXXXXX";}
            else {
                newOnEEdges[i].push_back(e1);newOnEEdges[i].push_back(e2);
                newonEdgesMap[i].push_back(onEdgesMap[i][j]);
            }
        }
        cout<<endl;
    }

    int newMergeOutEdgeIndacc = 0;
    int maxtmp = 0;
    vector<int>globalOutEVfaceind;
    for(int i=0;i<outEVacc;++i){
        if(mergeOutEdgeVind[i]==0){
            mergeOutEdgeVind[i]=newMergeOutEdgeIndacc++;
            for(int j=0;j<3;++j)globalVpos.push_back(outEVpos[i*3+j]);
            globalOutEVfaceind.push_back(outEVfaceind[i]);
        }
    }
    for(int i=0;i<outEVacc;++i)if(mergeOutEdgeVind[i]<0){
        mergeOutEdgeVind[i] = mergeOnEdgeVind[-mergeOutEdgeVind[i]-1]+newMergeOutEdgeIndacc;
        cout<<mergeOutEdgeVind[i]<<' ';
        maxtmp = max(maxtmp,mergeOutEdgeVind[i]);
    }
    if(isOnEdge)for(auto a:mergeOnEVpos)globalVpos.push_back(a);
    //for(int i=0;i<24;++i)globalVpos.push_back(mergeOnEVpos[i]);
    cout<<"newMergeOutEdgeIndacc: "<<newMergeOutEdgeIndacc<<' '<<maxtmp<<endl;


    for(int i=0;i<ssfacenum;++i){
        for(int j=0;j<outEdges[i].size()/2;++j){
            int aa = mergeOutEdgeVind[outEdges[i][j*2]];
            int bb = mergeOutEdgeVind[outEdges[i][j*2+1]];
            if(!(aa<0 && bb<0)){
                newOutEEdges[i].push_back(mergeOutEdgeVind[outEdges[i][j*2]]);
                newOutEEdges[i].push_back(mergeOutEdgeVind[outEdges[i][j*2+1]]);
            }
        }
        //for(auto a:outEdges[i])newOutEEdges[i].push_back(mergeOutEdgeVind[a]);
        for(int j=0;j<newOnEEdges[i].size();++j)
            if(newOnEEdges[i][j]<0){newOnEEdges[i][j] = mergeOnEdgeVind[-newOnEEdges[i][j]-1]+newMergeOutEdgeIndacc;}
            else newOnEEdges[i][j] = mergeOutEdgeVind[newOnEEdges[i][j]];
    }


    auto globalEdges = newOutEEdges;
    if(isOnEdge)for(int i=0;i<ssfacenum;++i)for(auto a:newOnEEdges[i])globalEdges[i].push_back(a);
    // auto globalEdges = newOnEEdges;
    auto globalEdgesMap = outEdgesMap;
    if(isOnEdge)for(int i=0;i<ssfacenum;++i)for(auto a:newonEdgesMap[i])globalEdgesMap[i].push_back(a);


    vector<vector<int> >globalEdgesMat(ssfacenum);
    for(int i=0;i<ssfacenum;++i){
        cout<<i<<": ";
        for(auto a:globalEdgesMap[i]){
            globalEdgesMat[i].push_back(ctrfedge[i][4*a+2]);
            globalEdgesMat[i].push_back(ctrfedge[i][4*a+3]);
        }
        for(auto a:globalEdgesMat[i])cout<<a<<' ';
        cout<<endl;
    }

    //for(int i=0;i<ssfacenum;++i){for(auto a:globalEdges[i])cout<<a<<' ';cout<<endl;}

    vector<int>faceedgeNum(ssfacenum);
    for(int i=0;i<ssfacenum;++i)faceedgeNum[i] = (globalEdges[i].size())/2;
    int totalEdgesnum = 0;
    for(auto a:faceedgeNum)totalEdgesnum+=a;
    int totalver = globalVpos.size()/3;


    auto  originfaceVind = faceVind;
    if(isOnEdge){

        auto inbetweenVertices = [](float *v1,float *v2,float *v3,float &v31len){

            float vecDir[3],vec31[3];
            MyUtility::minusVec(v2,v1,vecDir);
            MyUtility::minusVec(v3,v1,vec31);
            float aa = MyUtility::normVec(vecDir),bb = MyUtility::normVec(vec31);
            MyUtility::normalize(vecDir);
            MyUtility::normalize(vec31);

            if(abs(1- MyUtility::dot(vecDir,vec31))<1e-3)
                if(bb<aa){
                    v31len = bb;
                    return true;
                }



            return false;





        };

        vector<float>totalglobalVpos;
        for(auto &a:onEEdgesVInd)for(auto &b:a)b+=newMergeOutEdgeIndacc+ssvernum;
        for (int i = 0; i < ssvernum*3; ++i)totalglobalVpos.push_back(ssver[i]);
        for(auto v:globalVpos)totalglobalVpos.push_back(v);
        float* p_totalglobalVpos = totalglobalVpos.data();
        auto onestep = [p_totalglobalVpos,&onEEdgesVInd,&faceVind,inbetweenVertices](int i, int v1,int v2, int &j,vector<bool>&taken){
            vector<pair<float,int>>poslenInd;
            float tmplen;
            for(int k=0;k<onEEdgesVInd[i].size();++k)
                if(!taken[k] && inbetweenVertices(p_totalglobalVpos+v1*3,p_totalglobalVpos+v2*3,p_totalglobalVpos+onEEdgesVInd[i][k]*3,tmplen)){
                    taken[k] = true;
                    poslenInd.push_back(make_pair(tmplen,onEEdgesVInd[i][k]));
                    //cout<<i<<' '<<j<<' '<<onEEdgesVInd[i][k]<<endl;
                }
            if(poslenInd.size()!=0)sort(poslenInd.begin(),poslenInd.end());
            //for(auto a:poslenInd)cout<<a.first<<' ';cout<<endl;
            if(poslenInd.size()!=0)
                for(auto &a:poslenInd){faceVind[i].insert(faceVind[i].begin()+j+1,a.second);++j;}
            //if(j==9){break;}
        };
        for(int i = 0;i<ssfacenum;++i){
            //cout<<"onEEdgesVInd: "<<i<<": ";for(auto a:onEEdgesVInd[i])cout<<a<<' ';cout<<endl;
            //continue;
            int bb = faceVind[i].size();
            vector<bool>taken(onEEdgesVInd[i].size(),false);
            int j;
            for(j=0;j<faceVind[i].size()-1;++j){


                int v1 = faceVind[i][j],v2 = faceVind[i][j+1];

                onestep(i,v1,v2,j,taken);

            }

            onestep(i,faceVind[i][j],faceVind[i][0],j,taken);
            cout<<"faceVind: "<<i<<": ";for(auto a:faceVind[i])cout<<a<<' ';cout<<endl;
            cout<<"onEEdgesVInd: "<<i<<": ";for(auto a:onEEdgesVInd[i])cout<<a<<' ';cout<<endl;

            assert(bb+onEEdgesVInd[i].size()==faceVind[i].size());
        }

    }

    vector<int>globalEdges_cluster;
    for(auto &a:faceVind)for(int i=0;i<a.size()-1;++i){
        globalEdges_cluster.push_back(a[i]);
        globalEdges_cluster.push_back(a[i+1]);
    }
    vector<float>debugglobalVpos;
    for (int i = 0; i < ssvernum*3; ++i)debugglobalVpos.push_back(ssver[i]);
    for(auto a:globalVpos)debugglobalVpos.push_back(a);

    for(int i=0;i<ssfacenum;++i){for(auto a:globalEdges[i])globalEdges_cluster.push_back(a+ssvernum);}
    //outputCurve(string("../ConvertFolder/testmerge.curve"),debugglobalVpos,globalEdges_cluster);


    cout.rdbuf(orig_buf);
    //for(int i=0;i<globalEdges_cluster.size()/2;++i)
    // cout<<i<<' '<<(veclen(&(globalVpos[3*globalEdges_cluster[i*2]]),&(globalVpos[3*globalEdges_cluster[i*2+1]])))<<endl;
    //cout<<"XXXXXXXXXXXXX: "<<i<<endl;

    //  return;


    tetgenio in;
    tetgenio::facet* f;
    tetgenio::polygon* p;

    // PLC: pointlist.
    cout<<"PLC: pointlist. "<<ssvernum<<' '<<totaloutEv<<' '<<totalver<<endl;
    in.numberofpoints = ssvernum + totalver;
    //int tmpnum = 30;
    //in.numberofpoints = ssvernum + tmpnum;
    //in.numberofpoints = _frameVs.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    in.pointmarkerlist = new int[in.numberofpoints];
    //    for (int i = 0; i < in.numberofpoints; ++i) {
    //        in.pointmarkerlist[i] = -1;
    //    }
    for (int i = 0; i < ssvernum*3; ++i) {
        in.pointlist[i] = ssver[i];
    }
    for (int i = 0; i < totalver*3; ++i) {
        in.pointlist[ssvernum*3+i] = globalVpos[i];
    }

    for (int i = 0; i < ssvernum; ++i) {
        in.pointmarkerlist[i] = -1;
    }
    for(int i=0;i<newMergeOutEdgeIndacc;++i){
        in.pointmarkerlist[i+ssvernum] = -100-globalOutEVfaceind[i];
    }
    for (int i = 0; i < onEmergeind; ++i) {
        in.pointmarkerlist[i+ssvernum+newMergeOutEdgeIndacc] = -10000;
    }





    // PLC: facetlist.
    cout<<"PLC: facetlist."<<endl;
    in.numberoffacets = ssfacenum;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (int i = 0; i < in.numberoffacets; ++i) {
        f = &in.facetlist[i];
        //f->numberofpolygons = 1;
        f->numberofpolygons = 1+faceedgeNum[i];
        //f->numberofpolygons = 1+ctrfedgenum[i];

        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = faceVind[i].size();
        p->vertexlist = new int[p->numberofvertices];
        for (int j = 0; j < faceVind[i].size(); ++j) {
            p->vertexlist[j] = faceVind[i][j];
        }
        in.facetmarkerlist[i] = i+100;
        if(1)if(faceedgeNum[i]!=0){
            //int offset = ssvernum+ctrfverOutEnum_acc[i];
            int offset = ssvernum;
            for(int j=0;j<faceedgeNum[i];++j){
                p = &f->polygonlist[j+1];
                p->numberofvertices = 2;p->vertexlist = new int[p->numberofvertices];
                auto p_ev = &(globalEdges[i][j*2]);
                for (int k = 0; k < 2; ++k)p->vertexlist[k] = p_ev[k]+ offset;
            }
        }

    }

    in.numberofedges = totalEdgesnum;
    in.edgelist = new int[in.numberofedges*2];
    in.edgemarkerlist = new int[in.numberofedges];
    int tmpi = 0,tmpj=0;
    for (int i = 0; i < in.numberoffacets; ++i) {
        if(1)if(faceedgeNum[i]!=0){
            //int offset = ssvernum+ctrfverOutEnum_acc[i];
            int offset = ssvernum;
            for(int j=0;j<faceedgeNum[i];++j){
                auto p_ev = &(globalEdges[i][j*2]);
                in.edgemarkerlist[tmpj++] = 100+i;
                for (int k = 0; k < 2; ++k)in.edgelist[tmpi++] = p_ev[k]+ offset;
            }
        }
    }



    // PLC: region.
    cout<<"PLC: region."<<endl;
    in.numberofregions = ssspacenum;
    in.regionlist = new REAL[in.numberofregions * 5];
    for (int i = 0; i < in.numberofregions; ++i) {
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

    //float tetVolLimit = 50.0f;
    float tetVolLimit = 50.0f;

    //string switchesStr = "Qepq1.714Aa" + to_string(tetVolLimit);
    string switchesStr = "Qepq1.714Aa";//1.4.3
    //string switchesStr = "zQpq1.714Aa"+ to_string(tetVolLimit);//1.5

    std::vector<char> tmp(switchesStr.begin(), switchesStr.end());
    tmp.push_back('\0');
    char* switches = &tmp[0];

    //cout<<"out.numberofpoints "<<out.numberofpoints<<endl;

    //Utility::saveTetgenio((t_savedir+"myTetIn_"+".txt").c_str(),in);

    //in.save_nodes(const_cast<char *>((t_savedir+"test_in_").c_str()));
    //in.save_poly(const_cast<char *>((t_savedir+"test_in_").c_str()));


    tetrahedralize(switches, &in, &out);

    //char *tetout = (char*)(t_savedir.data());
    cout<<"out.numberofpoints "<<out.numberofpoints<<endl;
    //out.save_nodes(tetout);
    //out.save_elements(tetout);
    //out.save_faces(tetout);
    //out.save_poly(tetout);
    //out.save_edges(tetout);



    // _________________________Connection_________________________

    //vector<float>_frameVstack;
    //vector<vector<float>> _frameVs;
    //vector<vector<int>> _frameFs;
    //vector<vector<int>> _frameCells;
    //vector<vector<float>> _frameF2PlaneParas;

    cout<<"_________________________Connection_________________________"<<endl;

    _nCell = ssspacenum;
    _nCrossSection = planenum;
    for(int i=0;i<ssvernum*3;++i)_frameVstack.push_back(ssver[i]);
    //_frameFs = faceVind;
    _frameFs = originfaceVind;
    _frameCells.resize(ssspacenum);
    for(int i=0;i<ssspacenum;++i){
        for(int j=0;j<ssspacefacenum[i];++j){
            _frameCells[i].push_back( ssspace[i][j]);
        }
    }
    _Cs2PlaneParas.resize(_nCrossSection);
    _frameF2PlaneParas.resize(ssfacenum);
    for(int i=0;i<_nCrossSection;++i)for(int j=0;j<4;++j)_Cs2PlaneParas[i].push_back(pparam[i*4+j]);
    for(int i=0;i<ssfacenum;++i)_framef2Cs.push_back(ssface_planeindex[i]);
    for(int i=0;i<ssfacenum;++i)if(_framef2Cs[i]>=0)_frameF2PlaneParas[i] = _Cs2PlaneParas[_framef2Cs[i]];else _frameF2PlaneParas[i] = vector<float>({0,0,0,0});


    _globalVpos = globalVpos;
    _globalEdges = globalEdges;
    _globalEdgesMat = globalEdgesMat;

    for(auto &a:_globalEdges)for(auto &b:a)b+=ssvernum;

    cout<<"_nCell: "<<_nCell<<endl;
    vector<vector<int>>Cs2Frame(_nCrossSection);
    for(int i=0;i<_framef2Cs.size();++i)if(_framef2Cs[i]>=0)Cs2Frame[_framef2Cs[i]].push_back(i);

    vector<vector<int>>Cs2Ctr(_nCrossSection),Cs2CtrM(_nCrossSection);
    for(int i=0;i<_nCrossSection;++i)for(auto a:Cs2Frame[i]){
        for(auto b:_globalEdges[a])Cs2Ctr[i].push_back(b-ssvernum);
        for(auto b:_globalEdgesMat[a])Cs2CtrM[i].push_back(b);
    }

    //bool writeCtrGraphFile(string filename,const vector<float>&vertices,const vector<vector<int>>&edge2vertices,const vector<vector<int>>&edgeMat, const vector<vector<float>>&planepara);
    //writeCtrGraphFile("../ConvertFolder/CtrGraph",_globalVpos,Cs2Ctr,Cs2CtrM,_Cs2PlaneParas);



    return;

    /***************************************************************/
    /***************************************************************/




}
void outputCurve(string filename,vector<float>&vv,vector<int>&ee){

    ofstream fout(filename.data());
    if(fout.fail())cout<<"error"<<endl;
    fout<<vv.size()/3<<' '<<ee.size()/2<<endl;
    for(int i=0;i<vv.size()/3;++i){
        for(int j=0;j<3;++j)fout<<vv[i*3+j]<<' ';
        fout<<endl;
    }

    for(int i=0;i<ee.size()/2;++i){
        for(int j=0;j<2;++j)fout<<ee[i*2+j]<<' ';
        fout<<endl;
    }
}
void projectAllOnEdgeVertices(vector<float>&vpos,vector<vector<int>>&v2ps,float *pparam){


    auto projectErr = [pparam](float *pp,int pind){
        auto ppara = pparam+pind*4;
        float d = -ppara[3]/sqrt(ppara[0]*ppara[0]+ppara[1]*ppara[1]+ppara[2]*ppara[2]);
        return abs(pp[0]*ppara[0]+pp[1]*ppara[1]+pp[2]*ppara[2]+d);
    };
    auto project2p = [pparam](float *pp,int pind){
        auto ppara = pparam+pind*4;
        float d = -ppara[3]/sqrt(ppara[0]*ppara[0]+ppara[1]*ppara[1]+ppara[2]*ppara[2]);
        float err = (pp[0]*ppara[0]+pp[1]*ppara[1]+pp[2]*ppara[2]+d);
        for(int i=0;i<3;++i)pp[i]-=ppara[i]*err;
    };

    int npoints = v2ps.size();
    for(int i=0;i<npoints;++i){
        auto &pl = v2ps[i];
        float *pp = vpos.data()+i*3;
        float error = 100;
        while(error>1e-2){
            for(auto a:pl)project2p(pp,a);
            error = 0;
            for(auto a:pl)error+=projectErr(pp,a);
        }
    }


}

void CallTetgenForSingleSurface(vector<double>&vpos,  vector<uint>&fs,  vector<vector<float>>&out_vpos, vector<vector<int>>&out_tet){

    tetgenio in,out;
    tetgenio::facet* f;
    tetgenio::polygon* p;

    // PLC: pointlist.
    cout<<"PLC: pointlist. "<<endl;
    in.numberofpoints = vpos.size()/3;
    in.pointlist = new REAL[in.numberofpoints * 3];
    in.pointmarkerlist = new int[in.numberofpoints];

    for (int i = 0; i < in.numberofpoints*3; ++i) {
        in.pointlist[i] = vpos[i];
    }

    for (int i = 0; i < in.numberofpoints; ++i) {
        in.pointmarkerlist[i] = -1;
    }

    cout<<"PLC: facetlist."<<endl;
    in.numberoffacets = fs.size()/3;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (int i = 0; i < in.numberoffacets; ++i) {
        f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];
        for (int j = 0; j < 3; ++j) {
            p->vertexlist[j] = fs[i*3+j];
        }
        in.facetmarkerlist[i] = i;
    }


    // PLC: region.
    cout<<"PLC: region."<<endl;

    char tetin[] = "../tetgen1.4.3/c++tetin";
    // in.save_nodes(tetin);
    // in.save_poly(tetin);

    //float tetVolLimit = 50.0f;
    float tetVolLimit = 10000.0f;

    string switchesStr = "YQepq2.714Aa" + to_string(tetVolLimit);
    //string switchesStr = "Qepq1.714Aa";

    std::vector<char> tmp(switchesStr.begin(), switchesStr.end());
    tmp.push_back('\0');
    char* switches = &tmp[0];

    //cout<<"out.numberofpoints "<<out.numberofpoints<<endl;

    //Utility::saveTetgenio((t_savedir+"myTetIn_"+".txt").c_str(),in);

    //in.save_nodes(const_cast<char *>((t_savedir+"test_in_").c_str()));
    //in.save_poly(const_cast<char *>((t_savedir+"test_in_").c_str()));


    tetrahedralize(switches, &in, &out);

    //cout<<"out.numberofpoints "<<out.numberofpoints<<endl;

    cout<<"out.numberofpoints "<<out.numberofpoints<<endl;
    char tetout[] = {"../TetgenV/outtestgg"};
    out.save_nodes(tetout);
    out.save_elements(tetout);
    out.save_faces(tetout);
    out.save_poly(tetout);
    out.save_edges(tetout);

//    void MMLevelSet::LoadCell(const vector<vector<float>>& cellVs,
//                              const vector<int>&cellVMarkers,
//                              const vector<vector<int>>& cellTs,
//                              const vector<int>& cellInitLabel, int nMat)

    cout<<"Extract the output tets from tetgen data structure to our own data structure"<<endl;
    out_vpos.resize(out.numberofpoints, vector<float>(3));
//    _tetVMarkers.resize(out.numberofpoints);
//    _tetVMarkers_overlap.resize(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; ++i) {
        for (int j = 0; j < 3; ++j) {
            out_vpos[i][j] = out.pointlist[i * 3 + j];
        }
    }
//    for (int i = 0; i < out.numberofpoints; ++i) {
//        _tetVMarkers[i] = out.pointmarkerlist[i];
//    }
//    _frameF2tetFs.resize(_frameFs.size());
//    cout<<"_frameFs.size(): "<<_frameFs.size()<<endl;
    cout<<"step 1"<<endl;

    cout<<"step 2"<<endl;
    assert(out.numberoftetrahedronattributes == 1);
    vector<int> oneTetT(4);
    for (int i = 0; i < out.numberoftetrahedra; ++i) {
        for (int j = 0; j < 4; ++j) {
            oneTetT[j] = out.tetrahedronlist[i * 4 + j];

        }
        sort(oneTetT.begin(), oneTetT.end());
        out_tet.push_back(oneTetT);
    }






}
