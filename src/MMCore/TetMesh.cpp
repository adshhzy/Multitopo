#include "Mesh.h"
#include "UnionFind.h"
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <functional>
#include "gurobi_c++.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
using Eigen::Vector3f;

typedef unordered_map<vector<int>, int, container_hash<vector<int>>> efHash;
void TetMesh::GenerateIndicatorFunction(const int nMat, const vector<int>& labeling,vector<int>&mVs,
                                        vector<vector<float>>& v2f) {
    if (nMat < 2) return;
    int nV = Vs.size();
    int nBdV = mVs.size();
    vector<bool>b_isbV(nV,false);
    for(auto a:mVs)b_isbV[a] = true;
    int nXV = nV - nBdV;
    // By default, edge weight is set to be 1. Optionally, user could change the
    // weight with the density difference if a volume data is available.
    float weightIJ = 1.0f;
    // Divide the vertices into marked Vs (mVs) and unmarked Vs (uVs).
    vector<int> uVs;
    for (int vi = 0; vi < nV; ++vi) {
        if (!b_isbV[vi]) {
            uVs.push_back(vi);
        }
    }
    // Re-arrange vertices so that marked Vs are all in the begining.
    // Keep record of the mapping.
    vector<int> mapNewV2OrgV;
    mapNewV2OrgV.reserve(nV);
    mapNewV2OrgV.insert(mapNewV2OrgV.end(), mVs.begin(), mVs.end());
    mapNewV2OrgV.insert(mapNewV2OrgV.end(), uVs.begin(), uVs.end());

    vector<int> mapOrgV2NewV(nV);
    for (int vi = 0; vi < nV; ++vi) {
        mapOrgV2NewV[mapNewV2OrgV[vi]] = vi;
    }
    // Construct matrix L.
    // For more details, please refer to "Random Walks for Image Segmentation".
    // We only need partial of L matrix in the paper, i.e. Bt and Lu.
    // The matL here is nXV x nV in size.
    Eigen::SparseMatrix<float> matL(nXV, nV);
    vector<Eigen::Triplet<float>> tripletList;
    for (int nvi = nBdV; nvi < nV; ++nvi) {
        int vi = mapNewV2OrgV[nvi];
        float di = 0.0f;
        for (auto vj : Vs[vi].nbV) {
            float wij = weightIJ;
            int nvj = mapOrgV2NewV[vj];
            di += wij;
            tripletList.push_back(Eigen::Triplet<float>(nvi - nBdV, nvj, -wij));
        }
        tripletList.push_back(Eigen::Triplet<float>(nvi - nBdV, nvi, di));
    }
    matL.setFromTriplets(tripletList.begin(), tripletList.end());
    // Extract Bt and Lu matrix.
    Eigen::SparseMatrix<float> Bt = matL.leftCols(nBdV);
    Eigen::SparseMatrix<float> Lu = matL.rightCols(nXV);

    // Constrct the M matrix.
    // Since the probabilities at any node will sum to unity, only mMat − 1 sparse
    // linear systems must be solved. Here I skipped the background (mat 0).
    // Eigen::SparseMatrix<float> matM(nBdV, nMat - 1);
    // tripletList.clear();
    // for (int nvi = 0; nvi < nBdV; ++nvi) {
    //   int vi = mapNewV2OrgV[nvi];
    //   int mat = labeling[vi];
    //   if (mat != 0) {
    //     tripletList.push_back(Eigen::Triplet<float>(nvi, mat - 1, 1.0));
    //   }
    // }
    // matM.setFromTriplets(tripletList.begin(), tripletList.end());
    Eigen::MatrixXf matM = Eigen::MatrixXf::Zero(nBdV, nMat - 1);
    for (int nvi = 0; nvi < nBdV; ++nvi) {
        int vi = mapNewV2OrgV[nvi];
        int mat = labeling[vi];
        if (mat != 0) {
            matM(nvi, mat - 1) = 1.0f;
        }
    }

    // Solve the nMat-1 linear system.
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;
    solver.compute(Lu);
    if (solver.info() != Eigen::Success) {
        cout << "Error! Eigen LDLT factorization failed !!!! " << endl;
    }
    // auto matTest = -Bt * matM;
    // auto matX = solver.solve(matTest);
    Eigen::MatrixXf matX = solver.solve(-Bt * matM);
    if (solver.info() != Eigen::Success) {
        cout << "Error! Eigen solving failed !!!! " << endl;
    }

    // Save results to v2f.
    // cout << "nV = " << nV << endl;
    // cout << "nXV = " << nXV << endl;
    // cout << "nBdV = " << nBdV << endl;
    // cout << "nMat = " << nMat << endl;
    // cout << "nRow = " << matX.rows() << endl;
    // cout << "nCol = " << matX.cols() << endl;
    v2f.resize(nV, vector<float>(nMat, 0.0f));
    for (int vi = 0; vi < nV; ++vi) {
        if (b_isbV[vi]) {
            // if (!(vi >=0 && vi<nV && labeling[vi]>=0 && labeling[vi]<nMat)) {
            // cout << "vi = " << vi << endl;
            // cout << "label = " << labeling[vi] << endl;
            // }
            // cout << v2f[vi][labeling[vi]] << endl;
            v2f[vi][labeling[vi]] = 1.0f;
            continue;
        }
        int nvi = mapOrgV2NewV[vi] - nBdV;
        // cout << nvi << endl;
        float acc = 0.0f;
        for (int mi = 1; mi < nMat; ++mi) {
            v2f[vi][mi] = matX(nvi, mi - 1);
            acc += v2f[vi][mi];
        }
        v2f[vi][0] = 1.0f - acc;
    }
}

void TetMesh::GenerateIndicatorFunction(int nMat, const vector<int>& labeling,
                                        vector<vector<float>>& v2f) {
    if (nMat < 2) return;
    int nV = Vs.size();
    int nBdV = _bdVs.size();
    int nXV = nV - nBdV;
    // By default, edge weight is set to be 1. Optionally, user could change the
    // weight with the density difference if a volume data is available.
    float weightIJ = 1.0f;
    // Divide the vertices into marked Vs (mVs) and unmarked Vs (uVs).
    const auto& mVs = _bdVs;
    vector<int> uVs;
    for (int vi = 0; vi < nV; ++vi) {
        if (!Vs[vi].isBD) {
            uVs.push_back(vi);
        }
    }
    // Re-arrange vertices so that marked Vs are all in the begining.
    // Keep record of the mapping.
    vector<int> mapNewV2OrgV;
    mapNewV2OrgV.reserve(nV);
    mapNewV2OrgV.insert(mapNewV2OrgV.end(), mVs.begin(), mVs.end());
    mapNewV2OrgV.insert(mapNewV2OrgV.end(), uVs.begin(), uVs.end());

    vector<int> mapOrgV2NewV(nV);
    for (int vi = 0; vi < nV; ++vi) {
        mapOrgV2NewV[mapNewV2OrgV[vi]] = vi;
    }
    // Construct matrix L.
    // For more details, please refer to "Random Walks for Image Segmentation".
    // We only need partial of L matrix in the paper, i.e. Bt and Lu.
    // The matL here is nXV x nV in size.
    Eigen::SparseMatrix<float> matL(nXV, nV);
    vector<Eigen::Triplet<float>> tripletList;
    for (int nvi = nBdV; nvi < nV; ++nvi) {
        int vi = mapNewV2OrgV[nvi];
        float di = 0.0f;
        for (auto vj : Vs[vi].nbV) {
            float wij = weightIJ;
            int nvj = mapOrgV2NewV[vj];
            di += wij;
            tripletList.push_back(Eigen::Triplet<float>(nvi - nBdV, nvj, -wij));
        }
        tripletList.push_back(Eigen::Triplet<float>(nvi - nBdV, nvi, di));
    }
    matL.setFromTriplets(tripletList.begin(), tripletList.end());
    // Extract Bt and Lu matrix.
    Eigen::SparseMatrix<float> Bt = matL.leftCols(nBdV);
    Eigen::SparseMatrix<float> Lu = matL.rightCols(nXV);

    // Constrct the M matrix.
    // Since the probabilities at any node will sum to unity, only mMat − 1 sparse
    // linear systems must be solved. Here I skipped the background (mat 0).
    // Eigen::SparseMatrix<float> matM(nBdV, nMat - 1);
    // tripletList.clear();
    // for (int nvi = 0; nvi < nBdV; ++nvi) {
    //   int vi = mapNewV2OrgV[nvi];
    //   int mat = labeling[vi];
    //   if (mat != 0) {
    //     tripletList.push_back(Eigen::Triplet<float>(nvi, mat - 1, 1.0));
    //   }
    // }
    // matM.setFromTriplets(tripletList.begin(), tripletList.end());
    Eigen::MatrixXf matM = Eigen::MatrixXf::Zero(nBdV, nMat - 1);
    for (int nvi = 0; nvi < nBdV; ++nvi) {
        int vi = mapNewV2OrgV[nvi];
        int mat = labeling[vi];
        if (mat != 0) {
            matM(nvi, mat - 1) = 1.0f;
        }
    }

    // Solve the nMat-1 linear system.
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;
    solver.compute(Lu);
    if (solver.info() != Eigen::Success) {
        cout << "Error! Eigen LDLT factorization failed !!!! " << endl;
    }
    // auto matTest = -Bt * matM;
    // auto matX = solver.solve(matTest);
    Eigen::MatrixXf matX = solver.solve(-Bt * matM);
    if (solver.info() != Eigen::Success) {
        cout << "Error! Eigen solving failed !!!! " << endl;
    }

    // Save results to v2f.
    // cout << "nV = " << nV << endl;
    // cout << "nXV = " << nXV << endl;
    // cout << "nBdV = " << nBdV << endl;
    // cout << "nMat = " << nMat << endl;
    // cout << "nRow = " << matX.rows() << endl;
    // cout << "nCol = " << matX.cols() << endl;
    v2f.resize(nV, vector<float>(nMat, 0.0f));
    for (int vi = 0; vi < nV; ++vi) {
        if (Vs[vi].isBD) {
            // if (!(vi >=0 && vi<nV && labeling[vi]>=0 && labeling[vi]<nMat)) {
            // cout << "vi = " << vi << endl;
            // cout << "label = " << labeling[vi] << endl;
            // }
            // cout << v2f[vi][labeling[vi]] << endl;
            v2f[vi][labeling[vi]] = 1.0f;
            continue;
        }
        int nvi = mapOrgV2NewV[vi] - nBdV;
        // cout << nvi << endl;
        float acc = 0.0f;
        for (int mi = 1; mi < nMat; ++mi) {
            v2f[vi][mi] = matX(nvi, mi - 1);
            acc += v2f[vi][mi];
        }
        v2f[vi][0] = 1.0f - acc;
    }
}
void TetMesh::PreProcessBD(const vector<int>& label, int nLabel) {
    int nV = Vs.size();
    // Union-find data structure for bdComps.
    vector<int> parent(nV);
    MingUtility::fillVectorWithRange(0, 1, parent);

    // Find _bdComps through union-find _bdEs.
    vector<int> candFtEs;
    for (auto ei : _bdEs) {
        int v0 = Es[ei].v[0];
        int v1 = Es[ei].v[1];
        // Edge ei is completely inside of a bdComp, do union-fine.
        if (label[v0] == label[v1]) {
            MingUtility::unionfind_FindnUnion(v0, v1, parent);
        }
        // Edge ei lies on the frontier of two bdComps, save it for later.
        else {
            candFtEs.push_back(ei);
        }
    }
    unordered_map<int, int> bdCompPrV2bdCompID;
    for (auto vi : _bdVs) {
        int ri = MingUtility::unionfind_find(vi, parent);
        // New bdComp is found.
        if (bdCompPrV2bdCompID.find(ri) == bdCompPrV2bdCompID.end()) {
            bdCompPrV2bdCompID[ri] = _bdComps.size();
            _bdComps.push_back({vi});
        } else {
            _bdComps[bdCompPrV2bdCompID[ri]].push_back(vi);
        }
    }

    // Get _bdComp2label and _label2bdComps.
    int nBdComp = _bdComps.size();
    _bdComp2label.resize(nBdComp);
    _label2bdComps.resize(nLabel);
    for (int i = 0; i < _bdComps.size(); ++i) {
        int li = label[_bdComps[i][0]];
        _bdComp2label[i] = li;
        _label2bdComps[li].push_back(i);
    }

    // Find _bdComp2ftVs/Es/Fs.
    // vector<unordered_set<int>> bdComp2ftVsSets(nBdComp);
    // vector<unordered_set<int>> bdComp2ftFsSets(nBdComp);
    // _bdComp2ftVs.resize(nBdComp);
    _bdComp2ftEs.resize(nBdComp);
    // _bdComp2ftFs.resize(nBdComp);
    for (auto ei : candFtEs) {
        int v0 = Es[ei].v[0];
        int v1 = Es[ei].v[1];
        vector<int> nbBdFs;
        for (auto fi : Es[ei].nbF) {
            if (Fs[fi].isBD) {
                nbBdFs.push_back(fi);
            }
        }
        int c0 = bdCompPrV2bdCompID[MingUtility::unionfind_find(v0, parent)];
        int c1 = bdCompPrV2bdCompID[MingUtility::unionfind_find(v1, parent)];
        _bdComp2ftEs[c0].push_back(ei);
        _bdComp2ftEs[c1].push_back(ei);
        // bdComp2ftVsSets[c0].insert(v0);
        // bdComp2ftVsSets[c1].insert(v1);
        // bdComp2ftFsSets[c0].insert(nbBdFs[0]);
        // bdComp2ftFsSets[c0].insert(nbBdFs[1]);
        // bdComp2ftFsSets[c1].insert(nbBdFs[0]);
        // bdComp2ftFsSets[c1].insert(nbBdFs[1]);
    }
    // for (int i = 0; i < nBdComp; ++i) {
    //   _bdComp2ftVs[i].resize(bdComp2ftVsSets[i].size());
    //   _bdComp2ftFs[i].resize(bdComp2ftFsSets[i].size());
    //   copy(bdComp2ftVsSets[i].begin(), bdComp2ftVsSets[i].end(),
    //        _bdComp2ftVs[i].begin());
    //   copy(bdComp2ftFsSets[i].begin(), bdComp2ftFsSets[i].end(),
    //        _bdComp2ftFs[i].begin());
    // }

    // Find _bdLoops, _bdComp2Loops and _bdLoop2Comp.
    _bdComp2Loops.resize(nBdComp);
    UnionFind unionfind;
    vector<bool>EnF(Fs.size(),false);
    for (int ci = 0; ci < nBdComp; ++ci) {
        const auto& ftEs = _bdComp2ftEs[ci];
        unordered_set<int> ftEsSet(ftEs.begin(), ftEs.end());
        unionfind.SetElements(ftEs);
        for (auto ei : ftEs) {
            for(auto nf: Es[ei].nbF)if(Fs[nf].isBD)EnF[nf] = true;
            for (auto nbEi : Es[ei].nbE) {
                if (ftEsSet.find(nbEi) != ftEsSet.end()) {
                    //bool isn = false;
                    for(auto nf:Es[nbEi].nbF)if(Fs[nf].isBD)if(EnF[nf]){unionfind.Union(ei, nbEi);break;}
                    //unionfind.Union(ei, nbEi);
                }
            }
            for(auto nf: Es[ei].nbF)EnF[nf] = false;
        }
        vector<vector<int>> loops;
        unionfind.ExtractComponents(loops);
        for (const auto& loop : loops) {
            _bdComp2Loops[ci].push_back(_bdLoops.size());
            _bdLoops.push_back(loop);
            _bdLoop2Comp.push_back(ci);
        }
    }

    // Find _label2bdLoops.
    _label2bdLoops.resize(_label2bdComps.size());
    for (int li = 0; li < _label2bdComps.size(); ++li) {
        for (auto compID : _label2bdComps[li]) {
            for (auto loopID : _bdComp2Loops[compID]) {
                _label2bdLoops[li].push_back(loopID);
            }
        }
    }

    _bdLoopsV.resize(_bdLoops.size());
    for(int i=0;i<_bdLoops.size();++i){
        _bdLoopsV[i].clear();
        for(auto a:_bdLoops[i]){
            _bdLoopsV[i].push_back(Es[a].v[0]);
            _bdLoopsV[i].push_back(Es[a].v[1]);
        }
    }


    // Find boundary curves.
    GetActiveBDEsNFs(label);
    SmoothBoundaryCurves();
    AddbdEdgeMidpoint2bdCurve();
}
// For tet mesh, the vertex is NOT critical only if the vertices on its 1-ring
// neighborhood form exactly one single component for each of the two materials
// in interest. And neither of the two components have holes.
//
// TODO: Here I use union-find to analyze the connected component. An
// alternative way is to use DFS/BFS, which has better time complexity if an
// adjacent-list is used. But note that the DFS/BFS search should be limited to
// the 1-ring mesh, so the adjacent-list data-structure is not the one I
// currently stored in the Vs[vi].nbV data-structure, but should be a reduced
// set tailored for the 1-ring mesh for each distinct vertex. Without this key
// data-structure, BFS/DFS might be more time consuming, because more
// vertices/edges which are not on the 1-ring will be visited. One may
// pre-compute the data-structure and switch to BFS/DFS here. It might bring a
// speedup.
#include"../Utility/readers.h"
void testVisualizationFun(const vector<int>& voi,const vector<int>& foi,const vector<vector<float>> &VCoords,const vector<MeshFace3D>& Fs,const vector<int>& label){
    unordered_map<int,int>mappv;
    int ind = 0;
    vector<double>vp;
    vector<int>vMat;
    for(auto a:voi){
        mappv[a] = ind++;
        for(int i=0;i<3;++i)vp.push_back(VCoords[a][i]);
        vMat.push_back(label[a]);
    }
    vector<unsigned int>f2v;
    for(auto a:foi)for(int i=0;i<3;++i)f2v.push_back(mappv[Fs[a].v[i]]);


    string filename = string("/Users/Research/Geometry/MM/ConvertFolder/testVisual/")+to_string(100000);

    vector<int >facesMat(foi.size(),*max_element(vMat.begin(),vMat.end()));
    vector<int >faces2Cs(foi.size(),1);
    vector<int >faces2Cell(foi.size()*2,1);

    writeObjFile(filename,vp,f2v);
    writeVecFile(filename+"_f",facesMat);
    writeVecFile(filename+"_v",vMat);
    writeVecFile(filename+"_fc",faces2Cs);
    writeVecFile(filename+"_fce",faces2Cell);

}
extern bool isCJP;
bool TetMesh::isVCritical(int vi, const vector<int>& label, int mat0,
                          int mat1) const {
    // Vi's 1-ring neighbors.
    const vector<int>& voi = Vs[vi].nbV;
    const vector<int>& eoi = Vs[vi].rE;
    const vector<int>& foi = Vs[vi].rF;

    // The # of V/E/F that ENTIRELY belong to a material.
    // Used for detecting holes using Euler Characteristics.
    vector<vector<int>> mat2nVEF(2, vector<int>(3, 0));

    // Initialize UnionFind's parent data structure.
    vector<int> parent(voi.size());

    unordered_map<int, int> vi2localInd;
    for (int i = 0; i < voi.size(); ++i) {
        int v0 = voi[i];
        parent[i] = i;
        vi2localInd[v0] = i;
        if (label[v0] == mat0) {
            mat2nVEF[0][0]++;
        } else if (label[v0] == mat1) {
            mat2nVEF[1][0]++;
        }
    }

    // Special case: one of material has 0 component. It is the single-vertex
    // island situation, making vi critical.
    // TODO: In real world application, we don't like the single-vertex islands.
    // If you want to do something smart with those, here could be a good place.
    if (mat2nVEF[0][0] == 0 || mat2nVEF[1][0] == 0) return true;

    for (auto ei : eoi) {
        int v0 = Es[ei].v[0];
        int v1 = Es[ei].v[1];

        // Only union-find if the two vertices have the same label and the label is
        // one of the two labels we care about.
        if (label[v0] == label[v1] && (label[v0] == mat0 || label[v0] == mat1)) {
            // UnionFind : Find.
            int p1 = MingUtility::unionfind_find(vi2localInd[v0], parent);
            int p2 = MingUtility::unionfind_find(vi2localInd[v1], parent);
            // UnionFind : Union.
            parent[p1] = p2;
            // Count the #e that has two ends labeled with the same material.
            if (label[v0] == mat0) {
                mat2nVEF[0][1]++;
            } else {
                mat2nVEF[1][1]++;
            }
        }
    }
    vector<int> mats = {mat0, mat1};
    // Map a label to a representive vertex (root of union find) to check # of
    // components for each of the two materials.
    unordered_map<int, int> label2root;
    for (int i = 0; i < voi.size(); ++i) {
        int v0 = voi[i];
        if (label[v0] != mat0 && label[v0] != mat1) continue;
        int p1 = MingUtility::unionfind_find(i, parent);
        int l1 = label[v0];
        // If this label is unseen before, do nothing.
        if (label2root.find(l1) == label2root.end()) {
            label2root[l1] = p1;
        }
        // Otherwise, this label is seen before, check whether the root is the
        // same.
        else {
            if (label2root[l1] != p1) return true;
        }
    }
    // Now, we have exactly one component for each material, check holes.
    for (int mati = 0; mati < 2; ++mati) {
        int mat = mats[mati];
        for (auto fi : foi) {
            if (label[Fs[fi].v[0]] == mat && label[Fs[fi].v[1]] == mat &&
                    label[Fs[fi].v[2]] == mat) {
                mat2nVEF[mati][2]++;
            }
        }
        // If V-E+F!=1, this component has holes. It is a critical vertex.
        if (mat2nVEF[mati][0] - mat2nVEF[mati][1] + mat2nVEF[mati][2] != 1)
            return true;
    }

    int meme = 3;
    if(isCJP)meme = 2;
    if(meme == 1){
        set<int>allmat;
        for (int i = 0; i < voi.size(); ++i) {
            allmat.insert(label[voi[i]]);
        }
        UnionFind unifind;unifind.SetElements(voi);
        for (auto ei : eoi) {
            int v0 = Es[ei].v[0];
            int v1 = Es[ei].v[1];
            if(label[v0]==label[v1])unifind.Union(v0,v1);
        }
        vector<vector<int>>compnn;
        unifind.ExtractComponents(compnn);
        unordered_map<int, set<int>>v2connl;
        for (auto ei : eoi) {
            int v0 = Es[ei].v[0];
            int v1 = Es[ei].v[1];
            v2connl[v0].insert(label[v1]);
            v2connl[v1].insert(label[v0]);
        }
        for(auto &a:compnn)if(label[a[0]]!=mat0 && label[a[0]]!=mat1){
            set<int>containlable;
            for(auto b:a)for(auto c:v2connl[b])containlable.insert(c);
            if(containlable.find(mat0)==containlable.end() || containlable.find(mat1)==containlable.end()){
                cout<<"found connection critical! "<<vi<<' '<<mat0<<' '<<mat1<<endl;
                for(auto a:voi)cout<<label[a]<<' ';cout<<endl;

                //                for(auto b:containlable)cout<<b<<' ';cout<<endl;
                //                cout<<mat0<<' '<<mat1<<endl;
                //                cout<<a.size()<<' '<<compnn.size()<<' '<<label[a[0]]<<endl;
                //                for(auto &a:compnn)cout<<label[a[0]]<<' ';cout<<endl;
                //                cout<<voi.size()<<endl;

                //testVisualizationFun(voi,foi,VCoords,Fs,label);
                // exit(0);
                return true;
            }

        }


    }else if(meme == 2){
        int inc = 0,dec = 0;
        for (auto fi : foi) {
            unordered_set<int>fll;
            bool a1 = false,a0 = false;
            for(int i=0;i<3;++i) fll.insert( label[Fs[fi].v[i]] );
            if(fll.size()!=3)continue;
            a0 = fll.find(mat0)!=fll.end();
            a1 = fll.find(mat1)!=fll.end();
            if(a0 &&!a1)inc++;
            if(!a0 && a1)dec++;

        }
        if(inc!=dec){
            //cout<<"found junction critical! "<<vi<<' '<<mat0<<' '<<mat1<<' '<<inc<<' '<<dec<<endl;
            //for(auto a:voi)cout<<label[a]<<' ';cout<<endl;
            return true;
        }

    }

    return false;
}

bool TetMesh::ExtractLevelSets(int nLabel, const vector<int>& label,
                               vector<vector<int>>& comps,
                               vector<vector<int>>& conns) const {
    // Find all the active edges, and classify them to each label.
    // Note that each active edge will be classified to exactly two labels.
    vector<UnionFind*> label2activeEsUF(nLabel);
    for (int li = 0; li < nLabel; ++li) {
        label2activeEsUF[li] = new UnionFind();
    }
    unordered_set<int> activeEs;

    // cout << "bef |activeEs|" << activeEs.size() << endl;
    // cout << "bef |Es|" << Es.size() << endl;
    // for (int ti = 0; ti < nLabel; ++ti) {
    //   cout << ti << ":{" << label2activeEsUF[ti]->_elements.size() << ","
    //        << label2activeEsUF[ti]->_parent.size() << "}" << endl;
    // }

    for (int ei = 0; ei < Es.size(); ++ei) {
        int label0 = label[Es[ei].v[0]];
        int label1 = label[Es[ei].v[1]];
        // An active edge is found.
        if (label0 != label1) {
            activeEs.insert(ei);
            label2activeEsUF[label0]->AddElement(ei);
            label2activeEsUF[label1]->AddElement(ei);
        }
    }
    // cout << "------------" << endl;
    //
    // cout << "|activeEs|" << activeEs.size() << endl;
    // for (int ti = 0; ti < nLabel; ++ti) {
    //   cout << ti << ":{" << label2activeEsUF[ti]->_elements.size() << ","
    //        << label2activeEsUF[ti]->_parent.size() << "}" << endl;
    // }

    // return false;
    // For each active edge,
    // for (auto it = activeEs.begin(); it != activeEs.end(); ++it) {
    for (auto ei : activeEs) {
        vector<int> twoSideLabels = {label[Es[ei].v[0]], label[Es[ei].v[1]]};
        for (auto nbEi : Es[ei].nbE) {
            if (activeEs.find(nbEi) != activeEs.end()) {
                vector<int> shareSideLabel;
                for (auto labeli : twoSideLabels) {
                    if (labeli == label[Es[nbEi].v[0]] ||
                            labeli == label[Es[nbEi].v[1]]) {
                        shareSideLabel.push_back(labeli);
                    }
                }
                if (!(shareSideLabel.size() == 1 || shareSideLabel.size() == 2)) {
                    MingUtility::printVector(shareSideLabel);
                }
                // TODO remove after done debuging.
                assert(shareSideLabel.size() == 1 || shareSideLabel.size() == 2);
                for (auto labeli : shareSideLabel) {
                    label2activeEsUF[labeli]->Union(ei, nbEi);
                }
            }
        }
    }
    for (int i = 0; i < nLabel; ++i) {
        // Could do a GetNumOfComponents to quickly check whether there are surfaces
        // not using any of the _bdLoops. Todo for early termination.
        vector<vector<int>> labelcomps;
        vector<vector<int>> labelcomp2loopIDs;
        UnionFind* uf = label2activeEsUF[i];
        uf->ExtractComponents(labelcomps);
        unordered_map<int, int> root2labelcompID;

        for (int ci = 0; ci < labelcomps.size(); ++ci) {
            // cout << "-------------\n" << labelcomps.size() << endl;
            // MingUtility::printVector(labelcomps);
            root2labelcompID[uf->Find(labelcomps[ci][0])] = ci;
        }
        labelcomp2loopIDs.resize(labelcomps.size());
        for (auto loopID : _label2bdLoops[i]) {
            int firstEi = _bdLoops[loopID][0];
            int root = uf->Find(firstEi);
            int compi = root2labelcompID[root];
            labelcomp2loopIDs[compi].push_back(loopID);
        }
        for (const auto& ls : labelcomp2loopIDs) {
            // Some surface component is not connected to any of a bdLoops.
            // An island is found. Invalid results.
            if (ls.empty()) {
                // for (int ci = 0; ci < label2activeEsUF.size(); ++ci) {
                //   delete label2activeEsUF[ci];
                // }
                // cout << "------------" << endl;
                // for (int ti = 0; ti < nLabel; ++ti) {
                //   cout << "|activeEs|" << activeEs.size() << endl;
                //   cout << ti << ":{" << label2activeEsUF[ti]->_elements.size() << ","
                //        << label2activeEsUF[ti]->_parent.size() << "}" << endl;
                // }
                // MingUtility::printVector(labelcomp2loopIDs);
                return false;
            }
        }
        // TODO check genus for each surface.
        comps.reserve(comps.size() + labelcomp2loopIDs.size());
        comps.insert(comps.end(), labelcomp2loopIDs.begin(),
                     labelcomp2loopIDs.end());
    }
    // return false;
    // return false;
    // TODO remove the memory in label2activeEsUF, maybe using sharepointer
    // cout << label2activeEsUF.size() << endl;
    // for (int ci = 0; ci < label2activeEsUF.size(); ++ci) {
    //   delete label2activeEsUF[ci];
    // }
    return true;
}
bool TetMesh::ExtractLevelSetsBdComp(const vector<int>& label,
                                     vector<vector<int>>& comps,
                                     vector<vector<int>>& conns) const {
    // return false;
    // TODO ming!
    set<vector<int>> connsSet;
    vector<int> comp2label;
    vector<int> bdComp2comp(_bdComps.size());

    // unordered_set<int> cellActiveFs;
    unordered_map<int, vector<int>> cellActiveE2compPair;
    // unordered_set<int> cellActiveVs;
    unordered_set<int> processedBdComps;
    int compID = 0;
    for (int bdci = 0; bdci < _bdComps.size(); ++bdci) {
        // Case (2), this bdComp has been process in another comp, stop.
        if (processedBdComps.find(bdci) != processedBdComps.end()) {
            continue;
        }
        if (_bdComp2ftEs[bdci].empty()) {
            cout << "Check here. Unless there are bugs in my code, I only found one "
                    "material. If that's true, it's boring..." << endl;
            continue;
        }
        int firstEi = _bdComp2ftEs[bdci][0];
        int curBdCompLabel = _bdComp2label[bdci];
        // 3 cases for ei:
        // (0) Never seen this edge before, aka, 0 record in cellActiveE2compPair.
        //     = A new comp is found.
        // (1) Seen before, but none of its two end of comps has the same label of
        //     bdComp, aka, there is only 1 comp with diff label recorded for ei.
        //     = A new comp is found.
        // (2) Seen before, and at least one of its two end are marked with existing
        //     comps and has the same label of bdComp. (there could be 1 or 2 comps
        //     recorded for case 2)
        //     This case has been handled by processedBdComps.
        //     = This comp has been processed, skip.
        // Check which case it is for ei, by default, it is case (0).
        int whichCase = 0;
        if (cellActiveE2compPair.find(firstEi) != cellActiveE2compPair.end()) {
            whichCase = 1;
        }

        // Case (0) or (1), a new comp has been found.
        // Prepare date structures for level set extraction, using BFS.
        // comps.push_back({bdci});
        comp2label.push_back(curBdCompLabel);
        // bdComp2comp[bdci] = compID;
        unordered_set<int> visitedEs;
        queue<int> eFringe;
        for (auto ei : _bdComp2ftEs[bdci]) {
            eFringe.push(ei);
            visitedEs.insert(ei);
            if (whichCase == 1) {
                auto& compPair = cellActiveE2compPair[ei];
                compPair.push_back(compID);
                sort(compPair.begin(), compPair.end());
                connsSet.insert(compPair);
            } else {
                cellActiveE2compPair[ei] = {compID};
            }
        }
        while (!eFringe.empty()) {
            int ei = eFringe.front();
            eFringe.pop();
            for (auto nbEi : Es[ei].nbE) {
                // If the nb edge is visited in this comp before or it is not an active
                // edge at all, skip it. (An edge is active iff this edge has exactly
                // one vertex with label of interest.)
                //
                // MING: there is also a tradeoff can be made here. If detecting whether
                // an edge is active for curtain material is excuted too-many times, it
                // is probably better to pre-compute it in the beginning.
                // TODO: ming: actually, I need to pre-compute and find all activeE any
                // way, should I pre-compute whether an edge is active for all material?
                // Wasted storage though.
                int nEndWithCurLabel = 0;
                for (int i = 0; i < 2; ++i) {
                    if (label[Es[nbEi].v[i]] == curBdCompLabel) nEndWithCurLabel++;
                }
                if (nEndWithCurLabel != 1 || visitedEs.find(nbEi) != visitedEs.end())
                    continue;
                eFringe.push(nbEi);
                visitedEs.insert(nbEi);
                if (whichCase == 1) {
                    auto& compPair = cellActiveE2compPair[nbEi];
                    compPair.push_back(compID);
                    sort(compPair.begin(), compPair.end());
                    connsSet.insert(compPair);
                } else {
                    cellActiveE2compPair[nbEi] = {compID};
                }
            }
        }
        // A comp is found by now. Check all the bdComps to see which ones belong to
        // this comp. Use this information to mark processedBdComps and calculate
        // the genus of this comp.
        // TODO: need to skip boundary material.
        // TODO: MING: I am wrong!!!! we need to allow material to have >0 genus,
        // for example, the outter layer of a double cylinder.
        // I should probably flood from a single curve instead of a single material.
        // TODO: update: the genus=0 constraint still holds, but we should apply it
        // on a single "surface" not a component. A surface connects several
        // "curves". Each comp might have several surface, we want to enforce the
        // genus=0 property on each surface. But we can do it later
        //
        comps.push_back({});
        for (auto bdcj : _label2bdComps[curBdCompLabel]) {
            // if (bdcj == bdci) continue;
            if (_bdComp2ftEs[bdcj].empty()) {
                cout << "Check here. Unless there are bugs in my code, I only found "
                        "one material. If that's true, it's boring..." << endl;
                return false;
            }
            int firstEi = _bdComp2ftEs[bdcj][0];
            if (cellActiveE2compPair.find(firstEi) != cellActiveE2compPair.end()) {
                const auto& compPair = cellActiveE2compPair[firstEi];
                for (auto tCompID : compPair) {
                    if (tCompID == compID) {
                        comps[compID].push_back(bdcj);
                        bdComp2comp[bdcj] = compID;
                        processedBdComps.insert(bdcj);
                        break;
                    }
                }
            }
        }

        // TODO find other elements besides Es (dual Fs). Calculate genus for each
        // coponent.
        // return false if genus!=0.
        // unordered_set<int> activeFs;
        // unordered_set<int> activeTs;
        // for (const auto& e2compPair : cellActiveE2compPair) {
        //   int ei = e2compPair.first;
        //   for (auto fi : Es[ei].nbF) {
        //     activeFs.insert(fi);
        //   }
        //   for (auto ti : Es[ei].nbT) {
        //     activeTs.insert(ti);
        //   }
        // }
        // int nDualV = activeTs.size() +

        compID++;
    }
    // TODO check the active. Find if there are islands not calculated. if so
    // return false.
    int nActiveE = 0;
    for (int ei = 0; ei < Es.size(); ++ei) {
        if (label[Es[ei].v[0]] != label[Es[ei].v[1]]) {
            nActiveE++;
        }
    }
    return nActiveE == cellActiveE2compPair.size();
    // return false;
}
int TetMesh::ExtractLevelSets(bool allowHighGenus, const vector<int>& label,
                              vector<vector<int>>& comps,
                              vector<vector<int>>& conns) const {
    set<vector<int>> connsSet;
    vector<int> suf2label;
    unordered_map<int, vector<int>> cellActiveE2sufPair;
    int sufID = 0;
    for (int lpi = 0; lpi < _bdLoops.size(); ++lpi) {
        int firstEi = _bdLoops[lpi][0];
        int curSufLabel = _bdComp2label[_bdLoop2Comp[lpi]];
        // 3 cases for ei:
        // (0) Never seen this edge before, aka, 0 record in cellActiveE2sufPair.
        //     = A new suface is found.
        // (1) Seen before, there is only 1 suf recorded in cellActiveE2sufPair, and
        // the suf has a different label.
        //     = A new surface is found.
        // (2) Seen before, there is 1 or 2 sufs recorded in cellActiveE2sufPair,
        // and one of the sufs has the same label of curSufLabel.
        //     = This suf has been processed, skip.
        // Check which case it is, by default, it is case (0).
        int whichCase = 0;
        if (cellActiveE2sufPair.find(firstEi) != cellActiveE2sufPair.end()) {
            whichCase = 1;
            const auto& sufPair = cellActiveE2sufPair[firstEi];
            for (auto si : sufPair) {
                if (suf2label[si] == curSufLabel) {
                    whichCase = 2;
                    break;
                }
            }
        }
        if (whichCase == 2) continue;

        // Case (0) or (1), a new comp has been found.
        // Prepare date structures for level set extraction, using BFS.
        suf2label.push_back(curSufLabel);
        unordered_set<int> visitedEs;
        queue<int> eFringe;
        for (auto ei : _bdLoops[lpi]) {
            eFringe.push(ei);
            visitedEs.insert(ei);
            auto& sufPair = cellActiveE2sufPair[ei];
            if (sufPair.size() == 1) {
                sufPair.push_back(sufID);
                sort(sufPair.begin(), sufPair.end());
                // assert(sufPair.size() == 2);
                connsSet.insert(sufPair);
            } else {
                cellActiveE2sufPair[ei] = {sufID};
            }
        }
        while (!eFringe.empty()) {
            int ei = eFringe.front();
            eFringe.pop();
            for (auto nbEi : Es[ei].nbE) {
                // If the nb edge is visited in this comp before or it is not an active
                // edge at all, skip it. (An edge is active iff this edge has exactly
                // one vertex with label of interest.)
                //
                // MING: there is also a tradeoff can be made here. If detecting whether
                // an edge is active for curtain material is excuted too-many times, it
                // is probably better to pre-compute it in the beginning.
                // TODO: ming: actually, I need to pre-compute and find all activeE any
                // way, should I pre-compute whether an edge is active for all material?
                // Wasted storage though.
                int nEndWithCurLabel = 0;
                for (int i = 0; i < 2; ++i) {
                    if (label[Es[nbEi].v[i]] == curSufLabel) nEndWithCurLabel++;
                }
                if (nEndWithCurLabel != 1 || visitedEs.find(nbEi) != visitedEs.end())
                    continue;
                eFringe.push(nbEi);
                visitedEs.insert(nbEi);
                auto& sufPair = cellActiveE2sufPair[nbEi];
                if (sufPair.size() == 1) {
                    sufPair.push_back(sufID);
                    sort(sufPair.begin(), sufPair.end());
                    // assert(sufPair.size() == 2);
                    connsSet.insert(sufPair);
                } else {
                    cellActiveE2sufPair[nbEi] = {sufID};
                }
            }
        }
        // A suf is found by now. Check all the bdLoops to see which ones belong to
        // this suf. Use this information to get the new sufcomp and calculate
        // the genus of this suf.
        comps.push_back({});
        for (auto lpj : _label2bdLoops[curSufLabel]) {
            int firstEi = _bdLoops[lpj][0];
            if (visitedEs.find(firstEi) != visitedEs.end()) {
                comps[sufID].push_back(lpj);
            }
        }

        // Check whether V-E+F==2, where V = # of active tets + # of active faces    + # of active bdE
        // along the loops on the boundary; E = # of active Fs + # active edges + # of active bdE
        // along the loops on the boundary; F = # of active Es + # of loops.
        // Note that # of active faces along the loops on the boundary = # of active
        // edges along the loops on the boundary, because their dual in 2D (on the
        // boundary surface) are exactly the loops' vertices and edges.
        // So we just need to check:
        // # active tets - # active faces + # active edges + # loops == 2.
        // Terminate early if the genus requirement is not satisfied.
        unordered_set<int> activeFs;
        unordered_set<int> activeTs;
        for (auto ei : visitedEs) {
            for (auto fi : Es[ei].nbF) {
                activeFs.insert(fi);
            }
            for (auto ti : Es[ei].nbT) {
                activeTs.insert(ti);
            }
        }
        int nSufV = activeTs.size();
        int nSufE = activeFs.size();
        int nSufF = visitedEs.size() + comps[sufID].size();
        int eulerCharactor = nSufV - nSufE + nSufF;
        if (eulerCharactor != 2 && !allowHighGenus) return -1;
        // cout << "euler = " << eulerCharactor << endl;

        sufID++;
    }
    // TODO check the active. Find if there are islands not calculated. if so
    // return false.
    int nActiveE = 0;
    for (int ei = 0; ei < Es.size(); ++ei) {
        if (label[Es[ei].v[0]] != label[Es[ei].v[1]]) {
            nActiveE++;
        }
    }
    if (nActiveE != cellActiveE2sufPair.size() && !allowHighGenus) return -2;
    // Make the comps unique by sort this 2D vector.
    for (auto& comp : comps) {
        sort(comp.begin(), comp.end());
    }
    vector<int> compIDs(comps.size());
    iota(compIDs.begin(), compIDs.end(), 0);
    auto comparator = [&](int a, int b) { return comps[a] < comps[b]; };
    sort(compIDs.begin(), compIDs.end(), comparator);
    sort(comps.begin(), comps.end());
    vector<int> compIDsOld2New(comps.size());
    for (int i = 0; i < compIDs.size(); ++i) {
        compIDsOld2New[compIDs[i]] = i;
    }

    conns.clear();
    conns.resize(connsSet.size());
    copy(connsSet.begin(), connsSet.end(), conns.begin());
    for (int i = 0; i < conns.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
            conns[i][j] = compIDsOld2New[conns[i][j]];
        }
    }
    MingUtility::sort2Dvector(conns);
    return 0;
    // return nActiveE == cellActiveE2sufPair.size();
    // return false;
}


int TetMesh::ExtractLevelSets(bool allowHighGenus, const vector<int>& label,
                              vector<vector<int>>& comps,
                              vector<vector<int>>& conns,
                              vector<int> &comp2VEF, int &junctionPoints) const{



    vector<int>lableActiveEs(_label2bdComps.size(),0);
    set<vector<int>> connsSet;
    vector<int> suf2label;
    vector<int>comp2VEF_t;
    unordered_map<int, vector<int>> cellActiveE2sufPair;
    int sufID = 0;
    vector<bool> GlobalActiveTs(Ts.size(),false);
    for (int lpi = 0; lpi < _bdLoops.size(); ++lpi) {
        int firstEi = _bdLoops[lpi][0];
        int curSufLabel = _bdComp2label[_bdLoop2Comp[lpi]];
        // 3 cases for ei:
        // (0) Never seen this edge before, aka, 0 record in cellActiveE2sufPair.
        //     = A new suface is found.
        // (1) Seen before, there is only 1 suf recorded in cellActiveE2sufPair, and
        // the suf has a different label.
        //     = A new surface is found.
        // (2) Seen before, there is 1 or 2 sufs recorded in cellActiveE2sufPair,
        // and one of the sufs has the same label of curSufLabel.
        //     = This suf has been processed, skip.
        // Check which case it is, by default, it is case (0).
        int whichCase = 0;
        if (cellActiveE2sufPair.find(firstEi) != cellActiveE2sufPair.end()) {
            whichCase = 1;
            const auto& sufPair = cellActiveE2sufPair[firstEi];
            for (auto si : sufPair) {
                if (suf2label[si] == curSufLabel) {
                    whichCase = 2;
                    break;
                }
            }
        }
        if (whichCase == 2) continue;

        // Case (0) or (1), a new comp has been found.
        // Prepare date structures for level set extraction, using BFS.
        suf2label.push_back(curSufLabel);
        vector<bool> visitedEs(Es.size(),false);
        queue<int> eFringe;
        for (auto ei : _bdLoops[lpi]) {
            eFringe.push(ei);
            visitedEs[ei]=true;
            auto& sufPair = cellActiveE2sufPair[ei];
            lableActiveEs[curSufLabel]++;
            if (sufPair.size() == 1) {
                sufPair.push_back(sufID);
                sort(sufPair.begin(), sufPair.end());
                // assert(sufPair.size() == 2);
                connsSet.insert(sufPair);
            } else {
                cellActiveE2sufPair[ei] = {sufID};

            }
        }
        while (!eFringe.empty()) {
            int ei = eFringe.front();
            eFringe.pop();
            for (auto nbEi : Es[ei].nbE) {
                // If the nb edge is visited in this comp before or it is not an active
                // edge at all, skip it. (An edge is active iff this edge has exactly
                // one vertex with label of interest.)
                //
                // MING: there is also a tradeoff can be made here. If detecting whether
                // an edge is active for curtain material is excuted too-many times, it
                // is probably better to pre-compute it in the beginning.
                // TODO: ming: actually, I need to pre-compute and find all activeE any
                // way, should I pre-compute whether an edge is active for all material?
                // Wasted storage though.
                int nEndWithCurLabel = 0;
                for (int i = 0; i < 2; ++i) {
                    if (label[Es[nbEi].v[i]] == curSufLabel) nEndWithCurLabel++;
                }
                if (nEndWithCurLabel != 1 || visitedEs[nbEi])
                    continue;
                eFringe.push(nbEi);
                visitedEs[nbEi] = true;
                auto& sufPair = cellActiveE2sufPair[nbEi];
                lableActiveEs[curSufLabel]++;
                if (sufPair.size() == 1) {
                    sufPair.push_back(sufID);
                    sort(sufPair.begin(), sufPair.end());
                    // assert(sufPair.size() == 2);
                    connsSet.insert(sufPair);
                } else {
                    cellActiveE2sufPair[nbEi] = {sufID};
                }
            }
        }
        // A suf is found by now. Check all the bdLoops to see which ones belong to
        // this suf. Use this information to get the new sufcomp and calculate
        // the genus of this suf.
        comps.push_back({});
        for (auto lpj : _label2bdLoops[curSufLabel]) {
            int firstEi = _bdLoops[lpj][0];
            if (visitedEs[firstEi]) {
                comps[sufID].push_back(lpj);
            }
        }

        // Check whether V-E+F==2, where V = # of active tets + # of active faces    + # of active bdE
        // along the loops on the boundary; E = # of active Fs + # active edges + # of active bdE
        // along the loops on the boundary; F = # of active Es + # of loops.
        // Note that # of active faces along the loops on the boundary = # of active
        // edges along the loops on the boundary, because their dual in 2D (on the
        // boundary surface) are exactly the loops' vertices and edges.
        // So we just need to check:
        // # active tets - # active faces + # active edges + # loops == 2.
        // Terminate early if the genus requirement is not satisfied.
        //unordered_set<int> activeFs;
        //unordered_set<int> activeTs;
        int nFs = Fs.size(),nTs = Ts.size(),nEs = Es.size();
        int vEsize = 0;
        vector<bool> activeFs(nFs,false);
        vector<bool> activeTs(nTs,false);
        for(int ei = 0;ei<nEs;++ei)if(visitedEs[ei]){
            for (auto fi : Es[ei].nbF) {
                //activeFs.insert(fi);
                activeFs[fi] = true;
            }
            for (auto ti : Es[ei].nbT) {
                //activeTs.insert(ti);
                activeTs[ti] = true;
                GlobalActiveTs[ti] = true;
            }
            ++vEsize;
        }
        int nSufV = 0;
        int nSufE = 0;
        for(int mj=0;mj<nFs;++mj)if(activeFs[mj])++nSufE;
        for(int mj=0;mj<nTs;++mj)if(activeTs[mj])++nSufV;


        int nSufF = vEsize + comps[sufID].size();
        int eulerCharactor = nSufV - nSufE + nSufF;
        if (eulerCharactor != 2 && !allowHighGenus) return -1;
        // cout << "euler = " << eulerCharactor << endl;
        comp2VEF_t.push_back(eulerCharactor-comps[sufID].size());
        sufID++;
    }
    // TODO check the active. Find if there are islands not calculated. if so
    // return false.
    for(int curlable = 0;curlable<lableActiveEs.size();++curlable){
        int nActiveE = 0;
        for (int ei = 0; ei < Es.size(); ++ei) {
            if ( label[Es[ei].v[0]] != label[Es[ei].v[1]] && (label[Es[ei].v[0]] ==curlable || label[Es[ei].v[1]] ==curlable) ) {
                nActiveE++;
            }
        }
        if (nActiveE !=  lableActiveEs[curlable]&& !allowHighGenus) return -2;
    }
    // Make the comps unique by sort this 2D vector.
    for (auto& comp : comps) {
        sort(comp.begin(), comp.end());
    }
    vector<int> compIDs(comps.size());
    iota(compIDs.begin(), compIDs.end(), 0);
    auto comparator = [&](int a, int b) { return comps[a] < comps[b]; };
    sort(compIDs.begin(), compIDs.end(), comparator);
    sort(comps.begin(), comps.end());
    vector<int> compIDsOld2New(comps.size());
    for (int i = 0; i < compIDs.size(); ++i) {
        compIDsOld2New[compIDs[i]] = i;
    }

    conns.clear();
    conns.resize(connsSet.size());
    copy(connsSet.begin(), connsSet.end(), conns.begin());
    for (int i = 0; i < conns.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
            conns[i][j] = compIDsOld2New[conns[i][j]];
        }
    }
    MingUtility::sort2Dvector(conns);
    comp2VEF.resize(comps.size());
    for(int i=0;i<comps.size();++i)comp2VEF[i] = comp2VEF_t[compIDs[i]];

    vector<bool>kTjj(64,false);
    junctionPoints = 0;
    for(int ti=0;ti<Ts.size();++ti)if(GlobalActiveTs[ti]){
        //unordered_set<int>llll;
        bool isJ = true;
        auto &Tvv = Ts[ti].v;
        //for(int j=0;j<4;++j)llll.insert(label[Tvv[j]]);
        //if(llll.size()==4)junctionPoints++;
        for(int j=0;j<4;++j)if(kTjj[label[Tvv[j]]]){isJ = false;break;}else{kTjj[label[Tvv[j]]]=true;}
        if(isJ)junctionPoints++;
        for(int j=0;j<4;++j)kTjj[label[Tvv[j]]] = false;
    }


    return 0;
    // return nActiveE == cellActiveE2sufPair.size();
    // return false;
}
void TetMesh::ExportMeshToMathematica(const char* filename) const {
    ofstream ofs(filename, ofstream::out);
    ofs << "{";
    MingUtility::writeVector_wB_tbc(VCoords, ofs);
    // Vs
    ofs << "{";
    for (int i = 0; i < Vs.size(); ++i) {
        ofs << "{";
        Vs[i].write(ofs);
        if (i + 1 == Vs.size()) {
            ofs << "}";
        } else {
            ofs << "},";
        }
    }
    ofs << "},";
    // Es
    ofs << "{";
    for (int i = 0; i < Es.size(); ++i) {
        ofs << "{";
        Es[i].write(ofs);
        if (i + 1 == Es.size()) {
            ofs << "}";
        } else {
            ofs << "},";
        }
    }
    ofs << "},";
    // Fs
    ofs << "{";
    for (int i = 0; i < Fs.size(); ++i) {
        ofs << "{";
        Fs[i].write(ofs);
        if (i + 1 == Fs.size()) {
            ofs << "}";
        } else {
            ofs << "},";
        }
    }
    ofs << "},";
    // Ts
    ofs << "{";
    for (int i = 0; i < Ts.size(); ++i) {
        ofs << "{";
        Ts[i].write(ofs);
        if (i + 1 == Ts.size()) {
            ofs << "}";
        } else {
            ofs << "},";
        }
    }
    ofs << "},";
    // BD info
    ofs << "{";
    MingUtility::writeVector(_bdVs, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_bdEs, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_bdFs, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_bdComps, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_bdComp2label, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_label2bdComps, ofs);
    ofs << "},";
    // ofs << "{";
    // MingUtility::writeVector(_bdComp2ftVs, ofs);
    // ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_bdComp2ftEs, ofs);
    ofs << "},";
    // ofs << "{";
    // MingUtility::writeVector(_bdComp2ftFs, ofs);
    // ofs << "}";
    ofs << "{";
    MingUtility::writeVector(_bdLoops, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_bdLoop2Comp, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_bdComp2Loops, ofs);
    ofs << "},";
    ofs << "{";
    MingUtility::writeVector(_label2bdLoops, ofs);
    ofs << "}";
    ofs << "}";
    ofs.close();
}
void TetMesh::Load(const char* filename) {
    vector<vector<float>> inVs;
    vector<vector<int>> inTs;
    ifstream ifs(filename, ios::in);
    if (!ifs) {
        cout << "Error! " << filename << " cannot be found!" << endl;
        return;
    }
    int d = 0, nV = 0, nT = 0;
    if (ifs >> d >> nV >> nT) {
        inVs.resize(nV, vector<float>(3));
        inTs.resize(nT, vector<int>(4));
        int i = 0;
        while (i < nV && ifs >> inVs[i][0] >> inVs[i][1] >> inVs[i][2]) i++;
        i = 0;
        while (i < nT &&
               ifs >> inTs[i][0] >> inTs[i][1] >> inTs[i][2] >> inTs[i][3])
            i++;
    }
    ifs.close();
    Load(inVs, inTs);
}
void TetMesh::Load(const vector<vector<float>>& inVs,
                   const vector<int>& inVMarkers,
                   const vector<vector<int>>& inTs){

    VMarkers = inVMarkers;
    Load(inVs,inTs);
}
void TetMesh::Load(const vector<vector<float>>& inVs,
                   const vector<vector<int>>& inTs) {
    int nV = inVs.size();
    int nT = inTs.size();

    // xid: enumerators, eg. 1,2,3,4
    // xi: index of x in the mesh
    // x: the vertex list of the element x
    // Initialize VCoords, Vs, Es, Fs, Ts.
    VCoords = inVs;
    Vs.resize(nV);
    Ts.resize(nT);
    for (int ti = 0; ti < nT; ++ti) {
        auto& oneT = Ts[ti];
        for (int vid = 0; vid < 4; ++vid) {
            oneT.v[vid] = inTs[ti][vid];
        }
    }

    // Define hashes for edges and faces.
    efHash ehash;
    efHash fhash;

    // Combinations of taking triplets(face) out from 4 vertices; and pairs(edge)
    // out from 3 vertices.
    vector<vector<int>> t2fid;
    vector<vector<int>> f2eid;
    vector<vector<int>> t2eid;
    MingUtility::combination(4, 3, t2fid);
    MingUtility::combination(3, 2, f2eid);
    MingUtility::combination(4, 2, t2eid);
    // Candidate Fs/Es for each Tet/Face.
    vector<vector<int>> candFs(4, vector<int>(3));
    vector<vector<int>> candEs(3, vector<int>(2));
    vector<vector<int>> tetEs(6, vector<int>(2));

    // Read each tet to find edges and faces and fill the data structures.
    for (int ti = 0; ti < nT; ++ti) {
        const auto& t = inTs[ti];
        // cout << "ti=" << ti << endl;
        for (auto vi : t) {
            Vs[vi].nbT.push_back(ti);
        }

        // cout << "lala" << endl;
        // Gather candidate faces, sort each face's vi.
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 3; ++j) {
                // cout << t2fid[i][j] << endl;
                candFs[i][j] = t[t2fid[i][j]];
            }
            sort(candFs[i].begin(), candFs[i].end());
        }

        for (int fid = 0; fid < candFs.size(); ++fid) {
            // cout << "fid=" << fid << endl;
            const auto& f = candFs[fid];
            bool isNewFace = fhash.find(f) == fhash.end();
            int fi = -1;

            // A new face is found.
            if (isNewFace) {
                fi = Fs.size();
                MeshFace3D oneF;
                oneF.v[0] = f[0];
                oneF.v[1] = f[1];
                oneF.v[2] = f[2];
                oneF.nbT.push_back(ti);
                Fs.push_back(oneF);
                fhash[f] = fi;
                for (auto vi : f) {
                    Vs[vi].nbF.push_back(fi);
                }
            }
            // An existing face is met.
            else {
                fi = fhash[f];
                auto& f2nbT = Fs[fi].nbT;
                for (auto oti : f2nbT) {
                    Ts[ti].nbT.push_back(oti);
                }
                Fs[fi].nbT.push_back(ti);
            }
            Ts[ti].face[fid] = fi;

            // Gather candidate edges, sort each edge's vid.
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 2; ++j) {
                    // cout << f2eid[i][j] << endl;
                    candEs[i][j] = f[f2eid[i][j]];
                }
                sort(candEs[i].begin(), candEs[i].end());
            }

            for (int eid = 0; eid < candEs.size(); ++eid) {
                // cout << "eid=" << eid << endl;
                auto& e = candEs[eid];
                int v0 = e[0];
                int v1 = e[1];
                int ei = -1;
                // cout << "enen" << endl;

                // cout << "v0=" << v0 << " v1=" << v1 << endl;
                // A new edge is found.
                if (ehash.find(e) == ehash.end()) {
                    // cout << "enenIF" << endl;
                    ei = Es.size();
                    MeshEdge3D oneE;
                    oneE.v[0] = v0;
                    oneE.v[1] = v1;
                    oneE.nbF.push_back(fi);
                    // oneE.nbT.push_back(ti);  // MING: I will order the nbT later, so we
                    // can skip this step actually.
                    Es.push_back(oneE);
                    ehash[e] = ei;
                    Vs[v0].nbV.push_back(v1);
                    Vs[v1].nbV.push_back(v0);
                    for (auto vi : e) {
                        Vs[vi].nbE.push_back(ei);
                    }
                }
                // And existing edge is met.
                else {
                    // cout << "enenELSE" << endl;
                    ei = ehash[e];
                    if (isNewFace) {
                        Es[ei].nbF.push_back(fi);
                    }
                    // This is actually wrong.
                    // Es[ei].nbT.push_back(ti);  // MING: I will order the nbT later, so
                    // we
                    // can skip this step actually.
                }
                // cout << "keke" << endl;
                // Ts[ti].edge[eid] = ei;
                if (isNewFace) {
                    Fs[fi].edge[eid] = ei;
                }
                // cout << "lala" << endl;
            }
        }

        // Update the edges for tet ti.
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 2; ++j) {
                tetEs[i][j] = t[t2eid[i][j]];
            }
            sort(tetEs[i].begin(), tetEs[i].end());
            Ts[ti].edge[i] = ehash[tetEs[i]];
        }
    }
    for (int ei = 0; ei < Es.size(); ++ei) {
        for (auto fi : Es[ei].nbF) {
            for (int i = 0; i < 3; ++i) {
                if (Fs[fi].edge[i] != ei) {
                    Es[ei].nbE.push_back(Fs[fi].edge[i]);
                }
            }
        }
    }
    // Find the one ring neighborhood of each vertex.
    vector<int> oneF(3);
    for (int vi = 0; vi < nV; ++vi) {
        vector<int>& nbT = Vs[vi].nbT;
        vector<int>& rF = Vs[vi].rF;
        // For each its neighboring tet, find the single face that do not use vi.
        for (int ti : nbT) {
            int pos = 0;
            for (int vid = 0; vid < 4; ++vid) {
                int tvi = Ts[ti].v[vid];
                if (tvi == vi) continue;
                oneF[pos++] = tvi;
            }
            sort(oneF.begin(), oneF.end());
            rF.push_back(fhash[oneF]);
        }
        // For each face in the one ring neighborhood, add its edges to rE.
        unordered_set<int> rE;
        for (int fi : rF) {
            for (int eid = 0; eid < 3; ++eid) {
                rE.insert(Fs[fi].edge[eid]);
            }
        }
        Vs[vi].rE.resize(rE.size());
        copy(rE.begin(), rE.end(), Vs[vi].rE.begin());
    }
    // TODO isBD is not set yet because right now it is not necessary. Add the
    // code later when isBD is needed.

    // Pre-calculate all the centroid of each tet.
    _tVCoords.resize(Ts.size());
    for (int ti = 0; ti < Ts.size(); ++ti) {
        vector<float> oneV(3, 0.0f);
        for (int i = 0; i < 4; ++i) {
            MingUtility::SumUpTwoVectors(VCoords[Ts[ti].v[i]], oneV);
        }
        MingUtility::MultiplyVectorsByValue(0.25f, oneV);
        _tVCoords[ti] = oneV;
    }

    // Reorder the nbT of each E such that neighbor Ts are consecutive. The dual
    // forms a close loop.
    // Also adjust the orientation such that the normal of the 1-ring neighborhood
    // points to the first vertex of the edge.
    for (int ei = 0; ei < Es.size(); ++ei) {
        vector<vector<int>> conns;
        unordered_set<int> singleNbTs;
        for (auto fi : Es[ei].nbF) {
            if (Fs[fi].nbT.size() == 2) {
                conns.push_back({Fs[fi].nbT[0], Fs[fi].nbT[1]});
            } else {
                singleNbTs.insert(Fs[fi].nbT[0]);
            }
        }
        // If there is only one tet using ei, set Es[ei].nbT and continue.
        if (singleNbTs.size() == 1) {
            Es[ei].nbT.push_back(*singleNbTs.begin());
            continue;
        }
        // vector<int> tmpTs;
        // MingUtility::ExtractLineFromSegments(conns, tmpTs);
        // if (!(tmpTs.size() == Es[ei].nbT.size())) {
        //   MingUtility::printVector(tmpTs);
        //   MingUtility::printVector(Es[ei].nbT);
        //   MingUtility::printVector(conns);
        //
        // }
        // assert(tmpTs.size() == Es[ei].nbT.size());
        MingUtility::ExtractLineFromSegments(conns, Es[ei].nbT);
        // It is possible that conns is empty, when ei is used by only one tet.
        // In this case, there is only one tet, no need to change orientation.
        if (Es[ei].nbT.size() < 2) continue;
        int t0 = Es[ei].nbT[0];
        int t1 = Es[ei].nbT[1];
        int f1 = -1;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (Ts[t0].face[i] == Ts[t1].face[j]) {
                    f1 = Ts[t0].face[i];
                    break;
                }
            }
        }
        // if (shareFs.size() != 1) {
        //   cout << "shareFs " << shareFs.size() << endl;
        //   cout << t0 << " " << t1 << ", ei=" << ei << endl;
        //   MingUtility::printVector(shareFs);
        //   MingUtility::printVector(conns);
        // }
        // assert(shareFs.size() == 1);
        // f1 = shareFs[0];
        assert(f1 != -1);
        int vid1 =
                Fs[f1].v[0] + Fs[f1].v[1] + Fs[f1].v[2] - Es[ei].v[0] - Es[ei].v[1];
        int vid0 = Ts[t0].v[0] + Ts[t0].v[1] + Ts[t0].v[2] + Ts[t0].v[3] -
                Fs[f1].v[0] - Fs[f1].v[1] - Fs[f1].v[2];
        Vector3f v0(VCoords[vid0].data());
        Vector3f v1(VCoords[vid1].data());
        Vector3f ve0(VCoords[Es[ei].v[0]].data());
        Vector3f ve1(VCoords[Es[ei].v[1]].data());
        if (((v0 - ve0).cross(v1 - ve0)).dot(ve0 - ve1) < 0) {
            reverse(Es[ei].nbT.begin(), Es[ei].nbT.end());
        }
    }

    // Find _bdVs, _bdEs, _bdFs.
    // For each face, if its |nbTs| = 1, it is a bdF and its Es/Vs are bdEs/bdVs.
    unordered_set<int> bdVsSet;
    unordered_set<int> bdEsSet;
    for (int fi = 0; fi < Fs.size(); ++fi) {
        if (Fs[fi].nbT.size() == 1) {
            _bdFs.push_back(fi);
            Fs[fi].isBD = true;
            for (int i = 0; i < 3; ++i) {
                int vi = Fs[fi].v[i];
                int ei = Fs[fi].edge[i];
                bdVsSet.insert(vi);
                bdEsSet.insert(ei);
                Vs[vi].isBD = true;
                Es[ei].isBD = true;
            }
        }
    }
    _bdVs.resize(bdVsSet.size());
    _bdEs.resize(bdEsSet.size());
    copy(bdVsSet.begin(), bdVsSet.end(), _bdVs.begin());
    copy(bdEsSet.begin(), bdEsSet.end(), _bdEs.begin());

    // Compute nbWeight for each vertex.
    // nbWeight[x] = Sum[volume of x's nb tets]. This value is used for scoring.
    // Volume[tet] = 1/3! * |a.(bxc)| where a,b,c are edge vectors from a vertex.
    // Here for simplicity, constant 1/3! is skipped.
    vector<float> tetVols(Ts.size());
    for (int i = 0; i < Ts.size(); ++i) {
        Vector3f v0(VCoords[Ts[i].v[0]].data());
        Vector3f v1(VCoords[Ts[i].v[1]].data());
        Vector3f v2(VCoords[Ts[i].v[2]].data());
        Vector3f v3(VCoords[Ts[i].v[3]].data());
        tetVols[i] = abs((v1 - v0).dot((v2 - v0).cross(v3 - v0)));
    }
    _vNbWeights.resize(Vs.size());
    for (int i = 0; i < Vs.size(); ++i) {
        _vNbWeights[i] = 0.0f;
        for (auto ti : Vs[i].nbT) {
            _vNbWeights[i] += tetVols[ti];
        }
    }
}

void TetMesh::GetActiveBDEsNFs(const vector<int>& label) {
    unordered_set<int> activeEsSet;
    unordered_set<int> activeFsSet;
    _bdCurveEsFlag.clear();
    _bdEdgeMidPoints.clear();
    _bdEMInd.clear();
    _bdActiveE2Vs_fornonplcs.clear();


    bool isbdE = VMarkers.size()==VCoords.size();
    // isbdE = false;
    for (const auto& oneCompFtEs : _bdComp2ftEs) {
        copy(oneCompFtEs.begin(), oneCompFtEs.end(),
             inserter(activeEsSet, activeEsSet.end()));
    }
    int num = 0;
    for(auto a:VMarkers)if(a==-1){/*cout<<a<<' ';*/++num;}cout<<endl;
    cout<<num<<endl;
    for (auto ei : activeEsSet) {
        vector<int> fPair;
        for (auto fi : Es[ei].nbF) {
            if (Fs[fi].isBD) {
                activeFsSet.insert(fi);
                fPair.push_back(fi);
            }
        }
        _bdCurveEs.push_back(fPair);

        // add Zhiyang: to find the connect points for non-paralell crosssections
        if( isbdE && Es[ei].isBD && VMarkers[Es[ei].v[0]] ==-1 && VMarkers[Es[ei].v[1]] ==-1 ){
            _bdCurveEsFlag.push_back(_bdCurveEs.size()-1);
            vector<float>evc(3);
            for(int kk=0;kk<3;++kk)evc[kk]=(VCoords[Es[ei].v[0]][kk]+VCoords[Es[ei].v[1]][kk])/2;
            _bdEdgeMidPoints.push_back(evc);
            _bdActiveE2Vs_fornonplcs.push_back(Es[ei].v[0]);
            _bdActiveE2Vs_fornonplcs.push_back(Es[ei].v[1]);


        }
    }
    _bdActiveEs.resize(activeEsSet.size());
    _bdActiveFs.resize(activeFsSet.size());
    _bdActiveF2Vs.clear();


    copy(activeEsSet.begin(), activeEsSet.end(), _bdActiveEs.begin());
    copy(activeFsSet.begin(), activeFsSet.end(), _bdActiveFs.begin());
    unordered_map<int, int> activeF2crvVID;
    for (int vid = 0; vid < _bdActiveFs.size(); ++vid) {
        int fi = _bdActiveFs[vid];
        activeF2crvVID[fi] = vid;
    }
    for (int vid = 0; vid < _bdActiveFs.size(); ++vid) {
        int fi = _bdActiveFs[vid];
        for(int j=0;j<3;++j)_bdActiveF2Vs.push_back(Fs[fi].v[j]);
    }


    for (auto& oneE : _bdCurveEs) {
        for (auto& oneV : oneE) {
            oneV = activeF2crvVID[oneV];
        }
    }
    for (int i = 0; i < _bdActiveFs.size(); ++i) {
        _fi2bdActFInd[_bdActiveFs[i]] = i;
    }

    //for(int i=0;i<_bdEdgeMidPoints.size();++i)_bdEMInd.push_back(i+_bdActiveFs.size());

}

void TetMesh::AddbdEdgeMidpoint2bdCurve(){

    _bdCurveEsInd.clear();
    _bdCurveEsInd.resize(Es.size(),-1);
    for(int i=0;i<_bdEdgeMidPoints.size();++i)_bdEMInd.push_back(i+_bdCurveVs.size());

    for(auto &a:_bdEdgeMidPoints)_bdCurveVs.push_back(a);

    for(int i=0;i<_bdCurveEsFlag.size();++i){

        _bdCurveEsInd[_bdActiveEs[_bdCurveEsFlag[i]]] = _bdEMInd[i];
        auto &oriE1 = _bdCurveEs[_bdCurveEsFlag[i]];
        auto oriE2 = oriE1;
        oriE1[1] = _bdEMInd[i];oriE2[0] = _bdEMInd[i];
        //_bdCurveEs[_bdCurveEsFlag[i]] = oriE1;
        _bdCurveEs.push_back(oriE2);
    }

}
void TetMesh::ExtractBoundaryCurves(const vector<int>& label,
                                    vector<vector<float>>& crvVs,
                                    vector<vector<int>>& crvEs) const {
    crvEs = _bdCurveEs;
    for (int vid = 0; vid < _bdActiveFs.size(); ++vid) {
        int fi = _bdActiveFs[vid];
        vector<float> centroid(3, 0.0f);
        for (int i = 0; i < 3; ++i) {
            std::transform(centroid.begin(), centroid.end(),
                           VCoords[Fs[fi].v[i]].begin(), centroid.begin(),
                    std::plus<float>());
        }
        MingUtility::MultiplyVectorsByValue(1.0f / 3.0f, centroid);
        crvVs.push_back(centroid);
    }
}
void TetMesh::ExtractSmoothBoundaryCurves(const vector<int>& label,
                                          vector<vector<float>>& crvVs,
                                          vector<vector<int>>& crvEs) const {
    crvEs = _bdCurveEs;
    // Prepare data structure.
    unordered_map<int, vector<int>> crvV2nbCrvVs;
    for (const auto& oneE : crvEs) {
        crvV2nbCrvVs[oneE[0]].push_back(oneE[1]);
        crvV2nbCrvVs[oneE[1]].push_back(oneE[0]);
    }
    // Get original crvVs positions (centroid of each triangle).
    vector<vector<float>> orgCrvVs;
    for (int vid = 0; vid < _bdActiveFs.size(); ++vid) {
        int fi = _bdActiveFs[vid];
        vector<float> centroid(3, 0.0f);
        for (int i = 0; i < 3; ++i) {
            std::transform(centroid.begin(), centroid.end(),
                           VCoords[Fs[fi].v[i]].begin(), centroid.begin(),
                    std::plus<float>());
        }
        MingUtility::MultiplyVectorsByValue(1.0f / 3.0f, centroid);
        orgCrvVs.push_back(centroid);
    }
    // Triplets of adjacent crvVs.
    vector<vector<int>> tripletCrvVs;
    for (const auto& ele : crvV2nbCrvVs) {
        // MingUtility::printVector(ele.second);
        if (ele.second.size() == 2) {
            tripletCrvVs.push_back({ele.second[0], ele.first, ele.second[1]});
        } else {
            tripletCrvVs.push_back({ele.second[0], ele.first, ele.second[1]});
            tripletCrvVs.push_back({ele.second[0], ele.first, ele.second[2]});
            tripletCrvVs.push_back({ele.second[1], ele.first, ele.second[2]});
        }
    }
    int nCrvV = orgCrvVs.size();
    int nTriplets = tripletCrvVs.size();
    // Optimize.
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
        double lambda = 0.05;
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);
        // Quiet output from Gurobi.
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        vector<GRBVar> A(nCrvV);
        vector<GRBVar> B(nCrvV);
        vector<GRBVar> C(nCrvV);
        vector<GRBVar> X(nCrvV);
        vector<GRBVar> Y(nCrvV);
        vector<GRBVar> Z(nCrvV);

        // Create variables.
        for (int i = 0; i < nCrvV; ++i) {
            A[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            B[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            C[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            X[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            Y[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            Z[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }
        model.update();

        // Set objective, obj = smoothTerm + lambda*intactTerm.
        GRBQuadExpr obj = 0;
        GRBQuadExpr intactTerm = 0;
        for (int i = 0; i < nCrvV; ++i) {
            const vector<float>& centroid = orgCrvVs[i];
            intactTerm += (X[i] - centroid[0]) * (X[i] - centroid[0]) +
                    (Y[i] - centroid[1]) * (Y[i] - centroid[1]) +
                    (Z[i] - centroid[2]) * (Z[i] - centroid[2]);
        }
        GRBQuadExpr smoothTerm = 0;
        for (int j = 0; j < nTriplets; ++j) {
            const vector<int>& t = tripletCrvVs[j];
            GRBLinExpr tmp;
            tmp = (X[t[1]] - 0.5 * X[t[0]] - 0.5 * X[t[2]]);
            smoothTerm += tmp * tmp;
            tmp = (Y[t[1]] - 0.5 * Y[t[0]] - 0.5 * Y[t[2]]);
            smoothTerm += tmp * tmp;
            tmp = (Z[t[1]] - 0.5 * Z[t[0]] - 0.5 * Z[t[2]]);
            smoothTerm += tmp * tmp;
        }
        obj = smoothTerm + lambda * intactTerm;
        model.setObjective(obj, GRB_MINIMIZE);

        // Add linear constraint.
        for (int i = 0; i < nCrvV; ++i) {
            const vector<float>& v0 = VCoords[Fs[_bdActiveFs[i]].v[0]];
            const vector<float>& v1 = VCoords[Fs[_bdActiveFs[i]].v[1]];
            const vector<float>& v2 = VCoords[Fs[_bdActiveFs[i]].v[2]];
            model.addConstr(A[i] + B[i] + C[i] == 1.0);
            model.addConstr(X[i] - v0[0] * A[i] - v1[0] * B[i] - v2[0] * C[i] == 0.0);
            model.addConstr(Y[i] - v0[1] * A[i] - v1[1] * B[i] - v2[1] * C[i] == 0.0);
            model.addConstr(Z[i] - v0[2] * A[i] - v1[2] * B[i] - v2[2] * C[i] == 0.0);
        }

        // Optimize model.
        model.optimize();

        // Extract Results.
        crvVs.resize(nCrvV, vector<float>(3));
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            // cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
            for (int i = 0; i < nCrvV; ++i) {
                const vector<float>& v0 = VCoords[Fs[_bdActiveFs[i]].v[0]];
                const vector<float>& v1 = VCoords[Fs[_bdActiveFs[i]].v[1]];
                const vector<float>& v2 = VCoords[Fs[_bdActiveFs[i]].v[2]];
                crvVs[i][0] = X[i].get(GRB_DoubleAttr_X);
                crvVs[i][1] = Y[i].get(GRB_DoubleAttr_X);
                crvVs[i][2] = Z[i].get(GRB_DoubleAttr_X);
            }
        }
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
}
void TetMesh::SmoothBoundaryCurves() {
    // Prepare data structure.
    unordered_map<int, vector<int>> crvV2nbCrvVs;
    vector<bool>isCrossEdges(_bdCurveEs.size(),false);
    for(auto a:_bdCurveEsFlag)isCrossEdges[a] = true;
    for(int i=0;i<_bdCurveEs.size();++i)if(!isCrossEdges[i]){
        const auto& oneE  = _bdCurveEs[i];
        crvV2nbCrvVs[oneE[0]].push_back(oneE[1]);
        crvV2nbCrvVs[oneE[1]].push_back(oneE[0]);
    }
    //    for (const auto& oneE : _bdCurveEs) {
    //        crvV2nbCrvVs[oneE[0]].push_back(oneE[1]);
    //        crvV2nbCrvVs[oneE[1]].push_back(oneE[0]);
    //    }

    // Get original crvVs positions (centroid of each triangle).
    vector<vector<float>> orgCrvVs;
    for (int vid = 0; vid < _bdActiveFs.size(); ++vid) {
        int fi = _bdActiveFs[vid];
        vector<float> centroid(3, 0.0f);
        for (int i = 0; i < 3; ++i) {
            std::transform(centroid.begin(), centroid.end(),
                           VCoords[Fs[fi].v[i]].begin(), centroid.begin(),
                    std::plus<float>());
        }
        MingUtility::MultiplyVectorsByValue(1.0f / 3.0f, centroid);
        orgCrvVs.push_back(centroid);
    }
    // Triplets of adjacent crvVs.
    vector<vector<int>> tripletCrvVs;
    for (const auto& ele : crvV2nbCrvVs) {
        // MingUtility::printVector(ele.second);
        if (ele.second.size() == 2) {
            tripletCrvVs.push_back({ele.second[0], ele.first, ele.second[1]});
        } else if(ele.second.size() == 3){
            tripletCrvVs.push_back({ele.second[0], ele.first, ele.second[1]});
            tripletCrvVs.push_back({ele.second[0], ele.first, ele.second[2]});
            tripletCrvVs.push_back({ele.second[1], ele.first, ele.second[2]});
        }
    }
    int nCrvV = orgCrvVs.size();
    int nTriplets = tripletCrvVs.size();

    _bdCurveVs =  orgCrvVs;
    return;
    // Optimize.
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
        double lambda = 0.02;
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);
        // Quiet output from Gurobi.
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        vector<GRBVar> A(nCrvV);
        vector<GRBVar> B(nCrvV);
        vector<GRBVar> C(nCrvV);
        vector<GRBVar> X(nCrvV);
        vector<GRBVar> Y(nCrvV);
        vector<GRBVar> Z(nCrvV);

        // Create variables.
        for (int i = 0; i < nCrvV; ++i) {
            A[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            B[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            C[i] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            X[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            Y[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            Z[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }
        model.update();

        // Set objective, obj = smoothTerm + lambda*intactTerm.
        GRBQuadExpr obj = 0;
        GRBQuadExpr intactTerm = 0;
        for (int i = 0; i < nCrvV; ++i) {
            const vector<float>& centroid = orgCrvVs[i];
            intactTerm += (X[i] - centroid[0]) * (X[i] - centroid[0]) +
                    (Y[i] - centroid[1]) * (Y[i] - centroid[1]) +
                    (Z[i] - centroid[2]) * (Z[i] - centroid[2]);
        }
        GRBQuadExpr smoothTerm = 0;
        for (int j = 0; j < nTriplets; ++j) {
            const vector<int>& t = tripletCrvVs[j];
            GRBLinExpr tmp;
            tmp = (X[t[1]] - 0.5 * X[t[0]] - 0.5 * X[t[2]]);
            smoothTerm += tmp * tmp;
            tmp = (Y[t[1]] - 0.5 * Y[t[0]] - 0.5 * Y[t[2]]);
            smoothTerm += tmp * tmp;
            tmp = (Z[t[1]] - 0.5 * Z[t[0]] - 0.5 * Z[t[2]]);
            smoothTerm += tmp * tmp;
        }

        for(int j=0;j<_bdCurveEsFlag.size();++j){

            auto &be = _bdCurveEs[_bdCurveEsFlag[j]];
            auto &emidp = _bdEdgeMidPoints[j];
            GRBLinExpr tmp;
            for(int k=0;k<2;++k){
                auto ev = be[k];
                for(auto &vv : crvV2nbCrvVs[ev]){
                    tmp = (X[ev] - 0.5 * X[vv] - 0.5 * emidp[0]);
                    smoothTerm += tmp * tmp;
                    tmp = (Y[ev] - 0.5 * Y[vv] - 0.5 * emidp[1]);
                    smoothTerm += tmp * tmp;
                    tmp = (Z[ev] - 0.5 * Z[vv] - 0.5 * emidp[2]);
                    smoothTerm += tmp * tmp;
                }
            }





        }


        obj = smoothTerm + lambda * intactTerm;
        model.setObjective(obj, GRB_MINIMIZE);

        // Add linear constraint.
        for (int i = 0; i < nCrvV; ++i) {
            const vector<float>& v0 = VCoords[Fs[_bdActiveFs[i]].v[0]];
            const vector<float>& v1 = VCoords[Fs[_bdActiveFs[i]].v[1]];
            const vector<float>& v2 = VCoords[Fs[_bdActiveFs[i]].v[2]];
            model.addConstr(A[i] + B[i] + C[i] == 1.0);
            model.addConstr(X[i] - v0[0] * A[i] - v1[0] * B[i] - v2[0] * C[i] == 0.0);
            model.addConstr(Y[i] - v0[1] * A[i] - v1[1] * B[i] - v2[1] * C[i] == 0.0);
            model.addConstr(Z[i] - v0[2] * A[i] - v1[2] * B[i] - v2[2] * C[i] == 0.0);
        }

        // Optimize model.
        model.optimize();

        // Extract Results.
        _bdCurveVs.resize(nCrvV, vector<float>(3));
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            // cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
            for (int i = 0; i < nCrvV; ++i) {
                _bdCurveVs[i][0] = X[i].get(GRB_DoubleAttr_X);
                _bdCurveVs[i][1] = Y[i].get(GRB_DoubleAttr_X);
                _bdCurveVs[i][2] = Z[i].get(GRB_DoubleAttr_X);
            }
        }
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during optimization" << endl;
    }
    // MingUtility::writeVector(_bdCurveVs, "../crvVs.txt");
    // MingUtility::writeVector(_bdCurveEs, "../crvEs.txt");
}
void TetMesh::GenerateLevelSetFromLabeling(
        const vector<int>& label, vector<vector<float>>& sufVs,
        vector<vector<int>>& sufFs, vector<vector<int>>& sufFMats) const {
    unordered_map<int, int> t2vHash;
    vector<int> activeEs;
    sufVs = _bdCurveVs;
    vector<vector<int>> sufPs;  // Polygons.
    // Find activeEs, sufVs, sufPs.
    vector<int>bbE2V(Es.size(),-1);
    for (int ei = 0; ei < Es.size(); ++ei) {
        if (label[Es[ei].v[0]] == label[Es[ei].v[1]]) continue;
        activeEs.push_back(ei);
        vector<int> oneP;
        for (int ti : Es[ei].nbT) {
            if (t2vHash.find(ti) == t2vHash.end()) {
                t2vHash[ti] = sufVs.size();
                sufVs.push_back(_tVCoords[ti]);
            }
            oneP.push_back(t2vHash[ti]);
        }

        // Edges on the boudary will be processed differently: need to add dualv
        // on boundary Fs.
        if (Es[ei].isBD) {

            vector<int> nbBdFs;
            for (auto fi : Es[ei].nbF) {
                if (Fs[fi].isBD) {
                    nbBdFs.push_back(fi);
                }
            }
            assert(nbBdFs.size() == 2);
            bool inorder = false;
            int lastTi = Es[ei].nbT.back();
            // if (!(lastTi >=0 && lastTi <Ts.size())) {
            //   cout << "lastTi = " << lastTi << endl;
            //   MingUtility::printVector(Es[ei].nbT);
            // }
            // assert(lastTi >=0 && lastTi <Ts.size());





            if (Es[ei].nbT.size() > 1) {
                for (int i = 0; i < 4; ++i) {
                    if (Ts[lastTi].face[i] == nbBdFs[0]) {
                        inorder = true;
                        break;
                    }
                }
            }
            // If there is only one nbT, need to check the order by looking at the
            // normal.
            else {

                assert(nbBdFs.size() == 2);



                //bbE2V[ei] = sufVs.size();
                //vector<float>evc(3);
                //for(int kk=0;kk<3;++kk)evc[kk]=(VCoords[Es[ei].v[0]][kk]+VCoords[Es[ei].v[1]][kk])/2;
                //for(int kk=0;kk<3;++kk)evc[kk]=VCoords[Es[ei].v[0]][kk];
                //for(int kk=0;kk<3;++kk)evc[kk]=VCoords[Es[ei].v[1]][kk];
                //if (isflingeE)sufVs.push_back(evc);




                int f0 = nbBdFs[0];
                int f1 = nbBdFs[1];
                int vid0 =
                        Fs[f0].v[0] + Fs[f0].v[1] + Fs[f0].v[2] - Es[ei].v[0] - Es[ei].v[1];
                int vid1 =
                        Fs[f1].v[0] + Fs[f1].v[1] + Fs[f1].v[2] - Es[ei].v[0] - Es[ei].v[1];
                Vector3f v0(VCoords[vid0].data());
                Vector3f v1(VCoords[vid1].data());
                Vector3f ve0(VCoords[Es[ei].v[0]].data());
                Vector3f ve1(VCoords[Es[ei].v[1]].data());
                if (((v0 - ve0).cross(v1 - ve0)).dot(ve0 - ve1) < 0) {
                    inorder = true;
                }
            }

            bool isflingeE = _bdCurveEsInd[ei]!=-1;
            //assert(bbE2V[ei]>=0);
            if(_fi2bdActFInd.find(nbBdFs[0])==_fi2bdActFInd.end() || _fi2bdActFInd.find(nbBdFs[1])==_fi2bdActFInd.end()){
                cout<<"error"<<endl;
            }
            if (inorder) {
                oneP.push_back(_fi2bdActFInd.at(nbBdFs[0]));
                if (isflingeE) oneP.push_back(_bdCurveEsInd[ei]);
                oneP.push_back(_fi2bdActFInd.at(nbBdFs[1]));
            } else {
                oneP.push_back(_fi2bdActFInd.at(nbBdFs[1]));
                if (isflingeE) oneP.push_back(_bdCurveEsInd[ei]);
                oneP.push_back(_fi2bdActFInd.at(nbBdFs[0]));
            }
        }
        // if (oneP.size() < 3) {
        //   cout << oneP.size() << endl;
        // }
        assert(oneP.size() >= 3);
        sufPs.push_back(oneP);
    }
    // Break polygons into triangles and record the label of the two sides.
    for (int i = 0; i < activeEs.size(); ++i) {
        int ei = activeEs[i];
        const auto& oneP = sufPs[i];
        for(auto a:oneP)assert(a>-1);

        // Split.
        vector<int> mat = {label[Es[ei].v[0]], label[Es[ei].v[1]]};
        for (int j = 2; j < oneP.size(); ++j) {
            sufFs.push_back({oneP[0], oneP[j - 1], oneP[j]});
            sufFMats.push_back(mat);
        }
    }
}

//void TetMesh::GenerateLevelSetFromLabeling(
//    const vector<int>& label, vector<vector<float>>& sufVs,
//    vector<vector<int>>& sufFs, vector<vector<int>>& sufFMats) const {
//  unordered_map<int, int> t2vHash;
//  vector<int> activeEs;
//  sufVs = _bdCurveVs;
//  vector<vector<int>> sufPs;  // Polygons.
//  // Find activeEs, sufVs, sufPs.
//  for (int ei = 0; ei < Es.size(); ++ei) {
//    if (label[Es[ei].v[0]] == label[Es[ei].v[1]]) continue;
//    activeEs.push_back(ei);
//    vector<int> oneP;
//    for (int ti : Es[ei].nbT) {
//      if (t2vHash.find(ti) == t2vHash.end()) {
//        t2vHash[ti] = sufVs.size();
//        sufVs.push_back(_tVCoords[ti]);
//      }
//      oneP.push_back(t2vHash[ti]);
//    }

//    // Edges on the boudary will be processed differently: need to add dualv
//    // on boundary Fs.
//    if (Es[ei].isBD) {
//
//      vector<int> nbBdFs;
//      for (auto fi : Es[ei].nbF) {
//        if (Fs[fi].isBD) {
//          nbBdFs.push_back(fi);
//        }
//      }
//      assert(nbBdFs.size() == 2);
//      bool inorder = false;
//      int lastTi = Es[ei].nbT.back();
//      // if (!(lastTi >=0 && lastTi <Ts.size())) {
//      //   cout << "lastTi = " << lastTi << endl;
//      //   MingUtility::printVector(Es[ei].nbT);
//      // }
//      // assert(lastTi >=0 && lastTi <Ts.size());
//      if (Es[ei].nbT.size() > 1) {
//        for (int i = 0; i < 4; ++i) {
//          if (Ts[lastTi].face[i] == nbBdFs[0]) {
//            inorder = true;
//            break;
//          }
//        }
//      }
//      // If there is only one nbT, need to check the order by looking at the
//      // normal.
//      else {
//
//        assert(nbBdFs.size() == 2);
//        int f0 = nbBdFs[0];
//        int f1 = nbBdFs[1];
//        int vid0 =
//            Fs[f0].v[0] + Fs[f0].v[1] + Fs[f0].v[2] - Es[ei].v[0] - Es[ei].v[1];
//        int vid1 =
//            Fs[f1].v[0] + Fs[f1].v[1] + Fs[f1].v[2] - Es[ei].v[0] - Es[ei].v[1];
//        Vector3f v0(VCoords[vid0].data());
//        Vector3f v1(VCoords[vid1].data());
//        Vector3f ve0(VCoords[Es[ei].v[0]].data());
//        Vector3f ve1(VCoords[Es[ei].v[1]].data());
//        if (((v0 - ve0).cross(v1 - ve0)).dot(ve0 - ve1) < 0) {
//          inorder = true;
//        }
//      }
//      if (inorder) {
//        oneP.push_back(_fi2bdActFInd.at(nbBdFs[0]));

//        oneP.push_back(_fi2bdActFInd.at(nbBdFs[1]));
//      } else {
//        oneP.push_back(_fi2bdActFInd.at(nbBdFs[1]));

//        oneP.push_back(_fi2bdActFInd.at(nbBdFs[0]));
//      }
//    }
//    // if (oneP.size() < 3) {
//    //   cout << oneP.size() << endl;
//    // }
//    assert(oneP.size() >= 3);
//    sufPs.push_back(oneP);
//  }
//  // Break polygons into triangles and record the label of the two sides.
//  for (int i = 0; i < activeEs.size(); ++i) {
//    int ei = activeEs[i];
//    const auto& oneP = sufPs[i];

//    // Split.
//    vector<int> mat = {label[Es[ei].v[0]], label[Es[ei].v[1]]};
//    for (int j = 2; j < oneP.size(); ++j) {
//      sufFs.push_back({oneP[0], oneP[j - 1], oneP[j]});
//      sufFMats.push_back(mat);
//    }
//  }
//}

void TetMesh::GetCurveInfoForArrangement(vector<vector<int>>& activeFs,
                                         vector<vector<float>>& crvVs,
                                         vector<vector<int>>& crvEs) const {
    crvVs = _bdCurveVs;
    crvEs = _bdCurveEs;
    activeFs.resize(_bdActiveFs.size(), vector<int>(3, 0));
    for (int i = 0; i < activeFs.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            activeFs[i][j] = Fs[_bdActiveFs[i]].v[j];
        }
    }
}
void TetMesh::WriteLevelSet(const vector<vector<float>>& sufVs,
                            const vector<vector<int>>& sufFs,
                            const vector<vector<int>>& sufFMats,
                            const char* filename) const {
    cout<<"out suf: "<<filename<<' '<< sufVs.size() << " " << sufFs.size() << "\n";
    ofstream ofs(filename, ofstream::out);
    ofs << sufVs.size() << " " << sufFs.size() << "\n";
    // SufVs.
    for (const auto& oneV : sufVs) {
        ofs << oneV[0] << " " << oneV[1] << " " << oneV[2] << "\n";
    }
    // SufFs and SufFMats.
    for (int fi = 0; fi < sufFs.size(); ++fi) {
        ofs << sufFs[fi][0] << " " << sufFs[fi][1] << " " << sufFs[fi][2] << " "
                            << sufFMats[fi][0] << " " << sufFMats[fi][1] << "\n";
    }
    // Contour edges.
    ofs << _bdCurveEs.size() << "\n";
    for (const auto& oneE : _bdCurveEs) {
        ofs << oneE[0] << " " << oneE[1] << "\n";
    }
    ofs.close();
}

void TetMesh::WriteLevelSet(const vector<vector<float>>& sufVs,
                            const vector<vector<int>>& sufFs,
                            const vector<vector<int>>& sufFMats,
                            const vector<int>& MappingMat,
                            const char* filename) const {
    cout<<"out suf: "<<filename<<' '<< sufVs.size() << " " << sufFs.size() << "\n";
    ofstream ofs(filename, ofstream::out);
    ofs << sufVs.size() << " " << sufFs.size() << "\n";
    // SufVs.
    for (const auto& oneV : sufVs) {
        ofs << oneV[0] << " " << oneV[1] << " " << oneV[2] << "\n";
    }
    // SufFs and SufFMats.
    for (int fi = 0; fi < sufFs.size(); ++fi) {
        ofs << sufFs[fi][0] << " " << sufFs[fi][1] << " " << sufFs[fi][2] << " "
                            << MappingMat[sufFMats[fi][0]] << " " << MappingMat[sufFMats[fi][1]] << "\n";
    }
    // Contour edges.
    ofs << _bdCurveEs.size() << "\n";
    for (const auto& oneE : _bdCurveEs) {
        ofs << oneE[0] << " " << oneE[1] << "\n";
    }
    ofs.close();
}
void TetMesh::OutputCellBoundary(vector<double>&out_Vs, vector<uint>&out_Fs, vector<int>&Vori_Ind,vector<int>&outVmarkers){
    OutputCellBoundary(out_Vs,out_Fs,Vori_Ind);
    outVmarkers.resize(Vori_Ind.size());
    for(int i=0;i<Vori_Ind.size();++i)outVmarkers[i] = VMarkers[Vori_Ind[i]];

}
void TetMesh::OutputCellBoundaryOriInd(vector<int>&Vori_Ind){
    vector<double>out_Vs;vector<uint>out_Fs;
    OutputCellBoundary(out_Vs,out_Fs,Vori_Ind);
}
void TetMesh::OutputCellBoundary(vector<double>&out_Vs,vector<uint>&out_Fs,vector<int>&Vori_Ind){


    vector<int>isbV(Vs.size(),-1);
    int ind = 0;
    Vori_Ind.clear();
    out_Vs.clear();
    out_Fs.clear();
    //Vs[i].isBD  Fs[i].isBD
    for(int i=0;i<Vs.size();++i)if(Vs[i].isBD){
        isbV[i] = ind++;
        Vori_Ind.push_back(i);
        auto &pv = VCoords[i];
        for(int j=0;j<3;++j)out_Vs.push_back((pv[j]));
    }

    for(int i=0;i<Fs.size();++i)if(Fs[i].isBD){

        auto pv = Fs[i].v;
        for(int j=0;j<3;++j)out_Fs.push_back((isbV[pv[j]]));
    }

//    for(int i=0;i<Vs.size();++i){
//        isbV[i] = ind++;
//        Vori_Ind.push_back(i);
//        auto &pv = VCoords[i];
//        for(int j=0;j<3;++j)out_Vs.push_back((pv[j]));
//    }

//    for(int i=0;i<Fs.size();++i){

//        auto pv = Fs[i].v;
//        for(int j=0;j<3;++j)out_Fs.push_back((isbV[pv[j]]));
//    }
    for(auto a:out_Fs)assert(a<ind);
    //for(auto a:out_Fs)cout<<a<<' ';cout<<endl;




}

void TetMesh::CutTet(vector<double>&v2planeDist,vector<double>&outVpos,vector<uint>&outF2V){


    outVpos.clear();outF2V.clear();
//    outVpos.push_back(0);outVpos.push_back(0);outVpos.push_back(0);
//    outVpos.push_back(0);outVpos.push_back(1);outVpos.push_back(0);
//    outVpos.push_back(0);outVpos.push_back(0);outVpos.push_back(1);

//    outF2V.push_back(0);outF2V.push_back(1);outF2V.push_back(2);

    assert(v2planeDist.size()==Vs.size());
    int nV = v2planeDist.size();
    int nE = Es.size();
    int nF = Fs.size();
    int nT = Ts.size();


    vector<bool>VonOffPlane(nV);
    for(int i=0;i<nV;++i)VonOffPlane[i] = v2planeDist[i]>0;

    vector<int>isActiveE(nE,-1);
    vector<int>ActiveE;
    vector<bool>isActiveF(nF,false);
    vector<int>ActiveF;
    outVpos.clear();
    int vind = 0;
    double vpos[3];

    for(int i = 0;i<nE;++i)if(VonOffPlane[Es[i].v[0]]!=VonOffPlane[Es[i].v[1]]){
        isActiveE[i] = vind++;ActiveE.push_back(i);
        double w1 = fabs(v2planeDist[Es[i].v[0]]), w2 = fabs(v2planeDist[Es[i].v[1]]),ww = w1+w2;
        w1 /= ww;w2/=ww;
        auto &pv1 = VCoords[Es[i].v[0]],&pv2 =  VCoords[Es[i].v[1]];
        for(int j =0;j<3;++j)vpos[j] = w2*pv1[j] + w1*pv2[j];
        for(int j =0;j<3;++j)outVpos.push_back(vpos[j]);
    }


    vector<bool>EfaceTT(nF,false);
    for(int i = 0;i<nT;++i){
        auto &T = Ts[i];
        vector<int>aE;
        for(int j=0;j<6;++j)if(isActiveE[T.edge[j]]!=-1){aE.push_back(T.edge[j]);}
        if(aE.size()==0)continue;
        assert(aE.size()==3 || aE.size() ==4 );

        if(aE.size()==3){
            for(int j =0;j<3;++j)outF2V.push_back(isActiveE[aE[j]]);
        }else{

            for(auto a: Es[aE[0]].nbF)EfaceTT[a] = true;
            int diagonalE = -100;
            for(int j =1;j<4;++j){
                bool isdiagonalE = true;
                for(auto a: Es[aE[j]].nbF)if(EfaceTT[a]){isdiagonalE = false;break;}
                if(isdiagonalE){diagonalE = j;break;}
            }

            for(int j =1;j<4;++j)if(j!=diagonalE){
                outF2V.push_back(isActiveE[aE[0]]);
                outF2V.push_back(isActiveE[aE[diagonalE]]);
                outF2V.push_back(isActiveE[aE[j]]);
            }

            for(auto a: Es[aE[0]].nbF)EfaceTT[a] = false;



        }



    }









}
