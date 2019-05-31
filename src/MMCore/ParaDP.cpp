
#include "ParallelArrangement.h"
#include <Eigen/Dense>
#include "Util.h"
#include "UnionFind.h"
#include "tetgen.h"
#include <set>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iostream>
#include <string>

#include "../Utility/readers.h"
#include "../Utility/geo_sur.h"
static const int flagAdd = 100000;

void ParallelArrangement::FindLoopDP(int stepi, vector<int> &pickSeg,vector<vector<int>>&reloop){
    assert(stepi>0);
    assert(stepi<_nCell);

    reloop.clear();
    if(pickSeg.size()==0){
        reloop.push_back(vector<int>(0));
        return;
    }

    int nSeg = newCtrNet.mixCurve.nSeg;
    int nVer = newCtrNet.mixCurve.keyVer.size();
    auto& seg2Cells = newCtrNet.mixCurve.seg2Cells;

    auto& seg2keyV = newCtrNet.mixCurve.superedges;

    vector<int>preSetKey(nVer,0);
    vector<vector<int>>keyV2Seg(nVer);

    for(auto b:pickSeg)if(seg2keyV[b*2]!=-1){
        preSetKey[seg2keyV[b*2]]++;
        preSetKey[seg2keyV[b*2+1]]++;
        keyV2Seg[seg2keyV[b*2]].push_back(b);
        keyV2Seg[seg2keyV[b*2+1]].push_back(b);
    }



    vector<int>curSeg2Cell(nSeg,-1);
    vector<bool>ismergeC(_nCell,false);
    auto &mergeorder = newCtrNet.mergingorder;
    for(int i=0;i<=stepi;++i)ismergeC[mergeorder[i]]  = true;

    for(auto a:pickSeg){
        int c1 = seg2Cells[a*2],c2 = seg2Cells[a*2+1];
        assert(ismergeC[c1]!=ismergeC[c2]);
        if(ismergeC[c1])curSeg2Cell[a] = c1;
        else curSeg2Cell[a] = c2;
    }

    UnionFind unifind;
    unifind.SetElements(pickSeg);

    for(auto &a:keyV2Seg){
        assert(a.size()!=1);
        if(a.size()==2)unifind.Union(a[0],a[1]);
        else if(a.size()>2){
            vector<int>fff(a.size(),-1);
            int ind = 0;
            for(int i=0;i<a.size();++i)if(fff[i]==-1){
                int cc = curSeg2Cell[a[i]];
                fff[i] = cc;
                for(int j=i+1;j<a.size();++j)if(fff[j]==-1)if(curSeg2Cell[a[i]]==cc){
                    unifind.Union(a[i],a[j]);
                    fff[i] = cc;
                }
            }

        }

    }

    unifind.ExtractComponents(reloop);






}

int ParallelArrangement::CheckComponentNum(vector<vector<int>>&pickSet,vector<vector<int>>&preSetSeg){

    int nSeg = newCtrNet.mixCurve.nSeg;
    auto &SegNSeg = newCtrNet.mixCurve.SegNSeg;
    auto &superedges = newCtrNet.mixCurve.superedges;
    for(auto &a:pickSet)for(auto b:a)assert(b<nSeg);
    if(pickSet.size()==0)return 0;

    //    auto iscontain = [](int a,vector<int>&vec){
    //        for(auto b:vec)if(a==b)return true;
    //        return false;
    //    };
    //    UnionFind unifind;
    //    unifind.SetElements(pickSet);

    //    for(int i=0;i<pickSet.size();++i)for(int j=i+1;j<pickSet.size();++j)
    //        if(iscontain(pickSet[j],SegNSeg[pickSet[i]]))unifind.Union(pickSet[i],pickSet[j]);
    //    return unifind.GetNumOfComponents();
    int VEF = 0;
    vector<bool>ispreSetKey(newCtrNet.mixCurve.keyVer.size(),false);
    vector<bool>isCrossKey(newCtrNet.mixCurve.keyVer.size(),false);
    vector<bool>ispreSetSeg(nSeg,false);
    vector<int>preSetKey(newCtrNet.mixCurve.keyVer.size(),0);

    for(auto &a:preSetSeg){
        for(auto b:a){
            ispreSetSeg[b] = true;
            if(superedges[b*2]!=-1){
                preSetKey[superedges[b*2]]++;
                preSetKey[superedges[b*2+1]]++;
                ispreSetKey[superedges[b*2]] = true;
                ispreSetKey[superedges[b*2+1]] = true;
            }
        }
    }
    for(int i=0;i<preSetKey.size();++i)if(preSetKey[i]>2)isCrossKey[i] = true;

    for(auto &a:pickSet){
        vector<bool>ispreSetKeytmp(newCtrNet.mixCurve.keyVer.size(),false);
        for(auto b:a){
            if(ispreSetSeg[b])++VEF;
            if(superedges[b*2]!=-1){
                ispreSetKeytmp[superedges[b*2]]=true;
                ispreSetKeytmp[superedges[b*2+1]]=true;
            }
        }
        for(int i=0;i<newCtrNet.mixCurve.keyVer.size();++i){
            if(ispreSetKeytmp[i] && ispreSetKey[i])--VEF;
        }
        for(int i=0;i<newCtrNet.mixCurve.keyVer.size();++i){
            if(ispreSetKeytmp[i] && isCrossKey[i])--VEF;
        }
    }

    //    for(auto &a:preSetSeg){
    //        vector<bool>ispreSetKeytmp(newCtrNet.mixCurve.keyVer.size(),false);
    //        for(auto b:a){
    //            if(ispreSetSeg[b])++VEF;
    //            if(superedges[b*2]!=-1){
    //                ispreSetKeytmp[superedges[b*2]]=true;
    //                ispreSetKeytmp[superedges[b*2+1]]=true;
    //            }
    //        }
    //        for(int i=0;i<newCtrNet.mixCurve.keyVer.size();++i){
    //            if(ispreSetKeytmp[i] && ispreSetKey[i])--VEF;
    //        }
    //        for(auto b:a){
    //            ispreSetSeg[b] = true;
    //            if(superedges[b*2]!=-1){
    //                ispreSetKey[superedges[b*2]]=true;
    //                ispreSetKey[superedges[b*2+1]]=true;
    //            }
    //        }
    //    }
    //    for(auto &a:pickSet){
    //        vector<bool>ispreSetKeytmp(newCtrNet.mixCurve.keyVer.size(),false);
    //        for(auto b:a){
    //            if(ispreSetSeg[b])++VEF;
    //            if(superedges[b*2]!=-1){
    //                ispreSetKeytmp[superedges[b*2]]=true;
    //                ispreSetKeytmp[superedges[b*2+1]]=true;
    //            }
    //        }
    //        for(int i=0;i<newCtrNet.mixCurve.keyVer.size();++i){
    //            if(ispreSetKeytmp[i] && ispreSetKey[i])--VEF;
    //        }
    //    }


    return VEF;




}
void ParallelArrangement::LoopMergingDP(){



    //vector< vector< vector<int> > >loop2Segs;
    //vector<uint> mergingorder;

    //vector<vector<bool>>existSeg;
    //vector<vector<uint>>disSeg;
    //vector<vector<bool>>mergeSeg;
    //vector<vector<uint>>dismergeSeg;

    //_label2bdLoops
    auto iscontain = [](int a,vector<int>&vec){
        for(auto b:vec)if(a==b)return true;
        return false;
    };

    mergeLoop2Seg.clear();
    newloop2mergeLoop.clear();
    oldloop2mergeLoop.clear();
    label2mergeLoop.clear();
    label2mergeloopComp.clear();
    label2mergeLoopCompVE.clear();
    mergeLoop2label.clear();
    mergeLoop2VE.clear();
    mergeloopVoid.clear();

    int nSeg = newCtrNet.mixCurve.nSeg;
    mergeLoop2Seg.resize(_nCell);
    newloop2mergeLoop.resize(_nCell);
    oldloop2mergeLoop.resize(_nCell);
    label2mergeLoop.resize(_nCell);
    label2mergeloopComp.resize(_nCell);
    label2mergeLoopCompVE.resize(_nCell);
    mergeLoop2label.resize(_nCell);
    mergeLoop2VE.resize(_nCell);
    mergeloopVoid.resize(_nCell);

    maxExistStep.clear();
    maxExistStep.resize(_nMat,-1);

    for(int i=0;i<_nCell;++i){
        cout<<i<<endl;
        if(i== 36){
            cout<<i<<endl;
        }
        int curCell = newCtrNet.mergingorder[i];
        auto &curlabel2mergeLoop = label2mergeLoop[i];
        auto &curmergeLoop2Seg = mergeLoop2Seg[i];
        auto &MapNewloop2MergeLoop = newloop2mergeLoop[i];
        auto &MapOldloop2MergeLoop = oldloop2mergeLoop[i];
        auto &mergeLoop2VE_curC = mergeLoop2VE[i];
        if(i==0){
            mergeLoop2Seg[i] = loop2Segs[curCell];
            label2mergeLoop[i] = _label2bdLoops[curCell];
            mergeLoop2label[i] = _bdLoops2lable[curCell];
            MapNewloop2MergeLoop.resize(mergeLoop2Seg[i].size());
            for(int k = 0;k<mergeLoop2Seg[i].size();++k)MapNewloop2MergeLoop[k].push_back(k);
            MapOldloop2MergeLoop.clear();
            mergeLoop2VE_curC.resize(mergeLoop2Seg[i].size(),0);

            mergeloopVoid[i].resize(mergeLoop2Seg[i].size(),false);
            continue;
        }

        auto &bemergingLoops2Seg = loop2Segs[curCell];
        auto &prelabel2mergeLoop = label2mergeLoop[i-1];
        auto &premergeLoop2Seg = mergeLoop2Seg[i-1];
        auto &label2mergeloopComp_curC = label2mergeloopComp[i];
        auto &label2mergeLoopCompVE_curC = label2mergeLoopCompVE[i];
        auto &mergeLoop2label_curC = mergeLoop2label[i];
        curmergeLoop2Seg.clear();
        mergeLoop2VE_curC.clear();
        curlabel2mergeLoop.clear();
        curlabel2mergeLoop.resize(_nMat);

        label2mergeloopComp_curC.resize(_nMat);
        label2mergeLoopCompVE_curC.resize(_nMat);
        MapNewloop2MergeLoop.resize(bemergingLoops2Seg.size(),vector<int>(0));
        MapOldloop2MergeLoop.resize(premergeLoop2Seg.size(),vector<int>(0));

        if(bemergingLoops2Seg.size()==0){

            mergeLoop2Seg[i] = mergeLoop2Seg[i-1];
            newloop2mergeLoop[i] = newloop2mergeLoop[i-1];
            oldloop2mergeLoop[i] = oldloop2mergeLoop[i-1];
            label2mergeLoop[i] = label2mergeLoop[i-1];
            label2mergeloopComp[i] = label2mergeloopComp[i-1];
            label2mergeLoopCompVE[i] = label2mergeLoopCompVE[i-1];
            mergeLoop2label[i] = mergeLoop2label[i-1];
            mergeLoop2VE[i] = mergeLoop2VE[i-1];
            mergeloopVoid[i] = mergeloopVoid[i-1];
            continue;

        }


        for(int j=0;j<_nMat;++j){
            curlabel2mergeLoop[j].clear();
            auto& curMatbemergingloop = _label2bdLoops[curCell][j];
            if(curMatbemergingloop.size()==0){
                bool isFin = true;
                for(auto a:prelabel2mergeLoop[j])if(premergeLoop2Seg[a].size()!=0){
                    isFin = false;break;
                }
                //                for(auto a:prelabel2mergeLoop[j])if(premergeLoop2Seg[a].size()!=0){
                //                    MapOldloop2MergeLoop[a].push_back(curmergeLoop2Seg.size());
                //                    curlabel2mergeLoop[j].push_back(curmergeLoop2Seg.size());
                //                    curmergeLoop2Seg.push_back(premergeLoop2Seg[a]);
                //                    mergeLoop2label_curC.push_back(j);
                //                    mergeLoop2VE_curC.push_back(0);
                //                }
                if(isFin)continue;
            }

            maxExistStep[j] = i;
            //set<pair<int,int>>mmloop;
            UnionFind unifind;
            vector<int>element;
            for(int kk = 0;kk<prelabel2mergeLoop[j].size();++kk)if(premergeLoop2Seg[prelabel2mergeLoop[j][kk]].size()!=0)element.push_back(prelabel2mergeLoop[j][kk]);
            for(int kk = 0;kk<curMatbemergingloop.size();++kk)element.push_back(flagAdd+curMatbemergingloop[kk]);
            unifind.SetElements(element);
            vector<bool>prepool(prelabel2mergeLoop[j].size(),false);
            vector<bool>curpool(curMatbemergingloop.size(),false);
            for(auto a:newCtrNet.dismergeSeg[i]){
                vector<int>prePick,curPick;
                for(int k=0;k<prelabel2mergeLoop[j].size();++k)if(iscontain(a,premergeLoop2Seg[prelabel2mergeLoop[j][k]])){
                    prePick.push_back(prelabel2mergeLoop[j][k]);prepool[k]=true;
                }
                for(int k=0;k<curMatbemergingloop.size();++k)if(iscontain(a,bemergingLoops2Seg[curMatbemergingloop[k]])){
                    curPick.push_back(curMatbemergingloop[k]);curpool[k]=true;
                }
                //                for(auto b:prelabel2mergeLoop[j])if(iscontain(a,premergeLoop2Seg[b])){
                //                    prePick.push_back(b);
                //                }
                //                for(auto b:curlabel2mergeLoop[j])if(iscontain(a,curmergeLoop2Seg[b])){
                //                    curPick.push_back(b);
                //                }
                assert(prePick.size()<=1);
                assert(curPick.size()<=1);
                assert(prePick.size()==curPick.size());

                if(prePick.size()==1)unifind.Union(prePick[0],flagAdd+curPick[0]);

            }
            vector<vector<int>> &comp = label2mergeloopComp_curC[j];
            unifind.ExtractComponents(comp);

            label2mergeLoopCompVE_curC[j].clear();
            for(auto&oneConnect:comp){
                int n_bemerging = 0,n_premerge = 0,n_exsit = 0,n_mergeSeg = 0;
                int VEF=0;
                vector<bool>ploop(nSeg,false);
                vector<vector<int>>mergeSegSet;vector<vector<int>>involveSeg;
                for(auto&a:oneConnect)if(a<flagAdd){
                    involveSeg.push_back(premergeLoop2Seg[a]);
                    for(auto b:premergeLoop2Seg[a]){
                        ploop[b]=true;
                    }
                    ++n_premerge;
                }
                mergeSegSet.resize(1);mergeSegSet[0].clear();
                for(auto&a:oneConnect)if(a>flagAdd-1){
                    //mergeSegSet.push_back(bemergingLoops2Seg[a-flagAdd]);
                    for(auto b:bemergingLoops2Seg[a-flagAdd]){
                        //involveSeg.push_back(bemergingLoops2Seg[b]);
                        if(ploop[b]){ploop[b]=false;++n_mergeSeg;mergeSegSet[0].push_back(b);}
                        else ploop[b] = true;
                    }
                    ++n_bemerging;
                }
                vector<int>newloop;


                for(int k=0;k<nSeg;++k)if(ploop[k]){newloop.push_back(k);++n_exsit;}
                VEF = CheckComponentNum(mergeSegSet,involveSeg);
                //VEF=-n_bemerging*n_premerge;
                if(newloop.size()==0){
                    //label2mergeLoopCompVE_curC[j].push_back(0);
                    //continue;
                    VEF=0;
                }
                vector<vector<int>>splitloop;
                FindLoopDP(i,newloop,splitloop);
                int kda = true;
                for(auto &b:splitloop){
                    sort(b.begin(),b.end());
                    for(auto&a:oneConnect)if(a<flagAdd){MapOldloop2MergeLoop[a].push_back(curmergeLoop2Seg.size());}
                    for(auto&a:oneConnect)if(a>flagAdd-1){MapNewloop2MergeLoop[a-flagAdd].push_back(curmergeLoop2Seg.size());}
                    curlabel2mergeLoop[j].push_back(curmergeLoop2Seg.size());
                    curmergeLoop2Seg.push_back(b);
                    //assert(min(n_bemerging,n_premerge)<=1);
                    if(kda){label2mergeLoopCompVE_curC[j].push_back(VEF);mergeLoop2VE_curC.push_back(VEF);kda = false;}
                    else {label2mergeLoopCompVE_curC[j].push_back(0);mergeLoop2VE_curC.push_back(0);}
                    mergeLoop2label_curC.push_back(j);
                }



            }


        }

        vector<int>cnum(nSeg,0);
        for(auto &a:mergeLoop2Seg[i])for(auto b:a)++cnum[b];

        for(int k = 0;k<nSeg;++k){
            if(newCtrNet.existSeg[i][k])assert(cnum[k]==2);
            else assert(cnum[k]==0);
        }
        mergeloopVoid[i].clear();
        for(auto &a:mergeLoop2Seg[i])if(a.size()==0)mergeloopVoid[i].push_back(true);else mergeloopVoid[i].push_back(false);
        for(int k=0;k<mergeloopVoid[i].size();++k)if(mergeloopVoid[i][k])assert(mergeLoop2VE_curC[k]==0);
        if(i==_nCell-1)for(auto a: mergeloopVoid[i])assert(a);

    }

    cout<<"Mergeloop DP done!"<<endl;

}
void ParallelArrangement::Comp2label(vector<vector<int>>&comps,vector<int>&loop2label, vector<int>&relabel){

    relabel.resize(comps.size());
    for(int i=0;i<comps.size();++i){
        int labb = -1;
        for(auto a: comps[i])
            if(labb==-1){
                labb = relabel[i] = loop2label[a];
            }else{
                assert(labb == loop2label[a]);
            }
    }
}

bool ParallelArrangement::CheckValidation(int stepi,CellTopology &curCell){

    vector<int>comp2label(curCell._comps.size());
    Comp2label(curCell._comps,mergeLoop2label[stepi],comp2label);

    //    vector<int> comp2VEF;
    //    vector<int> comp2label;
    //    vector<int> isoComps2VEF;
    //    vector<int> isoComps2label;
    if(1){

        vector<int>comp2ge(curCell._comps.size());
        vector<vector<int>>label2openSurfGenus(_nMat);
        vector<vector<int>>label2closeSurfGenus(_nMat);
        for(int i=0;i<curCell._comps.size();++i){
            int ge = (2-curCell._comps[i].size()-curCell.comp2VEF[i]);
            assert(ge%2==0);
            assert(ge>=0);
            comp2ge[i] = ge/2;
            label2openSurfGenus[comp2label[i]].push_back(comp2ge[i]);
        }

        vector<int>isoComp2ge(curCell.isoComps2label.size());
        for(int i=0;i<curCell.isoComps2label.size();++i){
            int ge = (2-curCell.isoComps2VEF[i]);
            assert(ge%2==0);
            assert(ge>=0);
            isoComp2ge[i] = ge/2;
            label2closeSurfGenus[curCell.isoComps2label[i]].push_back(isoComp2ge[i]);
        }

        // first test: number of components


        for(int i=0;i<_nMat;++i){
            if(DP_label2NComp[i]==-1)continue;
            int nc = label2closeSurfGenus[i].size();
            if(nc>DP_label2NComp[i])return false;
            if(maxExistStep[i]>stepi && nc+1>DP_label2NComp[i])return false;
            if(stepi==_nCell-1)if(nc!=DP_label2NComp[i])return false;
        }



        // second test: combination possibility

        for(int i=0;i<_nMat;++i){
            if(DP_label2NComp[i]==-1)continue;
            auto &tge = DP_label2CGenus[i];
            vector<bool>mask(tge.size(),false);
            int curMax = -1,restSum = 0,openSum = 0;
            for(auto a:label2closeSurfGenus[i]){
                int j=0;
                for(;j<mask.size();++j)if(!mask[j])if(a == tge[j] || tge[j] > 100000){
                    mask[j] = true;break;
                }
                if(j==mask.size())return false;
            }
            for(int j=0;j<mask.size();++j)if(!mask[j]){curMax = max(curMax,tge[j]);restSum+=tge[j];}
            if(curMax>100000)continue;
            for(auto a:label2openSurfGenus[i]){
                //check max and accumulate
                if(a>curMax)return false;
                openSum+=a;
            }
            if(openSum>restSum)return false;
            //check max sum

        }

        //        for(int i=0;i<curCell.isoComps2label.size();++i){
        //            for(auto a:curCell.isoComps2allloops[i])if(DP_loops2groupping[a]!=-1){
        //                int ngroup = DP_loops2groupping[a];
        //                if(DP_grouppinggenus[ngroup ] <1000000 && DP_grouppinggenus[ngroup ]!=isoComp2ge[i])return false;
        //                vector<bool>lflag(_nLoops,false);
        //                for(auto b:curCell.isoComps2allloops[i])lflag[b] = true;
        //                for(auto b:DP_groupping[ngroup])if(!lflag[b])return false;
        //                break;
        //            }

        //        }

        return true;
    }

    vector<int>goalGenus(_nMat,0);//goalGenus[2]=1;//goalGenus[0]=1;
    //goalGenus[2]=1;goalGenus[3]=1;goalGenus[4]=3;
    goalGenus[2]=1;
    //goalGenus[2]=1;goalGenus[3]=6;
    //goalGenus[0]=0;


    if(stepi!=_nCell-1){
        //cout<<"stepi: "<<stepi<<endl;
        //if(curCell.isContainIso)return false;
        for(int i=0;i<curCell._comps.size();++i){
            int ge = (2-curCell._comps[i].size()-curCell.comp2VEF[i]);
            assert(ge%2==0);
            assert(ge>=0);
            int genus = ge/2;
            if(genus>goalGenus[comp2label[i]])return false;
        }
        if(curCell.isoComps2label.size()!=0){
            for(auto a : curCell.isoComps2label)if(maxExistStep[a]>stepi)return false;

            vector<int>numofSurpM(_nMat,0);
            for(auto a:curCell.isoComps2label)++numofSurpM[a];
            for(auto a:numofSurpM)if(a>1)return false;

            for(int i=0;i<curCell.isoComps2label.size();++i){
                int ge = (2-curCell.isoComps2VEF[i]);
                assert(ge%2==0);
                assert(ge>=0);
                int genus = ge/2;
                if(genus!=goalGenus[curCell.isoComps2label[i]])return false;
            }
        }
        return true;
    }else{

        //need change for multi material for early closing
        if(curCell.isoComps2label.size()!=_nMat){
            return false;
        }

        vector<int>numofSurpM(_nMat,0);
        for(auto a:curCell.isoComps2label)++numofSurpM[a];
        for(auto a:numofSurpM)if(a!=1)return false;


        for(int i=0;i<curCell.isoComps2label.size();++i){
            int ge = (2-curCell.isoComps2VEF[i]);
            assert(ge%2==0);
            assert(ge>=0);
            int genus = ge/2;
            if(genus!=goalGenus[curCell.isoComps2label[i]])return false;
        }

        return true;
    }


}
void ParallelArrangement::PrintMergeTopoInfo(int stepi,CellTopology &curCell){
    assert(stepi>=0);
    assert(stepi<_nCell);

    vector<int>TopoComp2label(curCell._comps.size());
    Comp2label(curCell._comps,mergeLoop2label[stepi],TopoComp2label);

    cout<<"step: "<<stepi<<endl;
    for(int i=0;i<curCell._comps.size();++i){
        cout<<i<<": ";
        for(auto a:curCell._comps[i])cout<<a<<' ';
        cout<<"lable: "<<TopoComp2label[i]<<endl;
    }
    cout<<"Iso Comp: "<<endl;
    for(auto a:curCell.isoComps2VEF)cout<<(2-a)/2<<' ';cout<<endl;
    for(auto a:curCell.isoComps2label)cout<<a<<' ';cout<<endl;
    cout<<"Junction Points: "<<curCell._junctionpoints<<endl;



}
void ParallelArrangement::PrintBeMergeTopoInfo(int stepi,CellTopology &curCell){
    assert(stepi>=0);
    assert(stepi<_nCell);

    int curCelli = newCtrNet.mergingorder[stepi];
    auto &_bdLoops2lable_curC = _bdLoops2lable[curCelli];
    vector<int>TopoComp2label(curCell._comps.size());
    Comp2label(curCell._comps,_bdLoops2lable_curC,TopoComp2label);

    cout<<"step: "<<stepi<<endl;
    for(int i=0;i<curCell._comps.size();++i){
        cout<<i<<": ";
        for(auto a:curCell._comps[i])cout<<a<<' ';
        cout<<"lable: "<<TopoComp2label[i]<<endl;
    }
    cout<<"Iso Comp: "<<endl;
    for(auto a:curCell.isoComps2VEF)cout<<(2-a)/2<<' ';cout<<endl;
    for(auto a:curCell.isoComps2label)cout<<a<<' ';cout<<endl;
    cout<<"Junction Points: "<<curCell._junctionpoints<<endl;



}

void ParallelArrangement::MapComponentLoops2CompactLoop(const int celli, const vector<vector<int>>&comps, vector<vector<int>> &comps2compactloop){


    int nloopwithinCell = loop2compactloopN[celli].size();
    auto &loop2compactloopNc = loop2compactloopN[celli];
    comps2compactloop = comps;
    for(auto &a:comps2compactloop)for(auto &b:a){
        assert(b<nloopwithinCell);
        b = loop2compactloopNc[b];
    }


}
void ParallelArrangement::MapComponentLoops2CompactLoop(const int celli, CellTopology &TPwithinASingleCell){
    MapComponentLoops2CompactLoop(celli, TPwithinASingleCell._comps,  TPwithinASingleCell.comps2allloops);
}

int ParallelArrangement::MergeTwoCellTopo(int stepi,CellTopology &preCell,CellTopology &curCell,CellTopology &mergeCell){

    assert(stepi>=0);
    assert(stepi<_nCell);


    mergeCell._score = preCell._score + curCell._score;
    mergeCell._junctionpoints = preCell._junctionpoints + curCell._junctionpoints;
    mergeCell.comp2VEF.clear();mergeCell._comps.clear();
    mergeCell.comps2allloops.clear();
    mergeCell.isoComps2allloops.clear();



    int curCelli = newCtrNet.mergingorder[stepi];


    auto &label2mergeLoop_curC = label2mergeLoop[stepi];
    auto &label2mergeloopComp_curC = label2mergeloopComp[stepi];
    auto &label2mergeLoopCompVE_curC = label2mergeLoopCompVE[stepi];
    auto &mergeloopVoid_curC = mergeloopVoid[stepi];


    auto &mergeLoop2label_preC = mergeLoop2label[stepi-1];
    auto &_bdLoops2lable_curC = _bdLoops2lable[curCelli];
    vector<vector<int>>tmpComp;

    vector<int>preTopoComp2label(preCell._comps.size());
    Comp2label(preCell._comps,mergeLoop2label_preC,preTopoComp2label);
    //    for(int i=0;i<preTopoComp2label.size();++i){
    //        int labb = -1;
    //        for(auto a: preCell._comps[i])
    //            if(labb==-1){
    //                labb = preTopoComp2label[i] = mergeLoop2label_preC[a];
    //            }else{
    //                assert(labb == mergeLoop2label_preC[a]);
    //            }
    //    }

    vector<int>curTopoComp2label(curCell._comps.size());
    Comp2label(curCell._comps,_bdLoops2lable_curC,curTopoComp2label);
    //    for(int i=0;i<curTopoComp2label.size();++i){
    //        int labb = -1;
    //        for(auto a: curCell._comps[i])
    //            if(labb==-1){
    //                labb = curTopoComp2label[i] = _bdLoops2lable_curC[a];
    //            }else{
    //                assert(labb == _bdLoops2lable_curC[a]);
    //            }
    //    }

    vector<vector<int>>curTopolabel2Comp(_nMat);
    for(int i=0;i<curTopoComp2label.size();++i)curTopolabel2Comp[curTopoComp2label[i]].push_back(i);


    auto &MapNewloop2MergeLoop = newloop2mergeLoop[stepi];
    auto &MapOldloop2MergeLoop = oldloop2mergeLoop[stepi];
    auto &mergeLoop2VE_curC = mergeLoop2VE[stepi];
    vector<set<int>>associatePreComp(mergeLoop2VE_curC.size());
    vector<set<int>>associateCurComp(mergeLoop2VE_curC.size());
    for(int i=0;i<_nMat;++i){

        UnionFind unifind;
        unifind.SetElements(label2mergeLoop_curC[i]);

        for(int j=0;j<preTopoComp2label.size();++j)if(preTopoComp2label[j]==i){
            auto &curcomp = preCell._comps[j];
            set<int>c2ml;
            for(auto a:curcomp)for(auto b:MapOldloop2MergeLoop[a]){
                c2ml.insert(b);
                associatePreComp[b].insert(j);
            }
            if(c2ml.size()>1){
                vector<int>tmpv;for(auto b:c2ml)tmpv.push_back(b);
                for(int k=1;k<tmpv.size();++k)unifind.Union(tmpv[k-1],tmpv[k]);
            }
        }

        for(int j=0;j<curTopoComp2label.size();++j)if(curTopoComp2label[j]==i){
            auto &curcomp = curCell._comps[j];
            set<int>c2ml;
            for(auto a:curcomp)for(auto b:MapNewloop2MergeLoop[a]){
                c2ml.insert(b);
                associateCurComp[b].insert(j);
            }
            if(c2ml.size()>1){
                vector<int>tmpv;for(auto b:c2ml)tmpv.push_back(b);
                for(int k=1;k<tmpv.size();++k)unifind.Union(tmpv[k-1],tmpv[k]);
            }
        }

        vector<vector<int>>unifindcomp;
        unifind.ExtractComponents(unifindcomp);
        for(auto &a:unifindcomp)sort(a.begin(),a.end());


        //        for(auto &a:unifindcomp)mergeCell._comps.push_back(a);
        for(auto &a:unifindcomp)tmpComp.push_back(a);



    }

    vector<int>nloop(mergeLoop2VE_curC.size(),0);
    for(auto &a:tmpComp)for(auto b:a)++nloop[b];
    for(auto a:nloop)assert(a==1);

    vector<int>tmpVEF;

    vector<vector<int>>tmp_comp2Allloop;




    tmpVEF.clear();
    tmpVEF.resize(tmpComp.size(),0);
    tmp_comp2Allloop.clear();
    tmp_comp2Allloop.resize(tmpComp.size());
    for(int i=0;i<tmpComp.size();++i){
        set<int>preComp,curComp;
        for(auto a:tmpComp[i])for(auto b:associatePreComp[a])preComp.insert(b);
        for(auto a:tmpComp[i])for(auto b:associateCurComp[a])curComp.insert(b);

        int VEF=0;
        for(auto a:preComp)VEF+=preCell.comp2VEF[a];
        for(auto a:curComp)VEF+=curCell.comp2VEF[a];
        for(auto a:tmpComp[i])VEF+=mergeLoop2VE_curC[a];
        tmpVEF[i] = VEF;

        set<int>mergeallloops;
        //vector<int>&mergeallloops = tmp_comp2Allloop[i];
        for(auto a:preComp)for(auto b:preCell.comps2allloops[a])mergeallloops.insert(b);
        for(auto a:curComp)for(auto b:curCell.comps2allloops[a])mergeallloops.insert(b);

        for(auto a:mergeallloops)tmp_comp2Allloop[i].push_back(a);
    }


    mergeCell.isContainIso = false;
    mergeCell.isoComps2label = preCell.isoComps2label;
    mergeCell.isoComps2VEF = preCell.isoComps2VEF;
    mergeCell.isoComps2allloops = preCell.isoComps2allloops;
    for(int i=0;i<tmpComp.size();++i){
        vector<int>tmpVec;
        auto &a = tmpComp[i];
        for(int j=0;j<a.size();++j)if(!mergeloopVoid_curC[a[j]]){tmpVec.push_back(a[j]);}
        if(tmpVec.size()==0){
            mergeCell.isContainIso = true;
            mergeCell.isoComps2label.push_back(mergeLoop2label[stepi][a[0]]);
            mergeCell.isoComps2VEF.push_back(tmpVEF[i]);

            mergeCell.isoComps2allloops.push_back(tmp_comp2Allloop[i]);
        }
        else {
            mergeCell._comps.push_back(tmpVec);
            mergeCell.comp2VEF.push_back(tmpVEF[i]);

            mergeCell.comps2allloops.push_back(tmp_comp2Allloop[i]);
        }
        //        if(stepi==_nCell-1){
        //            mergeCell._comps.push_back(a);
        //            mergeCell.comp2VEF.push_back(tmpVEF[i]);
        //        }


    }

    vector<int> compIDs(mergeCell._comps.size());
    for(int i=0;i<compIDs.size();++i)compIDs[i] = i;
    auto comparator = [&mergeCell](int a, int b) { return mergeCell._comps[a] < mergeCell._comps[b]; };
    sort(compIDs.begin(), compIDs.end(), comparator);
    sort(mergeCell._comps.begin(), mergeCell._comps.end());

    tmpVEF = mergeCell.comp2VEF;
    for(int i=0;i<mergeCell._comps.size();++i){
        mergeCell.comp2VEF[i] = tmpVEF[compIDs[i]];
    }

    tmp_comp2Allloop = mergeCell.comps2allloops;
    for(int i=0;i<mergeCell._comps.size();++i){
        mergeCell.comps2allloops[i] = tmp_comp2Allloop[compIDs[i]];
    }


    if(mergeCell.isoComps2label.size()>1){
        vector<int> isocompIDs(mergeCell.isoComps2label.size());
        for(int i=0;i<isocompIDs.size();++i)isocompIDs[i] = i;
        auto comparator = [&mergeCell](int a, int b) { return mergeCell.isoComps2VEF[a] < mergeCell.isoComps2VEF[b]; };
        sort(isocompIDs.begin(), isocompIDs.end(), comparator);
        sort(mergeCell.isoComps2VEF.begin(), mergeCell.isoComps2VEF.end());

        auto tmpIsoLabel = mergeCell.isoComps2label;
        for(int i=0;i<mergeCell.isoComps2VEF.size();++i){
            mergeCell.isoComps2label[i] = tmpIsoLabel[isocompIDs[i]];
        }

        auto tmpIsocomp2Allloop = mergeCell.isoComps2allloops;
        for(int i=0;i<mergeCell.isoComps2VEF.size();++i){
            mergeCell.isoComps2allloops[i] = tmpIsocomp2Allloop[isocompIDs[i]];
        }
    }



    mergeCell.mSeq.clear();
    for(auto a:preCell.mSeq)mergeCell.mSeq.push_back(a);
    for(auto a:curCell.mSeq)mergeCell.mSeq.push_back(a);


    return 0;




}


void ParallelArrangement::DynamicProgramming(vector<int>&pickTopo){

    //{0,1,-1,0,1,0,1,1,1,-1,0}
    //vector<int>pickTopo(_nCell,0);
    //for(int i=8;i<_nCell;++i)pickTopo[i]=1;

    //pickTopo.clear();
    //pickTopo.resize(_nCell,0);
    //for(int i=8;i<_nCell;++i)pickTopo[i]=1;
    //pickTopo[1]=1;
    vector<CellTopology>mergeTopo(_nCell+1);
    auto mergingorder = newCtrNet.mergingorder;
    mergeTopo[0]._score = 0;
    for(int i=0;i<_nCell;++i){
        if(i==35){
            //cout<<i<<endl;
        }
        auto newCellid = mergingorder[i];
        if(_cellTopo[newCellid].size()==0)
            mergeTopo[i+1] = mergeTopo[i];
        else{
            auto a = _cellTopo[newCellid][pickTopo[newCellid]];
            MapComponentLoops2CompactLoop(newCellid,a);
            MergeTwoCellTopo(i,mergeTopo[i],a,mergeTopo[i+1]);
        }

        PrintMergeTopoInfo(i,mergeTopo[i+1]);
        cout<<endl<<endl;
        cout<<"Be Merging Cell Topo"<<endl;
        if(_cellTopo[newCellid].size()!=0)PrintBeMergeTopoInfo(i,_cellTopo[newCellid][pickTopo[newCellid]]);
        else cout<<"Empty Cell"<<endl;
    }

    cout<<"DP Done!"<<endl;
    cout<<"Merging Order: ";for(auto a:mergingorder)cout<<a<<' ';cout<<endl;
    for(auto a:mergeTopo[_nCell].isoComps2VEF)cout<<(2-a)/2<<' ';cout<<endl;
    for(auto a:mergeTopo[_nCell].isoComps2label)cout<<a<<' ';cout<<endl;
    //for(int i=1;i<_nCell;++i)cout<<CheckValidation(i,mergeTopo[i+1])<<' ';cout<<endl;

    cout<<"Score: "<<mergeTopo[_nCell]._score<<endl;
    cout<<"Junction Points: "<<mergeTopo[_nCell]._junctionpoints<<endl;
    for(auto a:pickTopo)cout<<a<<' ';cout<<endl;

}


//static bool operator==(const CellTopology& lh, const CellTopology& rh) {
//    // return (lh._comps == rh._comps) && (lh._conns == rh._conns);
//    return (lh._comps == rh._comps) && (lh.comp2VEF == rh.comp2VEF);
//}
//static bool operator<(const CellTopology& lh, const CellTopology& rh) {
//    return (lh._score < rh._score );
//}
//static bool operator>(const CellTopology& lh, const CellTopology& rh) {
//    return (lh._score > rh._score );
//}
extern bool isCJP;
static bool operator==(const CellTopology& lh, const CellTopology& rh) {
    if(isCJP)return (lh._comps == rh._comps) && (lh.comp2VEF == rh.comp2VEF);
    return (lh._comps == rh._comps) && (lh.comp2VEF == rh.comp2VEF) && (lh._junctionpoints == rh._junctionpoints);
}
static bool operator<(const CellTopology& lh, const CellTopology& rh) {
    if(isCJP)if(lh._junctionpoints != rh._junctionpoints )return (lh._junctionpoints > rh._junctionpoints );
    return (lh._score < rh._score );
}
static bool operator>(const CellTopology& lh, const CellTopology& rh) {
    if(isCJP)if(lh._junctionpoints != rh._junctionpoints )return (lh._junctionpoints < rh._junctionpoints );
    return (lh._score > rh._score );
}
void ParallelArrangement::MergeSameTopo(vector<CellTopology>&topos,int stepi){

    int ntopo = topos.size();
    vector<int>groupind(ntopo,-1);
    vector<vector<vector<int>>>label2closeSurfGenus(ntopo,vector<vector<int>>(_nMat));
    for(int j=0;j<ntopo;++j){
        auto &label2closeSurfGenus_onetopo = label2closeSurfGenus[j];
        CellTopology &toposj = topos[j];
        for(int i=0;i<toposj.isoComps2label.size();++i)if(DP_label2NComp[toposj.isoComps2label[i]]!=-1){
            label2closeSurfGenus_onetopo[toposj.isoComps2label[i]].push_back(toposj.isoComps2VEF[i]);
        }
    }
    vector<vector<int>>openSurfGenus(ntopo);
    vector<vector<vector<int>>>openSurfComps(ntopo);
    for(int j=0;j<ntopo;++j){
        CellTopology &toposj = topos[j];
        vector<int>comp2label(toposj._comps.size());
        Comp2label(toposj._comps,mergeLoop2label[stepi],comp2label);
        auto &openSurfGenus_onetopo = openSurfGenus[j];
        auto &openSurfComps_onetopo = openSurfComps[j];
        for(int i=0;i<comp2label.size();++i)if(DP_label2NComp[comp2label[i]]!=-1){
            openSurfComps_onetopo.push_back(toposj._comps[i]);
            openSurfGenus_onetopo.push_back(toposj.comp2VEF[i]);
        }
    }
    vector<int>nJP(ntopo);
    for(int j=0;j<ntopo;++j)nJP[j] = topos[j]._junctionpoints;

    auto isequal = [this,&label2closeSurfGenus,&openSurfGenus,&openSurfComps,&nJP](int i,int j){
        if(openSurfComps[i]!=openSurfComps[j])return false;
        if(openSurfGenus[i]!=openSurfGenus[j])return false;
        if(label2closeSurfGenus[i]!=label2closeSurfGenus[j])return false;
        if(isCJP)if(nJP[i]!=nJP[j])return false;
        return true;

    };




    int ngroup = 0;
    for(int i=0;i<ntopo;++i){
        if(groupind[i]!=-1)continue;
        groupind[i] = ngroup;ngroup++;
        for(int j=i+1;j<ntopo;++j){
            if(groupind[j]!=-1)continue;
            if(isequal(i,j)/*topos[i]==topos[j]*/){
                groupind[j] = groupind[i];
            }
        }
    }
    vector<CellTopology>mergeTopo;
    for(int i=0;i<ngroup;++i){
        int maxSInd = -1;
        float maxS = -1e20;
        for(int j=0;j<ntopo;++j)if(groupind[j]==i)if(topos[j]._score>maxS){
            maxSInd = j;
            maxS = topos[j]._score;
        }
        mergeTopo.push_back(topos[maxSInd]);

    }

    topos = mergeTopo;




}
void triggerWarningbox(string s){
	cout<<s<<endl;
}
bool ParallelArrangement::DynamicProgramming(){

    vector<CellTopology>mergeTopoContainer1,mergeTopoContainer2;
    vector<CellTopology>*preContainer = &mergeTopoContainer1, *curContainer = &mergeTopoContainer2;

    auto mergingorder = newCtrNet.mergingorder;
    int beCellid = mergingorder[0];



    clock_t start = clock(), end,t1,t2,t3;


    double allowMinute = 6.5;
    if(1){
        string allowMinuteFname("../allowMinutes.txt");
        ifstream fin;
        fin.open(allowMinuteFname);
        if(fin.fail()){
            cout<<"cant open file: "<<allowMinuteFname<<endl;
        }else{
            fin>>allowMinute;
        }
    }

    cout<<"allowMinute: "<<allowMinute<<endl;


    int maxAllCandidate = 600000. * allowMinute;

    int maxCperCell = maxAllCandidate/(_nCell-4);

    CellTopology mergeTopo;

    vector<int>preTopo2beMInd;


    start = clock();

    for(auto &a:_cellTopo)for(int i=0;i<a.size();++i){
        a[i].mSeq.clear();a[i].mSeq.push_back(i);
    }
    int constraintTopo = DP_constraintTopos[beCellid];
    if(constraintTopo!=-1){
        (*preContainer).push_back(_cellTopo[beCellid][constraintTopo]);
    }else{
        *preContainer=_cellTopo[beCellid];
        if(preContainer->size()==0){
            preContainer->resize(1);
            (preContainer->at(0)).mSeq.push_back(-1);
        }
    }


    for(auto &a:*preContainer)MapComponentLoops2CompactLoop(beCellid,a);

    for(int stepi=1;stepi<_nCell;++stepi){

        t1 = clock();
        int newCellid = mergingorder[stepi];
        curContainer->clear();
        auto &bemergingTopo = _cellTopo[newCellid];

        int constraintTopo = DP_constraintTopos[newCellid];

        if(bemergingTopo.size()==0){
            for(auto &a:*preContainer)a.mSeq.push_back(-1);
            if(stepi==_nCell-1){
                for(int j=0;j<preContainer->size();++j ){
                    if(CheckValidation(stepi,preContainer->at(j))){
                        curContainer->push_back(preContainer->at(j));
                    }
                }
                swap(preContainer,curContainer);
            }

            continue;
        }

        for(auto &a:bemergingTopo)MapComponentLoops2CompactLoop(newCellid,a);

        int nAllc = preContainer->size()*bemergingTopo.size();
        vector<pair<double,int>>pickPreCandidates(preContainer->size());
        for(int ddd = 0;ddd<pickPreCandidates.size();++ddd)pickPreCandidates[ddd] = make_pair(preContainer->at(ddd)._score,ddd);
        if(nAllc>maxCperCell){
            //random_shuffle(pickPreCandidates.begin(),pickPreCandidates.end());
            sort(pickPreCandidates.begin(),pickPreCandidates.end(),[](pair<double,int>&x, pair<double,int>&y){ return x.first > y.first; });
            int cutend = double(maxCperCell)/double(nAllc) * double(pickPreCandidates.size());
            pickPreCandidates.resize(cutend);
            cout<<"cut: "<<preContainer->size()<<" -> "<<pickPreCandidates.size()<<endl;
        }


        for(int i=0;i<bemergingTopo.size();++i){
            if(constraintTopo!=-1 && i!= constraintTopo)continue;


            for(int j=0;j<pickPreCandidates.size();++j ){
                MergeTwoCellTopo(stepi,preContainer->at(pickPreCandidates[j].second),bemergingTopo[i],mergeTopo);
                if(CheckValidation(stepi,mergeTopo)){
                    curContainer->push_back(mergeTopo);
                }
            }

        }
        t2 = clock();
        auto nValid = curContainer->size();
        auto nAll = pickPreCandidates.size()*bemergingTopo.size();
        cout<<"DP: "<<stepi<<": "<<nAll<<" -> "<<nValid<<" -> ";
        MergeSameTopo(*curContainer,stepi);
        t3 = clock();
        cout<< curContainer->size() <<endl;
        cout<< double(t2 - t1) / (CLOCKS_PER_SEC) << "\t"<< double(t3 - t2) / (CLOCKS_PER_SEC) <<endl;
        cout<< double(t2 - t1) / (0.001*CLOCKS_PER_SEC)/(nAll) <<"\t"
            << double(t3 - t2) / (0.001*CLOCKS_PER_SEC)/nValid<<"\t"
            <<  double(t3 - t1) / (0.001*CLOCKS_PER_SEC)/(nAll) <<endl;
        //cout<<"candidates: "<<stepi<<' '<<curContainer->size()<<endl;
        //        if(stepi==12){
        //            cout<<"final: ";
        //            for(int i=0;i<curContainer->size();++i)cout<<CheckValidation(stepi,curContainer->at(i))<<' ';
        //            cout<<endl;
        //        }
        swap(preContainer,curContainer);
    }
    pickTopoInd.clear();
    pickTopoInd.resize(_nCell,0);

    cout<<"DP Done! ================================================================================================================"<<endl;
    end = clock();
    double time = double(end - start) / (0.001*CLOCKS_PER_SEC);
    cout << "Stage 1 time: "<<time_stage1<<" s"<<endl;
    cout << "DP time: " << time << " ms"<<endl;

    cout<<"Merging Order: ";for(auto a:mergingorder)cout<<a<<' ';cout<<endl;
    if(preContainer->size()!=0){
        auto &ReturnTopo = *max_element(preContainer->begin(),preContainer->end());

        for(int i=0;i<_nCell;++i)pickTopoInd[mergingorder[i]] =  ReturnTopo.mSeq[i];

        cout<<"Pick Topo: "<<CheckValidation(_nCell-1,ReturnTopo)<<endl;
        for(auto a:ReturnTopo.comp2VEF)cout<<(2-a)/2<<' ';cout<<endl;
        for(auto a:pickTopoInd)cout<<a<<' ';cout<<endl;
        cout<<"components: "<<endl;
        for(auto a:ReturnTopo.isoComps2VEF)cout<<(2-a)/2<<' ';cout<<endl;
        for(auto a:ReturnTopo.isoComps2label)cout<<a<<' ';cout<<endl;
        cout<<"Score: "<<ReturnTopo._score<<endl;
        cout<<"Junction Points: "<<ReturnTopo._junctionpoints<<endl;



        OutputDPInfo(outfoldername_meta,pickTopoInd,ReturnTopo._score,time,ReturnTopo.isoComps2VEF,ReturnTopo.isoComps2label);

        return true;
    }else{

        cout<<"Not found!"<<endl;
        triggerWarningbox("Not found!");
        vector<int>fake(_nMat,-2);
        OutputDPInfo(outfoldername_meta,pickTopoInd,0,-1.,fake,fake);
        return false;
    }


}

void ParallelArrangement::OutputDPInfo(string filepath, vector<int>pickTopoIndex, double totalscore, double time,vector<int>isoComp2VEF, vector<int>isoComps2label){

    string filename = filepath + "pickTopo";
    writeVecFile(filename,pickTopoIndex);


    vector<float>scoreAndtime;
    scoreAndtime.push_back(totalscore);scoreAndtime.push_back(time);


    filename = filepath + "scoreAndtime";
    writeVecFile(filename,scoreAndtime);

    for(auto &a:isoComp2VEF)a=(2-a)/2;
    filename = filepath + "isoComp2genus";
    writeVecFile(filename,isoComp2VEF);


    filename = filepath + "isoComps2label";
    writeVecFile(filename,isoComps2label);

    vector<int>label_constrains;
    for(int i=0;i<DP_label2CGenus.size();++i)label_constrains.push_back(DP_label2CGenus[i][0]);

    for(auto &a:label_constrains)if(a>10000000)a=-1;

    filename = filepath + "protocol";
    writeVecFile(filename,label_constrains);

    vector<float>stage_onetime;
    stage_onetime.push_back(time_stage1*1000);

    filename = filepath + "timeStage1";
    writeVecFile(filename,stage_onetime);
}


void ParallelArrangement::DynamicProgramming(string protocalname, vector<int> &constraintTopos){

    //    vector<int>DP_label2NComp;
    //    vector<vector<int>>DP_label2CGenus;




    ifstream fin;
    fin.open(protocalname);
    if(fin.fail()){
        cout<<"cant open file: "<<protocalname<<endl;
        return;
    }

    if(1){
        DP_label2NComp.clear();DP_label2CGenus.clear();
        DP_label2NComp.resize(_nMat);
        for(int i=0;i<_nMat;++i){
            fin>>DP_label2NComp[i];
        }
        DP_label2CGenus.resize(_nMat);
        int inNum;
        for(int i=0;i<_nMat;++i){
            if(DP_label2NComp[i]==-1){
                fin>>inNum;
                if(inNum!=-1){
                    cout<<"Wrong protocal input!"<<endl;
                    return;
                }
                DP_label2CGenus[i].push_back(inNum);
            }else{
                for(int j=0;j<DP_label2NComp[i];++j){
                    fin>>inNum;
                    DP_label2CGenus[i].push_back(inNum);
                }
            }
        }

        DP_groupping.clear();DP_grouppinggenus.clear();DP_loops2groupping.clear();
        DP_loops2groupping.resize(_nLoops,-1);
        //string aaa = '@';
        //string oneline;
        //      while(true){
        //            getline( fin, oneline );
        //            fin>>aaa;
        //           if(aaa=="#")break;
        //        vector<vector<int>>DP_groupping;
        //        vector<int>DP_grouppinggenus;
        //        vector<int>DP_loops2groupping;
        //            while( getline( fin, oneline ) ){
        //                stringstream strs( oneline );
        //                string prefix;

        //                strs >> prefix;


        //            vector<int>grouptest({3,11});
        //            DP_groupping.push_back(grouptest);
        //            DP_grouppinggenus.push_back(0);

        //            DP_loops2groupping.resize(_nLoops,-1);
        //            for(int i = 0;i<DP_groupping.size();++i)for(auto a:DP_groupping[i])DP_loops2groupping[a] = i;
        //        }

    }else{
        vector<int>a = {-1,2,2};
        vector<vector<int>>b = {{-1},{0,0},{-1,-1}};

        DP_label2NComp = a;
        DP_label2CGenus = b;
    }

    for(int i =0;i<DP_label2NComp.size();++i)if(DP_label2NComp[i]!=-1)if(DP_label2NComp[i]!=DP_label2CGenus[i].size()){
        cout<<"Illegal topology input! "<<endl;
        return;
    }
    int maxFlag = numeric_limits<int>::max();
    for(auto &c:DP_label2CGenus)for(auto &d:c)if(d==-1)d = maxFlag;
    for(auto &c:DP_label2CGenus)sort(c.begin(),c.end());

    for(auto &c:DP_label2CGenus){for(auto &d:c)cout<<d<<" ";cout<<endl;}



    DP_constraintTopos = constraintTopos;
    cout<<"DP_constraintTopos: ";for(auto a:DP_constraintTopos)cout<<a<<' ';cout<<endl;
    for(int i=0;i<_nCell;++i){
        if(DP_constraintTopos[i] < -1 || DP_constraintTopos[i]>=int(_cellTopo[i].size())){
            cout<<"Wrong Constraint Topos "<<i<<' '<< DP_constraintTopos[i] <<' '<<_cellTopo[i].size()<<endl;
            return;
        }
    }



    DynamicProgramming();


}



/******************************old backups************************************/
/*
extern bool isCJP;
static bool operator==(const CellTopology& lh, const CellTopology& rh) {
    if(isCJP)return (lh._comps == rh._comps) && (lh.comp2VEF == rh.comp2VEF);
    return (lh._comps == rh._comps) && (lh.comp2VEF == rh.comp2VEF) && (lh._junctionpoints == rh._junctionpoints);
}
static bool operator<(const CellTopology& lh, const CellTopology& rh) {
    if(isCJP)if(lh._junctionpoints != rh._junctionpoints )return (lh._junctionpoints > rh._junctionpoints );
    return (lh._score < rh._score );
}
static bool operator>(const CellTopology& lh, const CellTopology& rh) {
    if(isCJP)if(lh._junctionpoints != rh._junctionpoints )return (lh._junctionpoints < rh._junctionpoints );
    return (lh._score > rh._score );
}
void MergeSameTopo(vector<CellTopology>&topos){

    int ntopo = topos.size();
    vector<int>groupind(ntopo,-1);
    int ngroup = 0;
    for(int i=0;i<ntopo;++i){
        if(groupind[i]!=-1)continue;
        groupind[i] = ngroup;ngroup++;
        for(int j=i+1;j<ntopo;++j){
            if(groupind[j]!=-1)continue;
            if(topos[i]==topos[j]){
                groupind[j] = groupind[i];
            }
        }
    }
    vector<CellTopology>mergeTopo;
    for(int i=0;i<ngroup;++i){
        int maxSInd = -1;
        float maxS = -1e20;
        for(int j=0;j<ntopo;++j)if(groupind[j]==i)if(topos[j]._score>maxS){
            maxSInd = j;
            maxS = topos[j]._score;
        }
        mergeTopo.push_back(topos[maxSInd]);

    }

    topos = mergeTopo;




}

bool ParallelArrangement::CheckValidation(int stepi,CellTopology &curCell){

    vector<int>goalGenus(_nMat,0);//goalGenus[2]=1;//goalGenus[0]=1;
    //goalGenus[2]=1;goalGenus[3]=1;goalGenus[4]=3;
    goalGenus[2]=1;
    //goalGenus[2]=1;goalGenus[3]=6;
    //goalGenus[0]=0;


    vector<int>TopoComp2label(curCell._comps.size());
    Comp2label(curCell._comps,mergeLoop2label[stepi],TopoComp2label);


    if(stepi!=_nCell-1){
        //cout<<"stepi: "<<stepi<<endl;
        //if(curCell.isContainIso)return false;
        for(int i=0;i<curCell._comps.size();++i){
            int ge = (2-curCell._comps[i].size()-curCell.comp2VEF[i]);
            assert(ge%2==0);
            assert(ge>=0);
            int genus = ge/2;
            if(genus>goalGenus[TopoComp2label[i]])return false;
        }
        if(curCell.isoComps2label.size()!=0){
            for(auto a : curCell.isoComps2label)if(maxExistStep[a]>stepi)return false;

            vector<int>numofSurpM(_nMat,0);
            for(auto a:curCell.isoComps2label)++numofSurpM[a];
            for(auto a:numofSurpM)if(a>1)return false;

            for(int i=0;i<curCell.isoComps2label.size();++i){
                int ge = (2-curCell.isoComps2VEF[i]);
                assert(ge%2==0);
                assert(ge>=0);
                int genus = ge/2;
                if(genus!=goalGenus[curCell.isoComps2label[i]])return false;
            }
        }
        return true;
    }else{

        //need change for multi material for early closing
        if(curCell.isoComps2label.size()!=_nMat){
            return false;
        }

        vector<int>numofSurpM(_nMat,0);
        for(auto a:curCell.isoComps2label)++numofSurpM[a];
        for(auto a:numofSurpM)if(a!=1)return false;


        for(int i=0;i<curCell.isoComps2label.size();++i){
            int ge = (2-curCell.isoComps2VEF[i]);
            assert(ge%2==0);
            assert(ge>=0);
            int genus = ge/2;
            if(genus!=goalGenus[curCell.isoComps2label[i]])return false;
        }

        return true;
    }


}



*/
/******************************old backups************************************/
