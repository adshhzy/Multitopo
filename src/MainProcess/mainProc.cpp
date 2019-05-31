#include "GroupingManager/GroupingManager.h"
#include <time.h>
#include "Utility/Utility.h"
#include "Utility/geo_sur.h"
#include "Utility/geo_curv.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <math.h>
#include "MMLevelSet.h"
#include "LP_solver.h"
#include "Cell.h"
#include "ParallelArrangement.h"
#include <iterator>
#include <string>
#include <dirent.h>
#include <sys/stat.h>

ParallelArrangement pArr;
vector<MMLevelSet>pLs;

void makemethodfolder(string foldername){
    struct stat sb;
    if(!(stat(foldername.data(), &sb) == 0 && S_ISDIR(sb.st_mode)))mkdir(foldername.data(),0777);
    return;
}


int az2(string incontour, string outdir, int MMprop_f2vmethod,
        int endstep, int density, vector<int>&numofTopo, vector<int>&cellSeq, int &nMat, int computeSaveLoad, string iofile){

    string cs_folder = outdir+string("cross-section/");
    string suf_folder = outdir+string("suf/");
    string meta_folder = outdir+string("meta/");


    makemethodfolder(cs_folder);
    makemethodfolder(suf_folder);
    makemethodfolder(meta_folder);

    //set output path for metadata
    pArr.SetMetaFolder(meta_folder);

    GroupingManager gm;
    clock_t start, end;
    vector<float> times;


    float upbound = 1.0f;
    float segLenLowBound = 0.0f;


    bool outputPerCell = true;
    bool allowHighGenus = false;



    cout << "\n[1] Preprocessing data & Exploring topologies in cells ..." << endl;
    start = clock();
    gm.setPara(outdir.data(), -1,-1,-1,-1,-1);
    gm.Contour2ArrangementMM(incontour.data(), NULL,NULL);
    end = clock();
    times.push_back( (float)(end - start) / CLOCKS_PER_SEC );
    cout << "Time: "<< times[times.size()-1] << "s \n\n\n";



    //return 0;


    pArr.loadInfoFromSpaceDivision(gm.ar,MMprop_f2vmethod,cs_folder);
    cellSeq = pArr.newCtrNet.mergingorder;

    nMat = pArr._nMat;

    if(endstep == 0)return pArr._frameFs.size();

    numofTopo.clear();


    cout << "------------ Process Each Cell ------------" << endl;
    clock_t totalstart = clock(), totalend;
    totalstart = clock();




    ifstream ifs;ofstream ofs;
    if(computeSaveLoad==2 || computeSaveLoad ==3){
        ifs.open(iofile,ios::in | ios::binary);
        if(ifs.fail()){
            cout<<"cant open file: "<<iofile<<endl;
            return 0;
        }
    }else if(computeSaveLoad==1){
        ofs.open(iofile,ios::out | ios::binary);
        if(ofs.fail()){
            cout<<"cant open file: "<<iofile<<endl;
            return 0;
        }
    }



    int nCell = pArr._nCell;
    pLs.clear();

    pLs.resize(nCell);


    pArr._cellNexplored2Topogroup.resize(nCell);
    for (int ci = 0; ci < nCell; ++ci) {
        cout << "\n++ cell[" << ci << "] ++" << endl;

        if(ci==nCell-1){
            cout<<ci<<endl;
        }

        vector<vector<float>> cellVs;
        vector<vector<int>> cellTs;
        vector<int> cellInitLabel,cellVMarkers;
        int cellNMat;
        pArr.GetCell(ci, cellVs, cellVMarkers,cellTs, cellInitLabel, cellNMat);
        cout << "nV = " << cellVs.size() << " nT = " << cellTs.size()
             << " nMat = " << cellNMat << endl;


        //nMat = max(nMat,cellNMat);
        if(cellNMat==1){numofTopo.push_back(0);pArr._cellNexplored2Topogroup[ci]=vector<int>(0);continue;}
        MMLevelSet &ls = pLs[ci];


        ls.LoadCell(cellVs, cellVMarkers,cellTs, cellInitLabel, cellNMat);



        pArr.MapLoopToSegs(ci,ls._mesh->_bdLoopsV);
        pArr.MapActiveFCP(ci,ls._mesh->_bdActiveF2Vs,ls._mesh->_bdCurveVs,ls._mesh->_bdActiveE2Vs_fornonplcs);



        OffsetVector seedOffset = vector<float>(cellNMat, 0.0);

        Labeling seedLabeling;
        ls.OffsetVector2Labeling(seedOffset, seedLabeling);
        //resCells[ci]._activeFs = activeFs;
        if (false) {
            string outfile = outdir  + to_string(ci) + ".contour";
            pArr.ExportCrossSectionCurves(ci, outfile.c_str());
        }


        ls.setVarUpbound(upbound);
        ls.PreProcess();//for exhaustive search
        // Output the 0-offset interface.

        auto &mappingMat = pArr._cell2labels[ci];
        if (outputPerCell) {
            cout<<"outputPerCell"<<endl;
            vector<vector<float>> sufVs;
            vector<vector<int>> sufFs;
            vector<vector<int>> sufFMats;
            ls._mesh->GenerateLevelSetFromLabeling(seedLabeling, sufVs, sufFs,
                                                   sufFMats);//isosurface
            //string outfile = outdir + "/suf_" + to_string(ci) + "_0offset.suf";
            string outfile = suf_folder + "suf_" + to_string(ci) + "_0.suf";
            ls._mesh->WriteLevelSet(sufVs, sufFs, sufFMats,mappingMat, outfile.c_str());
        }




        clock_t start = clock(), end;
        start = clock();

        if(computeSaveLoad==2 || computeSaveLoad ==3)ls.ReadTopologyFromFile(ifs);
        else ls.ExploreTopologiesAlongRays(allowHighGenus, density, segLenLowBound);//ray shooting

        if(computeSaveLoad==1)ls.WriteTopologyToFile(ofs);

        end = clock();
        float time = float(end - start) / CLOCKS_PER_SEC;
        cout << "time: " << time << endl;

        int nTopo = ls._cellTopologies.size();

        pArr.GetCellTopo(ci,ls._cellTopologies,ls._mesh->_label2bdLoops,ls._cellNexplored2Topogroup);

        numofTopo.push_back(max(1,nTopo));
        if (nTopo < 1) {
            cout << "No topology found in cell[" << ci
                 << "], check the tetrahedralization again" << endl;
            //return 0;
            continue;
        }

        //resCells[ci]._sufs.resize(nTopo);

        if(nTopo>500){
            cout<<"more than 500 topo!"<<endl;
           // exit(500);
        }
        Labeling topolables;
        if(computeSaveLoad!=3)for (int i = 0; i < nTopo; ++i) {

            cout<<i<<' '<<nTopo<<' '<<ls._cellTopologies[i]._junctionpoints<<' '<<ls._cellTopologies[i]._score<<endl;
            if (outputPerCell && true/*i < NMAXSUFOUTPUT*/) {
                vector<vector<float>> sufVs;
                vector<vector<int>> sufFs;
                vector<vector<int>> sufFMats;
                ls.OffsetVector2Labeling(ls._cellTopologies[i]._offset, topolables);
                ls._mesh->GenerateLevelSetFromLabeling(topolables,
                                                       sufVs, sufFs, sufFMats);
                string outfile =
                        suf_folder + "suf_" + to_string(ci) + "_" + to_string(i) + ".suf";
                ls._mesh->WriteLevelSet(sufVs, sufFs, sufFMats,mappingMat, outfile.c_str());
            }
        }



    }

    totalend = clock();
    float time = float(totalend - totalstart) / CLOCKS_PER_SEC;
    cout << "Total time: " << time << endl<< endl<< endl<< endl;

    pArr.time_stage1 = time;
    //pArr.WriteToMappingInfo();


    pArr.LoopMergingDP();
    pArr.MakeupLoop2SegsLists();

    //cout << "------------ Save Curves ------------" << endl;
    //pArr.ExportCrossSectionCurvesNonParallel("../OutputCrossSec.contour");//all curves

    cout << "Done" << endl;


    if(computeSaveLoad==2 || computeSaveLoad==3){
        ifs.close();
    }else if(computeSaveLoad==1){
        ofs.close();
    }




    return pArr._frameFs.size();

//    cout << "------------ Save Complete Surface ------------" << endl;
//    vector<int> topoIDs(nCell, 0);
//    vector<vector<float>> sufVs;
//    vector<vector<int>> sufFs;
//    vector<vector<int>> sufFMats;
//    vector<vector<int>> segEs;
//    pArr.GatherAllActiveFAndCrvVE();
//    pArr.CombineInterfaces(resCells, topoIDs, sufVs, sufFs, sufFMats, segEs);
//    MingUtility::WriteLevelSet(sufVs, sufFs, sufFMats, segEs, "../OutputSuf.suf");//all surfaces(materials)
    cout << "Done" << endl;




}

vector<int> runDynamicProgramming(string protocalname, vector<int>&TopoConstraintInd){

      pArr.DynamicProgramming(protocalname,TopoConstraintInd);
      return pArr.pickTopoInd;

}


void outputCombineSurface(string outfolder, vector<int>&Cell2NTopo, vector<int>&pickTopoInd){

    vector<vector<double>>TopoVs; vector<vector<uint>>TopoFs;
    vector<vector<int> > TopoFMs;vector<vector<uint>>TopoCtrs;

    n_rf::SufStructure pSuf;


    n_rf::MultiCellTopo displayMMCT;
    n_rf::MultiSurface displaySurface_multi;

    displayMMCT.ReadCellTopo(outfolder+string("suf"),Cell2NTopo,false);
    displayMMCT.GlobalSmoothing();

    displayMMCT.GetAllCellTopos(pickTopoInd, TopoVs, TopoFs,TopoFMs,TopoCtrs);
    pArr.CreateTotalSurf(TopoVs, TopoFs,TopoFMs,TopoCtrs,true,&pSuf);
    displaySurface_multi.ImportFromSuf(pSuf,true,false,true);
    displaySurface_multi.SmoothSurfNet_JuFair(50);

    displaySurface_multi.WriteObj(outfolder + "outcombine");
    displaySurface_multi.WriteSuf(outfolder + "outcombine");


}
