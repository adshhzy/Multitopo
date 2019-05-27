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

using namespace std;

int adsmain222(int argc, char** argv) {
  bool argExploreAll = true;
  int density = 2;
  float upbound = 1.0f;
  float segLenLowBound = 0.0f;
  float spacing = 100.0f;
  float tetVolLimit = 10000.0f;
  string csListFile = "../CrossSectionFiles.txt";
  string outDir = "../CrossSectionFiles.txt";
  bool outputPerCell = false;
  bool allowHighGenus = false;
  int c;
  //for(int i=0;i<argc;++i)cout<<argv[i]<<endl;
  optind=1;
  while ((c = getopt(argc, argv, "lgd:b:u:s:t:f:o:")) != -1) {
    switch (c) {
      case 'l':
        argExploreAll = false;
        cout << "Do not run the exhaustive search." << endl;
        break;
      case 'd':
        density = atoi(optarg);
        break;
      case 'u':
        upbound = atof(optarg);
        break;
      case 'b':
        segLenLowBound = atof(optarg);
        break;
      case 's':
        spacing = atof(optarg);
        break;
      case 't':
        tetVolLimit = atof(optarg);//v4.3
        break;
      case 'f':
        csListFile = optarg;
        break;
      case 'o':
        outDir = optarg;
        outputPerCell = true;
        break;
      case 'g':
        allowHighGenus = true;
        break;
      case '?':
        cout << "Bad argument setting!" << endl;
        break;
    }
  }
  cout << "Ray density = " << density << endl;
  cout << "Si upper bound = " << upbound << endl;
  cout << "Segment length lower bound = " << segLenLowBound << endl;
  cout << "Spacing between cross sections = " << spacing << endl;
  cout << "Tet vol limit = " << tetVolLimit << endl;
  cout << "Allow high genus = " << allowHighGenus << endl;
  if (outputPerCell)
    cout << "Output dir for per-cell interfaces = " << outDir << endl;
  cout << "input file: "<<csListFile<<endl;
  vector<string> csFiles;
  MingUtility::LoadFileToVectorOfStrings(csListFile.c_str(), csFiles);
  if (csFiles.size() < 2) return 0;
  cout << "load cross sections from:";
  MingUtility::printVector(csFiles);

  cout << "------------ Build Up Arrangement ------------" << endl;
  ParallelArrangement pArr;
  pArr.loadCrossSections(csFiles, spacing, tetVolLimit);
  cout << "Done" << endl;

  // pArr.ExportToMathematica("../ParallelArrangement.txt");

  cout << "------------ Process Each Cell ------------" << endl;
  int nCell = pArr._nCell;
  vector<Cell> resCells(nCell);
  for (int ci = 0; ci < nCell; ++ci) {
    cout << "\n++ cell[" << ci << "] ++" << endl;


    vector<vector<float>> cellVs;
    vector<vector<int>> cellTs;
    vector<int> cellInitLabel;
    int cellNMat;
    pArr.GetCell(ci, cellVs, cellTs, cellInitLabel, cellNMat);
    cout << "nV = " << cellVs.size() << " nT = " << cellTs.size()
         << " nMat = " << cellNMat << endl;
    MMLevelSet ls;


    vector<vector<int>> activeFs;
    vector<vector<float>> crvVs;
    vector<vector<int>> crvEs;
    ls.LoadCell(cellVs, cellTs, cellInitLabel, cellNMat);




    ls._mesh->GetCurveInfoForArrangement(activeFs, crvVs, crvEs);

    pArr.ReportCrossSectionCurve(ci, activeFs, crvVs, crvEs);
    OffsetVector seedOffset = vector<float>(cellNMat, 0.0);
    Labeling seedLabeling;
    ls.OffsetVector2Labeling(seedOffset, seedLabeling);
    resCells[ci]._activeFs = activeFs;
    if (outputPerCell) {
      string outfile = outDir + "/" + to_string(ci) + ".contour";
      pArr.ExportCrossSectionCurves(ci, outfile.c_str());
    }


    ls.setVarUpbound(upbound);
    ls.PreProcess();//for exhaustive search
    // Output the 0-offset interface.
    if (outputPerCell) {
        cout<<"outputPerCell"<<endl;
      vector<vector<float>> sufVs;
      vector<vector<int>> sufFs;
      vector<vector<int>> sufFMats;
      ls._mesh->GenerateLevelSetFromLabeling(seedLabeling, sufVs, sufFs,
                                             sufFMats);//isosurface
      string outfile = outDir + "/suf_" + to_string(ci) + "_0offset.suf";
      ls._mesh->WriteLevelSet(sufVs, sufFs, sufFMats, outfile.c_str());
    }
    clock_t start = clock(), end;
    start = clock();
    ls.ExploreTopologiesAlongRays(allowHighGenus, density, segLenLowBound);//ray shooting
    end = clock();
    float time = float(end - start) / CLOCKS_PER_SEC;
    cout << "time: " << time << endl;

    int nTopo = ls._cellTopologies.size();

    if (nTopo < 1) {
      cout << "No topology found in cell[" << ci
           << "], check the tetrahedralization again" << endl;
      return 0;
    }

    resCells[ci]._sufs.resize(nTopo);

    for (int i = 0; i < nTopo; ++i) {

      vector<vector<float>> sufVs;
      vector<vector<int>> sufFs;
      vector<vector<int>> sufFMats;
      ls._mesh->GenerateLevelSetFromLabeling(ls._cellTopologies[i]._label,
                                             sufVs, sufFs, sufFMats);
      resCells[ci]._sufs[i] = Interface(sufVs, sufFs, sufFMats);
      cout<<i<<' '<<nTopo<<endl;
      if (outputPerCell && i < NMAXSUFOUTPUT) {
        string outfile =
            outDir + "/suf_" + to_string(ci) + "_" + to_string(i) + ".suf";
        ls._mesh->WriteLevelSet(sufVs, sufFs, sufFMats, outfile.c_str());
      }
    }

  }
  cout << "Done" << endl;

  cout << "------------ Save Curves ------------" << endl;
  pArr.ExportCrossSectionCurves("../InputCurves.contour");//all curves
  cout << "Done" << endl;

  cout << "------------ Save Complete Surface ------------" << endl;
  vector<int> topoIDs(nCell, 0);
  vector<vector<float>> sufVs;
  vector<vector<int>> sufFs;
  vector<vector<int>> sufFMats;
  vector<vector<int>> segEs;
  pArr.GatherAllActiveFAndCrvVE();
  pArr.CombineInterfaces(resCells, topoIDs, sufVs, sufFs, sufFMats, segEs);
  MingUtility::WriteLevelSet(sufVs, sufFs, sufFMats, segEs, "../OutputSuf.suf");//all surfaces(materials)
  cout << "Done" << endl;
  return 0;
}




int createCrosssections(string inputfile) {
  bool argExploreAll = true;
  int density = 3;
  float upbound = 1.0f;
  float segLenLowBound = 0.0f;
  float spacing = 100.0f;
  float tetVolLimit = 1000.0f;
  //string csListFile = "/Users/Research/Geometry/MM/data/toy_3col/inputall.txt";
  string csListFile = inputfile;
  string outDir = "/Users/Research/Geometry/MM/result/toy_3col/";
  bool outputPerCell = false;
  bool allowHighGenus = false;

  cout << "Ray density = " << density << endl;
  cout << "Si upper bound = " << upbound << endl;
  cout << "Segment length lower bound = " << segLenLowBound << endl;
  cout << "Spacing between cross sections = " << spacing << endl;
  cout << "Tet vol limit = " << tetVolLimit << endl;
  cout << "Allow high genus = " << allowHighGenus << endl;
  if (outputPerCell)
    cout << "Output dir for per-cell interfaces = " << outDir << endl;
  cout << "input file: "<<csListFile<<endl;
  vector<string> csFiles;
  MingUtility::LoadFileToVectorOfStrings(csListFile.c_str(), csFiles);
  if (csFiles.size() < 2) return 0;
  cout << "load cross sections from:";
  MingUtility::printVector(csFiles);

  cout << "------------ Build Up Arrangement ------------" << endl;
  ParallelArrangement pArr;
  pArr.loadCrossSections(csFiles, spacing, tetVolLimit);
  cout << "Done" << endl;

  cout << "------------ Process Each Cell ------------" << endl;
  int nCell = pArr._nCell;
  vector<Cell> resCells(nCell);
  for (int ci = 0; ci < nCell; ++ci) {
    cout << "\n++ cell[" << ci << "] ++" << endl;

    vector<vector<float>> cellVs;
    vector<vector<int>> cellTs;
    vector<int> cellInitLabel;
    int cellNMat;
    pArr.GetCell(ci, cellVs, cellTs, cellInitLabel, cellNMat);
    cout << "nV = " << cellVs.size() << " nT = " << cellTs.size()
         << " nMat = " << cellNMat << endl;
    MMLevelSet ls;

    vector<vector<int>> activeFs;
    vector<vector<float>> crvVs;
    vector<vector<int>> crvEs;
    ls.LoadCell(cellVs, cellTs, cellInitLabel, cellNMat);


    ls._mesh->GetCurveInfoForArrangement(activeFs, crvVs, crvEs);

    pArr.ReportCrossSectionCurve(ci, activeFs, crvVs, crvEs);
  }
  cout << "Done" << endl;

  cout << "------------ Save Curves ------------" << endl;
  pArr.ExportCrossSectionCurves("../ConvertFolder/ConvertCurves.contour");//all curves
  cout << "Done" << endl;


  return 0;
}
