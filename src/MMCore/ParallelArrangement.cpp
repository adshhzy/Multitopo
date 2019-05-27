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
// using Eigen::Vector3f;
// using namespace cimg_library;

typedef vector<vector<rgb>> imagecols;

void ParallelArrangement::ExportToMathematica(const char* filename) {
    ofstream ofs(filename, ofstream::out);
    ofs << "{";

    MingUtility::writeVector_wB_tbc(_frameVs, ofs);
    MingUtility::writeVector_wB_tbc(_frameFs, ofs);
    MingUtility::writeVector_wB_tbc(_frameCells, ofs);
    MingUtility::writeVector_wB_tbc(_frameF2PlaneParas, ofs);

    MingUtility::writeVector_wB_tbc(_tetVs, ofs);
    MingUtility::writeVector_wB_tbc(_frameF2tetFs, ofs);
    MingUtility::writeVector_wB_tbc(_cell2tetTs, ofs);
    MingUtility::writeVector_wB_tbc(_cell2tetVIds, ofs);

    MingUtility::writeVector_wB_tbc(_constTetV2Label, ofs);

    MingUtility::writeVector_wB(_cell2labeling, ofs);

    ofs << "}";
    ofs.close();
}

void ParallelArrangement::loadCrossSections(const vector<string>& filenames,
                                            float spacing, float tetVolLimit) {
    if (filenames.empty()) return;
    vector<float> spacings(filenames.size() - 1, spacing);
    loadCrossSections(filenames, spacings, tetVolLimit);
}

void ParallelArrangement::loadCrossSections(const vector<string>& filenames,
                                            vector<float> spacings,
                                            float tetVolLimit) {
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
    // cout << "_nMat=" << _nMat << endl;

    //--------------------- Construct the frame ---------------------
    // <1> frameVs.
    // Two addtional layers are added on the top and bottom. The spaceing of each
    // addtional layer to it adjacent layer is set to be the minimum spacing among
    // all the spacings.
    // |spacings| = _nCrossSection-1; |extSpacings| = _nCrossSection+2 = nLayer.
    float bdSpacing = *min_element(spacings.begin(), spacings.end());
    vector<float> extSpacings = spacings;
    // For the second layer (first true layer).
    extSpacings.insert(extSpacings.begin(), bdSpacing);
    // For the first layer.
    extSpacings.insert(extSpacings.begin(), 0.0);
    // For the last layer.
    extSpacings.push_back(bdSpacing);
    int nLayer = _nCrossSection + 2;
    _nCell = _nCrossSection + 1;
    for (int i = 1; i < extSpacings.size(); ++i) {
        extSpacings[i] += extSpacings[i - 1];
    }
    _frameVs.reserve(nLayer * 4);
    for (int i = 0; i < nLayer; ++i) {
        _frameVs.push_back({0, 0, extSpacings[i]});
        _frameVs.push_back({1.0f * width, 0., extSpacings[i]});
        _frameVs.push_back({1.0f * width, 1.0f * height, extSpacings[i]});
        _frameVs.push_back({0, 1.0f * height, extSpacings[i]});
    }
    // <2> frameFs.
    // frameFs = shared Fs and side Fs.
    // Shared Fs.
    for (int i = 0; i < nLayer; ++i) {
        _frameFs.push_back({i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 3});
    }
    // Side Fs.
    for (int i = 0; i < _nCell; ++i) {
        int pre = i * 4;
        _frameFs.push_back({pre + 0, pre + 1, pre + 5, pre + 4});
        _frameFs.push_back({pre + 1, pre + 2, pre + 6, pre + 5});
        _frameFs.push_back({pre + 2, pre + 3, pre + 7, pre + 6});
        _frameFs.push_back({pre + 3, pre + 0, pre + 4, pre + 7});
    }
    // <3> frameCells.
    for (int i = 0; i < _nCell; ++i) {
        int pre = i * 4 + nLayer;
        _frameCells.push_back({i, i + 1, pre, pre + 1, pre + 2, pre + 3});
    }
    // <4> frameF2PlaneParas.
    // Four parameters used to define each plane {a,b,c,d} in the form of
    // ax+by+cz+d=0.
    for (int i = 0; i < nLayer; ++i) {
        _frameF2PlaneParas.push_back({0, 0, 1, -1.0f * extSpacings[i]});
    }

    // ----------------------- Generate tets ------------------------
    // Find the centroid of each cell for classifying tets in to cells.
    vector<vector<float>> cellCenter(_nCell, vector<float>(3));
    for (int i = 0; i < _nCell; ++i) {
        auto& res = cellCenter[i];
        for (int p = i * 4; p < i * 4 + 8; ++p) {
            std::transform(res.begin(), res.end(), _frameVs[p].begin(), res.begin(),
                           std::plus<float>());
        }
        MingUtility::MultiplyVectorsByValue(0.125f, res);
    }
    // Prepare the data structure for tetgen.
    tetgenio in, out;
    tetgenio::facet* f;
    tetgenio::polygon* p;
    // PLC: pointlist.
    in.numberofpoints = _frameVs.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < _frameVs.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            in.pointlist[i * 3 + j] = _frameVs[i][j];
        }
    }
    // PLC: facetlist.
    in.numberoffacets = _frameFs.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (int i = 0; i < _frameFs.size(); ++i) {
        f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = _frameFs[i].size();
        p->vertexlist = new int[p->numberofvertices];
        for (int j = 0; j < _frameFs[i].size(); ++j) {
            p->vertexlist[j] = _frameFs[i][j];
        }
        in.facetmarkerlist[i] = i;
    }

    // PLC: region.
    in.numberofregions = _frameCells.size();
    in.regionlist = new REAL[in.numberofregions * 5];
    for (int i = 0; i < _frameCells.size(); ++i) {
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

    // Tetrahedralize the PLC.
    // Switches are chosen to read PLC (p), do quality mesh generation (q) with a
    // specified quality bound (1.414), apply a maximum volume constraint
    // (a10000) and omit terminal output except errors (Q).
    // -Qpq1.4a10000
    string switchesStr = "Qpq1.414Aa" + to_string(tetVolLimit);
    // char switches[switchesS.size()] = switchesS.c_str();
    // char switches[] = "Qpq1.414Aa10000";
    std::vector<char> tmp(switchesStr.begin(), switchesStr.end());
    tmp.push_back('\0');
    char* switches = &tmp[0];
    // char switches[] = "Qpq1.414Aa10";
    tetrahedralize(switches, &in, &out);

    char tetout[] = "../tetgen1.4.3/c++tetout";
    // out.save_nodes(tetout);
    // out.save_elements(tetout);
    // out.save_faces(tetout);

    // Extract the output tets from tetgen data structure to our own data
    // structure.
    _tetVs.resize(out.numberofpoints, vector<float>(3));
    for (int i = 0; i < out.numberofpoints; ++i) {
        for (int j = 0; j < 3; ++j) {
            _tetVs[i][j] = out.pointlist[i * 3 + j];
        }
    }
    _frameF2tetFs.resize(_frameFs.size());
    vector<int> oneTetF(3);
    for (int i = 0; i < out.numberoftrifaces; ++i) {
        int marker = out.trifacemarkerlist[i];
        for (int j = 0; j < 3; ++j) {
            oneTetF[j] = out.trifacelist[i * 3 + j];
        }
        sort(oneTetF.begin(), oneTetF.end());
        _frameF2tetFs[marker].push_back(oneTetF);
    }

    _cell2tetTs.resize(_frameCells.size());
    vector<unordered_set<int>> cellVsSet(_frameCells.size());
    assert(out.numberoftetrahedronattributes == 1);
    vector<int> oneTetT(4);
    for (int i = 0; i < out.numberoftetrahedra; ++i) {
        int marker = -out.tetrahedronattributelist[i] - 1;
        for (int j = 0; j < 4; ++j) {
            oneTetT[j] = out.tetrahedronlist[i * 4 + j];
            cellVsSet[marker].insert(oneTetT[j]);
        }
        sort(oneTetT.begin(), oneTetT.end());
        _cell2tetTs[marker].push_back(oneTetT);
    }
    _cell2tetVIds.resize(_frameCells.size());
    for (int i = 0; i < _frameCells.size(); ++i) {
        _cell2tetVIds[i].resize(cellVsSet[i].size());
        copy(cellVsSet[i].begin(), cellVsSet[i].end(), _cell2tetVIds[i].begin());
    }
    _cell2MapTetV2CellV.resize(_nCell);
    for (int i = 0; i < _nCell; ++i) {
        auto& vMap = _cell2MapTetV2CellV[i];
        for (int j = 0; j < _cell2tetVIds[i].size(); ++j) {
            int tetVi = _cell2tetVIds[i][j];
            vMap[tetVi] = j;
        }
    }

    // ----------------------- Label constrainted tetVs -----------------------
    // By default, everything is labeled as background: label 0.
    _constTetV2Label.resize(_tetVs.size(), 0);
    for (int i = 0; i < _nCrossSection; ++i) {
        int fi = i + 1;
        unordered_set<int> faceVs;
        for (const auto& oneF : _frameF2tetFs[fi]) {
            for (int j = 0; j < 3; ++j) {
                faceVs.insert(oneF[j]);
            }
        }
        for (auto vi : faceVs) {
            int x = min(max((int)_tetVs[vi][0], 0), width - 1);
            int y = min(max((int)_tetVs[vi][1], 0), height - 1);
            int label = rgb2label[image2cols[i][x][y]];
            _constTetV2Label[vi] = label;
        }
    }
    // MingUtility::printVector(_constTetV2Label);

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
        for (auto vi : _cell2tetVIds[i]) {
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
    _frameCell2csFrameIDs.push_back({1});  // first cell.
    for (int i = 1; i < _nCell - 1; ++i) {
        _frameCell2csFrameIDs.push_back({i, i + 1});
    }
    _frameCell2csFrameIDs.push_back({_nCrossSection});  // last cell.
    _isFrameFReported.resize(nFrameF, true);
    // The cross sections are frameF {1,2,...,_nCrossSection}.
    vector<vector<int>> ePair;
    MingUtility::combination(3, 2, ePair);
    for (int fi = 1; fi <= _nCrossSection; ++fi) {
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
            assert(oneE.second.size() == 2);
            _frameF2CrvEs[fi].push_back(oneE.second);
        }
    }
}
void ParallelArrangement::GetCell(int celli, vector<vector<float>>& cellVs,
                                  vector<vector<int>>& cellTs,
                                  vector<int>& cellInitLabel, int& cellNMat) {
    if (celli >= _nCell) return;
    cellTs = _cell2tetTs[celli];
    cellInitLabel = _cell2labeling[celli];
    cellNMat = _cell2labels[celli].size();
    int nV = _cell2tetVIds[celli].size();
    cellVs.resize(nV);
    for (int i = 0; i < nV; ++i) {
        cellVs[i] = _tetVs[_cell2tetVIds[celli][i]];
    }
    for (int ti = 0; ti < cellTs.size(); ++ti) {
        for (int i = 0; i < 4; ++i) {
            cellTs[ti][i] = _cell2MapTetV2CellV[celli][cellTs[ti][i]];
        }
    }
}

void ParallelArrangement::ReportCrossSectionCurve(int celli,
                                                  const fGroup& activeFs,
                                                  const vGroup& crvVs,
                                                  const eGroup& crvEs) {
    const auto& csFFIds = _frameCell2csFrameIDs[celli];
    // Check in the beginning whether all related cs are all reported.
    bool allReported = true;
    for (auto ffi : csFFIds) {
        if (!_isFrameFReported[ffi]) {
            allReported = false;
            break;
        }
    }
    if (allReported) return;
    // Use a map to remember the mapping from activeF (3 vertices) to crvV id.
    map<vector<int>, int> mapActiveF2Vid;
    vector<int> oneF(3);
    const auto& mapCellV2TetV = _cell2tetVIds[celli];
    assert(activeFs.size() == crvVs.size());
    for (int cvi = 0; cvi < activeFs.size(); ++cvi) {
        const auto& actF = activeFs[cvi];
        for (int i = 0; i < 3; ++i) {
            // Convert per cell index of a vertex to the arrangement index.
            oneF[i] = mapCellV2TetV[actF[i]];
        }
        sort(oneF.begin(), oneF.end());
        mapActiveF2Vid[oneF] = cvi;
    }
    // For each unreported cross section, find all its activeFs. Save crvV & E.
    for (auto ffi : csFFIds) {
        if (_isFrameFReported[ffi]) continue;
        // cout << "reporting face: " << ffi << endl;
        auto& activeFs = _frameF2ActiveFs[ffi];
        auto& curveVs = _frameF2CrvVs[ffi];
        auto& curveEs = _frameF2CrvEs[ffi];
        // Save the crvVs.
        // cout << "activeFs.size() = " << activeFs.size() << endl;
        for (const auto& oneTri : activeFs) {
            if (mapActiveF2Vid.find(oneTri) == mapActiveF2Vid.end()) {
                cout << "Error! An active face is not reported!" << endl;
                MingUtility::printVector(oneTri);
                return;
            } else {
                curveVs.push_back(crvVs[mapActiveF2Vid[oneTri]]);
            }
        }
        // Find the two side material of the crvEs.
        // cout << "curveEs.size() = " << curveEs.size() << endl;
        for (int ei = 0; ei < curveEs.size(); ++ei) {
            auto& oneE = curveEs[ei];
            // cout << "oneE.size() = " << oneE.size() << endl;
            // MingUtility::printVector(oneE);
            const auto& f0 = activeFs[oneE[0]];
            const auto& f1 = activeFs[oneE[1]];
            unordered_map<int, int> v2n;
            for (auto vi : f0) {
                v2n[vi]++;
            }
            for (auto vi : f1) {
                v2n[vi]++;
            }
            vector<int> sharedVs;
            for (const auto& oneV : v2n) {
                // cout << "oneV: " << oneV.first << "," << oneV.second << endl;
                if (oneV.second == 2) {
                    sharedVs.push_back(oneV.first);
                }
            }
            assert(sharedVs.size() == 2);

            Eigen::Vector3f ve0(curveVs[oneE[0]].data());
            Eigen::Vector3f ve1(curveVs[oneE[1]].data());
            Eigen::Vector3f vm0(_tetVs[sharedVs[0]].data());
            Eigen::Vector3f vm1(_tetVs[sharedVs[1]].data());
            Eigen::Vector3f norm(0.0f, 0.0f, 1.0f);
            if (((ve1 - ve0).cross(vm1 - vm0)).dot(norm) > 0) {
                reverse(sharedVs.begin(), sharedVs.end());
            }
            oneE.push_back(_constTetV2Label[sharedVs[0]]);
            oneE.push_back(_constTetV2Label[sharedVs[1]]);
        }
        _isFrameFReported[ffi] = true;
    }
}
void ParallelArrangement::ExportCrossSectionCurves(const char* filename) {
    int nPlane = 0;
    for (int fi = 1; fi <= _nCrossSection; ++fi) {
        if (_isFrameFReported[fi]) {
            nPlane++;
        }
    }
    if (nPlane < 1) {cout<<"no report!"<<endl;return;}
    ofstream ofs(filename, ofstream::out);
    ofs << nPlane << "\n";
    for (int fi = 1; fi <= _nCrossSection; ++fi) {
        if (!_isFrameFReported[fi]) continue;
        const auto& paras = _frameF2PlaneParas[fi];
        const auto& crvVs = _frameF2CrvVs[fi];
        const auto& crvEs = _frameF2CrvEs[fi];
        ofs << paras[0] << " " << paras[1] << " " << paras[2] << " " << paras[3]
                        << "\n";
        ofs << crvVs.size() << " " << crvEs.size() << "\n";
        for (const auto& oneV : crvVs) {
            ofs << oneV[0] << " " << oneV[1] << " " << oneV[2] << "\n";
        }
        for (const auto& oneE : crvEs) {
            ofs << oneE[0] << " " << oneE[1] << " " << oneE[2] << " " << oneE[3]
                           << "\n";
        }
    }
    ofs.close();
}

void ParallelArrangement::ExportCrossSectionCurvesNonParallel(const char* filename) {
    int nPlane = 0;
    for (int fi = 0; fi < _frameFs.size(); ++fi) {
        if (_framefisCs[fi] && _isFrameFReported[fi]) {
            nPlane++;
        }
    }
    if (nPlane < 1) {cout<<"no report!"<<endl;return;}
    ofstream ofs(filename, ofstream::out);
    ofs << nPlane << "\n";
    for (int fi = 0; fi < _frameFs.size(); ++fi) {
        if (!_framefisCs[fi] || !_isFrameFReported[fi] ) continue;
        const auto& paras = _frameF2PlaneParas[fi];
        const auto& crvVs = _frameF2CrvVs[fi];
        const auto& crvEs = _frameF2CrvEs[fi];
        ofs << paras[0] << " " << paras[1] << " " << paras[2] << " " << paras[3]
                        << "\n";
        ofs << crvVs.size() << " " << crvEs.size() << "\n";
        for (const auto& oneV : crvVs) {
            ofs << oneV[0] << " " << oneV[1] << " " << oneV[2] << "\n";
        }
        for (const auto& oneE : crvEs) {
            ofs << oneE[0] << " " << oneE[1] << " " << oneE[2] << " " << oneE[3]
                           << "\n";
        }
    }
    ofs.close();
}

// I assume all the related frame Fs are reported.
void ParallelArrangement::ExportCrossSectionCurves(int celli,
                                                   const char* filename) {
    const auto& csFrameIDs = _frameCell2csFrameIDs[celli];
    ofstream ofs(filename, ofstream::out);
    if (!ofs) {
        cout << "Error! " << filename << " cannot be found!" << endl;
        return;
    }
    ofs << csFrameIDs.size() << "\n";
    for (auto fi : csFrameIDs) {
        const auto& paras = _frameF2PlaneParas[fi];
        const auto& crvVs = _frameF2CrvVs[fi];
        const auto& crvEs = _frameF2CrvEs[fi];
        ofs << paras[0] << " " << paras[1] << " " << paras[2] << " " << paras[3]
                        << "\n";
        ofs << crvVs.size() << " " << crvEs.size() << "\n";
        for (const auto& oneV : crvVs) {
            ofs << oneV[0] << " " << oneV[1] << " " << oneV[2] << "\n";
        }
        for (const auto& oneE : crvEs) {
            ofs << oneE[0] << " " << oneE[1] << " " << oneE[2] << " " << oneE[3]
                           << "\n";
        }
    }
    ofs.close();
}
void ParallelArrangement::GatherAllActiveFAndCrvVE() {
    int cnt = 0;
    // fGroup activeFs;
    for (int fi = 1; fi <= _nCrossSection; ++fi) {
        auto& actFs = _frameF2ActiveFs[fi];
        auto& crvVs = _frameF2CrvVs[fi];
        auto& crvEs = _frameF2CrvEs[fi];
        for (int i = 0; i < actFs.size(); ++i) {
            _mapActF2CrvVi[actFs[i]] = _allCrvVs.size();
            // cout << _mapActF2CrvVi[actFs[i]] << endl;
            // activeFs.push_back(actFs[i]);
            _allCrvVs.push_back(crvVs[i]);
        }
        for (const auto& oneE : crvEs) {
            _allCrvEs.push_back({oneE[0] + cnt, oneE[1] + cnt});
        }
        cnt += actFs.size();
    }
}
void ParallelArrangement::CombineInterfaces(const vector<Cell>& cells,
                                            const vector<int>& topoIDs,
                                            vector<vector<float>>& sufVs,
                                            vector<vector<int>>& sufFs,
                                            vector<vector<int>>& sufFMats,
                                            vector<vector<int>>& segEs) {
    sufVs = _allCrvVs;
    segEs = _allCrvEs;
    assert(cells.size() == _nCell && topoIDs.size() == _nCell);
    for (int ci = 0; ci < _nCell; ++ci) {
        int topoID = topoIDs[ci];
        const auto& actFs = cells[ci]._activeFs;
        const auto& interface = cells[ci]._sufs[topoID];
        const auto& cSufVs = interface._sufVs;
        const auto& cSufFs = interface._sufFs;
        const auto& cSufFMats = interface._sufFMats;
        int nActFs = actFs.size();
        int nV = cSufVs.size();
        int cnt = sufVs.size();

        vector<int> locV2globalV(cSufVs.size());
        for (int i = 0; i < nActFs; ++i) {
            auto oneF = actFs[i];
            for (int j = 0; j < 3; ++j) {
                oneF[j] = _cell2tetVIds[ci][oneF[j]];
            }
            sort(oneF.begin(), oneF.end());
            locV2globalV[i] = _mapActF2CrvVi[oneF];
        }
        for (int i = nActFs; i < nV; ++i) {
            locV2globalV[i] = cnt + i - nActFs;
        }

        sufVs.reserve(sufVs.size() + nV - nActFs);
        sufVs.insert(sufVs.end(), cSufVs.begin() + nActFs, cSufVs.end());
        for (int fi = 0; fi < cSufFs.size(); ++fi) {
            auto oneF = cSufFs[fi];
            for (int j = 0; j < 3; ++j) {
                oneF[j] = locV2globalV[oneF[j]];
            }
            sufFs.push_back(oneF);
            auto oneMPair = cSufFMats[fi];
            for (int j = 0; j < 2; ++j) {
                oneMPair[j] = _cell2labels[ci][oneMPair[j]];
            }
            // MingUtility::printVector(_cell2MapTetL2CellL[ci]);
            sufFMats.push_back(oneMPair);
        }
    }
}


