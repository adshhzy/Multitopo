#include <iostream>
#include <vector>
#include <queue>
#include "MMLevelSet.h"
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iomanip>
#include "tetgen.h"
#include "ParallelArrangement.h"

using namespace std;

void MMLevelSet::LoadCell(const vector<vector<float>>& cellVs,
                          const vector<vector<int>>& cellTs,
                          const vector<int>& cellInitLabel, int nMat) {
    _N = cellVs.size();
    _K = nMat;
    //TetMesh a;
    //_mesh = &a;
    if(_mesh!=NULL)delete _mesh;
    _mesh = new TetMesh();
    _mesh->setDimension(3);
    _mesh->Load(cellVs, cellTs);
    // TODO note that, right now, the indicator function is calculated for each
    // cell independently. If more than one cell shall be considered when
    // calculating the function for a specific cell (e.g. using tri-laplacian),
    // GenerateIndicatorFunction should be implemented in the ParallelArrangement
    // class, instead of in class Mesh.
    _initLabel = cellInitLabel;
    _mesh->GenerateIndicatorFunction(nMat, cellInitLabel, _V2F);
    PreProcessMeshBD();
}
void MMLevelSet::LoadCell(const vector<vector<float>>& cellVs,
                          const vector<int>&cellVMarkers,
                          const vector<vector<int>>& cellTs,
                          const vector<int>& cellInitLabel, int nMat) {
    _N = cellVs.size();
    _K = nMat;
    //TetMesh a;
    //_mesh = &a;
    if(_mesh!=NULL)delete _mesh;
    _mesh = new TetMesh();
    _mesh->setDimension(3);
    _mesh->Load(cellVs,cellVMarkers, cellTs);
    // TODO note that, right now, the indicator function is calculated for each
    // cell independently. If more than one cell shall be considered when
    // calculating the function for a specific cell (e.g. using tri-laplacian),
    // GenerateIndicatorFunction should be implemented in the ParallelArrangement
    // class, instead of in class Mesh.
    _initLabel = cellInitLabel;
    //_mesh->GenerateIndicatorFunction(nMat, cellInitLabel, _V2F);
    _mesh->GenerateIndicatorFunction(nMat, cellInitLabel,_mesh->_bdVs, _V2F);
    PreProcessMeshBD();

    GetCellBoundaryInfo();

}
void MMLevelSet::LoadMesh(const char* filename) {
    ifstream ifs(filename, ios::in);
    if (!ifs) {
        cout << "Error! " << filename << " cannot be found!" << endl;
        return;
    }
    int dim;
    if (ifs >> dim) {
        if (dim == 2) {
            //_mesh = new TriMesh();
        } else {
            _mesh = new TetMesh();
        }
    }
    ifs.close();
    _mesh->setDimension(dim);
    _mesh->Load(filename);
    PreProcessMeshBD();
}

void MMLevelSet::PreProcessMeshBD() {
    // V2F (on the boundary) must be set before calling PreProcessBD of mesh.
    Labeling seedLabeling;
    OffsetVector seedOffset(_K, 0.0);
    OffsetVector2Labeling(seedOffset, seedLabeling);
    // TODO _K should be nLabel instead or is it? Anyway, make it consistent later
    // when there are more than 1 cell and the cell does not use all the
    // materials.
    _mesh->PreProcessBD(seedLabeling, _K);
}
void MMLevelSet::GetCellBoundaryInfo(){

    _mesh->OutputCellBoundary(CellBoundary_Vs,CellBoundary_Fs,CellBoundary_VoriInd,CellBoundary_Vmarkers);

    CellBoundary_OriLabel = _initLabel;
}
// O(Nlog(N) * K^2)
void MMLevelSet::PreProcess() {
    bool db = false;
    // db = true;
#if db1
    db = true;
#endif

    // 1. Find the boundary vertices. (_bdVs)
    // Find all the vertices with one (and only one) of its scalar value to be 1.
    // Those are the boundary vertices with boundary conditions set.
    // Pre-calculate log(fi) for each V, for score calculation. (_V2LogF)
    _bdVs.resize(_N, false);
    _V2LogF.resize(_N, vector<float>(_K, 0.0));
    for (int p = 0; p < _N; ++p) {
        for (int k = 0; k < _K; ++k) {
            _V2LogF[p][k] = log(_V2F[p][k]);
            if (_V2F[p][k] == 1.0) {
                // _bdVs.insert(p);
                _bdVs[p] = true;
                // break;
            }
        }
    }

#if DOPERTURBE
    cout << "Perturbation on vertices" << endl;
    // 2. Add perturbation to avoid overlapping hyperplane. (_V2F)
    // Avoid overlapping hyper plane or v, by perturbing the V2F.
    // Do not perturbe boundary vertices.
    float diff, r, V2Fi, V2Fj;
    srand(time(NULL));
    bool perturbed = true;
    int nPerturb = -1;
    while (perturbed) {
        perturbed = false;
        nPerturb++;
        for (int i = 0; i < _K; ++i) {
            for (int j = i + 1; j < _K; ++j) {
                // unordered_set<float> exist;
                unordered_set<float> exist;
                for (int p = 0; p < _N; ++p) {
                    // avoid boundary points.
                    if (_bdVs[p]) continue;
                    diff = _V2F[p][i] - _V2F[p][j];
                    diff = round(diff * ROUNDFACTOR) / ROUNDFACTOR;
                    // TODO if do it this way, change V2Fi & j might make previous
                    // hyperplane to overlap. For now I have a while loop in the outside
                    // to prevent this. But I am concerned that if #V becomes huge, the
                    // while loop could run a huge number of times because my perturbation
                    // percision is fixed, will it?
                    while (exist.find(diff) != exist.end()) {
                        if (db) {
                            // if (i==0 && j==2) {
                            cout << "here! (i,j,p) = (" << i << "," << j << "," << p << ")"
                                 << endl;
                        }
                        perturbed = true;
                        // random number = 0~0.999 * PERTURBATION - 0.5*PERTURBATION
                        // 0 ~ 0.000999 - 0.0005 = -0.0005 ~ 0.000499
                        r = PERTURBATION * 0.001 * (rand() % 1000) - 0.5 * PERTURBATION;
                        V2Fi = _V2F[p][i] + r;
                        V2Fj = _V2F[p][j] - r;
                        if (V2Fi >= 0 && V2Fi <= 1 && V2Fj >= 0 && V2Fj <= 1) {
                            _V2F[p][i] = V2Fi;
                            _V2F[p][j] = V2Fj;
                            diff = _V2F[p][i] - _V2F[p][j];
                            diff = round(diff * ROUNDFACTOR) / ROUNDFACTOR;
                        }
                    }
                    exist.insert(diff);
                }
            }
        }
    }
    // cout << "nPerturb = " << nPerturb << endl;
    cout << nPerturb << endl;
    if (db) {
        cout << "V2F:" << endl;
        MingUtility::printVector(_V2F);
    }
    // MingUtility::writeVector(_V2F, "../new_v2f.txt");
#endif

    // 3. Find the ordering of vertices for each pair of material. (Orders)
    _Orders.resize(_K, vector<Ordering>(_K));
    for (int i = 0; i < _K; ++i) {
        for (int j = i + 1; j < _K; ++j) {
            vector<pair<float, int>> order(_N);
            // Ordering order(_N);
            for (int p = 0; p < _N; ++p) {
                order[p].first = _V2F[p][i] - _V2F[p][j];
                order[p].second = p;
            }
            // Sort vertices by decreasing order of fi(v) - fj(v).
            sort(order.begin(), order.end(), greater<pair<float, int>>());
            _Orders[i][j].resize(_N);
            for (int p = 0; p < _N; ++p) {
                // _Orders is only used in Labeling2Range, only (i,j) is needed
                // no need to store the symmetric (j,i)
                _Orders[i][j][p] = order[p].second;
            }
#if db2
            if (db) {
                cout << "_Orders[" << i << "][" << j << "] = ";
                MingUtility::printVector(_Orders[i][j]);
                vector<float> diffs;
                for (int p = 0; p < _N; ++p) {
                    diffs.push_back(_V2F[_Orders[i][j][p]][i] -
                            _V2F[_Orders[i][j][p]][j]);
                }
                MingUtility::printVector(diffs);
                cout << endl;
            }
#endif
        }
    }

    // 4. Initialize the LP solver. (_solver)
    _solver.setVarUpbound(_upbound);
    _solver.addVariables(_K);
#if db2
    cout << "boundaryVs: ";
    MingUtility::printVector(_bdVs);
    for (int p = 0; p < _N; ++p) {
        if (_bdVs[p]) {
            cout << p << " ";
        }
    }
    cout << endl << endl;
#endif
}

// O(N * Klog(K))
void MMLevelSet::OffsetVector2Labeling(const OffsetVector& offsetvector,
                                       Labeling& label) const {
    label.resize(_N);
    float maxVal, tmpVal;
    int maxK;
    for (int p = 0; p < _N; ++p) {
        maxVal = FLT_MIN;
        maxK = -1;
        for (int i = 0; i < _K; ++i) {
            tmpVal = _V2F[p][i] + offsetvector[i];
            if (tmpVal > maxVal) {
                maxVal = tmpVal;
                maxK = i;
            }
        }
        label[p] = maxK;
    }
}

void MMLevelSet::OffsetVector2LabelingForVs(const OffsetVector& offsetvector,
                                            const vector<int>& Vs,
                                            Labeling& label) const {
    label.resize(_N);
    float maxVal, tmpVal;
    int maxK;
    for (auto p : Vs) {
        maxVal = FLT_MIN;
        maxK = -1;
        for (int i = 0; i < _K; ++i) {
            tmpVal = _V2F[p][i] + offsetvector[i];
            if (tmpVal > maxVal) {
                maxVal = tmpVal;
                maxK = i;
            }
        }
        label[p] = maxK;
    }
}

// O(N * K^2)
void MMLevelSet::Labeling2Range(const Labeling& label, Range& range) const {
    // Range has at most K*(K-1) boundary, each of which is defined by a vertex.
    // For simplicity, here use K*K vector of vertices to represent a Range.
    bool db = false;
#if db1
    // db = true;
    // if (_score == -5) db = true;
#endif
    range.resize(_K * _K, -1);
    for (int i = 0; i < _K; ++i) {
        for (int j = i + 1; j < _K; ++j) {
            int lastI = -1;
            int firstJ = -1;
            for (int p = 0; p < _N; ++p) {
                int vi = _Orders[i][j][p];
                if (label[vi] == i) {
                    lastI = vi;
                } else if (label[vi] == j) {
                    firstJ = vi;
                    // Once firstJ is found, there is no i afterwards.
                    break;
                }
            }
            range[i * _K + j] = lastI;
            range[j * _K + i] = firstJ;
            if (db) {
                cout << "v[" << i + 1 << "][" << j + 1 << "] = " << lastI + 1 << endl;
                cout << "v[" << j + 1 << "][" << i + 1 << "] = " << firstJ + 1 << endl;
            }
        }
    }
}
// O(K^3 + LP)
// Fomulate the Range checking / Range2Offset to a LP problem and solve it.
bool MMLevelSet::SolveRangeWithLP(const Range& range) {
    _solver.clearConstraints();
    for (int i = 0; i < _K; ++i) {
        for (int j = 0; j < _K; ++j) {
            if (i == j) continue;
            int vij = range[i * _K + j];
            if (vij != -1) {
                // // new constraint found
                // // V2F[vij][i] + si > V2F[vij][j] + sj
                // // => si - sj > V2F[vij][j] - V2F[vij][i]
                _solver.addConstraint(i, j, _V2F[vij][j] - _V2F[vij][i]);
            }
        }
    }
    return _solver.solve();
}
// O(K^3 + LP)
// Fomulate the Range checking / Range2Offset to a LP problem and solve it.
bool MMLevelSet::SolveRangeWithLP(const Range& range, LPSolver& solver) {
    solver.addVariables(_K);
    for (int i = 0; i < _K; ++i) {
        for (int j = 0; j < _K; ++j) {
            if (i == j) continue;
            int vij = range[i * _K + j];
            if (vij != -1) {
                _solver.addConstraint(i, j, _V2F[vij][j] - _V2F[vij][i]);
            }
        }
    }
    return solver.solve();
}
// O(K^3 + LP)
// Convert a range to an offset, using LP.
void MMLevelSet::Range2Offset(const Range& range, OffsetVector& offset) {
#if REUSESOLVER
    if (!SolveRangeWithLP(range)) {
        cout << "DEBUG! This should not happen!" << endl;
        return;
    }
    _solver.getResults(offset);
#else
    LPSolver solver;
    if (!SolveRangeWithLP(range, solver)) {
        cout << "DEBUG! This should not happen!" << endl;
        return;
    }
    solver.getResults(offset);
#endif
}

// O(K^3 + LP) + O(N * Klog(K))
// Convert a range to a labeling.
void MMLevelSet::Range2Labeling(const Range& range, Labeling& label) {
    OffsetVector curOffset;
    Range2Offset(range, curOffset);
    if (curOffset.size() != 0) {
        OffsetVector2Labeling(curOffset, label);
    }
}

// O(K^3 + LP)
// Use a LP library to check this.
bool MMLevelSet::IsValidRange(const Range& range) {
#if REUSESOLVER
    return SolveRangeWithLP(range);
#else
    LPSolver solver;
    return SolveRangeWithLP(range, solver);
#endif
}
bool MMLevelSet::IsBoundaryValid(const Range& range, int mat0, int mat1) {
    // Debug
    // if (visited[range].first == 6 && mat0 == 1 && mat1 == 2) {
    //   cout << "here" << endl;
    //   MingUtility::printVector(range);
    //   for (int i = 0; i < _K; ++i) {
    //     for (int j = 0; j < _K; ++j) {
    //       if (i == j) continue;
    //       int vij = range[i * _K + j];
    //
    //       if (vij != -1) {
    //         if (i == mat0 && j == mat1) {
    //           cout << "s"<< i << " - s" << j << " == " << _V2F[vij][j] -
    //           _V2F[vij][i] << endl;
    //           cout << "vij = " << vij << endl;
    //         } else {
    //           cout << "s"<< i << " - s" << j << " >= " << _V2F[vij][j] -
    //           _V2F[vij][i] << endl;
    //           cout << "vij = " << vij << endl;
    //         }
    //       }
    //     }
    //   }
    // }
    _solver.clearConstraints();
    for (int i = 0; i < _K; ++i) {
        for (int j = 0; j < _K; ++j) {
            if (i == j) continue;
            int vij = range[i * _K + j];

            if (vij != -1) {
                if (i == mat0 && j == mat1) {
                    _solver.addEqualConstraint(i, j, _V2F[vij][j] - _V2F[vij][i]);
                } else {
                    _solver.addConstraint(i, j, _V2F[vij][j] - _V2F[vij][i]);
                }
            }
        }
    }
    return _solver.solveBoundary();
}
// O(K^3 + LP) + O(N * Klog(K)) + O(K^2 * ((N * K^2) + (K^3 + LP)))
void MMLevelSet::Range2NeighborRanges(const Range& range, vector<Range>& nbs) {
    // Labeling label0;
    // Range2Labeling(range, label0);
    Labeling label0 = visited[range].second;
    // For the K*(K-1) neighbors, flip the label of corresponding vertex, find the
    // range according to the new labeling, and check the emptiness of the range.
    bool db = false;
#if db1
    // db = true;
    // if (_score == -1) db = true;
#endif
    if (db) {
        // cout << "\nscore = " << _score << endl;
        cout << "vid: ";
        MingUtility::printRange(1, _N, 1);
        cout << "lab: ";
        MingUtility::printVector(label0);
    }
    Labeling org = label0;
    // For critical edge analysis.
    int rangeID1 = visited[range].first;
    for (int i = 0; i < _K; ++i) {
        for (int j = 0; j < _K; ++j) {
            if (i == j) continue;
            int vij = range[i * _K + j];
            if (vij != -1) {
                // For boundary vertices, do not change the label.
                // if (_bdVs.find(vij) != _bdVs.end()) {
                if (_bdVs[vij]) {
                    continue;
                }
                // // If boundary i,j doesn't exist, continue.
                if (!IsBoundaryValid(range, i, j)) {
                    continue;
                }
                // Flip.
                label0[vij] = j;
                Range nbRange;
                if (db) {
                    cout << "-\n";
                }
                Labeling2Range(label0, nbRange);
                // Critical Edge analysis.
                if (!IsValidRange(nbRange)) {
                    // Flip back.
                    label0[vij] = i;
                    continue;
                }
                // TODO debug here, check whether the label is valid
                // if (rangeID1==6 && i==1 && j==2) {
                //   cout << "range = ";
                //   MingUtility::printVector(range);
                //
                //   cout << "IsBoundaryValid? " << IsBoundaryValid(range, i, j) <<
                //   endl;
                //   vector<float> tmpRes;
                //   _solver.getResults(tmpRes);
                //   cout << "res=";
                //   MingUtility::printVector(tmpRes);
                //
                //         for (int ii = 0; ii < _K; ++ii) {
                //           for (int jj = ii + 1; jj < _K; ++jj) {
                //             // cout << "order_new[" << ii << "][" << jj << "]  = ";
                //             bool matii = false, matjj = false;
                //             for (int pp = 0; pp < _N; ++pp) {
                //               int vi = _Orders[ii][jj][pp];
                //               if (label0[vi] == ii) matii = true;
                //               if (label0[vi] == jj) matjj = true;
                //               if (matii && matjj && label0[vi] == ii) {
                //                 cout << rangeID1 << ", " << i << ", " << j << endl;
                //               }
                //             }
                //           }
                //         }
                //
                //             cout << "org label = ";
                //             MingUtility::printVector(org);
                //             for (int ii = 0; ii < _K; ++ii) {
                //               for (int jj = ii + 1; jj < _K; ++jj) {
                //                 cout << "order_org[" << ii << "][" << jj << "]  = ";
                //                 for (int pp = 0; pp < _N; ++pp) {
                //                   int vi = _Orders[ii][jj][pp];
                //                   cout << org[vi] << " ";
                //                 }
                //                 cout << endl;
                //               }
                //             }
                //
                //             cout << "  label0  = ";
                //             MingUtility::printVector(label0);
                //             for (int ii = 0; ii < _K; ++ii) {
                //               for (int jj = ii + 1; jj < _K; ++jj) {
                //                 cout << "order_lb0[" << ii << "][" << jj << "]  = ";
                //                 for (int pp = 0; pp < _N; ++pp) {
                //                   int vi = _Orders[ii][jj][pp];
                //                   cout << label0[vi] << " ";
                //                 }
                //                 cout << endl;
                //               }
                //             }
                //
                // }

                bool isOld = true;
                // If nbRange was never seen before, and is valid
                if (visited.find(nbRange) ==
                        visited.end() /*  && IsValidRange(nbRange) */) {
                    isOld = false;
                    visited[nbRange] = make_pair(int(allRanges.size()), label0);
                    allRanges.push_back(nbRange);
                    nbs.push_back(nbRange);
                    if (db) {
                        cout << "i=" << i << ", j=" << j << ", vij=" << vij + 1 << endl;
                        cout << "vid=";
                        MingUtility::printRange(1, _N, 1);
                        cout << "lab=";
                        MingUtility::printVector(label0);

                        for (int ii = 0; ii < _K; ++ii) {
                            for (int jj = 0; jj < _K; ++jj) {
                                cout << "<" << ii + 1 << "," << jj + 1
                                     << "> :" << nbRange[ii * _K + jj] + 1 << endl;
                            }
                        }
                    }
                }
                // Critical Edge analysis.
                int rangeID2 = visited[nbRange].first;
                int id1 = min(rangeID1, rangeID2);
                int id2 = max(rangeID1, rangeID2);
                pair<int, int> curE = make_pair(id1, id2);
                // New edge is found.
                if (edge2id.find(curE) == edge2id.end()) {
                    edge2id[curE] = allEdges.size();
                    allEdges.push_back(curE);
                }

                // if (criticalEdgeAnalysis(label0, vij, i, j)) {
                if (_mesh->isVCritical(vij, label0, i, j)) {
                    criticalEandV.insert(make_pair(edge2id[curE], vij));
                    // #if db2
                    // if (edge2id[curE] == 109) {
                    // if (rangeID1 == 6) {
                    if (false) {
                        // Debug
                        int nPos = 0;
                        Labeling lb1 = visited[allRanges[rangeID1]].second;
                        // Labeling lb2 = visited[allRanges[rangeID2]].second;
                        Labeling lb2 = visited[nbRange].second;
                        for (int i = 0; i < lb1.size(); ++i) {
                            if (lb1[i] != lb2[i]) nPos++;
                        }
                        cout << "nPos = " << nPos << endl;
                        cout << "isOld = " << isOld << endl;
                        cout << "{mat0, mat1} = {" << i << ", " << j << "}" << endl;
                        cout << "vij = " << vij << endl;
                        if (nPos != 1) {
                            cout << "{id1, id2} = " << rangeID1 << ", " << rangeID2 << endl;
                            // cout << "max allRanges id = " << allRanges.size() - 1
                            //      << ", isOld = " << isOld << endl;
                            cout << "label1 = ";
                            MingUtility::printVector(lb1);
                            cout << "label2 = ";
                            MingUtility::printVector(lb2);
                            cout << "lb0(2) = ";
                            MingUtility::printVector(label0);

                            Labeling oldLabel = lb2;
                            Range oldR;
                            Labeling newLabel = label0;
                            Range newR;
                            Labeling2Range(oldLabel, oldR);

                            cout << "                   ";
                            for (int pp = 0; pp < _N; ++pp) {
                                cout << pp % 10 << " ";
                            }
                            cout << endl;

                            for (int ii = 0; ii < _K; ++ii) {
                                for (int jj = ii + 1; jj < _K; ++jj) {
                                    cout << "order_old[" << ii << "][" << jj << "]  = ";
                                    for (int pp = 0; pp < _N; ++pp) {
                                        int vi = _Orders[ii][jj][pp];
                                        cout << oldLabel[vi] << " ";
                                    }
                                    cout << endl;
                                }
                            }
                            Labeling2Range(newLabel, newR);
                            cout << "oldLabel = ";
                            MingUtility::printVector(oldLabel);
                            cout << "oldR = ";
                            MingUtility::printVector(oldR);

                            cout << "                   ";
                            for (int pp = 0; pp < _N; ++pp) {
                                cout << pp % 10 << " ";
                            }
                            cout << endl;
                            for (int ii = 0; ii < _K; ++ii) {
                                for (int jj = ii + 1; jj < _K; ++jj) {
                                    cout << "order_new[" << ii << "][" << jj << "]  = ";
                                    for (int pp = 0; pp < _N; ++pp) {
                                        int vi = _Orders[ii][jj][pp];
                                        cout << newLabel[vi] << " ";
                                    }
                                    cout << endl;
                                }
                            }

                            cout << "newLabel = ";
                            MingUtility::printVector(newLabel);
                            cout << "newR = ";
                            MingUtility::printVector(newR);
                            cout << endl;

                            cout << "                   ";
                            for (int pp = 0; pp < _N; ++pp) {
                                cout << pp % 10 << " ";
                            }
                            cout << endl;
                            for (int ii = 0; ii < _K; ++ii) {
                                for (int jj = ii + 1; jj < _K; ++jj) {
                                    cout << "order_rg1[" << ii << "][" << jj << "]  = ";
                                    for (int pp = 0; pp < _N; ++pp) {
                                        int vi = _Orders[ii][jj][pp];
                                        cout << lb1[vi] << " ";
                                    }
                                    cout << endl;
                                }
                            }

                            // cout << "newLabel = ";
                            // MingUtility::printVector(newLabel);
                            // cout << "newR = ";
                            // MingUtility::printVector(newR);
                            // cout << endl;
                        }
                    }
                    // #endif
                }

                if (db) {
                    cout << "--\n";
                }
                // Flip back.
                label0[vij] = i;
            }
        }
    }
}

// O(N * Klog(K)) + O(N * K^2)
void MMLevelSet::OffsetVector2Range(const OffsetVector& offsetvector,
                                    Range& range) const {
    Labeling label;
    OffsetVector2Labeling(offsetvector, label);
    Labeling2Range(label, range);
}

// TODO: given a range, calculate the score and also the topology here together.
// score = sum( Log( V2F[p][label of p] ) );
float MMLevelSet::Range2Score(const Range& range) {
    // _score -= 1.0;
    // return _score
    const Labeling& label = visited[range].second;
    return Labeling2Score(label);
}
float MMLevelSet::OffsetVector2Score(const OffsetVector& offset) {
    Labeling label;
    OffsetVector2Labeling(offset, label);
    return Labeling2Score(label);
}

float MMLevelSet::Labeling2Score(const Labeling& label) {
    float score = 0.0;
    for (int p = 0; p < _N; ++p) {
        // score += _V2LogF[p][label[p]];
        score += _V2LogF[p][label[p]] * _mesh->_vNbWeights[p];
    }
    return score;
}

void MMLevelSet::ExploreTopologies(const OffsetVector& seedOffsetVector,
                                   int maxRangeToExplore, int maxTopoToFind,
                                   vector<Range>& topologies) {
    // // Pre processing to generate the helper data structure, e.g. Orders.
    // PreProcess();
    // cout << "1" << endl;

    // The ranges on the fringe of exploration, sorted by the score of the range.
    priority_queue<pair<float, Range>> pq;

    // Highest-Score-First Search.
    Range seedRange;
    OffsetVector2Range(seedOffsetVector, seedRange);
    Labeling seedLabeling;
    Range2Labeling(seedRange, seedLabeling);
    visited[seedRange] = make_pair(int(allRanges.size()), seedLabeling);
    allRanges.push_back(seedRange);
    pq.push(make_pair(Range2Score(seedRange), seedRange));
    // cout << "2" << endl;
    // TODO: explore topology, how to represent a unique topo?
    // topologies.push_back(seedRange);
    int counter = 0;

    // Debug
    // cout << "  seed label  = ";
    // MingUtility::printVector(seedLabeling);
    // for (int ii = 0; ii < _K; ++ii) {
    //   for (int jj = ii + 1; jj < _K; ++jj) {
    //     cout << "order_seed[" << ii << "][" << jj << "]  = ";
    //     for (int pp = 0; pp < _N; ++pp) {
    //       int vi = _Orders[ii][jj][pp];
    //       cout << seedLabeling[vi] << " ";
    //     }
    //     cout << endl;
    //   }
    // }

    while (!pq.empty() && counter < maxRangeToExplore &&
           topologies.size() < maxTopoToFind) {
        ++counter;
        pair<float, Range> cur = pq.top();
        pq.pop();
        vector<Range> nbs;
        Range2NeighborRanges(cur.second, nbs);
        // topologies.push_back(cur.second);
        for (const auto& nb : nbs) {
            pq.push(make_pair(Range2Score(nb), nb));
        }
    }
    // cout << "3" << endl;
}
void MMLevelSet::clearupVectorSpace(){



    _V2F.clear();
    _V2LogF.clear();

    _Orders.clear();
    _bdVs.clear();
    _rays.clear();

    // Critical points on each ray.
    _criticalPtsPerRay.clear();

    // Topology points on each ray.
    _topologyPtsPerRay.clear();

    // Unique cell topologies.
    _cellTopologies.clear();
    // For debug purpose.
    _cellTopo2topoPtID.clear();
    _cellTopo2AllTopoPtID.clear();
    _badCellTopoPtID.clear();

    _cellNexplored2Topogroup.clear();

    _V2nbVs.clear();

    allRanges.clear();
    allEdges.clear();
    edge2id.clear();
    criticalEandV.clear();

    _initLabel.clear();



}

/*********************************************************************************************/
/*********************************************************************************************/


void MMLevelSet::OutputCellBoundary(vector<double>&Vs,vector<uint>&Fs,vector<int>&VMats){

    vector<int>Vori_Ind;
    _mesh->OutputCellBoundary(Vs,Fs,Vori_Ind);

    VMats.clear();
    for(auto a:Vori_Ind)VMats.push_back(_initLabel[a]);


}

void MMLevelSet::CutTet(vector<double>&v2planeDist,vector<double>&outVpos,vector<uint>&outF2V){


    _mesh->CutTet(v2planeDist,outVpos,outF2V);



}

extern ParallelArrangement pArr;
void CallTetgenForSingleSurface(vector<double>&vpos,  vector<uint>&fs,  vector<vector<float>>&out_vpos, vector<vector<int>>&out_tet);
int MMLevelSet::AddConstraintPointInsideCell(vector<vector<float>> &aVs,vector<int>&aVMs,string outdir,int ci){


    vector<int> Vori_Ind = CellBoundary_VoriInd;
    vector<double>Vs = CellBoundary_Vs;
    vector<uint>Fs = CellBoundary_Fs;
    vector<int>outVmarkers = CellBoundary_Vmarkers;


    vector<vector<float>>out_vpos; vector<vector<int>>out_tet;

    //9.215060397484454  13.10298581614663  30.01930571426029

//    for(int i=0;i<20;++i){
//        Vs.push_back(19);Vs.push_back(13);Vs.push_back(10+i);
//    }
//    for(int i=0;i<20;++i){
//        Vs.push_back(0);Vs.push_back(-3+0.3*i);Vs.push_back(0);
//    }
    //aVMs.clear();
    //aVMs.resize(20,1);
    for(auto &a:aVs)for(auto b:a)Vs.push_back(b);

    CallTetgenForSingleSurface(Vs, Fs, out_vpos,out_tet);

    int obsize = Vori_Ind.size();

    _N = out_vpos.size();
    auto orilable = CellBoundary_OriLabel;


    outVmarkers.resize(_N);
    for(int i=Vori_Ind.size();i<_N;++i)outVmarkers[i] = 0;

    clearupVectorSpace();//there is something not welly clearup before use again, so I just clear all vector spaces.
    if(_mesh!=NULL)delete _mesh;
    _mesh = new TetMesh();
    _mesh->setDimension(3);
    _mesh->Load(out_vpos,outVmarkers, out_tet);




    _initLabel.resize(_N);
    for(int i=0;i<Vori_Ind.size();++i)_initLabel[i] = orilable[Vori_Ind[i]];
    for(int i=0;i<aVMs.size();++i)_initLabel[i+obsize] = aVMs[i];
    vector<int>fixv;
    for(int i=0;i<Vori_Ind.size();++i)fixv.push_back(i);
    for(int i=0;i<aVMs.size();++i)fixv.push_back(i+obsize);
    //_mesh->GenerateIndicatorFunction(nMat, cellInitLabel, _V2F);
    cout<<"1"<<endl;
    _mesh->GenerateIndicatorFunction(_K, _initLabel,fixv, _V2F);
    cout<<"2"<<endl;
    PreProcessMeshBD();
    cout<<"3"<<endl;


    PreProcess();

    cout<<"4"<<endl;

    clock_t start = clock(), end;
    start = clock();
    ExploreTopologiesAlongRays(false, 2, 0.0f);//ray shooting
    end = clock();
    float time = float(end - start) / CLOCKS_PER_SEC;
    cout << "time: " << time << endl;

    cout<<"5"<<endl;

    int nTopo = _cellTopologies.size();

    if(nTopo>500){
        cout<<"more than 500 topo!"<<endl;
        exit(500);
    }
    Labeling topolables;

    vector<vector<int>>new_bdLoopE2V(_mesh->_bdLoopsV.size());
    for(int i=0;i<_mesh->_bdLoopsV.size();++i){
        auto &nb = new_bdLoopE2V[i];
        for(auto a:_mesh->_bdLoopsV[i])nb.push_back(Vori_Ind[a]);
    }


    pArr.MapLoopToSegs(ci,new_bdLoopE2V);

    vector<int>new_bdActiveF2Vs;
    for(auto a:_mesh->_bdActiveF2Vs)new_bdActiveF2Vs.push_back(Vori_Ind[a]);
    vector<int>new_bdActiveE2Vs_fornonplcs;
    for(auto a:_mesh->_bdActiveE2Vs_fornonplcs)new_bdActiveE2Vs_fornonplcs.push_back(Vori_Ind[a]);


    pArr.MapActiveFCP(ci,new_bdActiveF2Vs,_mesh->_bdCurveVs,new_bdActiveE2Vs_fornonplcs,false);
    pArr.GetCellTopo(ci,_cellTopologies,_mesh->_label2bdLoops,_cellNexplored2Topogroup);


    auto &mappingMat = pArr._cell2labels[ci];
    for (int i = 0; i < nTopo; ++i) {

        vector<vector<float>> sufVs;
        vector<vector<int>> sufFs;
        vector<vector<int>> sufFMats;
        OffsetVector2Labeling(_cellTopologies[i]._offset, topolables);
        _mesh->GenerateLevelSetFromLabeling(topolables,
                                            sufVs, sufFs, sufFMats);

        cout<<i<<' '<<nTopo<<' '<<_cellTopologies[i]._junctionpoints<<' '<<_cellTopologies[i]._score<<endl;
        if (true) {
            string outfile = outdir + "/suf_" + to_string(ci) + "_" + to_string(i) + ".suf";
            _mesh->WriteLevelSet(sufVs, sufFs, sufFMats,mappingMat, outfile.c_str());
        }
    }



    pArr.LoopMergingDP();
    //vector<int>pickTopo({0,0,1});
    //pArr.DynamicProgramming();

    return nTopo;



}

void MMLevelSet::ReadTopologyFromFile(ifstream &ifs){

    auto readVector = [&ifs](vector<int>&vec){
        int num;
        ifs.read((char*)(&num),sizeof(int));
        vec.resize(num);
        ifs.read((char*)(vec.data()),sizeof(int)*vec.size());

    };

    readVector(_cellNexplored2Topogroup);
    int num;
    ifs.read((char*)(&num),sizeof(int));
    _cellTopologies.resize(num);
    for(int i=0;i<num;++i)_cellTopologies[i].readbinary(ifs);

}
void MMLevelSet::WriteTopologyToFile(ofstream &ofs){

    auto writeVector = [&ofs](vector<int>&vec){
        int num;
        num = vec.size();
        ofs.write((char*)(&num),sizeof(int));

        ofs.write((char*)(vec.data()),sizeof(int)*vec.size());

    };
    writeVector(_cellNexplored2Topogroup);
    int num = _cellTopologies.size();
    ofs.write((char*)(&num),sizeof(int));
    for(int i=0;i<num;++i)_cellTopologies[i].writebinary(ofs);


}

