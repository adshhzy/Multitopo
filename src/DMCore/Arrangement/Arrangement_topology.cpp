#include "Arrangement.h"

void Arrangement::explorTopologyInCell(int ci, bitIntHash& allSets_b){

	unordered_set<int> cellV, segV;
	unordered_map<int, int> segV2cyci;
	unordered_map<int, int> segV2segi;
	unordered_map<int, int> cellV2tetV, segV2tetV;
	vector<int> labelInVs, labelOutVs, labelXVs;
	vector<vector<int> > tetV2nbTet;
	vector<int> tetV2cyci;
	vector<queue<int> > isoSufFringe; 
	vector<float> VPro;
	vector<float> VProVol;
	UnDirGraph graph;
	vector<float> criticalIsoVal;
	/* tetrahedralize each cell with tetgen lib */
	int v1, v2;
	tetgenio tetin, tetout;
	tetgenio tetaddin;
	for (int fi : cell2F[ci]){
		for (int ei : cellFs[fi]){
			v1 = cellEs[ei][0];
			v2 = cellEs[ei][1];
			cellV.insert(v1);
			cellV.insert(v2);
		}
	}
	
	for (int cyci = 0; cyci < cell2cyc[ci].size(); ++cyci){
		for (int si : cyc2seg[cell2cyc[ci][cyci]]){
			for (int vi : seg2V[si]){
				segV.insert(vi);
				if(segV2segi[vi]==0){
					segV2segi[vi] = TagMapping::getSegiTag(si);
				} else {
					segV2segi[vi] = (TagMapping::getSegiTag(si) + segV2segi[vi])/2;
				}
			}
		}
	}

	// (1) Tetdrahedralize cell ci with TETGEN
	tetCell_tetgen(ci, tetin, tetout, tetaddin, cellV, segV, segV2segi, cellV2tetV, segV2tetV);

	// (2) Generate data structures for the tets for future use
	tetCell_genCellTet(ci, tetout);

	// (3) Find components & split edges & calc volume & indentify in and out among the boundary Vs
	tetCell_findComponent(ci, tetout, segV2tetV, graph, tetV2nbTet, labelInVs, labelOutVs, labelXVs, isoSufFringe, tetV2cyci);
	tetCell_volWeight(ci, VProVol, labelXVs);

	// (4) Calculate the implicite function with Random Walk
#if VAXMAN
	tetCell_Vaxman(ci, labelInVs, labelOutVs, labelXVs, VPro);
#else
	tetCell_randomWalk(graph, labelInVs, labelOutVs, labelXVs, VPro, ci);
#endif

	// (5) Find critical iso values, which covers all distinct topologies
	tetCell_criticalPoint(tetV2nbTet, VPro, tetout, labelXVs, criticalIsoVal, ci);

	// (6) Find unique groupings at each critical point
	tetCell_findGroupings(ci, VPro, VProVol, labelXVs, graph, criticalIsoVal, isoSufFringe, tetV2cyci, allSets_b);
}

void Arrangement::tetCell_tetgen(int ci, tetgenio &tetin, tetgenio &tetout, tetgenio &tetaddin, unordered_set<int> &cellV,  unordered_set<int> &segV,unordered_map<int, int> &segV2segi, unordered_map<int, int> &cellV2tetV,  unordered_map<int, int> &segV2tetV){

	int v1, v2, v=0;
	tetgenio::facet *f;
	tetgenio::polygon *p;

	// prepare for tetgen datastructure
	tetin.deinitialize(); tetin.initialize();
	tetout.deinitialize(); tetout.initialize();
	tetaddin.deinitialize(); tetaddin.initialize();
	tetin.firstnumber = 0; tetout.firstnumber = 0; tetaddin.firstnumber = 0; 

	// point coordinates
	tetin.numberofpoints = cellV.size() + segV.size();
	tetin.pointlist = new REAL[tetin.numberofpoints * 3];
	tetin.pointmarkerlist = new int[tetin.numberofpoints]; 
	v = 0;

	for (int vi : cellV){
		tetin.pointlist[v * 3] = cellVs[vi][0];
		tetin.pointlist[v * 3 + 1] = cellVs[vi][1];
		tetin.pointlist[v * 3 + 2] = cellVs[vi][2];
		tetin.pointmarkerlist[v] = TagMapping::tetCellVTag(); 
		cellV2tetV[vi] = v++;

	}

	for (int vi : segV){
		tetin.pointlist[v * 3] = Vs[vi][0];
		tetin.pointlist[v * 3 + 1] = Vs[vi][1];
		tetin.pointlist[v * 3 + 2] = Vs[vi][2];
		tetin.pointmarkerlist[v] = TagMapping::tetSegVTag(segV2segi[vi]); 
		segV2tetV[vi] = v++;
	}

	// facet list
	tetin.numberoffacets = cell2F[ci].size();
	tetin.facetlist = new tetgenio::facet[tetin.numberoffacets];
	tetin.facetmarkerlist = new int[tetin.numberoffacets];
	for (int fi = 0; fi < tetin.numberoffacets; ++fi){
		tetin.facetmarkerlist[fi] = TagMapping::tetFaceTag(fi);
		f = &tetin.facetlist[fi];
		f->numberofholes = 0;
		f->holelist = NULL;

		f->numberofpolygons = 1;
		for (int si : F2seg[cell2F[ci][fi]]){
			f->numberofpolygons += (seg2V[si].size() - 1);
		}
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		int pi = 0;
		// for seg polys
		for (int si : F2seg[cell2F[ci][fi]]){
			for(int vi = 0; vi < seg2V[si].size() - 1; ++vi){
				p = &f->polygonlist[pi++];
				p->numberofvertices = 2;
				p->vertexlist = new int[p->numberofvertices];
				p->vertexlist[0] = segV2tetV[seg2V[si][vi]];
				p->vertexlist[1] = segV2tetV[seg2V[si][vi+1]];
			}
		}
		// for face poly
		int lastendv = -1;
		vector<int> tmpverlist;
		for (int ei : cellFs[cell2F[ci][fi]]){
			v1 = cellEs[ei][0];
			v2 = cellEs[ei][1];
			if (v2 != lastendv){
				tmpverlist.push_back(cellV2tetV[v1]);
				for (auto it = sorted_verIdOnCellE[ei].begin(); it != sorted_verIdOnCellE[ei].end(); ++it){
					tmpverlist.push_back(segV2tetV[*it]);
				}
				lastendv = v2;
			}
			//reverse
			else{
				tmpverlist.push_back(cellV2tetV[v2]);
				for (auto it = sorted_verIdOnCellE[ei].rbegin(); it != sorted_verIdOnCellE[ei].rend(); ++it){
					tmpverlist.push_back(segV2tetV[*it]);
				}
				lastendv = v1;
			}
		}
		p = &f->polygonlist[pi++];
		p->numberofvertices = tmpverlist.size();
		p->vertexlist = new int[p->numberofvertices];
		for (int vi = 0; vi < p->numberofvertices; ++vi){
			p->vertexlist[vi] = tmpverlist[vi];
		}
	}
	
#if _DEBUG_MING_TETGEN
	Utility::saveTetgenio((t_savedir+"myTetIn_"+Utility::toStr(ci)+".txt").c_str(),tetin);

	tetin.save_nodes(const_cast<char *>((t_savedir+"test_in_"+Utility::toStr(ci)).c_str()));
	tetin.save_poly(const_cast<char *>((t_savedir+"test_in_"+Utility::toStr(ci)).c_str()));
#endif 

	// Q: quiet | pq1.1: CDT | a0.05: volumn constraint
	string teta = Utility::toStr(t_teta);
	string tetpq = "QnT1e-5pq1.4ea"; 
#if _CG_LOG
	cout << "Tetgen = " << tetpq + teta << endl;
	cout << "cell2cyc[ci].size() = " << cell2cyc[ci].size() << endl;
#endif
	cout << "-> Calling Tetgen in cell[" << ci << "]... ";
	tetrahedralize(const_cast<char *>((tetpq+teta).c_str()), &tetin, &tetout);
	cout << "DONE. " << endl;

#if _DEBUG_MING_TETGEN
	tetout.save_elements(const_cast<char *>((t_savedir+"test_out_"+Utility::toStr(ci)).c_str()));
	tetout.save_faces(const_cast<char *>((t_savedir+"test_out_"+Utility::toStr(ci)).c_str()));
	tetout.save_edges(const_cast<char *>((t_savedir+"test_out_"+Utility::toStr(ci)).c_str()));
	tetout.save_nodes(const_cast<char *>((t_savedir+"test_out_"+Utility::toStr(ci)).c_str()));//maybe overlapped by tetCell_findComponent()
	tetout.save_neighbors(const_cast<char *>((t_savedir+"test_out_"+Utility::toStr(ci)).c_str()));
#endif
}

void Arrangement::tetCell_genCellTet(int ci, tetgenio &tetout){
	unordered_map<vector<int>, int, myVecHash> faceHash;
	unordered_map<vector<int>, int, myVecHash> edgeHash;
	int cnTet=0, cnFace=0, cnEdge=0, cnV=0;
	cnTet = tetout.numberoftetrahedra;
	cnV = tetout.numberofpoints;
	cnFace = 2*cnTet + tetout.numberoftrifaces/2;
	allCellTets.cells[ci].init_tet(cnTet);
	allCellTets.cells[ci].init_v(cnV, tetout.pointlist);
	allCellTets.cells[ci].init_face(cnFace);

	vector<int> oneF2E(3);
	int faceVID [4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };
	int edgeVID [3][2] = { {0,1}, {0,2}, {1,2} };
	vector<int> oneFace(3);
	vector<int> oneEdge(2);

	int faceID = 0;
	vector<vector<int> > EnbF;
	// find cnFace & cnEdge
	for (int ti=0; ti<cnTet; ++ti){
		for (int j=0; j<4; ++j){
			// cellVertices.nbTet & cellTets & cellFace
			allCellTets.cells[ci].cellVertices[tetout.tetrahedronlist[ti*4+j]].nbTet.push_back(ti);
			allCellTets.cells[ci].cellTets[ti].v[j] = (tetout.tetrahedronlist[ti*4+j]);
			allCellTets.cells[ci].cellTets[ti].nbTet[j] = (tetout.neighborlist[ti*4+j]);
			oneFace[0] = tetout.tetrahedronlist[ti*4+faceVID[j][0]];
			oneFace[1] = tetout.tetrahedronlist[ti*4+faceVID[j][1]];
			oneFace[2] = tetout.tetrahedronlist[ti*4+faceVID[j][2]];
			sort(oneFace.begin(), oneFace.end());
			// a new face is found
			int curFaceID = -1;
			if (faceHash.find(oneFace)==faceHash.end()){
				curFaceID = faceID;
				allCellTets.cells[ci].cellTets[ti].face[j] = curFaceID; 
				allCellTets.cells[ci].cellFaces[faceID].v[0] = oneFace[0];
				allCellTets.cells[ci].cellFaces[faceID].v[1] = oneFace[1];
				allCellTets.cells[ci].cellFaces[faceID].v[2] = oneFace[2];
				for (int k=0; k<3; ++k){
					oneEdge[0]=oneFace[edgeVID[k][0]]; 
					oneEdge[1]=oneFace[edgeVID[k][1]];
					int curEdgeID = -1;
					if (edgeHash.find(oneEdge)==edgeHash.end()){
						// a new edge is found
						curEdgeID = cnEdge;
						vector<int> oneEnbF;
						EnbF.push_back(oneEnbF);
						edgeHash[oneEdge] = cnEdge++;
					} else {
						curEdgeID = edgeHash[oneEdge];
					}
					oneF2E[k] = curEdgeID;
					EnbF[curEdgeID].push_back(curFaceID);
				}
				allCellTets.cells[ci].cellFaces[faceID].edge[0] = oneF2E[0];
				allCellTets.cells[ci].cellFaces[faceID].edge[1] = oneF2E[1];
				allCellTets.cells[ci].cellFaces[faceID].edge[2] = oneF2E[2];
				faceHash[oneFace] = faceID++;
			} else {
				curFaceID = faceHash[oneFace];
				allCellTets.cells[ci].cellTets[ti].face[j] = curFaceID;
			}
			allCellTets.cells[ci].cellFaces[curFaceID].nbTet.push_back(ti);

		}
	}

	allCellTets.cells[ci].init_edge(cnEdge);
	for(auto it=edgeHash.begin(); it!=edgeHash.end(); ++it){
		allCellTets.cells[ci].cellEdges[it->second].v[0] = (it->first)[0];
		allCellTets.cells[ci].cellEdges[it->second].v[1] = (it->first)[1];
		allCellTets.cells[ci].cellEdges[it->second].nbFace = EnbF[it->second];
		allCellTets.cells[ci].cellVertices[it->first[0]].nbV.push_back(it->first[1]);
		allCellTets.cells[ci].cellVertices[it->first[0]].nbE.push_back(it->second);
		allCellTets.cells[ci].cellVertices[it->first[1]].nbV.push_back(it->first[0]);
		allCellTets.cells[ci].cellVertices[it->first[1]].nbE.push_back(it->second);
	}

	for (int vi=0; vi<tetout.numberofpoints; ++vi){
		allCellTets.cells[ci].cellVertices[vi].xyz[0] = (float)(tetout.pointlist[vi*3]);
		allCellTets.cells[ci].cellVertices[vi].xyz[1] = (float)(tetout.pointlist[vi*3+1]);
		allCellTets.cells[ci].cellVertices[vi].xyz[2] = (float)(tetout.pointlist[vi*3+2]);
		allCellTets.cells[ci].cellVertices[vi].VtetMark = tetout.pointmarkerlist[vi];
	}

	allCellTets.cells[ci].cellBBFaces.reserve(tetout.numberoftrifaces);
	BigInt BBEdge_bi(allCellTets.cells[ci].cellEdges.size());
	for (int fi=0; fi<tetout.numberoftrifaces; ++fi){
		oneFace[0] = tetout.trifacelist[fi*3];
		oneFace[1] = tetout.trifacelist[fi*3+1];
		oneFace[2] = tetout.trifacelist[fi*3+2];
		sort(oneFace.begin(), oneFace.end());
		int curFaceID = faceHash[oneFace];
		allCellTets.cells[ci].cellBBFaces.push_back(curFaceID);
		allCellTets.cells[ci].cellFaces[curFaceID].FtetMark = tetout.trifacemarkerlist[fi];
		for (int eind=0; eind<3; ++eind){
			int curEdgeID = allCellTets.cells[ci].cellFaces[curFaceID].edge[eind]; 
			allCellTets.cells[ci].cellEdges[curEdgeID].EtetMark = TagMapping::E_FACE;
			BBEdge_bi.set(curEdgeID);
		}
	}
	allCellTets.cells[ci].cellBBEdges = BBEdge_bi.getOnesInd();
	for (int ei=0; ei<tetout.numberofedges; ++ei){
		oneEdge[0] = tetout.edgelist[ei*2];
		oneEdge[1] = tetout.edgelist[ei*2+1];
		sort(oneEdge.begin(), oneEdge.end());
		int curEdgeID = edgeHash[oneEdge];
		allCellTets.cells[ci].cellEdges[curEdgeID].EtetMark = TagMapping::E_BD; 
	}
}

void Arrangement::tetCell_findComponent(int ci, tetgenio &tetout, unordered_map<int, int> &segV2tetV, UnDirGraph &graph, vector<vector<int> > &tetV2nbTet, vector<int> &labelInVs, vector<int> &labelOutVs, vector<int> &labelXVs, vector<queue<int> > &isoSufFringe, vector<int> &tetV2cyci){
	tetV[ci] = vector<double>(tetout.pointlist, tetout.pointlist + tetout.numberofpoints * 3);
	tetQ[ci] = vector<int>(tetout.tetrahedronlist, tetout.tetrahedronlist + tetout.numberoftetrahedra * 4);

	graph.resize(allCellTets.cells[ci].cellVertices.size());

	vector<int> nxtSegVID(allCellTets.cells[ci].cellVertices.size(), TagMapping::NXTSEGV_NAN);
	isoSufFringe.resize(cell2cyc[ci].size());
	tetV2cyci.resize(allCellTets.cells[ci].cellVertices.size(), -1);
	for (int cyci = 0; cyci < cell2cyc[ci].size(); ++cyci){
		for (int si : cyc2seg[cell2cyc[ci][cyci]]){
			// find second point along the seg direction
			int v_start = segV2tetV[seg2V[si][0]];
			int v_dir = segV2tetV[seg2V[si][1]];
			int v_end = segV2tetV[seg2V[si].back()];

			vector<int> oneSeg2TetV;

			int v_nxt = v_dir; float maxDir = -INFINITY;
			int e_nxt = -1;
Eigen::Vector3f v_s(allCellTets.cells[ci].cellVertices[v_start].xyz[0],allCellTets.cells[ci].cellVertices[v_start].xyz[1],allCellTets.cells[ci].cellVertices[v_start].xyz[2]);
			Eigen::Vector3f v_d(allCellTets.cells[ci].cellVertices[v_dir].xyz[0],allCellTets.cells[ci].cellVertices[v_dir].xyz[1],allCellTets.cells[ci].cellVertices[v_dir].xyz[2]);

			// find an edge steming from v_start that follows the direction of the seg, which is important because we recorded the left/right material with fixed seg orientation
			for (int nbe : allCellTets.cells[ci].cellVertices[v_start].nbE){
				int nbv = allCellTets.cells[ci].cellEdges[nbe].getOtherV(v_start);

				if (nbv == v_dir){
					v_nxt = nbv;
					e_nxt = nbe;
					break;
				}
				Eigen::Vector3f v_nb(allCellTets.cells[ci].cellVertices[nbv].xyz[0],allCellTets.cells[ci].cellVertices[nbv].xyz[1],allCellTets.cells[ci].cellVertices[nbv].xyz[2]);
				float dir = (v_d-v_s).dot((v_nb-v_s).normalized());

				if (dir > maxDir){
					v_nxt = nbv;
					e_nxt = nbe;
					maxDir = dir;
				}
			}

			//v_start:		the 1st pt of the seg; 
			//v_dir:		the 2nd pt of the seg; 
			//v_nxt:		the steiner pt after v_start along v_dir direction
			// trace out the seg, change v on the seg to be -10-cyci
			int v_cur = v_start, v_pre;
			tetV2cyci[v_start] = cyci;

			// NOTE that, a start V on a seg might be shared by multiple segs, so, we can't set nxtSegVID[v_start] = v_nxt, I use SHARE to didentify that. Otherwise, labeling in/out for the first edge of the segment might went wrong
			nxtSegVID[v_start] = TagMapping::NXTSEGV_SHARE ;
			allCellTets.cells[ci].cellEdges[e_nxt].EtetMark = TagMapping::tetSegETag(si);
			oneSeg2TetV.push_back(v_start);
			while (v_nxt != v_end){
				oneSeg2TetV.push_back(v_nxt);
				tetout.pointmarkerlist[v_nxt] = TagMapping::tetSegVTag(TagMapping::getSegiTag(si));
				allCellTets.cells[ci].cellVertices[v_nxt].VtetMark = TagMapping::tetSegVTag(TagMapping::getSegiTag(si));
				v_pre = v_cur;
				v_cur = v_nxt;

				// there are two and exactly two bd/seg edge on point v_cur now
				e_nxt = -e_nxt;
				for (int nbe : allCellTets.cells[ci].cellVertices[v_cur].nbE){
					if ((allCellTets.cells[ci].cellEdges[nbe].EtetMark!=TagMapping::E_INTERIOR)&&(allCellTets.cells[ci].cellEdges[nbe].EtetMark!=TagMapping::E_FACE)){
						e_nxt += nbe;
					}
				}
				v_nxt = allCellTets.cells[ci].cellEdges[e_nxt].getOtherV(v_cur);
				nxtSegVID[v_cur] = v_nxt;
				allCellTets.cells[ci].cellEdges[e_nxt].EtetMark = TagMapping::tetSegETag(si);
				tetV2cyci[v_cur] = cyci;
			}
			oneSeg2TetV.push_back(v_end);
			tetV2cyci[v_end] = cyci;
			isoSufFringe[cyci].push(v_start);

			allCellTets.cells[ci].cellSeg2TetV.push_back(oneSeg2TetV);
		}
	}

	tetVmarker[ci] = vector<int>(tetout.pointmarkerlist, tetout.pointmarkerlist + tetout.numberofpoints);
	allCellTets.cells[ci].init_Vmark(tetout.numberofpoints, tetout.pointmarkerlist);


	// Step2: find the dual graph for BBox faces and edges
	//
	// dural graph
	UnDirGraph bbGraph;
	bbGraph.resize(allCellTets.cells[ci].cellBBFaces.size());
	// by default, the mat is set to be -1, later they will all be replaced by in/out mat mark
	vector<int> boundaryMats(bbGraph.size(), TagMapping::MAT_NAN);
	// triangle mapping, to reduce search within bbFaces
	vector<int> F2BBFInd(allCellTets.cells[ci].cellFaces.size(), -1);
	for (size_t fi=0; fi<allCellTets.cells[ci].cellBBFaces.size(); ++fi){
		F2BBFInd[allCellTets.cells[ci].cellBBFaces[fi]] = fi;
	}
	vector<int> nbBBFace;
	for (auto ei : allCellTets.cells[ci].cellBBEdges){
		nbBBFace.clear();
		// find two nb BBFace, if there are more than 2, error!
		for (auto fi : allCellTets.cells[ci].cellEdges[ei].nbFace){
			// found a BBFace
			if(TagMapping::isFaceOnBB(allCellTets.cells[ci].cellFaces[fi].FtetMark)){
				// debug
				if(nbBBFace.size()==2){cout << "MING! BUG: in nbBBFace! too many..." << endl;}
				nbBBFace.push_back(fi);
				if(nbBBFace.size()==2){
					break;
				}
			}
		}
		// debug
		if(nbBBFace.size()!=2){cout << "MING! BUG: in nbBBFace! too few..." << endl;}


		// for edges on seg, do not add dual edge into bbGraph, but label nb 2 triangle in/out
		if (TagMapping::isEdgeOnSeg(allCellTets.cells[ci].cellEdges[ei].EtetMark)){
			int mat0, mat1;
			mat0 = boundaryMats[F2BBFInd[nbBBFace[0]]];
			mat1 = boundaryMats[F2BBFInd[nbBBFace[1]]];
			if (mat0==TagMapping::MAT_NAN){
				if (mat1==TagMapping::MAT_NAN){ 
					mat0 = tetCell_getMatOfFaceNextToSeg(ci, nbBBFace[0], ei, nxtSegVID);
					boundaryMats[F2BBFInd[nbBBFace[0]]] = mat0;
					mat1 = TagMapping::getOtherMat(mat0, matMarks);
					boundaryMats[F2BBFInd[nbBBFace[1]]] = mat1;
				} else {
					mat0 = TagMapping::getOtherMat(mat1, matMarks);
					boundaryMats[F2BBFInd[nbBBFace[0]]] = mat0;
				}
			} else {
				if (mat1==TagMapping::MAT_NAN){
					mat1 = TagMapping::getOtherMat(mat0, matMarks);
					boundaryMats[F2BBFInd[nbBBFace[1]]] = mat1;
				}
			}
		}
		// otherwise, edges on bbox, but not on input segments, add dual edge
		else{
			bbGraph.addEdge(F2BBFInd[nbBBFace[0]],F2BBFInd[nbBBFace[1]]);
		}
		

	}
#if _DEBUG_MING_MORE
	Utility::writeOutOneVector((t_savedir+"boundaryMats_"+Utility::toStr(ci)+".txt").c_str(), boundaryMats);
#endif

	// step3, label faces on the cell boundary to be in/out
	BigInt visited(bbGraph.size());
	vector<int> F2compID(bbGraph.size(),-1);
	vector<int> compID2mat;
	int compID = -1;
	while (!visited.isAllOne()){
		compID++;
		queue<int> fringe;
		bool compLabelSet = false;
		int curF = visited.find1st(0);
		visited.set(curF);
		fringe.push(curF);
		F2compID[curF] = compID;
		if (boundaryMats[curF]!=TagMapping::MAT_NAN && compLabelSet==false){
			compID2mat.push_back(boundaryMats[curF]);
			compLabelSet = true;
		}

		while (!fringe.empty()){
			curF = fringe.front();
			fringe.pop();

			for(auto nbF : bbGraph.adj[curF]){
				if (!visited.get(nbF)){
					fringe.push(nbF);
					visited.set(nbF);
					F2compID[nbF] = compID;
					if (boundaryMats[nbF]!=TagMapping::MAT_NAN && compLabelSet==false){
						compID2mat.push_back(boundaryMats[nbF]);
						compLabelSet = true;
					}
				}
			}
		}
	}
	//label face/edge/vertices to be in/out/onSeg/interior
	for (int bbfi = 0; bbfi<bbGraph.size(); ++bbfi){
		int curFaceID = allCellTets.cells[ci].cellBBFaces[bbfi];
		int curInOutMark = compID2mat[F2compID[bbfi]];
		allCellTets.cells[ci].cellFaces[curFaceID].FinoutMark = curInOutMark; 
		for(auto ei : allCellTets.cells[ci].cellFaces[curFaceID].edge){
			if (allCellTets.cells[ci].cellEdges[ei].EinoutMark==TagMapping::MAT_NAN){
				allCellTets.cells[ci].cellEdges[ei].EinoutMark = curInOutMark;
			} else if(allCellTets.cells[ci].cellEdges[ei].EinoutMark != curInOutMark){
				allCellTets.cells[ci].cellEdges[ei].EinoutMark = TagMapping::MAT_SEG;
			}

		}
		for(auto vi : allCellTets.cells[ci].cellFaces[curFaceID].v){
			if (allCellTets.cells[ci].cellVertices[vi].VinoutMark==TagMapping::MAT_NAN){
				allCellTets.cells[ci].cellVertices[vi].VinoutMark = curInOutMark;
			} else if(allCellTets.cells[ci].cellVertices[vi].VinoutMark != curInOutMark){
				allCellTets.cells[ci].cellVertices[vi].VinoutMark = TagMapping::MAT_SEG;
			}
		}

	}

	vector<int> edgeToSplit;
	for(size_t ei=0; ei<allCellTets.cells[ci].cellEdges.size(); ++ei){
		int eInOutMark = allCellTets.cells[ci].cellEdges[ei].EinoutMark;
		int v0InOutMark = allCellTets.cells[ci].cellVertices[allCellTets.cells[ci].cellEdges[ei].v[0]].VinoutMark;
		int v1InOutMark = allCellTets.cells[ci].cellVertices[allCellTets.cells[ci].cellEdges[ei].v[1]].VinoutMark;
		
		// edges in the interior
		if (eInOutMark==TagMapping::MAT_NAN){
			if ((v0InOutMark==TagMapping::MAT_IN && v1InOutMark==TagMapping::MAT_IN)||(v0InOutMark==TagMapping::MAT_OUT && v1InOutMark==TagMapping::MAT_OUT)
				||(v0InOutMark==TagMapping::MAT_SEG && v1InOutMark==TagMapping::MAT_SEG)||(v0InOutMark==TagMapping::MAT_SEG && v1InOutMark==TagMapping::MAT_IN)||(v0InOutMark==TagMapping::MAT_IN && v1InOutMark==TagMapping::MAT_SEG)){
				edgeToSplit.push_back(ei);
			}
		} 
		// edges in the "outside", don't need to do anything if edge is in the "inside" (which is not possible to have its two ens to be outside) or on the input seg
		else if (eInOutMark==TagMapping::MAT_OUT){
			if ((v0InOutMark==TagMapping::MAT_IN||v0InOutMark==TagMapping::MAT_SEG) && (v1InOutMark==TagMapping::MAT_IN||v1InOutMark==TagMapping::MAT_SEG)){
				edgeToSplit.push_back(ei);
			}

		}
		else if (eInOutMark==TagMapping::MAT_IN){
			if ((v0InOutMark==TagMapping::MAT_OUT||v0InOutMark==TagMapping::MAT_SEG) && (v1InOutMark==TagMapping::MAT_OUT||v1InOutMark==TagMapping::MAT_SEG)){
				edgeToSplit.push_back(ei);
			}

		}
	}

	bool ep = false;
	if (ep) cout << "      # of edges to split: " << edgeToSplit.size() << endl;

	int ecount = 0;

	for(auto curSei : edgeToSplit){
		ecount ++;
		if(ep) cout << "cur ecount = "<< ecount << ", eToSplit = " << curSei << endl;
		vector<int> cenNBv;
		
		vector<int> topNBe;
		vector<int> cenNBe;
		vector<int> btmNBe;

		vector<int> topNBf;
		vector<int> cenNBf;
		vector<int> btmNBf;

		vector<int> topNBt;
		vector<int> cenNBt;
		vector<int> btmNBt;

		bool isOpenCircle = false;

		int v_top = allCellTets.cells[ci].cellEdges[curSei].v[0];
		int v_btm = allCellTets.cells[ci].cellEdges[curSei].v[1];

		vector<int> orgNBf = allCellTets.cells[ci].cellEdges[curSei].nbFace;
		vector<vector<int> > in_nbTPairs;
		int in_startPairID = 0;
		int in_startEleID = allCellTets.cells[ci].cellFaces[orgNBf[0]].nbTet[0];

		if (ep) cout << "stop 1..." << endl;
		for(int find=0; find<orgNBf.size(); ++find){
			vector<int> onePair = allCellTets.cells[ci].cellFaces[orgNBf[find]].nbTet;
			if (onePair.size()==1){
				in_startPairID = find;
				in_startEleID = onePair[0];
				onePair.push_back(-1);
				isOpenCircle = true;
			}
			in_nbTPairs.push_back(onePair);
		}
		if (ep) cout << "stop 2..." << endl;
		
		Utility::findOrderedPairs(in_nbTPairs, orgNBf, in_startPairID, in_startEleID, cenNBf, cenNBt);
		if (ep) cout << "stop 3..." << endl;

		if (isOpenCircle){
			cenNBt.pop_back();
		}
		if (ep) cout << "stop 4..." << endl;

		for(auto fi : cenNBf){
			// cenNBv
			int v_cen = allCellTets.cells[ci].cellFaces[fi].v[0] + allCellTets.cells[ci].cellFaces[fi].v[1] + allCellTets.cells[ci].cellFaces[fi].v[2] - v_top - v_btm;
			cenNBv.push_back(v_cen);

			// topNBe, btmNBe
			int e_top, e_btm;
			for(int eind=0; eind<3; ++eind){
				int e_tmp = allCellTets.cells[ci].cellFaces[fi].edge[eind];
				if (e_tmp!=curSei){
					if (allCellTets.cells[ci].cellEdges[e_tmp].v[0]==v_top || allCellTets.cells[ci].cellEdges[e_tmp].v[1]==v_top){
						e_top = e_tmp;
					} else {
						e_btm = e_tmp;
					}
				}
			}
			topNBe.push_back(e_top);
			btmNBe.push_back(e_btm);
		}
		if (ep) cout << "stop 5..." << endl;
		for (int tind=0; tind<cenNBt.size(); ++tind){
			int curTetID = cenNBt[tind];
			int nbFi0 = tind;
			int nbFi1 = tind+1;
			if ((!isOpenCircle) && tind == cenNBt.size()-1){
				nbFi1 = 0;
			}

			// topNBf, btmNBf
			int f_0 = cenNBf[nbFi0];
			int f_1 = cenNBf[nbFi1];
			int f_top, f_btm;
			for(int find=0; find<4; ++find){
				int f_tmp = allCellTets.cells[ci].cellTets[curTetID].face[find];
				if (f_tmp!=f_0 && f_tmp!=f_1){
					if (allCellTets.cells[ci].cellFaces[f_tmp].v[0]==v_top || allCellTets.cells[ci].cellFaces[f_tmp].v[1]==v_top ||allCellTets.cells[ci].cellFaces[f_tmp].v[2]==v_top){
						f_top = f_tmp;
					} else {
						f_btm = f_tmp;
					}
				}
			}
			topNBf.push_back(f_top);
			btmNBf.push_back(f_btm);

			// cenNBe
			int e_cen = allCellTets.cells[ci].cellFaces[f_top].edge[0] + allCellTets.cells[ci].cellFaces[f_top].edge[1] + allCellTets.cells[ci].cellFaces[f_top].edge[2] - topNBe[nbFi0]- topNBe[nbFi1];
			cenNBe.push_back(e_cen);

			// topNBt, btmNBt
			int t_top = -1; 
			int t_btm = -1; 
			if (allCellTets.cells[ci].cellFaces[f_top].nbTet.size()>1){
				t_top=allCellTets.cells[ci].cellFaces[f_top].nbTet[0] + allCellTets.cells[ci].cellFaces[f_top].nbTet[1] - curTetID;
			}
			if (allCellTets.cells[ci].cellFaces[f_btm].nbTet.size()>1){
				t_btm=allCellTets.cells[ci].cellFaces[f_btm].nbTet[0] + allCellTets.cells[ci].cellFaces[f_btm].nbTet[1] - curTetID;
			}
			topNBt.push_back(t_top);
			btmNBt.push_back(t_btm);
		}
		if (ep) cout << "stop 6..." << endl;

		// finished constructing the helping datastructures, for new cell elements
		
		// prepare indexes
		int curNedge = allCellTets.cells[ci].cellEdges.size();
		int curNface = allCellTets.cells[ci].cellFaces.size();
		int curNtet = allCellTets.cells[ci].cellTets.size();
		int addNedge = cenNBf.size()+2;
		int addNface = 2*cenNBf.size()+cenNBt.size();
		int addNtet = 2*cenNBt.size();

		// prepare datastructure for additional cell elements
		if (ep) cout << "stop 7..." << endl;
		CellVertex addV;
		vector<CellEdge> addEs; addEs.reserve(addNedge);
		vector<CellFace> addFs; addFs.reserve(addNface);
		vector<CellTet> addTets; addTets.reserve(addNtet);
		CellTet oneTet;
		CellFace oneFace;
		CellEdge oneEdge;

		// addV
		int v_x = allCellTets.cells[ci].cellVertices.size();
		if (ep) cout << "stop 8..." << endl;
		for (int i=0; i<3; ++i){
			addV.xyz[i] = (allCellTets.cells[ci].cellVertices[v_top].xyz[i]+allCellTets.cells[ci].cellVertices[v_btm].xyz[i])/2.0;
		}
		if (ep) cout << "stop 9..." << endl;
		addV.nbV = cenNBv;
		addV.nbV.push_back(v_top);
		addV.nbV.push_back(v_btm);
		addV.nbE.reserve(addV.nbV.size());
		addV.nbE.push_back(curSei);

		for(int i=curNedge; i<=curNedge+cenNBf.size(); ++i){
			addV.nbE.push_back(i);
		}
		addV.nbTet = cenNBt;
		for(int i=curNtet; i<curNtet+cenNBt.size(); ++i){
			addV.nbTet.push_back(i);
		}
		addV.VinoutMark = allCellTets.cells[ci].cellEdges[curSei].EinoutMark;
		if (ep) cout << "stop 10..." << endl;
		// addEs
		// (1) (v_top, v_x)
		oneEdge.v[0] = v_top;
		oneEdge.v[1] = v_x;
		oneEdge.nbFace = cenNBf;
		oneEdge.EinoutMark = allCellTets.cells[ci].cellEdges[curSei].EinoutMark;
		addEs.push_back(oneEdge);

		// (2) (v_btm, v_x)
		oneEdge.v[0] = v_btm;
		for(int i=0; i<cenNBf.size(); ++i){
			oneEdge.nbFace[i] = curNface + i;
		}
		addEs.push_back(oneEdge);

		// (3~n) (cenNBv[i], v_x)
		for(int i=0; i<cenNBf.size(); ++i){
			oneEdge.v[0] = cenNBv[i];
			oneEdge.nbFace.clear();
			oneEdge.nbFace.push_back(cenNBf[i]);
			oneEdge.nbFace.push_back(curNface+i);
			if (i==0){
				oneEdge.nbFace.push_back(curNface+cenNBf.size()+i);
				if(!isOpenCircle){
					oneEdge.nbFace.push_back(curNface+cenNBf.size()+cenNBt.size()-1);
				}
			}
			else if (i==cenNBf.size()-1){
				oneEdge.nbFace.push_back(curNface+cenNBf.size()+i-1);
				if(!isOpenCircle){
					oneEdge.nbFace.push_back(curNface+cenNBf.size()+i);
				}
			} else {
				oneEdge.nbFace.push_back(curNface+cenNBf.size()+i-1);
				oneEdge.nbFace.push_back(curNface+cenNBf.size()+i);
			}
			oneEdge.EinoutMark = allCellTets.cells[ci].cellFaces[cenNBf[i]].FinoutMark;
			addEs.push_back(oneEdge);
		}
		if (ep) cout << "stop 11..." << endl;
		// addFs
		// (1) top on cenNBf
		oneFace.v[0]=v_x;
		oneFace.v[1]=v_top;
		oneFace.edge[0]=curSei;
		for(int i=0; i<cenNBf.size(); ++i){
			oneFace.v[2]=cenNBv[i];
			oneFace.edge[1]=topNBe[i];
			oneFace.edge[2]=curNedge+1+i;
			oneFace.nbTet.clear();
			if (i==0){
				oneFace.nbTet.push_back(cenNBt[i]);
				if(!isOpenCircle){
					oneFace.nbTet.push_back(cenNBt[cenNBt.size()-1]);
				}
			}
			else if (i==cenNBf.size()-1){
				oneFace.nbTet.push_back(cenNBt[i-1]);
				if(!isOpenCircle){
					oneFace.nbTet.push_back(cenNBt[i]);
				}
			} else {
				oneFace.nbTet.push_back(cenNBt[i]);
				oneFace.nbTet.push_back(cenNBt[i-1]);
			}
			oneFace.FinoutMark = allCellTets.cells[ci].cellFaces[cenNBf[i]].FinoutMark;
			addFs.push_back(oneFace);
		}
		// (2) btm on cenNBf
		oneFace.v[0]=v_x;
		oneFace.v[1]=v_btm;
		oneFace.edge[0]=curNedge;
		for(int i=0; i<cenNBf.size(); ++i){
			oneFace.v[2]=cenNBv[i];
			oneFace.edge[1]=btmNBe[i];
			oneFace.edge[2]=curNedge+1+i;
			oneFace.nbTet.clear();
			if (i==0){
				oneFace.nbTet.push_back(curNtet+i);
				if(!isOpenCircle){
					oneFace.nbTet.push_back(curNtet+cenNBt.size()-1);
				}
			}
			else if (i==cenNBf.size()-1){
				oneFace.nbTet.push_back(curNtet+i-1);
				if(!isOpenCircle){
					oneFace.nbTet.push_back(curNtet+i);
				}
			} else {
				oneFace.nbTet.push_back(curNtet+i);
				oneFace.nbTet.push_back(curNtet+i-1);
			}
			oneFace.FinoutMark = allCellTets.cells[ci].cellFaces[cenNBf[i]].FinoutMark;
			addFs.push_back(oneFace);
		}
		// (3) faces in cenNBt
		oneFace.v[0]=v_x;
		oneFace.nbTet.resize(2);
		for (int i=0; i<cenNBt.size(); ++i){
			int ind_0 = i;
			int ind_1 = i+1;
			if((i==cenNBt.size()-1) && (!isOpenCircle)){
				ind_1 = 0;
			}
			oneFace.v[1] = cenNBv[ind_0];
			oneFace.v[2] = cenNBv[ind_1];
			oneFace.edge[0]=cenNBe[i];
			oneFace.edge[1]=curNedge+1+ind_0;
			oneFace.edge[2]=curNedge+1+ind_1;
			oneFace.nbTet[0]=cenNBt[i];
			oneFace.nbTet[1]=curNtet+i;
			oneFace.FinoutMark = TagMapping::MAT_NAN;
			addFs.push_back(oneFace);
		}
		if (ep) cout << "stop 12..." << endl;
		//addTets
		// (1) tets on top
		oneTet.v[0]=v_top;
		oneTet.v[1]=v_x;
		for (int i=0; i<cenNBt.size(); ++i){
			int ind_0 = i;
			int ind_1 = i+1;
			if((i==cenNBt.size()-1) && (!isOpenCircle)){
				ind_1 = 0;
			}
			oneTet.v[2]=cenNBv[ind_0];
			oneTet.v[3]=cenNBv[ind_1];
			oneTet.face[0]=curNface+cenNBf.size()+i;
			oneTet.face[1]=topNBf[i];
			oneTet.face[2]=cenNBf[ind_1];
			oneTet.face[3]=cenNBf[ind_0];
			oneTet.nbTet[0]=curNtet+i;
			oneTet.nbTet[1]=topNBt[i];
			if(i==0){
				oneTet.nbTet[2]=cenNBt[i+1];
				if(isOpenCircle){
					oneTet.nbTet[3]=-1;
				}else{
					oneTet.nbTet[3]=cenNBt[cenNBt.size()-1];
				}
			}
			else if(i==cenNBt.size()-1){
				oneTet.nbTet[3]=cenNBt[i-1];
				if(isOpenCircle){
					oneTet.nbTet[2]=-1;
				}else{
					oneTet.nbTet[2]=cenNBt[0];
				}
			}else{
				oneTet.nbTet[2]=cenNBt[i+1];
				oneTet.nbTet[3]=cenNBt[i-1];
			}
			addTets.push_back(oneTet);
		}
			
		// (2) tets on btm 
		oneTet.v[0]=v_btm;
		oneTet.v[1]=v_x;
		for (int i=0; i<cenNBt.size(); ++i){
			int ind_0 = i;
			int ind_1 = i+1;
			if((i==cenNBt.size()-1) && (!isOpenCircle)){
				ind_1 = 0;
			}
			oneTet.v[2]=cenNBv[ind_0];
			oneTet.v[3]=cenNBv[ind_1];
			oneTet.face[0]=curNface+cenNBf.size()+i;
			oneTet.face[1]=btmNBf[i];
			oneTet.face[2]=curNface+ind_1;
			oneTet.face[3]=curNface+ind_0;
			oneTet.nbTet[0]=cenNBt[i];
			oneTet.nbTet[1]=btmNBt[i];
			if(i==0){
				oneTet.nbTet[2]=curNtet+i+1;
				if(isOpenCircle){
					oneTet.nbTet[3]=-1;
				}else{
					oneTet.nbTet[3]=curNtet+cenNBt.size()-1;
				}
			}
			else if(i==cenNBt.size()-1){
				oneTet.nbTet[3]=curNtet+i-1;
				if(isOpenCircle){
					oneTet.nbTet[2]=-1;
				}else{
					oneTet.nbTet[2]=curNtet;
				}
			}else{
				oneTet.nbTet[2]=curNtet+i+1;
				oneTet.nbTet[3]=curNtet+i-1;
			}
			addTets.push_back(oneTet);
		}
		if (ep) cout << "stop 13..." << endl;
		//----------------------------------------
		//now change the nbinfo of existing V/E/F/T which were impacted by split
		//V
		//(1) v_top
		for(int i=0; i<allCellTets.cells[ci].cellVertices[v_top].nbV.size(); ++i){
			if(allCellTets.cells[ci].cellVertices[v_top].nbV[i]==v_btm){
				allCellTets.cells[ci].cellVertices[v_top].nbV[i]=v_x;
				break;
			}
		}
		//don't need to change nbE,nbTet for v_top, because indexes are the same
		//(2) v_btm
		for(int i=0; i<allCellTets.cells[ci].cellVertices[v_btm].nbV.size(); ++i){
			if(allCellTets.cells[ci].cellVertices[v_btm].nbV[i]==v_top){
				allCellTets.cells[ci].cellVertices[v_btm].nbV[i]=v_x;
				break;
			}
		}
		for(int i=0; i<allCellTets.cells[ci].cellVertices[v_btm].nbE.size(); ++i){
			if(allCellTets.cells[ci].cellVertices[v_btm].nbE[i]==curSei){
				allCellTets.cells[ci].cellVertices[v_btm].nbE[i]=curNedge;
				break;
			}
		}
		vector<int> btm_nbTet = allCellTets.cells[ci].cellVertices[v_btm].nbTet;
		for(auto ti : cenNBt){
			btm_nbTet.erase(std::remove(btm_nbTet.begin(),btm_nbTet.end(),ti),btm_nbTet.end());
		}
		for(int i=0; i<cenNBt.size(); ++i){
			btm_nbTet.push_back(curNtet+i);
		}
		allCellTets.cells[ci].cellVertices[v_btm].nbTet = btm_nbTet;

		// (3) v of cenNBv
		for(int i=0; i<cenNBv.size(); ++i){
			int curVID = cenNBv[i];
			allCellTets.cells[ci].cellVertices[curVID].nbV.push_back(v_x);
			allCellTets.cells[ci].cellVertices[curVID].nbE.push_back(curNedge+1+i);
			vector<int> cen_nbTet = allCellTets.cells[ci].cellVertices[curVID].nbTet;
			if(i==0){
				cen_nbTet.push_back(curNtet+i);
				if(!isOpenCircle){
					cen_nbTet.push_back(curNtet+cenNBt.size()-1);
				}
			}else if (i==cenNBv.size()-1){
				cen_nbTet.push_back(curNtet+i-1);
				if(!isOpenCircle){
					cen_nbTet.push_back(curNtet+i);
				}
			}else{
				cen_nbTet.push_back(curNtet+i);
				cen_nbTet.push_back(curNtet+i-1);
			}
			allCellTets.cells[ci].cellVertices[curVID].nbTet = cen_nbTet;
		}
		if (ep) cout << "stop 14..." << endl;
		//E
		//topNBe doesn't change
		//(1) btmNBe
		for(int i=0; i<btmNBe.size(); ++i){
			int curEID = btmNBe[i];
			std::replace(allCellTets.cells[ci].cellEdges[curEID].nbFace.begin(), allCellTets.cells[ci].cellEdges[curEID].nbFace.end(), cenNBf[i], curNface+i);
		}

		//(2) cenNBe
		for(int i=0; i<cenNBe.size(); ++i){
			int curEID = cenNBe[i];
			allCellTets.cells[ci].cellEdges[curEID].nbFace.push_back(curNface+cenNBf.size()+i);
		}
		if (ep) cout << "stop 15..." << endl;
		//F
		//topNBf doesn't change, cenNBf doesn't exist anymore
		//(1) btmNBf
		for(int i=0; i<btmNBf.size(); ++i){
			int curFID = btmNBf[i];
			std::replace(allCellTets.cells[ci].cellFaces[curFID].nbTet.begin(), allCellTets.cells[ci].cellFaces[curFID].nbTet.end(), cenNBt[i], curNtet+i);
		}
		if (ep) cout << "stop 16..." << endl;
		//Tet
		//topNBt doesn't change, cenNBt doesn't exist anymore
		//(1) btmNBt
		for(int i=0; i<btmNBt.size(); ++i){
			int curTID = btmNBt[i];
			if(curTID==-1){
				continue;
			}
			for(int j=0; j<4; ++j){
				if(allCellTets.cells[ci].cellTets[curTID].nbTet[j]==cenNBt[i]){
					allCellTets.cells[ci].cellTets[curTID].nbTet[j] = curNtet+i;
					break;
				}
			}
		}
		if (ep) cout << "stop 17..." << endl;
		// now add the additional V/E/F/T into allCellTets
		// V
		allCellTets.cells[ci].cellVertices.push_back(addV);
		// E
		allCellTets.cells[ci].cellEdges[curSei] = addEs[0];
		allCellTets.cells[ci].cellEdges.insert(allCellTets.cells[ci].cellEdges.end(), addEs.begin()+1, addEs.end());
		// F
		for (int i=0; i<cenNBf.size(); ++i){
			int curFID = cenNBf[i];
			allCellTets.cells[ci].cellFaces[curFID] = addFs[i];
		}
		allCellTets.cells[ci].cellFaces.insert(allCellTets.cells[ci].cellFaces.end(), addFs.begin()+cenNBf.size(), addFs.end());
		// T
		for (int i=0; i<cenNBt.size(); ++i){
			int curTID = cenNBt[i];
			allCellTets.cells[ci].cellTets[curTID] = addTets[i];
		}
		allCellTets.cells[ci].cellTets.insert(allCellTets.cells[ci].cellTets.end(), addTets.begin()+cenNBt.size(), addTets.end());
		if (ep) cout << "stop 18..." << endl;

		if (ep) cout << "+++++++++++++++++++" << endl;
	}
	
	if (ep) cout << "Edge split is DONE!" << endl;

	// split the Vs into inside/outside/unknown
	tetV2nbTet.reserve(allCellTets.cells[ci].cellVertices.size());
	for(size_t i=0; i<allCellTets.cells[ci].cellVertices.size(); ++i){
		int VinoutMark = allCellTets.cells[ci].cellVertices[i].VinoutMark;
		if(VinoutMark==TagMapping::MAT_NAN){
			labelXVs.push_back(i);
		}else if(VinoutMark==TagMapping::MAT_OUT){
			labelOutVs.push_back(i);
		}
		// MAT_IN or MAT_SEG are all considered as in
		else{
			labelInVs.push_back(i);
		}
		tetV2nbTet.push_back(allCellTets.cells[ci].cellVertices[i].nbTet);
	}

	if (ep) cout << "In/Out split of Vs is DONE!" << endl;

}

void Arrangement::tetCell_Vaxman(const int ci, const vector<int> &labelInVs, const vector<int> &labelOutVs, const vector<int> &labelXVs, vector<float> &VPro){
	
	cout << " We are in VAXMAN!!! "<< endl;
	
	vector<vector<float> > Verts;
	vector<vector<int> > Tris;
	vector<bool> Flags;

	vector<int> tmpTs_in;
	vector<int> tmpTs_out;
	vector<int> tmpVs;

	tmpVs.reserve(allCellTets.cells[ci].cellVertices.size()/2);
	tmpTs_in.reserve(allCellTets.cells[ci].cellFaces.size()/4);
	tmpTs_out.reserve(allCellTets.cells[ci].cellFaces.size()/4);

	unordered_map<int, int> V2Vhash;

	Verts.reserve(allCellTets.cells[ci].cellVertices.size()/2);
	Tris.reserve(allCellTets.cells[ci].cellFaces.size()/4);
	Flags.reserve(allCellTets.cells[ci].cellFaces.size()/4);

	int fmark;
	vector<int> oneF(3);
	vector<float> oneV(3);
	for(size_t fi=0; fi<allCellTets.cells[ci].cellFaces.size(); ++fi){
		fmark = allCellTets.cells[ci].cellFaces[fi].FinoutMark;
		if((fmark == TagMapping::MAT_IN) || (fmark == TagMapping::MAT_OUT)){
			for(int i=0; i<3; ++i){
				oneF[i] = allCellTets.cells[ci].cellFaces[fi].v[i];
				
				//MINGDEBUG
				tmpVs.push_back(oneF[i]);

				// hashed already
				if(V2Vhash.find(oneF[i])!=V2Vhash.end()){
					oneF[i] = V2Vhash[oneF[i]];
				}
				// new vertex to hash
				else{
					oneV[0] = allCellTets.cells[ci].cellVertices[oneF[i]].xyz[0];
					oneV[1] = allCellTets.cells[ci].cellVertices[oneF[i]].xyz[1];
					oneV[2] = allCellTets.cells[ci].cellVertices[oneF[i]].xyz[2];
					V2Vhash[oneF[i]] = (int) (Verts.size());
					oneF[i] = (int) (Verts.size());
					Verts.push_back(oneV);
				}
			}
			Tris.push_back(oneF);
			if(fmark == TagMapping::MAT_IN){
				Flags.push_back(true);

				//MINGDEBUG
				tmpTs_in.push_back(fi);
			} else {
				Flags.push_back(false);

				//MINGDEBUG
				tmpTs_out.push_back(fi);
			}
		}
	}

	 Utility::forceSufOriented(Tris);

	 MVC weightMVC(Verts, Tris, Flags);

	VPro.resize(labelInVs.size()+labelOutVs.size()+labelXVs.size(), INFINITY);
	for (auto v : labelInVs){
		VPro[v] = 1;
	}
	for (auto v : labelOutVs){
		VPro[v] = 0;
	}
	for (auto v : labelXVs){
		VPro[v] = weightMVC.getValue(allCellTets.cells[ci].cellVertices[v].xyz);
	}
	allCellTets.cells[ci].VPro = VPro;

}

void Arrangement::tetCell_randomWalk(UnDirGraph &graph,vector<int> &labelInVs, vector<int> &labelOutVs, vector<int> &labelXVs, vector<float> &VPro, int tmpCi){

	int vCount=0;
	int nVIn = labelInVs.size();
	int nVOut = labelOutVs.size();
	int nVX = labelXVs.size();
	int nV= nVIn + nVOut + nVX;
	vector<int> tetV2orderedV(nV);
	for (auto v : labelInVs){
		tetV2orderedV[v] = vCount++;
	}
	for (auto v : labelOutVs){
		tetV2orderedV[v] = vCount++;
	}
	for (auto v : labelXVs){
		tetV2orderedV[v] = vCount++;
	}
	vector<int> orderedV2tetV = labelInVs;
	orderedV2tetV.reserve(nV);
	orderedV2tetV.insert(orderedV2tetV.end(), labelOutVs.begin(), labelOutVs.end());
	orderedV2tetV.insert(orderedV2tetV.end(), labelXVs.begin(), labelXVs.end());

#if _DEBUG_MING_MORE
	if(!volInfo.empty()){
		vector<float> V2VolVal(nV);
		for(int i=0; i<nV; ++i){
			V2VolVal[i] = getVolValAtPos(allCellTets.cells[tmpCi].cellVertices[i].xyz);
		}
		Utility::writeOutOneVector((t_savedir+"V2VolVal_"+Utility::toStr(tmpCi)+".txt").c_str(), V2VolVal);
	}
#endif

	Eigen::SparseMatrix<float> matL(nVX, nV);
	vector<Eigen::Triplet<float> >tripletList;
	tripletList.reserve(nVX*6);
	for (int i=nVIn+nVOut; i<nV; ++i){
		
		float di = 0.0f;
		for (auto j : allCellTets.cells[tmpCi].cellVertices[orderedV2tetV[i]].nbV){
			
			int ii = orderedV2tetV[i];
			int jj = j;
			float wij = getVolWeight(allCellTets.cells[tmpCi].cellVertices[ii].xyz,allCellTets.cells[tmpCi].cellVertices[jj].xyz);
			di += wij;
			tripletList.push_back(Eigen::Triplet<float>(i-nVIn-nVOut,tetV2orderedV[j],-wij));
		}
		
		tripletList.push_back(Eigen::Triplet<float>(i-nVIn-nVOut,i,di));
	}
	matL.setFromTriplets(tripletList.begin(), tripletList.end());
	Eigen::SparseMatrix<float> Bt = matL.leftCols(nVIn+nVOut);
	Eigen::SparseMatrix<float> Lu = matL.rightCols(nVX);

	Eigen::VectorXf b = Eigen::VectorXf::Zero(nVIn+nVOut);
	for (int i=0; i<nVIn; ++i){
		b[i] = 1;
	}

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver;
	solver.compute(Lu);
	if (solver.info()!=Eigen::Success){
		cout << "Eigen LDLT factorization failed !!!! " << endl;
	}
	Eigen::VectorXf x = solver.solve(-Bt*b);
	if (solver.info()!=Eigen::Success){
		cout << "Eigen solving failed !!!! " << endl;
	}

	VPro.resize(nV, INFINITY);
	for (auto v : labelInVs){
		VPro[v] = 1;
	}
	for (auto v : labelOutVs){
		VPro[v] = 0;
	}
	for (int i=0; i<nVX; ++i){
		VPro[orderedV2tetV[nVIn+nVOut+i]] = x[i];
	}
	allCellTets.cells[tmpCi].VPro = VPro;
}

void Arrangement::tetCell_volWeight(int ci, vector<float> &VProVol, vector<int> &labelXVs){
	vector<float> tetVol(allCellTets.cells[ci].cellTets.size());
	
	for (int ti=0; ti<allCellTets.cells[ci].cellTets.size(); ++ti){
		vector<Eigen::Vector3f> tetVs;
		tetVs.reserve(4);
		for (int vi=0; vi<4; ++vi){
			int vind = allCellTets.cells[ci].cellTets[ti].v[vi];
			tetVs.push_back(Eigen::Vector3f(allCellTets.cells[ci].cellVertices[vind].xyz[0], allCellTets.cells[ci].cellVertices[vind].xyz[1], allCellTets.cells[ci].cellVertices[vind].xyz[2]));
		}
		tetVol[ti] = abs((tetVs[1]-tetVs[0]).dot((tetVs[2]-tetVs[0]).cross(tetVs[3]-tetVs[0])))/6.0;
	}

	VProVol.resize(allCellTets.cells[ci].cellVertices.size());

	for (auto vi : labelXVs){
		float nbTetVolSum = 0.0f;
		for (auto ti : allCellTets.cells[ci].cellVertices[vi].nbTet){
			if(ti==-1) continue;
			nbTetVolSum += tetVol[ti];
		}
		VProVol[vi] = nbTetVolSum; 
	}

}

int Arrangement::tetCell_criticalPoint_nComp(int ci, int vID, vector<int> &nbTets, vector<float> &VPro, tetgenio &tetout){
	float cenVal = VPro[vID];
	UnDirGraph lGraph(VPro.size());
	for (auto ti : nbTets){
		for (int i = 0; i < 3; ++i){
			for ( int j = i+1; j < 4; ++j){
				int v1 = allCellTets.cells[ci].cellTets[ti].v[i];
				int v2 = allCellTets.cells[ci].cellTets[ti].v[j];
				lGraph.addEdge(v1, v2);
			}
		}
	}
	lGraph.removeDuplicates();
	vector<int> &nbVs = lGraph.adj[vID];

	BigInt unVisit(allCellTets.cells[ci].cellVertices.size());
	for (auto vi : nbVs){
		unVisit.set(vi);
	}
	int curV = unVisit.find1st(1);

	int nComp = 0;
	bool sign = false;
	while(curV!=-1){
		nComp++;
		if (nComp > 2) break;
		queue<int> fringe;
		fringe.push(curV);
		unVisit.unset(curV);
		sign = VPro[curV] >= cenVal;
		while(!fringe.empty()){
			curV = fringe.front();
			fringe.pop();
			for (auto v : lGraph.adj[curV]){
				if ((unVisit.get(v)==1)&&(sign == (VPro[v]>=cenVal))){
					unVisit.unset(v);
					fringe.push(v);
				}
			}
		}
		curV = unVisit.find1st(1);
	}
	return nComp;
}

void Arrangement::tetCell_criticalPoint(vector<vector<int> > &tetV2nbTet, vector<float> &VPro, tetgenio &tetout, vector<int> &labelXVs, vector<float> &criticalIsoVal, int tmpCi){
#if VAXMAN
	criticalIsoVal.push_back(0.5);
#else
	vector<float> criticalIsoValTrue;
	for (auto vi : labelXVs){
		if ((tetCell_criticalPoint_nComp(tmpCi, vi, tetV2nbTet[vi], VPro, tetout) != 2) && (VPro[vi]!=0.0) && (VPro[vi]!=1.0)){
			criticalIsoValTrue.push_back(VPro[vi]);
		}
	}
	criticalIsoValTrue.push_back(0.0f);
	criticalIsoValTrue.push_back(1.0f);
	sort(criticalIsoValTrue.begin(), criticalIsoValTrue.end());
	for (size_t i=0; i<criticalIsoValTrue.size()-1; ++i){

		float trueIso;
		if(criticalIsoValTrue[i]<=0.5 && criticalIsoValTrue[i+1]>= 0.5){
			trueIso = 0.5;
		} else if (criticalIsoValTrue[i]>0.5){
			trueIso = criticalIsoValTrue[i] + (criticalIsoValTrue[i+1]-criticalIsoValTrue[i])*0.1;
		} else{
			trueIso = criticalIsoValTrue[i+1] - (criticalIsoValTrue[i+1]-criticalIsoValTrue[i])*0.1;
		}
		criticalIsoVal.push_back(trueIso);
	}
	criticalIsoVal.erase(unique(criticalIsoVal.begin(), criticalIsoVal.end()), criticalIsoVal.end());
#endif
}

float Arrangement::tetCell_findGroupings_GroupScore(vector<float> &VPro, vector<float> &VProVol, vector<int> &labelXVs, float iso){
	float res = 0.0f;
	for (auto vi : labelXVs){
		if (VPro[vi] > iso){
			res += VProVol[vi]*logf(VPro[vi]);
		} else {
			float val = 1 - VPro[vi];
			// numerical issue
			if(val<=0){
				val=0.000001;
			}
			res += VProVol[vi]*logf(val);
		}
	}
	return res;
}
float Arrangement::tetCell_findGroupings_mix_GroupScore(int ci, const vector<int> &curMixIsoGrouping, const vector<int> &cycFirstV, const vector<BigInt> &loc_set2actEdge, int nLocSet, const vector<float> &VPro, const vector<float> &VProVol, const vector<int> &labelXVs){
	BigInt IsEdgeOnSuf(allCellTets.cells[ci].cellEdges.size());
	int curEdge = -1;
	for(auto loc_subGroupID : curMixIsoGrouping){
		BigInt actEdges = loc_set2actEdge[loc_subGroupID];
		while(!actEdges.isZero()){
			curEdge = actEdges.popOne();
			IsEdgeOnSuf.set(curEdge);
		}
	}
	BigInt visited(allCellTets.cells[ci].cellVertices.size());
	// flood from each cycle, until all cycle Vs are visited
	for(auto srcV : cycFirstV){
		// if srcV is visited, then this V has been flooded by other srcV, skip
		if(visited.get(srcV)) continue;
		queue<int> fringe;
		fringe.push(srcV);
		visited.set(srcV);
		int curV;
		int nbV;
		while(!fringe.empty()){
			curV = fringe.front();
			fringe.pop();
			for(int i=0; i<allCellTets.cells[ci].cellVertices[curV].nbE.size(); ++i){
				// if the nb edge is on surface, stop flood
				if(IsEdgeOnSuf.get(allCellTets.cells[ci].cellVertices[curV].nbE[i])) continue;
				nbV = allCellTets.cells[ci].cellVertices[curV].nbV[i];
				if(!visited.get(nbV)){
					fringe.push(nbV);
					visited.set(nbV);
				}
			}
		}
	}

	BigInt &IsVInside = visited;

	float res = 0.0f;
	for (auto vi : labelXVs){
		if (IsVInside.get(vi)){
			res += VProVol[vi]*logf(VPro[vi]);
		} else {
			float val = 1 - VPro[vi];
			// numerical issue
			if(val<=0){
				val=0.000001;
			}
			res += VProVol[vi]*logf(val);
		}
	}
	return res;
}

float Arrangement::tetCell_findGroupings_GroupScore_noVol(vector<float> &VPro, vector<float> &VProVol, vector<int> &labelXVs, float iso){
	float res = 0.0f;
	for (auto vi : labelXVs){
		if (VPro[vi] > iso){
			res += logf(VPro[vi]);
		} else {
			res += logf(1 - VPro[vi]);
		}
	}
	return res;
}


void Arrangement::tetCell_findGroupings(int ci, vector<float> &VPro, vector<float> &VProVol,  vector<int> &labelXVs, UnDirGraph &graph, vector<float> &criticalIsoVal, vector<queue<int> > &isoSufFringe, vector<int> &tetV2cyci, bitIntHash& allSets_b){
#if VAXMAN
	allGroupingIso[ci] = criticalIsoVal;

#else

	int nCyc = isoSufFringe.size();
	unordered_map<vector<int>, vector<float>, myVecHash > cellGroupings;

	BigInt actEdge(allCellTets.cells[ci].cellEdges.size());
	BigInt actFace(allCellTets.cells[ci].cellFaces.size());
	BigInt actTet(allCellTets.cells[ci].cellTets.size());
	bitIntHash loc_allSets_b;
	vector<vector<int> > cellEdge2loc_subGroupIDs(allCellTets.cells[ci].cellEdges.size());
	vector<vector<float> > cellEdge2loc_subGroupIsos(allCellTets.cells[ci].cellEdges.size());
	vector<vector<float> > loc_set2isos;
	vector<float> & loc_set2iso = cell2Loc_setIso[ci];
	vector<int> & loc_set2global_set = cell2Loc_Global_setMap[ci];
	int loc_nSets=0;
	for (size_t isoi=0; isoi<criticalIsoVal.size(); ++isoi){		
        UnionFind_o UF(nCyc);
		float isoVal = criticalIsoVal[isoi];
        //if(isoVal<0.15 || isoVal > 0.85)continue;
		tetCell_tagActiveEle(ci, isoVal, VPro, actTet, actFace, actEdge);

		BigInt toVisit = actTet;
		BigInt groupedCyc(nCyc);

		vector<int> tet2SetID(allCellTets.cells[ci].cellTets.size(), -1);
		int setID = -1;
		vector<int> floodOrder;
		while(!toVisit.isZero()){
			setID++;
			int curTet = toVisit.find1st(1);
			queue<int> fringe;
			fringe.push(curTet);
			toVisit.unset(curTet);
			tet2SetID[curTet] = setID;
			
			while (!fringe.empty()){
				curTet = fringe.front();
				fringe.pop();
				floodOrder.push_back(curTet);
				for (size_t i=0; i<4; ++i){
					int nxtTet = allCellTets.cells[ci].cellTets[curTet].nbTet[i];
					if (nxtTet==-1) continue;
					int nxtFace = allCellTets.cells[ci].cellTets[curTet].face[i];
					// if next tet is an unvisited tet and the shared face is active, then flood into nxtTet
					if(toVisit.get(nxtTet) && actFace.get(nxtFace)){
						toVisit.unset(nxtTet);
						tet2SetID[nxtTet] = setID;
						fringe.push(nxtTet);
					}
				}
			}
		}

		float groupingGenus = tetCell_calcGenus(setID+1, nCyc, actEdge.getOnesNum(), actFace.getOnesNum(), actTet.getOnesNum());

		if (groupingGenus!=0.0f){
			continue;
		}

		vector<BigInt> curSets(setID+1, BigInt(maxNCyc));
		for(size_t i=0; i<isoSufFringe.size(); ++i){
			int curV = isoSufFringe[i].front();

			for (auto nbTet : allCellTets.cells[ci].cellVertices[curV].nbTet){
				if (tet2SetID[nbTet]!=-1){
					curSets[tet2SetID[nbTet]].set(i);
					break;
				}
			}
		}

		bool toContinue = false;
		for(int i=0; i<curSets.size(); ++i){
			if(curSets[i].isZero()){
                //cout << "Error! Iso = " << isoVal << " ci = " << ci << endl;
				toContinue = true;
			}
		}
		if(toContinue){
			continue;
		}

		vector<int> curGrouping;
		float curGroupingScore = tetCell_findGroupings_GroupScore(VPro, VProVol, labelXVs, isoVal);//0.0f;
		
		for (size_t i=0; i<curSets.size(); ++i){
			if(allSets_b.find(curSets[i]) == allSets_b.end()){
				allSets_b[curSets[i]] = nSets++;
			}
			//for mixing sets
			if(loc_allSets_b.find(curSets[i]) == loc_allSets_b.end()){
				loc_allSets_b[curSets[i]] = loc_nSets++;
				vector<float> one_loc_sets2isos(1,isoVal);
				loc_set2isos.push_back(one_loc_sets2isos);
				loc_set2global_set.push_back(allSets_b[curSets[i]]);
			}else{
				loc_set2isos[loc_allSets_b[curSets[i]]].push_back(isoVal);
			}
			//end for
			curGrouping.push_back(allSets_b[curSets[i]]);
		}

		sort(curGrouping.begin(), curGrouping.end());
		vector<float> curGroupInfo;
		curGroupInfo.push_back(isoVal);
		curGroupInfo.push_back(curGroupingScore);
		if (cellGroupings.find(curGrouping)==cellGroupings.end()){
			cellGroupings[curGrouping] = curGroupInfo;
		} else if (cellGroupings[curGrouping][TagMapping::GI_SCORE] < curGroupInfo[TagMapping::GI_SCORE]){
			cellGroupings[curGrouping] = curGroupInfo;
		}

		//MING: prepare datastructure of edge->loc sub-grouping id , which helps to calculate mix-iso-grouping scores and gen mix-iso-surfaces
		// only visit edges of actEdge, find one of its nbTet(nbFace[0]->nbTet[0]), add the loc_setID of that tet to the current edge
		// tet2SetID <-> curSets, loc_allSets_b[curSet[i]] = loc sub-grouping id
		BigInt actEdgeCopy(actEdge);
		while(!actEdgeCopy.isZero()){
			int curActEdge = actEdgeCopy.popOne();
			int curActFace = allCellTets.cells[ci].cellEdges[curActEdge].nbFace[0];
			int curActTet = allCellTets.cells[ci].cellFaces[curActFace].nbTet[0] == -1? allCellTets.cells[ci].cellFaces[curActFace].nbTet[1] : allCellTets.cells[ci].cellFaces[curActFace].nbTet[0];
			cellEdge2loc_subGroupIDs[curActEdge].push_back(loc_allSets_b[curSets[tet2SetID[curActTet]]]);
			cellEdge2loc_subGroupIsos[curActEdge].push_back(isoVal);
		}
	}

	int nLocSet = loc_set2isos.size();
	// for each local sets, find a unique iso value to represent it (the iso value that is closest to 0.5)
	loc_set2iso.reserve(loc_set2isos.size());
	for(int i=0; i<loc_set2isos.size(); ++i){
		float repIso=-1;
		float curDiff=0.0f;
		float minDiff=INFINITY;
		for(auto iso : loc_set2isos[i]){
			curDiff = abs(iso-0.5);
			if (curDiff < minDiff){
				repIso = iso;
				minDiff = curDiff;
			}
		}
		loc_set2iso.push_back(repIso);
	}

	// map each loc set to its active edges, not this active edge is not the same of previous activeEdge with one single iso value, 
	// this active edge is usually partial of the activeEdge generated with single iso
	vector<BigInt> & loc_set2actEdge = cell2Loc_setActEdge[ci];
	loc_set2actEdge.resize(loc_set2iso.size(), BigInt(allCellTets.cells[ci].cellEdges.size()));
	for(int i=0; i<allCellTets.cells[ci].cellEdges.size(); ++i){
		allCellTets.cells[ci].cellEdges[i].loc_SubGroupIDs.resize(loc_set2iso.size());
		for(int j=0; j<cellEdge2loc_subGroupIDs[i].size(); ++j){
			int curLocSubGroupID = cellEdge2loc_subGroupIDs[i][j];
			// if true, this loc_subgroup choose this iso value and this edge as its representive 
			// otherwise, skip processing
			if(cellEdge2loc_subGroupIsos[i][j] == loc_set2iso[curLocSubGroupID]){
				allCellTets.cells[ci].cellEdges[i].loc_SubGroupIDs.set(curLocSubGroupID);
				loc_set2actEdge[curLocSubGroupID].set(i);
			}
		}
	}

	// Step 1: find sets of the cycles in cell, and their weights & hard collisions
	vector<vector<int> > loc_sets; // local index of sub-grouping 2 the VECTOR format of this sub-grouping
	vector<BigInt> loc_sets2setsb; // local index of sub-grouping 2 the BIGINT format of this sub-grouping
	vector<float> weights;
	vector<vector<int> > conflicts;
	loc_sets.resize(loc_allSets_b.size());
	loc_sets2setsb.resize(loc_allSets_b.size());
	for (auto set = loc_allSets_b.begin(); set!=loc_allSets_b.end(); ++set){
		loc_sets[set->second]=(set->first.getOnesInd());
		loc_sets2setsb[set->second]=set->first;
	}
	weights.resize(loc_sets.size());
	
	vector<int> cycFirstV(nCyc);
	for(int i=0; i<nCyc; ++i){
		cycFirstV[i] = isoSufFringe[i].front();
	}

	// Step 2: find the top set covers through exact set cover (groupings)tetCell_findGroupings_mix_GroupScore
	vector<vector<int> > mix_groupings;
	vector<float> mix_groupingCosts;
    cout << "start exact set cover: "<<endl;
        AlgoX::AlgoXExactSetCover(nCyc, (int)loc_sets.size(),
            loc_sets, weights, conflicts, mix_groupings, mix_groupingCosts);

        cout << "finish exact set cover: "<< mix_groupings.size() <<endl;

    int n_mix = 0;
	for(int mgi=0; mgi<mix_groupings.size(); ++mgi){
        if(n_mix > 10)continue;
		vector<int> curMixIsoGrouping_loc = mix_groupings[mgi];
		vector<int> curMixIsoGrouping_global = curMixIsoGrouping_loc;

		for(int mgj=0; mgj<curMixIsoGrouping_global.size(); ++mgj){
			// mix_groupings holds the mix-iso-grouping with local indexes
			// map local indexes to global sub-grouping indexes, and save it to mix-curGrouping
			curMixIsoGrouping_global[mgj] = allSets_b[loc_sets2setsb[curMixIsoGrouping_global[mgj]]];
		}
		sort(curMixIsoGrouping_global.begin(), curMixIsoGrouping_global.end());

		// if this mix-iso-grouping is new, store it in the cellGrouping hash
		if (cellGroupings.find(curMixIsoGrouping_global)==cellGroupings.end()){
			vector<float> curGroupInfo;
			curGroupInfo.push_back(MIXISOFAKEVALUE); // invalid iso value, to identify that this iso-grouping is a mix-iso-grouping
			float mix_score = tetCell_findGroupings_mix_GroupScore(ci, curMixIsoGrouping_loc, cycFirstV, loc_set2actEdge, nLocSet, VPro, VProVol, labelXVs);
			curGroupInfo.push_back(mix_score);

			cellGroupings[curMixIsoGrouping_global] = curGroupInfo;
            ++n_mix;
		}
	}


	for (auto it = cellGroupings.begin(); it!= cellGroupings.end(); ++it){
		if(it->first.size()==0) continue;
		allGroupings[ci].push_back(it->first);
		allGroupingScore[ci].push_back(it->second[TagMapping::GI_SCORE]);
		allGroupingIso[ci].push_back(it->second[TagMapping::GI_ISO]);
	}
    cout << "n_mix: "<<n_mix<<endl;

#if _SAVEDATA_MING
	Utility::writeVector(allGroupings[ci], (t_savedir+"allGroupings_"+Utility::toStr(ci)+".txt").c_str());
	Utility::writeVector(allGroupingScore[ci], (t_savedir+"allGroupingScore_"+Utility::toStr(ci)+".txt").c_str());
	Utility::writeVector(allGroupingIso[ci], (t_savedir+"allGroupingIso_"+Utility::toStr(ci)+".txt").c_str());
#endif

#endif //end vaxmanIF
}


void Arrangement::tetCell_tagActiveEdge(int ci, float iso, BigInt &actEdge){
	for (int i=0; i<actEdge.size(); ++i){
		if ((allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellEdges[i].v[0]]<=iso)^(allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellEdges[i].v[1]]<=iso)){
			actEdge.set(i);
		} else {
			actEdge.unset(i);
		}
	}
}

void Arrangement::tetCell_tagActiveEle(int ci, float iso, vector<float> &VPro, BigInt &actTet, BigInt &actFace, BigInt &actEdge){
	for (int i=0; i<actEdge.size(); ++i){
		if ((allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellEdges[i].v[0]]<=iso)^(allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellEdges[i].v[1]]<=iso)){
			actEdge.set(i);
		} else {
			actEdge.unset(i);
		}
	}
	for (int i=0; i<actFace.size(); ++i){
		if (((allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellFaces[i].v[0]]<=iso)^(allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellFaces[i].v[1]]<=iso))||
			((allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellFaces[i].v[0]]<=iso)^(allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellFaces[i].v[2]]<=iso))||
			((allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellFaces[i].v[1]]<=iso)^(allCellTets.cells[ci].VPro[allCellTets.cells[ci].cellFaces[i].v[2]]<=iso))
			){
			actFace.set(i);
		}else{
			actFace.unset(i);
		}
	}
	vector<int> curTetVs(4);
	for (int i=0; i<actTet.size(); ++i){
		curTetVs[0] = allCellTets.cells[ci].cellTets[i].v[0];
		curTetVs[1] = allCellTets.cells[ci].cellTets[i].v[1];
		curTetVs[2] = allCellTets.cells[ci].cellTets[i].v[2];
		curTetVs[3] = allCellTets.cells[ci].cellTets[i].v[3];
		if(((allCellTets.cells[ci].VPro[curTetVs[0]]<=iso)&&(allCellTets.cells[ci].VPro[curTetVs[1]]<=iso)&&(allCellTets.cells[ci].VPro[curTetVs[2]]<=iso)&&(allCellTets.cells[ci].VPro[curTetVs[3]]<=iso))||
			((allCellTets.cells[ci].VPro[curTetVs[0]]>iso)&&(allCellTets.cells[ci].VPro[curTetVs[1]]>iso)&&(allCellTets.cells[ci].VPro[curTetVs[2]]>iso)&&(allCellTets.cells[ci].VPro[curTetVs[3]]>iso))){
				actTet.unset(i);
		}else{
			actTet.set(i);
		}
	}
} 
int Arrangement::tetCell_getMatOfFaceNextToSeg(int ci, int fi, int ei, vector<int> &nxtSegVID){
	int eV1, eV2, fV3;
	eV1 = allCellTets.cells[ci].cellEdges[ei].v[0];
	eV2 = allCellTets.cells[ci].cellEdges[ei].v[1];
	fV3 = allCellTets.cells[ci].cellFaces[fi].v[0] 
		+ allCellTets.cells[ci].cellFaces[fi].v[1]
		+ allCellTets.cells[ci].cellFaces[fi].v[2]
		- eV1 - eV2;
	int curSegi, curFaceti, curPi;
	curSegi = TagMapping::getSegiFromEdgemark(allCellTets.cells[ci].cellEdges[ei].EtetMark);
	curFaceti = seg2F[curSegi];

	curPi = F2PlaneID[curFaceti];
	Eigen::Vector3f n(planeParas[curPi][0], planeParas[curPi][1], planeParas[curPi][2]);

	// here need to chekc whether eV1, or eV2 is the end of a seg, because if so, they maybe shared by multiple segs, so, cannot be evaluated by nxtSegVID at this point, apply on the other end point instead
	// if both of them are on cell edge (shared by other segs)
	if ((nxtSegVID[eV1]==TagMapping::NXTSEGV_SHARE) && (nxtSegVID[eV2]==TagMapping::NXTSEGV_SHARE)){
		cout << "MING: should never happen, we handled this !!" << endl;
	}
	else if(((nxtSegVID[eV1]!=TagMapping::NXTSEGV_SHARE)&&(nxtSegVID[eV1]!=eV2))||((nxtSegVID[eV2]!=TagMapping::NXTSEGV_SHARE)&&(nxtSegVID[eV2]==eV1))){
		int tmp_eV = eV1;
		eV1 = eV2;
		eV2 = tmp_eV;
	}

	Eigen::Vector3f v3(allCellTets.cells[ci].cellVertices[fV3].xyz[0],allCellTets.cells[ci].cellVertices[fV3].xyz[1],allCellTets.cells[ci].cellVertices[fV3].xyz[2]);
	Eigen::Vector3f v1(allCellTets.cells[ci].cellVertices[eV1].xyz[0],allCellTets.cells[ci].cellVertices[eV1].xyz[1],allCellTets.cells[ci].cellVertices[eV1].xyz[2]);
	Eigen::Vector3f v2(allCellTets.cells[ci].cellVertices[eV2].xyz[0],allCellTets.cells[ci].cellVertices[eV2].xyz[1],allCellTets.cells[ci].cellVertices[eV2].xyz[2]);

	if(isVOnLeftSideOfSeg(v3, v1, v2, n)){
		return seg2lmat[curSegi];
	} else {
		return TagMapping::getOtherMat(seg2lmat[curSegi], matMarks);
	}
}

float Arrangement::getVolWeight(float xyzi [], float xyzj []){
	//if there is no vol info, use 1 as the weight
	if(volInfo.empty()||t_randwb==0.0f){
		return 1.0f; 
	}else{
		float gi = getVolValAtPos(xyzi);
		float gj = getVolValAtPos(xyzj);
		return exp(-1.0f*t_randwb*pow(gi-gj,2));
	}
}
float Arrangement::getVolValAtPos(float xyz[]){

	float xyz2 [3];
	xyz2[0] = xyz[0];
	xyz2[1] = xyz[1];
	xyz2[2] = xyz[2];

	int xyz0 [3];
	int xyz1 [3];
	float xyzD1 [3];
	float xyzD0 [3];
	for(int i=0; i<3; ++i){
		float val = 1.0f*(xyz2[i]-volInfo.getLowerCorner(i))/volInfo.getUnit(i);
		val = val < 0.0f ? 0.0f : val; //need this, because of numerical issue
		xyz0[i] = floor(val);
		xyz1[i] = xyz0[i] == (volInfo.getGridSize(i)-1) ? (volInfo.getGridSize(i)-1) : xyz0[i]+1;
		xyzD1[i] = val - xyz0[i];
		xyzD0[i] = 1 - xyzD1[i];
	}
	float val = 
		volInfo.getDataAt(xyz0[0],xyz0[1],xyz0[2])*xyzD0[0]*xyzD0[1]*xyzD0[2]
	+   volInfo.getDataAt(xyz1[0],xyz0[1],xyz0[2])*xyzD1[0]*xyzD0[1]*xyzD0[2]
	+   volInfo.getDataAt(xyz1[0],xyz1[1],xyz0[2])*xyzD1[0]*xyzD1[1]*xyzD0[2]
	+   volInfo.getDataAt(xyz0[0],xyz1[1],xyz0[2])*xyzD0[0]*xyzD1[1]*xyzD0[2]
	+   volInfo.getDataAt(xyz0[0],xyz0[1],xyz1[2])*xyzD0[0]*xyzD0[1]*xyzD1[2]
	+   volInfo.getDataAt(xyz1[0],xyz0[1],xyz1[2])*xyzD1[0]*xyzD0[1]*xyzD1[2]
	+   volInfo.getDataAt(xyz1[0],xyz1[1],xyz1[2])*xyzD1[0]*xyzD1[1]*xyzD1[2]
	+   volInfo.getDataAt(xyz0[0],xyz1[1],xyz1[2])*xyzD0[0]*xyzD1[1]*xyzD1[2];

	return val;
}

void Arrangement::tetCell_tetgen_scaleUp(tetgenio &tetin, float scale, float translate []){
	for(int i=0; i<tetin.numberofpoints; ++i){
		tetin.pointlist[3*i] = (tetin.pointlist[3*i] + translate[0])*scale;
		tetin.pointlist[3*i+1] = (tetin.pointlist[3*i+1] + translate[1])*scale;
		tetin.pointlist[3*i+2] = (tetin.pointlist[3*i+2] + translate[2])*scale;
	}
}
void Arrangement::tetCell_tetgen_scaleBack(tetgenio &tetout, float scale, float translate []){
	for(int i=0; i<tetout.numberofpoints; ++i){
		tetout.pointlist[3*i] = tetout.pointlist[3*i]/scale - translate[0];
		tetout.pointlist[3*i+1] = tetout.pointlist[3*i+1]/scale - translate[1];
		tetout.pointlist[3*i+2] = tetout.pointlist[3*i+2]/scale - translate[2];
	}
}
