#include "GroupingManager.h"

GroupingManager::GroupingManager() : goalGenus(0), maxGenus(15) {}
GroupingManager::~GroupingManager(){}

// main grouping algorithm
void GroupingManager::Grouping(){
#if VAXMAN
	for (int fti = 0; fti < ar.nCell; ++fti){
		if (ar.cellEmpty[Order[fti]]){
			continue;
		}
		resGrouping.push_back(0);
	}
	solGrouping.push_back(vector<int>(1));
#else
	vector<GState> f0States;
	vector<GState> f1States;
	GStateHash f0StateHash;
	GStateHash f1StateHash;
	vector<vector<int> > f0SolGrouping;
	vector<vector<int> > f1SolGrouping;
	vector<float> f0SolGroupingCost;
	vector<float> f1SolGroupingCost;

	solGrouping.clear();

	maxGenus_log = (int)ceil(log(1.0f + maxGenus)/log(2.0));
	nfrontcycs_log.resize(ar.nCell, 1);
	for (int fti = 0; fti < ar.nCell; ++fti){
		if (!ar.front2cyc[fti].empty()){
			nfrontcycs_log[fti] = max(1, (int)ceil(log((float)ar.front2cyc[fti].size())/log(2.0)));
		}
	}
	ncellcycs_log.resize(ar.nCell, 1);
	for (int celli = 0; celli < ar.nCell; ++celli){
		if (ar.cell2ncyc[celli] != 0){
			ncellcycs_log[celli] = max(1, (int)ceil(log((float)ar.cell2ncyc[celli])/log(2.0)));
		}
	}

	// no valid representation for initState and lastState
	// they are specifically made up in the program
	initState.resize(0);
	lastState.resize(0);

	int initGrouping = -1;
	f1SolGrouping.resize(1);
	f1SolGrouping[0].push_back(initGrouping);
	f1SolGroupingCost.push_back(0);
	f1StateHash[initState] = 0;
	f1States.push_back(initState);
	for (int fti = 0; fti < ar.nCell; ++fti){
		if (ar.cellEmpty[Order[fti]]){
			// cout << "empty cell, skip..." << endl;
			continue;
		}
		swap(f0StateHash, f1StateHash);
		swap(f0States, f1States);
		swap(f0SolGrouping, f1SolGrouping);
		swap(f0SolGroupingCost, f1SolGroupingCost);
		f1StateHash.clear();
		f1States.clear();
		f1SolGrouping.resize(1);	// the first sol is left as invalid
		f1SolGroupingCost.resize(1);
		for (int sti = 0; sti < f0States.size(); ++sti){
			Propagate(f0States[sti], fti, f1States,
				f0StateHash, f1StateHash,
				f0SolGrouping, f1SolGrouping,
				f0SolGroupingCost, f1SolGroupingCost);
		}
        cout << f0States.size() << " --> " << f0States.size()*ar.allGroupings[Order[fti]].size() << " --> " << f1States.size()<<endl;
	}
	// finished serching
	// must left with one solution (with a fake one)
	if (f1SolGroupingCost.size() == 2){
		resGrouping = f1SolGrouping[1];
		resGrouping.erase(resGrouping.begin());
		for (int gi = 1; gi < f1SolGrouping.size(); ++gi){
			for (int ggi = 1; ggi < f1SolGrouping[gi].size(); ++ggi){
				int celli = ar.nonEmpty_Order[ggi - 1];

				vector<int> grouping = ar.allGroupings[celli][f1SolGrouping[gi][ggi]];
				for (auto seti : grouping){
					vector<int> oneSet = ar.allSets[seti];
					for (size_t osi=0; osi<oneSet.size(); ++osi){
						oneSet[osi] = ar.cell2cyc[celli][oneSet[osi]];
					}
					solGrouping.push_back(oneSet);
				}
			}
		}
	}
	else{
		cout << "No Solution" << endl;
	}

#endif

	#if _SAVEDATA_MING
		writeSolGrouping((outDir+"solGrouping_"+Utility::toStr(goalGenus)+".txt").c_str());
	#endif

}

void GroupingManager::GenSurface(){
#if USETRIANGULATION
    //ar.sufCell_genSurface_Triangulation(resGrouping);
#else
    //ar.sufCell_genSurface(resGrouping);
#endif
}

void GroupingManager::SmoothSurface(){
    //ar.sufCell_smoothSurface();
}
void GroupingManager::SaveResBefSmooth(){
    //ar.sufCell_saveFinalSufForSmooth(0);
}

// a GState = {k, S1,...,Sn, g1,...,gk, 000}
void GroupingManager::PackGState(int fronti, vector<int> &finalSufs, vector<int> &finalSufs_genus, vector<BigInt> &sufBoundaries, GState &nxtState){

	int n = (int)ar.front2cyc[fronti].size();
	int k = (int)finalSufs.size();
	int nbit_k = nfrontcycs_log[fronti] + 1;
	int nbit_S = nfrontcycs_log[fronti];
	int nbit_g = maxGenus_log;
	nxtState.resize(nbit_k + n*(nbit_S + nbit_g));

	vector<int> Ss(n);
	vector<int> gs(k);
	vector<int> Cs;
	for (int si = 0; si < finalSufs.size(); ++si){
		Cs = sufBoundaries[finalSufs[si]].getOnesInd();
		for (int ci : Cs){
			Ss[ci] = si;
		}
	}
	// need to keep a unique order of surface for s1...sn, to avoid duplicated states
	// e.g. from small to big
	vector<int> umap;
	FindUniqueOrderMap(k, Ss, umap);
	for (int i = 0; i < n; ++i){
		Ss[i] = umap[Ss[i]];
	}
	for (int i = 0; i < k; ++i){
		gs[i] = finalSufs_genus[umap[i]];
	}

	nxtState.pack(nbit_k, nbit_S, nbit_g, n, Ss, k, gs);
}

void GroupingManager::UnPackGState(int fronti, GState &state, vector<vector<int> > &suf2cycid, vector<int> &suf2genus){

	if (state.nBit == 0){ return; }

	int n = (int)ar.front2cyc[fronti - 1].size();
	int k = -1;
	int nbit_k = nfrontcycs_log[fronti - 1] + 1;
	int nbit_S = nfrontcycs_log[fronti - 1];
	int nbit_g = maxGenus_log;
	vector<int> Ss;
	vector<int> gs;
	state.unpack(nbit_k, nbit_S, nbit_g, n, Ss, k, suf2genus);
	suf2cycid.resize(k);
	for (int ci = 0; ci < n; ++ci){
		suf2cycid[Ss[ci]].push_back(ci);
	}
}

void GroupingManager::PackGGrouping(int celli, vector<int> &curGrouping, vector<vector<int> > &sets, GGrouping &ggroup){
	int n = ar.cell2ncyc[celli];
	ggroup.resize(n);
	vector<int> Gs(n);
	for (int gi = 0; gi < curGrouping.size(); ++gi){
		for (int ci : sets[curGrouping[gi]]){
			Gs[ci] = gi;
		}
	}
	ggroup.pack(ncellcycs_log[celli], Gs);
}

void GroupingManager::UnPackGGrouping(int celli, GGrouping &ggroup, vector<vector<int> > &grouping){

	if (ggroup.nBit == 0){ return; }

	int n = ar.cell2ncyc[celli];
	int nbit_G = ncellcycs_log[celli];
	vector<int> Gs;

	ggroup.unpack(nbit_G, n, Gs);
	int k = -1;
	for (int gi : Gs){
		k = gi > k ? gi : k;
	}
	k++;
	grouping.resize(k);
	for (int gi = 0; gi < n; ++gi){
		grouping[Gs[gi]].push_back(gi);
	}
}



void GroupingManager::FindUniqueOrderMap(int maxSi, vector<int> &Ss, vector<int> &map){
	map.resize(maxSi, -1);
	int id = 0;
	for (int si : Ss){
		if (map[si] == -1){
			map[si] = id;
			id++;
		}
	}
}

void GroupingManager::GenerateAllSingleDoubleTripleCurveSets(int nCyc, vector<vector<int> > &sets){
	int nSet = nCyc * (nCyc*nCyc + 5) / 6;
	sets.resize(nSet);
	// single
	for (int i = 0; i < nCyc; ++i){
		sets[i].push_back(i);
	}
	// double
	int id = nCyc;
	for (int i = 0; i < nCyc; ++i){
		for (int j = i + 1; j < nCyc; ++j){
			sets[id].push_back(i);
			sets[id].push_back(j);
			id++;
		}
	}
	// triple
	for (int i = 0; i < nCyc; ++i){
		for (int j = i + 1; j < nCyc; ++j){
			for (int k = j + 1; k < nCyc; ++k){
				sets[id].push_back(i);
				sets[id].push_back(j);
				sets[id].push_back(k);
				id++;
			}
		}
	}
}

void GroupingManager::CalculateSetCost(int celli, vector<vector<int> > &sets, vector<float> &weights){
	weights.resize(sets.size());
	//TODO, calculate weights. For now, the the weights are the same, 0.0f;
}

void GroupingManager::DetectConflicts(int celli, vector<vector<int> > &sets, vector<vector<int> > &conflicts){
	//TODO: find hard colisions. For now, no conflicts
}

void GroupingManager::Contour2Arrangement(const char* ctrfilename, const char* volfilename, const char* bboxfilename){
	ar.setPara(outDir, t_teta, t_randwb, t_nSmBefLoop, t_nLoop, t_nSmInLoop, goalGenus);
	ar.setFilename(ctrfilename);
	Ctr2ArHandler::readContour(ctrfilename, volfilename, bboxfilename, ar);
}

void GroupingManager::Contour2ArrangementMM(const char* ctrfilename, const char* volfilename, const char* bboxfilename){
    ar.setPara(outDir, t_teta, t_randwb, t_nSmBefLoop, t_nLoop, t_nSmInLoop, goalGenus);
    ar.setFilename(ctrfilename);
    Ctr2ArHandler::processContourMM(ctrfilename, ar);
}
// find the "best" search order among the cell graph
// which minimized the max number of cycles on each front.
// for now, it is just a BFS traversal of the graph.
// TODO: implement the search for the real "best" order 
// (minimize width of each layer)
void GroupingManager::CalcSearchOrder(){
	vector<bool> visited(ar.nCell, false);
	queue<int> fringe;
	int curCell;

	fringe.push(0);
	visited[0] = true;
	while (!fringe.empty()){
		curCell = fringe.front();
		fringe.pop();
		Order.push_back(curCell);
		for (int nxtCell : ar.cellGraph[curCell]){
			if (!visited[nxtCell]){
				visited[nxtCell] = true;
				fringe.push(nxtCell);
			}
		}
	}
	ar.findCloseCyclesInFronts(Order);
#if _SAVEDATA_MING
	writeSearchOrder((outDir+"order.txt").c_str());
	ar.writeArrangement((outDir+"arrangement.txt").c_str());
#endif
}

void GroupingManager::writeSearchOrder(const char* filename){
	ofstream ofs(filename, ofstream::out);
	ofs << "{"; Utility::writeVector(Order, ofs); ofs << "}";
	ofs.close();
}

void GroupingManager::writeSolGrouping(const char* filename){
#if VAXMAN
#else
	ofstream ofs(filename, ofstream::out);
	ofs << "{"; Utility::writeVector(solGrouping, ofs); ofs << "}";
	ofs.close();
#endif
}

void GroupingManager::writeOutputs(const char* filename, vector<float> &data){
	ofstream ofs(filename, ofstream::out);
	ofs << "{"; Utility::writeVector(data, ofs); ofs << "}";
	ofs.close();
}
