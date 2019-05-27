#include "GroupingManager.h"

void GroupingManager::Propagate(
	GState &f0State, int fronti,
	vector<GState> &f1States,
	GStateHash &f0StateHash,
	GStateHash &f1StateHash,
	vector<vector<int> > &f0SolGrouping,
	vector<vector<int> > &f1SolGrouping,
	vector<float> &f0SolGroupingCost,
	vector<float> &f1SolGroupingCost
	){

	int celli = Order[fronti];
	int lastFrontID = ar.nCell - 1; 
	// find the last front id
	while (ar.cellEmpty[Order[lastFrontID]]){
		lastFrontID--;
	}
	vector<int> cellcyc = ar.cell2cyc[celli];
	vector<int> f0cyc;
	if (fronti != 0){ // for init state, no f0cyc
		f0cyc = ar.front2cyc[fronti - 1];
	}
	int nf0cyc = (int)f0cyc.size();
	int nf1cyc = (int)ar.front2cyc[fronti].size();
	int ncellcyc = (int)cellcyc.size();

	vector<vector<int> > &sets = ar.allSets;
	vector<vector<int> > &groupings = ar.allGroupings[celli];
	vector<float> &groupingCosts = ar.allGroupingScore[celli];
	// Step 3: for each grouping, generate the next state
	// 3.1 | Parse state from bit representation to format easier to process
	vector<vector<int> > ssuf2cyc;
	vector<int> ssuf2genus;
	UnPackGState(fronti, f0State, ssuf2cyc, ssuf2genus);
	int nSSuf = (int)ssuf2genus.size();		// # of existing surfaces in f0 state
	int nGSuf = 0;						// # of grouped surface in the new grouping of celli
	int nSuf = 0;						// = nSSuf + nGSuf

    // 3.2 | init UnionFind_o/suf attributes/suf Boudary, for f0 state info
	vector<int> cyc2suf_i(nf0cyc + ncellcyc);
    UnionFind_o UF_i(nf0cyc + ncellcyc);
	vector<SufAttribute> sufAttributes_i(nSSuf);
	vector<BigInt> sufBoundaries_i(nSSuf, BigInt(nf1cyc));
	for (int si = 0; si < nSSuf; ++si){
		for (int ci : ssuf2cyc[si]){
			cyc2suf_i[ci] = si;
			sufBoundaries_i[si] |= ar.front_B_f0cyc2f1cyc[fronti][ci];
		}
		for (int ci = 1; ci < ssuf2cyc[si].size(); ++ci){
			UF_i.Union(ssuf2cyc[si][0], ssuf2cyc[si][ci]);
		}
		sufAttributes_i[si].k = 1;
		sufAttributes_i[si].sum_Bi = (int)ssuf2cyc[si].size();
		sufAttributes_i[si].sum_gi = ssuf2genus[si];
	}

	// 3.3 | for each grouping (exact set cover), generate the nxt state
	vector<int> curGrouping;
	float curCost;
	for (int gi = 0; gi < groupings.size(); ++gi){
		curGrouping = groupings[gi]; // seti
		curCost = groupingCosts[gi];
		nGSuf = (int)curGrouping.size();
		nSuf = nSSuf + nGSuf; 

        // 3.3.1 | Copy the init UnionFind_o/Attributs/Boundary info
		// add in the info from gsuf
		vector<int> cyc2suf = cyc2suf_i;
        UnionFind_o UF = UF_i;
		vector<SufAttribute> sufAttributes = sufAttributes_i;
		vector<BigInt> sufBoundaries = sufBoundaries_i;
		sufAttributes.resize(nSuf);
		sufBoundaries.resize(nSuf);

		for (int si = 0; si < nGSuf; ++si){
			BigInt f1Boundary(nf1cyc);
			for (int ci : sets[curGrouping[si]]){
				cyc2suf[ci + nf0cyc] = si + nSSuf;
				f1Boundary |= ar.front_B_cellcyc2f1cyc[fronti][ci];
			}

			for (int ci = 1; ci < sets[curGrouping[si]].size(); ++ci){
				UF.Union(sets[curGrouping[si]][0] + nf0cyc, sets[curGrouping[si]][ci] + nf0cyc);
			}
			sufAttributes[si + nSSuf].k = 1;
			sufAttributes[si + nSSuf].sum_Bi = (int)sets[curGrouping[si]].size();
			sufBoundaries[si + nSSuf] = f1Boundary;
		}

		// 3.3.2 | 
		// 1) For each gsuf, find the f0 cycles that share seg with it
		// 2) UF.find the corresponding ssufs; Union all ssufs and this gsuf
		// 3) Calculate the "sum" of attributes of those ssufs and this gsuf
		BigInt sharef0cycs_b(nf0cyc);
		BigInt sufsToGroup_b(nSuf);
		vector<int> sharef0cycs;
		for (int si = 0; si < nGSuf; ++si){

			sharef0cycs_b.clear();
			sufsToGroup_b.clear();
			int sumL = 0;
			// 1) 
			int tmpCount = 0;
			for (int ci : sets[curGrouping[si]]){
				tmpCount++;
				sharef0cycs_b |= ar.front_B_cellcyc2f0cyc[fronti][ci];
				sumL += ar.front_L_cellcyc2nf0cos[fronti][ci];
			}
			sharef0cycs = sharef0cycs_b.getOnesInd();
			// by default, the unioned surface is current gsuf
			// if exist ssufs (|sharef0cycs| > 0), union them with gsuf
			int unionedCi = UF.Find(si + nSSuf);
			int unionedSuf;
			if (!sharef0cycs.empty()){
				// 2)
				sufsToGroup_b.set(cyc2suf[UF.Find(sharef0cycs[0])]);
				for (int ci = 1; ci < sharef0cycs.size(); ++ci){
					sufsToGroup_b.set(cyc2suf[UF.Find(sharef0cycs[ci])]);
				}
				for (int ci = 1; ci < sharef0cycs.size(); ++ci){
					UF.Union(sharef0cycs[0], sharef0cycs[ci]);
				}
				sufsToGroup_b.set(cyc2suf[UF.Find(sets[curGrouping[si]][0] + nf0cyc)]);
				unionedCi = UF.Union(sharef0cycs[0], sets[curGrouping[si]][0] + nf0cyc);
				unionedSuf = cyc2suf[unionedCi];
				sufsToGroup_b.unset(unionedSuf);
				vector<int> ssufsToGroup = sufsToGroup_b.getOnesInd();
				// 3)
				// between ssufs, directly add all the attributes and boundaries
				// between ssufs and gsuf, "add" the additional attribute L
				for (int i : ssufsToGroup){
					sufAttributes[unionedSuf] += sufAttributes[i];
					sufBoundaries[unionedSuf] |= sufBoundaries[i];
				}
				sufAttributes[unionedSuf].L += sumL;
			}
		}
		// 3.3.3 |
		// get the unioned sufs; calculate its genus; 
		// if any of sufs is invalid: 1) uses none of f1 cyc; or 2) genus already exceeds the bound
		// otherwise, they are valid, create one nxtState
		BigInt finalSufs_b(nSuf);
		vector<int> finalSufs;
		vector<int> finalSufs_genus;
		bool isValid = true;
		for (int ci = 0; ci < nf0cyc + ncellcyc; ++ci){
			finalSufs_b.set(cyc2suf[UF.Find(ci)]);
		}
		finalSufs = finalSufs_b.getOnesInd();
		finalSufs_genus.reserve(finalSufs.size());
		if (fronti != lastFrontID){
			for (int sufi : finalSufs){
				// invalid: 1) isolated surface that do not use f1cycs
				if (sufBoundaries[sufi].isZero()){
					isValid = false;
					break;
				}
				int attr_B = sufBoundaries[sufi].getOnesNum();
				finalSufs_genus.push_back(sufAttributes[sufi].calcGenus(attr_B));
				// invalid: 2) genus exceeds the bound
				if (finalSufs_genus.back() > goalGenus){
					isValid = false;
					break;
				}
			}
		}
		// for last front
		else{
			// invalid: 1) left with multiple surfaces
			if (finalSufs.size() != 1){
				isValid = false;
				continue;
			}
			int attr_B = sufBoundaries[0].getOnesNum();
			finalSufs_genus.push_back(sufAttributes[0].calcGenus(attr_B));
			// invalid: 2) genus is no equavalent to goalGenus
			if (finalSufs_genus.back() != goalGenus){
				isValid = false;
				continue;
			}
		}
		
		if (isValid){
			GState nxtState;
			if (fronti != lastFrontID){
				PackGState(fronti, finalSufs, finalSufs_genus, sufBoundaries, nxtState);
			}
			else{
				nxtState = lastState;
			}

			float f1cost = f0SolGroupingCost[f0StateHash[f0State]] + curCost;
			if (f1StateHash[nxtState] == 0){
				f1States.push_back(nxtState);
				f1StateHash[nxtState] = (int)f1SolGroupingCost.size();
				GGrouping ggroup;
				PackGGrouping(celli, curGrouping, sets, ggroup);
				f1SolGrouping.push_back(f0SolGrouping[f0StateHash[f0State]]);
				f1SolGrouping.back().push_back(gi);
				f1SolGroupingCost.push_back(f1cost);
			}
			else {
				if (f1cost > f1SolGroupingCost[f1StateHash[nxtState]]){
					f1SolGroupingCost[f1StateHash[nxtState]] = f1cost;
					GGrouping ggroup;
					PackGGrouping(celli, curGrouping, sets, ggroup);
					f1SolGrouping[f1StateHash[nxtState]] = f0SolGrouping[f0StateHash[f0State]];
					f1SolGrouping[f1StateHash[nxtState]].push_back(gi);
				}
			}

		}
	}
}
