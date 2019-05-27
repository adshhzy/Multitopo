#include "Arrangement.h"

Arrangement::Arrangement() : inputVol(false),inputBBox(false){}
Arrangement::~Arrangement(){}

void Arrangement::findCloseCyclesInCells(){

	seg2cyc.resize(nSeg, -1);
	vector<int> oneCyc;
	BigInt oneCyc_b(nSeg);

	// create cycles for closed segs
	for (int seg = 0; seg < nSeg; ++seg){
		if (isSegClosed(seg)){
			oneCyc.clear(); oneCyc_b.clear();
			oneCyc.push_back(seg);
			oneCyc_b.set(seg);

			seg2cyc[seg] = (int)cyc2seg.size();
			cyc2seg.push_back(oneCyc);
			cyc2seg_b.push_back(oneCyc_b);
		}
	}

	// go through all cells, create cycles that are formed by multiple segs
	// and store the list of cycles of the cell into cell2cyc[celli]
	cell2cyc.resize(nCell);
	cell2ncyc.resize(nCell);
	vector<bool> visited(nSeg, false);
	BigInt segs_b, nxtSeg_b;
	nCycs = 0;
	for (int celli = 0; celli < nCell; ++celli){
		segs_b = cell2seg_b[celli];
		traceCycsFromSegs(segs_b, cell2cyc[celli]);
		cell2ncyc[celli] = (int)(cell2cyc[celli].size());
	}
	nCycs = (int)cyc2seg.size();
}

void Arrangement::findCloseCyclesInFronts(vector<int>& _order){
	Order = _order;
	for (int oi : Order){
		if (!cellEmpty[oi]){
			nonEmpty_Order.push_back(oi);
		}
	}
	/*the first front is the first cell in order*/
	front2cyc.resize(nCell);
	front2seg_b.resize(nCell, BigInt(nSeg));

	BigInt preSegs_b = cell2seg_b[Order[0]];
	front2cyc[0] = cell2cyc[Order[0]];
	front2seg_b[0] = preSegs_b;

	vector<int> curCycs;
	BigInt curSegs_b, cellSegs_b, cmnSegs_b;

	// for each front, trace segs to form cycles, store that to cyc2seg
	// and store the list of cycles into Fronts[fronti]
	// similar to findCloseCyclesInCell()
	int celli;
	for (size_t i = 1; i < Order.size(); ++i){
		celli = Order[i];
		// if cell is empty, copy the result from the last valid front
		if (cellEmpty[celli]){
			front2cyc[i] = front2cyc[i - 1];
			front2seg_b[i] = front2seg_b[i - 1];
		}
		else{
			cellSegs_b = cell2seg_b[celli];
			cmnSegs_b = preSegs_b & cellSegs_b;
			curSegs_b = preSegs_b ^ cellSegs_b;
			traceCycsFromSegs(curSegs_b, front2cyc[i]);
			front2seg_b[i] = curSegs_b;

			preSegs_b = curSegs_b;
		}
	}
	nCycs = (int)cyc2seg.size();
	preCalcFront_L_B();
	// this hash has fulfilled its duty, time to retire
	cycsegb2cycid.clear();
}

// reorder the cell edges in each cell face (i.e. change the order of the vector of cellFs[fi])
// make the edges conseg
void Arrangement::reOrderCellEinCellF(){
	vector<int> faceE;
	unordered_map<int, int> faceV2E;
	for (int ci = 0; ci < nCell; ++ci){
		for (int fi : cell2F[ci]){
			faceV2E.clear();
			faceE = cellFs[fi];
			for (int ei : faceE){
				faceV2E[cellEs[ei][0]] += ei;
				faceV2E[cellEs[ei][1]] += ei;
			}
			int lasteind = cellFs[fi][0];
			int lastendv = cellEs[lasteind][1];
			for (int ei = 1; ei < faceE.size(); ++ei){
				cellFs[fi][ei] = faceV2E[lastendv] - lasteind;
				lasteind = cellFs[fi][ei];
				lastendv = cellEs[lasteind][0] + cellEs[lasteind][1] - lastendv;
			}
		}
	}
}

void Arrangement::preCalcFront_L_B(){
	int celli, nf1cycs, nf0cycs, cellcyc, f0cyc, f1cyc;
	vector<int> cellcycs, f0cycs, f1cycs;
	BigInt cmnSegs_b;
	front_L_cellcyc2nf0cos.resize(nCell);
	front_B_cellcyc2f1cyc.resize(nCell);
	front_B_cellcyc2f0cyc.resize(nCell);
	front_B_f0cyc2f1cyc.resize(nCell);
	// for each front
	for (int fti = 0; fti < nCell; ++fti){

		celli = Order[fti];
		cellcycs = cell2cyc[celli];
		if (fti > 0){
			f0cycs = front2cyc[fti - 1];
		}
		f1cycs = front2cyc[fti];
		nf1cycs = (int)f1cycs.size();
		nf0cycs = (int)f0cycs.size();
		front_B_cellcyc2f1cyc[fti].resize(cellcycs.size(), BigInt(nf1cycs));
		front_B_cellcyc2f0cyc[fti].resize(cellcycs.size(), BigInt(nf0cycs));
		front_B_f0cyc2f1cyc[fti].resize(f0cycs.size(), BigInt(nf1cycs));
		front_L_cellcyc2nf0cos[fti].resize(cellcycs.size(), 0);

		// calc front_B_f0cyc2f1cyc[fti]
		// for each f0 cyc, find the cycs on f1 that share seg with it
		for (int i = 0; i < f0cycs.size(); ++i){
			f0cyc = f0cycs[i];
			for (int j = 0; j < f1cycs.size(); ++j){
				f1cyc = f1cycs[j];
				if (cycShareSeg(f0cyc, f1cyc)){
					front_B_f0cyc2f1cyc[fti][i].set(j);
				}
			}
		}

		// calc front_B_cellcyc2f1cyc[fti] & front_B_cellcyc2f0cyc[fti]
		// for each cell cyc, find the cycs on f1/f0 that share seg with it
		for (int i = 0; i < cellcycs.size(); ++i){
			cellcyc = cellcycs[i];
			for (int j = 0; j < f1cycs.size(); ++j){
				f1cyc = f1cycs[j];
				if (cycShareSeg(cellcyc, f1cyc)){
					front_B_cellcyc2f1cyc[fti][i].set(j);
				}
			}
			for (int j = 0; j < f0cycs.size(); ++j){
				f0cyc = f0cycs[j];
				if (cycShareSeg(cellcyc, f0cyc)){
					front_B_cellcyc2f0cyc[fti][i].set(j);
				}
			}
		}

		// calc front_L_cellcyc2snum[fti]
		// for each cell cycle, count # of continous open seg shared with f0 
		for (int i = 0; i < cellcycs.size(); ++i){
			cellcyc = cellcycs[i];
			for (int j = 0; j < f0cycs.size(); ++j){
				f0cyc = f0cycs[j];
				cmnSegs_b = cyc2seg_b[cellcyc] & cyc2seg_b[f0cyc];
				front_L_cellcyc2nf0cos[fti][i] += countContinousOpenSeg(cmnSegs_b);
			}
		}
	}
}

// # of continous oepn seg = # of open ends / 2
// assume the data is "good", meaning there is no tree branch in the segments
int Arrangement::countContinousOpenSeg(BigInt &cmnSegs_b){
	int nOpenEnds = 0;
	vector<int> cmnSegs = cmnSegs_b.getOnesInd();
	unordered_map<int, int> tmpmap;
	for (int i = 0; i < cmnSegs.size(); ++i){
		tmpmap[seg2V[cmnSegs[i]].front()]++;
		tmpmap[seg2V[cmnSegs[i]].back()]++;
	}
	for (auto& end : tmpmap){
		if (end.second == 1){
			nOpenEnds++;
		}
	}
	return nOpenEnds / 2;
}

// allVs stores all vertices in a plain array {x1,y1,z1,x2,y2,z2,...}
// this function return the position of oneV in allVs
int Arrangement::findMatchVer(vector<float>& allVs, vector<float>& oneV){
	int size = (int)allVs.size() / 3;
	for (int i = 0; i < size; ++i){
		if (allVs[3 * i + 0] == oneV[0] && allVs[3 * i + 1] == oneV[1] && allVs[3 * i + 2] == oneV[2]){
			return i;
		}
	}
	return -1;
}

void Arrangement::dealDegeneratedOneSegWithTwoV(vector<int> & oneSeg){
	if(oneSeg.size()==2){
		vector<float> oneV(3);
		for(int i=0; i<3; ++i){
			oneV[i] = (Vs[oneSeg[0]][i]+Vs[oneSeg[1]][i])/2;
		}
		oneSeg.insert(oneSeg.begin()+1, Vs.size());
		Vs.push_back(oneV);
	}
}

// handle the bug where the newly added intersection points are too close to its nb, so tetgen will crash later
// remove the close nbV(s)
void Arrangement::cleanTwoEndsOfOneSeg(vector<int> & oneSeg, vector<int> & toDelSegVs){
	if(oneSeg.size() == 2) {
		return;
	}
	vector<Eigen::Vector3f> SegVCoords;
	vector<float> Dists;
	SegVCoords.reserve(oneSeg.size());
	Dists.reserve((oneSeg.size()-1));
	float minDist = INFINITY;
	SegVCoords.push_back(Eigen::Vector3f (Vs[oneSeg[0]][0],Vs[oneSeg[0]][1],Vs[oneSeg[0]][2]));
	for(size_t i=1; i<oneSeg.size(); ++i){
		SegVCoords.push_back(Eigen::Vector3f (Vs[oneSeg[i]][0],Vs[oneSeg[i]][1],Vs[oneSeg[i]][2]));
		Dists.push_back((SegVCoords[i]-SegVCoords[i-1]).norm());
	}
	for(size_t i=1; i<Dists.size()-1; ++i){
		minDist = Dists[i] < minDist ? Dists[i] : minDist;
	}
	minDist = minDist*CLEANENDTHRESHOLDRATIO;

	if(Dists[0] < minDist){
		toDelSegVs.push_back(oneSeg[1]);
		oneSeg.erase(oneSeg.begin()+1);
	}
	if(oneSeg.size() == 2) {
		return;
	}
	if(Dists.back() < minDist){
		toDelSegVs.push_back(oneSeg[oneSeg.size()-2]);
		oneSeg.erase(oneSeg.end()-2);
	}
	if(oneSeg.size() == 2) {
		return;
	}
}

// from an end of an edge, trace out a segment, return the seg with this input seg "end" as the first vertices
void Arrangement::traceOneSeg(vector<vector<int> > &nbVs,vector<vector<int> > &nbVsLeftMaterial, vector<int> &ind2ind, vector<bool> &unvisited, vector<int> &oneSeg, int end, int &leftMaterial){
	int preV, curV, nxtV;
	oneSeg.clear();
	preV = end;
	curV = nbVs[end].back();
	nbVs[end].pop_back();
	leftMaterial = nbVsLeftMaterial[end].back();
	nbVsLeftMaterial[end].pop_back();

	oneSeg.push_back(ind2ind[end]);
	oneSeg.push_back(ind2ind[curV]);
	unvisited[end] = false;
	unvisited[curV] = false;
	// start trace until meet another end
	while (nbVs[curV].size() > 1){
		nxtV = nbVs[curV][0] == preV ? nbVs[curV][1] : nbVs[curV][0];
		oneSeg.push_back(ind2ind[nxtV]);
		unvisited[nxtV] = false;
		preV = curV;
		curV = nxtV;
	}
	nbVs[curV].pop_back();
	nbVsLeftMaterial[curV].pop_back();
}

// given segs, trace out the cycs form by those segs, store them in cyc2seg and cyc2seg
// this function must run after findCycleInCells().
// assumption: segs must be able to form closed cycles, no open seg left alone in the input.
// cycsegb2cycid is a hash used to avoid storing duplicated cycs multiple times (easy to appear on fronts)
// it will be cleaned after finding all cycles are found (after processing all fronts)
void Arrangement::traceCycsFromSegs(BigInt &segs_b, vector<int>& resCycs){
	resCycs.clear();

	vector<bool> visited(nSeg, false);
	int start, end, curEnd, nxtSeg, curSeg;
	vector<int> oneCyc;
	BigInt oneCyc_b(nSeg);
	BigInt nxtSeg_b;
	vector<int> canSegs = segs_b.getOnesInd();
	for (int seg : canSegs){
		if (visited[seg]) continue;
		visited[seg] = true;
		if (isSegClosed(seg)){
			resCycs.push_back(seg2cyc[seg]);
		}
		else{
			start = seg2V[seg].front();
			end = seg2V[seg].back();
			curEnd = start;
			oneCyc.clear(); oneCyc_b.clear();
			oneCyc.push_back(seg); oneCyc_b.set(seg);
			curSeg = seg;
			// start tracing one closed cycle
			// assumption here: segments in a cell must form closed cycles, and an intersection point is used twice exactly
			while (curEnd != end){
				nxtSeg_b = segs_b & segGraph_b[V2seggV[curEnd]];
				nxtSeg_b.unset(curSeg);
				nxtSeg = nxtSeg_b.getOnesInd()[0];
				oneCyc.push_back(nxtSeg); oneCyc_b.set(nxtSeg);
				visited[nxtSeg] = true;
				curEnd = seg2V[nxtSeg].front() == curEnd ? seg2V[nxtSeg].back() : seg2V[nxtSeg].front();
				curSeg = nxtSeg;
			}
			int cycid = cycsegb2cycid[oneCyc_b];
			if (cycid == 0){
				resCycs.push_back((int)(cyc2seg.size()));
				cyc2seg.push_back(oneCyc);
				cyc2seg_b.push_back(oneCyc_b);
				cycsegb2cycid[oneCyc_b] = (int)(cyc2seg.size());
			}
			// this cyc exists 
			else{
				resCycs.push_back(cycid - 1);
				cycsegb2cycid[oneCyc_b] = cycid;
			}
		}
	}
}

void Arrangement::loadVol(const char * filename, const vector<float> & transVec, const float scaleRatio){
	//if beta=0, do not use vol, then do not read it
	if(t_randwb==0){
		return;
	}
	if(Utility::isFileReadable(filename)){
		string strFullname(filename);
		int lastdotindex = strFullname.find_last_of("."); 
		string extension = strFullname.substr(lastdotindex+1);
        //transform(extension.begin(), extension.end(), extension.begin(), tolower);
        transform(extension.begin(), extension.end(), extension.begin(), [](int i){ return std::tolower(i); });

		if (extension == "mrc"){
			// if it is a .MRC file
			MRCReader mrcReader(filename);
			mrcReader.getVolume(volInfo);
			float scale = STDBBOXMAXSCALE*POLYMENDERENLARGERATIO;
			float lowerCorner = -scale/2;
			float unit = scale/VOLCELLN; 
			volInfo.setLowerCornerXYZ(lowerCorner, lowerCorner, lowerCorner);
			volInfo.setUnitXYZ(unit, unit, unit);
			inputVol = true;
		} else if (extension == "hdr"){
			// else if it is a mediacal vol .HDR
			string datfile = strFullname.substr(0, lastdotindex) + ".dat"; 

			VolReader volReader(filename, datfile.c_str());
			volReader.getVolume(volInfo);
			volInfo.setLowerCornerXYZ(scaleRatio * (volReader.getLowCorner(0) + transVec[0]), scaleRatio * (volReader.getLowCorner(1) + transVec[1]), scaleRatio * (volReader.getLowCorner(2) + transVec[2]));
			volInfo.setUnitXYZ(scaleRatio * volReader.getUnit(0), scaleRatio * volReader.getUnit(1), scaleRatio * volReader.getUnit(2));
			volInfo.normalize(-1, 0);

			inputVol = true; 
		}
	}
}

void Arrangement::loadBBox(const char * filename){
	ifstream fin ( filename );
	string line;
	vector<string> lines;
	if (fin.is_open()) {
		bbox_org.reserve(6);
		while (getline(fin, line)) {
			bbox_org.push_back(atof(line.c_str()));
		}
		fin.close();
		inputBBox = true;
	}
	
}

// store the core data-structures of arrangement, for debug/verification use
// file is stored in the format that Mathematica can read directly
void Arrangement::writeArrangement(const char* filename){
	ofstream ofs(filename, ofstream::out);
	ofs << "{";

	// write Vs
	ofs << "{"; Utility::writeVector(Vs, ofs); ofs << "}";

	// write Segs
	ofs << ",{"; Utility::writeVector(seg2V, ofs); ofs << "}";

	// write Fs
	ofs << ",{"; Utility::writeVector(F2seg, ofs); ofs << "}";

	// write Cells
	ofs << ",{"; Utility::writeVector(cell2F, ofs); ofs << "}";

	// write cell2Fside
	ofs << ",{"; Utility::writeVector(cell2Fside, ofs); ofs << "}";

	// write F2PlaneID
	ofs << ",{"; Utility::writeVector(F2PlaneID, ofs); ofs << "}";

	// write planeParas
	ofs << ",{"; Utility::writeVector(planeParas, ofs); ofs << "}";

	// write cellGraph
	ofs << ",{"; Utility::writeVector(cellGraph, ofs); ofs << "}";

	// write cellEmpty
	ofs << ",{"; Utility::writeVector(cellEmpty, ofs); ofs << "}";

	// write cell frame
	ofs << ",{"; Utility::writeVector(cellVs, ofs); ofs << "}";
	ofs << ",{"; Utility::writeVector(cellEs, ofs); ofs << "}";
	ofs << ",{"; Utility::writeVector(cellFs, ofs); ofs << "}";

	// write closed cycles of each cell
	ofs << ",{"; Utility::writeVector(cyc2seg, ofs); ofs << "}";
	ofs << ",{"; Utility::writeVector(cell2cyc, ofs); ofs << "}";

	// write closed cycles of each front
	ofs << ",{"; Utility::writeVector(front2cyc, ofs); ofs << "}";
	//ofs << ",{"; Utility::writeVector(Order, ofs); ofs << "}";

	// write pre-calculated L values
	ofs << ",{"; Utility::writeVector(front_L_cellcyc2nf0cos, ofs); ofs << "}";

	ofs << "}";

	ofs.close();
}



void Arrangement::connectParallelArrangement(vector<int>_framef2Cs,
                                vector<float>_frameVstack,
                                vector<vector<float>> _frameVs,
                                vector<vector<int>> _frameFs,
                                vector<vector<int>> _frameCells,
                                vector<vector<float>> _frameF2PlaneParas,
                                int _nCell,
                                int _nCrossSection){

    _framef2Cs = this->_framef2Cs;
    _frameVstack = this->_frameVstack;
    _frameVs = this->_frameVs;
    _frameFs = this->_frameFs;
    _frameCells = this->_frameCells;
    _frameF2PlaneParas = this->_frameF2PlaneParas;
    _nCell = this->_nCell;
    _nCrossSection = this->_nCrossSection;

}
