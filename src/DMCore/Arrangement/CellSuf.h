#ifndef _CELLSUF_H_
#define _CELLSUF_H_

#include "../Utility/Utility.h"
#include <stdint.h>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <fstream>

using namespace std;

class OneCellSuf{
public:
	vector<vector<float> > sufV;
	vector<vector<int> > sufF;
	vector<vector<int> > sufSeg;
	vector<vector<int> > sufSeg2F;

	OneCellSuf(){} 
	~OneCellSuf(){}
	void writeSufInfo(const char* filename){
		ofstream ofs(filename, ofstream::out);
		writeSufInfo(ofs);
		ofs.close();
	}
	void writeSufInfo(ofstream &ofs){
		ofs << "{";

		// sufV
		ofs << "{";
		Utility::writeVector(sufV, ofs);
		ofs << "},";

		// sufF
		ofs << "{";
		Utility::writeVector(sufF, ofs);
		ofs << "},";

		// sufSeg
		ofs << "{";
		Utility::writeVector(sufSeg, ofs);
		ofs << "},";

		// sufSeg2F
		ofs << "{";
		Utility::writeVector(sufSeg2F, ofs);
		ofs << "}";

		ofs << "}";
	}
	void saveObj(const char * filename){
		Utility::saveTilingObj(filename, sufV, sufF);
	}
	void saveSeg(const char * filename){
		Utility::writeVectorForm(sufSeg, filename);
	}
	void saveSuf(const char * filename){
		vector<vector<int> > tmpSegs;
		tmpSegs.reserve(sufSeg.size()*20);
		vector<int> oneSeg(2);
		for(int i=0; i<sufSeg.size(); ++i){
			for(int j=0; j<sufSeg[i].size()-1; ++j){
				oneSeg[0] = sufSeg[i][j];
				oneSeg[1] = sufSeg[i][j+1];
				tmpSegs.push_back(oneSeg);
			}
		}
		Utility::saveLuSuf(filename,sufV, sufF, tmpSegs);
	}
};


class AllCellSuf{
public:
	int nCell;
	bool initCell;
	OneCellSuf * cells;
	AllCellSuf(){
		initCell = false;
	}
	void init(int _nCell){
		nCell = _nCell;
		cells = new OneCellSuf [nCell];
		initCell = true;
	}
	~AllCellSuf(){
		if (initCell){
			delete [] cells;
		}
	}

	void writeAllSufInfo(const char* filename){
		ofstream ofs(filename, ofstream::out);
		ofs << "{";

		for(int i=0; i<nCell; ++i){
			if(i!=0){
				ofs << ",";
			}
			cells[i].writeSufInfo(ofs);
		}

		ofs << "}";
		ofs.close();
	}
};

#endif //_CELLSUF_H_
