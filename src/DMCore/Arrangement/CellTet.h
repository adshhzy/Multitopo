#ifndef _CELLTET_H_
#define _CELLTET_H_


#include "../Utility/Utility.h"
#include <stdint.h>
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

class CellTet{
public:
	int v [4];
	int face [4];

	int nbTet [4];

	CellTet(){}
	void write(ofstream& ofs){
		ofs << "{" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << "},";
		ofs << "{" << face[0] << "," << face[1] << "," << face[2] << "," << face[3] << "},";
		ofs << "{" << nbTet[0] << "," << nbTet[1] << "," << nbTet[2] << "," << nbTet[3] << "}";
	}
};

class CellFace{
public:
	int v [3];
	int edge [3];

	int FtetMark; // interior, facei
	int FinoutMark; // in, out, unknown

	vector<int> nbTet;

	CellFace():FtetMark(0),FinoutMark(-1){}
	void write(ofstream& ofs){
		ofs << "{" << v[0] << "," << v[1] << "," << v[2] << "},";
		ofs << "{" << edge[0] << "," << edge[1] << "," << edge[2] << "},";
		ofs << "{" << FtetMark << "," << FinoutMark << "},";
		ofs << "{"; Utility::writeVector(nbTet, ofs); ofs << "}";
	}
};

class CellEdge{
public:
	int v [2];

	int EtetMark; // interior, onsegi, onbd, inface
	int EinoutMark; // in, out, unknown, onseg

	vector<int> nbFace;

	BigInt loc_SubGroupIDs;

	CellEdge():EtetMark(0),EinoutMark(-1){}
	void write(ofstream& ofs){
		ofs << "{" << v[0] << "," << v[1] << "},";
		ofs << "{" << EtetMark << "," << EinoutMark << "},";
		ofs << "{"; Utility::writeVector(nbFace, ofs); ofs << "},";
        //ofs << "{"; Utility::writeVector(loc_SubGroupIDs.getOnesInd(), ofs); ofs << "}";
	}

	inline int getOtherV(int v1){
		return v[0] + v[1] - v1;
	}
};
class CellVertex{
public:

	float xyz[3];
	float pro;
	float proVol;

	int VtetMark; // interior, segi, bd, facei
	int VinoutMark; // in, out, unknown

	vector<int> nbV;
	vector<int> nbE;
	vector<int> nbTet;
	CellVertex():pro(0),proVol(0),VtetMark(0),VinoutMark(-1){}
	void write(ofstream& ofs){
		ofs << "{" << xyz[0] << "," << xyz[1] << "," << xyz[2] << "},";
		ofs << "{" << pro << "," << proVol << "},";
		ofs << "{" << VtetMark << "," << VinoutMark << "},";
		ofs << "{"; Utility::writeVector(nbV, ofs); ofs << "},";
		ofs << "{"; Utility::writeVector(nbE, ofs); ofs << "},";
		ofs << "{"; Utility::writeVector(nbTet, ofs); ofs << "}";
	}
};

class OneCellTet{
public:
	vector<int> cellBBFaces;
	vector<int> cellBBEdges;

	double * VCoord;
	vector<float> VPro;
	int * VMark;

	vector<CellTet> cellTets;
	vector<CellFace> cellFaces;
	vector<CellEdge> cellEdges;
	vector<CellVertex> cellVertices;
	bool initTet;
	bool initFace;
	bool initEdge;
	bool initV;

	vector<vector<int> > cellSeg2TetV;

	OneCellTet(){
		initTet = initFace = initEdge = initV = false;
	}
	void init_tet(int _nTet){
		cellTets.resize(_nTet);
		initTet = true;
	}
	void init_face(int _nFace){
		cellFaces.resize(_nFace);
		initFace = true;
	}
	void init_edge(int _nEdge){
		cellEdges.resize(_nEdge);
		initEdge = true;
	}
	void init_v(int _nV, double * _VCoord){
		cellVertices.resize(_nV);
		VCoord = new double [_nV*3];
		VMark = new int [_nV];
		memcpy(VCoord, _VCoord, _nV*3*sizeof(double));
		initV = true;
	}
	void init_Vmark(int _nV, int * _VMark){
		memcpy(VMark, _VMark, _nV*sizeof(int));
	}
	~OneCellTet(){
		if (initV){
			delete [] VCoord;
			delete [] VMark;
		}
	}
	void writeTetInfo(const char* filename){
		ofstream ofs(filename, ofstream::out);
		ofs << "{";

		// CellVertex
		ofs << "{";
		for (int i=0; i<cellVertices.size(); ++i){
			if(i==0) ofs << "{";
			else ofs << ",{";
			cellVertices[i].write(ofs);
			ofs << "}";
		}
		ofs << "},";

		// CellEdge
		ofs << "{";
		for (int i=0; i<cellEdges.size(); ++i){
			if(i==0) ofs << "{";
			else ofs << ",{";
			cellEdges[i].write(ofs);
			ofs << "}";
		}
		ofs << "},";

		// CellFace
		ofs << "{";
		for (int i=0; i<cellFaces.size(); ++i){
			if(i==0) ofs << "{";
			else ofs << ",{";
			cellFaces[i].write(ofs);
			ofs << "}";
		}
		ofs << "},";

		// CellTet
		ofs << "{";
		for (int i=0; i<cellTets.size(); ++i){
			if(i==0) ofs << "{";
			else ofs << ",{";
			cellTets[i].write(ofs);
			ofs << "}";
		}
		ofs << "},";

		// cellBBFaces
		ofs << "{";
		Utility::writeVector(cellBBFaces, ofs);
		ofs << "},";

		// cellBBEdges
		ofs << "{";
		Utility::writeVector(cellBBEdges, ofs);
		ofs << "},";

		// VPro
		ofs << "{";
		Utility::writeVector(VPro, ofs);
		ofs << "}";

		ofs << "}";
		ofs.close();
	}
};

class AllCellTet{
public:
	int nCell;
	bool initCell;
	OneCellTet * cells;
	AllCellTet(){
		initCell = false;
	}
	void init(int _nCell){
		nCell = _nCell;
		cells = new OneCellTet [nCell];
		initCell = true;
	}
	~AllCellTet(){
		if (initCell){
			delete [] cells;
		}
	}
};

#endif //_CELLTET_H_
