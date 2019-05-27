#ifndef _UTILITY_H_
#define _UTILITY_H_

#ifndef _DO_PARALLEL_
#define _DO_PARALLEL_ 1
#endif

#ifndef _SAVERES_MING
#define _SAVERES_MING 0
#endif

#ifndef _CG_LOG 
#define _CG_LOG 1
#endif

#ifndef _SAVEDATA_MING
#define _SAVEDATA_MING 0
#endif

#ifndef _DEBUG_MING_MORE
#define _DEBUG_MING_MORE 0
#endif

#ifndef _DEBUG_MING_TETGEN
#define _DEBUG_MING_TETGEN 0
#endif

#ifndef VAXMAN
#define VAXMAN 0
#endif

#ifndef USETRIANGULATION
#define USETRIANGULATION 0
#endif

#ifndef VOLCELLN
#define VOLCELLN 256
#endif

#ifndef STDBBOXMINSCALE 
#define STDBBOXMAXSCALE 80
#endif

#ifndef POLYMENDERENLARGERATIO
#define POLYMENDERENLARGERATIO 1.25
#endif

#ifndef MIXISOFAKEVALUE
#define MIXISOFAKEVALUE -0.5f
#endif

#ifndef CLEANENDTHRESHOLD
#define CLEANENDTHRESHOLD 0.0005
#endif


#ifndef CLEANENDTHRESHOLDRATIO
#define CLEANENDTHRESHOLDRATIO 0.3
#endif

#ifndef INFINITY
#define INFINITY 99999999999999
#endif

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "BigInt.h"
#include "Graph.h"
#include "tetgen.h"

extern int db_celli;

using namespace std;

class myVecHash{
public:
	inline size_t operator() (const vector<int> &x) const {

		int nHash = 0;
		for (size_t i = 0; i < x.size(); i++) {
			nHash += x[i];
			nHash += (nHash << 10);
			nHash ^= (nHash >> 6);
		}
		nHash += (nHash << 3);
		nHash ^= (nHash >> 11);
		nHash += (nHash << 15);
		return (size_t)nHash;
	}
};

class Utility{
public:
	enum INOUTSIDE {OUTSIDE = 0, INSIDE = 1, UNKNOWN = -2};
	enum INPARA {EXE = 0, INCONTOUR, INVOL, INBBOX, OUTDIR, TETLIMIT, RANDWBETA, NBEFLOOP, NLOOP, NINLOOP, FIRSTGENUS};

	template <class T>
	static void writeVector(vector<vector<T>>& vec, ofstream& ofs);

	template <class T>
	static void writeVector(vector<T>& vec, ofstream& ofs);

	template <class T>
	static void writeVector(T& vec, const char* filename);

	template <class T>
	static void writeVectorForm(vector<vector<T>>& vec, ofstream& ofs);

	template <class T>
	static void writeVectorForm(vector<T>& vec, ofstream& ofs);

	template <class T>
	static void writeVectorForm(T& vec, const char* filename);

	template <class T>
	static void printVector(vector<vector<T> >& in);

	template <class T>
	static void printVector(vector<T>& in);

	template <class T>
	static string toStr(T& in);

	template <class T>
	static void writeAnArray(const char* filename, const T * arr, const int num);

	template <class T>
	static void flatten2DVectors(const vector<vector<T>>& vec, vector<T> & res);

	template <class T>
	static void flatten2DVectors(const vector<vector<T>>& vec, T *res);


	template <class T>
	static inline int findNEleSmallerInSortedList(const vector<T> &sortedList, const T val){
		int index;
		for(index=0; (index<sortedList.size())&&(sortedList[index]<val); ++index);
		return index;
	}

	template<class T>
	static inline void saveTriangulationCurveFile(const char* filename, const vector<vector<vector<T> > > & Vs){
		int nV = 0;
		for(size_t i=0; i<Vs.size(); ++i){
			nV += Vs[i].size();
		}
		ofstream ofs(filename, ofstream::out);
		ofs << Vs.size() <<"\n";
		ofs << nV <<"\n";
		for(size_t i=0; i<Vs.size(); ++i){
			ofs << Vs[i].size() << " 0 1\n" ;
			for(size_t j=0; j<Vs[i].size(); ++j){
				ofs << Vs[i][j][0] << " " << Vs[i][j][1] << " " << Vs[i][j][2] << "\n";
			}
		}

		ofs.close();
	}

    static void flipSufOrientation(vector<vector<int> > & sufF){
		int tmpV = -1;
		for(int fi=0; fi<sufF.size(); ++fi){
			tmpV = sufF[fi][0];
			sufF[fi][0] = sufF[fi][1];
			sufF[fi][1] = tmpV;
		}
	}

    static void forceSufOriented(vector<vector<int> > & sufF){
		// find the edges E2V and edge-face connections E2nbF
		unordered_map<vector<int>, int, myVecHash> eHash;
		vector<vector<int> > E2V;
		vector<vector<int> > E2nbF;
		vector<int> oneE2V(2);
		for(int fi=0; fi<sufF.size(); ++fi){
			for(int i=0; i<3; ++i){
				oneE2V[0] = sufF[fi][i];
				oneE2V[1] = i == 2? sufF[fi][0] : sufF[fi][i+1];
				sort(oneE2V.begin(), oneE2V.end());
				//new edge is found
				if(eHash.find(oneE2V)==eHash.end()){
					eHash[oneE2V] = E2V.size();
					E2V.push_back(oneE2V);
					vector<int> oneE2nbF(1,fi);
					E2nbF.push_back(oneE2nbF);
				}else{
					E2nbF[eHash[oneE2V]].push_back(fi);
				}
			}
		}

		// construct the dual graph, and record the two ends of non-manifold edges
		UnDirGraph dGraph(sufF.size());
		unordered_set<int> nonManifoldVs;
		unordered_set<int> nonManifoldEs;
		for(int ei=0; ei<E2V.size(); ++ei){
			if(E2nbF[ei].size()==2){
				dGraph.addEdge(E2nbF[ei][0],E2nbF[ei][1]);
			}
		}

		BigInt visited(sufF.size());

		vector<int> commonVPos_curF(2);
		vector<int> commonVInd_curF(2);
		int otherVPos_curF;
		int otherVInd_nbF;

		queue<int> fringe;
		int curF = 0;
		fringe.push(curF);
		visited.set(curF);
		while(!fringe.empty()){
			curF = fringe.front();
			fringe.pop();
			for(auto nbF : dGraph.adj[curF]){
				if(!visited.get(nbF)){
					int count = 0;
					otherVPos_curF = 3;
					otherVInd_nbF = sufF[nbF][0] + sufF[nbF][1] +sufF[nbF][2];
					for(int i=0; i<3; ++i){
						for(int j=0; j<3; ++j){
							if(sufF[curF][i] == sufF[nbF][j]){
								commonVPos_curF[count] = i;
								commonVInd_curF[count] = sufF[curF][i];
								otherVPos_curF -= i;
								otherVInd_nbF -= sufF[curF][i];

								count++;

							}
						}
					}
					sufF[nbF][otherVPos_curF] = otherVInd_nbF;
					sufF[nbF][commonVPos_curF[0]] = commonVInd_curF[1];
					sufF[nbF][commonVPos_curF[1]] = commonVInd_curF[0];

					fringe.push(nbF);
					visited.set(nbF);
				}
			}
		}

		if(!visited.isAllOne()){
			cout << "Error! the input sufF is not a single connected component! " << endl;
		}


	}

	template <class T>
	static void writeOutOneVector(const char* filename, T &vec); 

	inline static bool isFileReadable (const char* filename) {
		if (FILE *file = fopen(filename, "r")) {
			fclose(file);
			return true;
		} else {
			return false;
		} 
	}


	static void findOrderedPairs(vector<vector<int> > &in_Pairs, vector<int> &in_PairsInd, int in_startPairID, int in_startEleID, vector<int>& out_PairIDs, vector<int>&out_EleIDs){
		out_PairIDs.push_back(in_PairsInd[in_startPairID]);
		out_EleIDs.push_back(in_startEleID);

		vector<vector<int> > leftPairs = in_Pairs;
		vector<int> leftPairInds = in_PairsInd;
		leftPairs.erase(leftPairs.begin()+in_startPairID);
		leftPairInds.erase(leftPairInds.begin()+in_startPairID);
		int preEleID = in_startEleID;

		while(leftPairs.size()!=0){ 
			int curPairID = 0;
			for(curPairID=0; curPairID<leftPairs.size(); curPairID++){
				if (leftPairs[curPairID][0]==preEleID){
					preEleID = leftPairs[curPairID][1];
					break;
				}
				else if (leftPairs[curPairID][1]==preEleID){
					preEleID = leftPairs[curPairID][0];
					break;
				}

			}
			out_PairIDs.push_back(leftPairInds[curPairID]);
			out_EleIDs.push_back(preEleID);
			leftPairs.erase(leftPairs.begin()+curPairID);
			leftPairInds.erase(leftPairInds.begin()+curPairID);
		}
	}

	static void findMergingOfTwoLists(vector<float> &list1, vector<float> &list2, vector<vector<int> > &merge1, vector<vector<int> > &merge2, float err){
		int n1 = list1.size();
		int n2 = list2.size();
		int i1=0, i2=0;
		merge1.resize(n1+1);
		merge2.resize(n2+1);
		while (i1<n1 && i2<n2){
			if (abs(list1[i1]-list2[i2]) < err){
				i1++; i2++;
			} else if(list1[i1]<list2[i2]){
				merge2[i2].push_back(i1);
				i1++;
			} else if(list1[i1]>list2[i2]){
				merge1[i1].push_back(i2);
				i2++;
			} 
		}
		if (i1==n1){
			for (int i=i2; i<n2; ++i){
				merge1[n1].push_back(i);
			}
		}else if (i2==n2){
			for (int i=i1; i<n1; ++i){
				merge2[n2].push_back(i);
			}
		}
	}

	static inline void readSkipNWord(ifstream &ifs, int n){
		string word;
		for(int i=0; i<n; ++i){
			ifs >> word;
		}
	}
	template <class T>
	static void read2DArrayVector(const char* filename, vector<vector<T> > &data){
		ifstream ifs(filename, ofstream::in);
		if (!ifs.good()) {
			cout << "Can not read file " << filename << endl;
			return;
		}

		int num;
		ifs >> num;
		data.resize(num);

		// read vertices
		for(int i=0; i<data.size(); i++){
			ifs >> num;
			data[i].resize(num);
			for(int j=0; j<num; ++j){
				ifs >> data[i][j];
			}
		}
		ifs.close();
	}

	static void readTilingObj(const char* tilefile, vector<vector<float> > &sufV, vector<vector<int> > &sufF){
		ifstream ifs(tilefile, ofstream::in);
		if (!ifs.good()) {
			cout << "Can not read OBJ file " << tilefile << endl;
			return;
		}
		string line;
		string line1;
		int nV, nF;
		readSkipNWord(ifs, 6);
		// #V #F
		ifs >> line >> line1 >> nV;
		cout << line << endl;
		cout << line1 << endl;
		ifs >> line >> line1 >> nF;
		cout << line << endl;
		cout << line1 << endl;
		sufV.resize(nV, vector<float>(3));
		sufF.resize(nF, vector<int>(3));

		// read vertices
		for(int i=0; i<nV; i++){
			ifs >> line >> sufV[i][0] >>  sufV[i][1] >> sufV[i][2];
		}
		// read faces
		int v1, v2, v3;
		for(int i=0; i<nF; i++){
			ifs >> line >> v1 >>  v2 >> v3;
			sufF[i][0] = v1-1;
			sufF[i][1] = v2-1;
			sufF[i][2] = v3-1;
		}
		ifs.close();
	}

	static void saveTilingObj(const char* tilefile, vector<vector<float> > &sufV, vector<vector<int> > &sufF){
		ofstream writer(tilefile, ofstream::out);
		if (!writer.good()) {
			cout << "Can not write OBJ file " << tilefile << endl;
			return;
		}
		writer << "# OBJ File Generated by CycleGrouping\n";
		writer << "# Vertices: " << sufV.size() << "\n";
		writer << "# Faces: " << sufF.size() << "\n";
		// write vertices
		for(int i=0; i<sufV.size(); i++){
			writer << "v "<< sufV[i][0] << " " << sufV[i][1] << " " << sufV[i][2] <<"\n";
		}
		// write faces
		for (int i=0; i<sufF.size(); i++){
			writer << "f "<< sufF[i][0]+1 << " " << sufF[i][1]+1 << " " << sufF[i][2]+1 <<"\n";
		}
		writer.close();

	}

	static void saveLuSuf(const char* suffile, vector<vector<float> > &sufV, vector<vector<int> > &sufF, vector<vector<int> > &sufSeg){
		ofstream writer(suffile, ofstream::out);
		if (!writer.good()) {
			cout << "Can not write OBJ file " << suffile << endl;
			return;
		}
		writer << sufV.size() << " " << sufF.size() << "\n";
		// write vertices
		for(int i=0; i<sufV.size(); i++){
			writer << sufV[i][0] << " " << sufV[i][1] << " " << sufV[i][2] <<"\n";
		}
		// write faces
		for (int i=0; i<sufF.size(); i++){
			writer << sufF[i][0] << " " << sufF[i][1] << " " << sufF[i][2] << " 0 1" <<"\n";
		}
		writer << sufSeg.size() << "\n";
		// write segs
		for (int i=0; i<sufSeg.size(); i++){
			writer << sufSeg[i][0] << " " << sufSeg[i][1] <<"\n";
		}
		writer.close();
	}

	static void saveTilingObj(const char* tilefile, vector<float> &sufV, vector<int> &sufF){
		ofstream writer(tilefile, ofstream::out);
		if (!writer.good()) {
			cout << "Can not write OBJ file " << tilefile << endl;
			return;
		}
		writer << "# OBJ File Generated by CycleGrouping\n";
		writer << "# Vertices: " << sufV.size()/3 << "\n";
		writer << "# Faces: " << sufF.size()/3 << "\n";
		// write vertices
		for(int i=0; i<sufV.size()/3; i++){
			writer << "v "<< sufV[i*3+0] << " " << sufV[i*3+1] << " " << sufV[i*3+2] <<"\n";
		}
		// write faces
		for (int i=0; i<sufF.size()/3; i++){
			writer << "f "<< sufF[i*3+0]+1 << " " << sufF[i*3+1]+1 << " " << sufF[i*3+2]+1 <<"\n";
		}
		writer.close();
	}
	static void saveTilingObj(const char* tilefile, const int nV, const float * sufV, const int nF, const int * sufF){
		ofstream writer(tilefile, ofstream::out);
		if (!writer.good()) {
			cout << "Can not write OBJ file " << tilefile << endl;
			return;
		}
		writer << "# OBJ File Generated by CycleGrouping\n";
		writer << "# Vertices: " << nV << "\n";
		writer << "# Faces: " << nF << "\n";
		// write vertices
		for(int i=0; i<nV; i++){
			writer << "v "<< sufV[i*3+0] << " " << sufV[i*3+1] << " " << sufV[i*3+2] <<"\n";
		}
		// write faces
		for (int i=0; i<nF; i++){
			writer << "f "<< sufF[i*3+0]+1 << " " << sufF[i*3+1]+1 << " " << sufF[i*3+2]+1 <<"\n";
		}
		writer.close();
	}
	static void scaleUpVers(vector<vector<float> > & ctrvers, const float scale, const vector<float> & translate){
		for(int i=0; i<ctrvers.size(); ++i){
			for(int j=0; j<ctrvers[i].size()/3; ++j){
				ctrvers[i][3*j] = (ctrvers[i][3*j] + translate[0])*scale;
				ctrvers[i][3*j+1] = (ctrvers[i][3*j+1] + translate[1])*scale;
				ctrvers[i][3*j+2] = (ctrvers[i][3*j+2] + translate[2])*scale;
			}
		}
	}
	static void scaleDnVers(vector<vector<float> > & ctrvers, const float scale, const vector<float> & translate){
		for(int i=0; i<ctrvers.size(); ++i){
			for(int j=0; j<ctrvers[i].size()/3; ++j){
				ctrvers[i][3*j] = ctrvers[i][3*j]/scale - translate[0];
				ctrvers[i][3*j+1] = ctrvers[i][3*j+1]/scale - translate[1];
				ctrvers[i][3*j+2] = ctrvers[i][3*j+2]/scale - translate[2];
			}
		}
	}
	static void scaleDnVers(vector<float> & ctrvers, const float scale, const vector<float> & translate){
		for(int j=0; j<ctrvers.size()/3; ++j){
			ctrvers[3*j] = ctrvers[3*j]/scale - translate[0];
			ctrvers[3*j+1] = ctrvers[3*j+1]/scale - translate[1];
			ctrvers[3*j+2] = ctrvers[3*j+2]/scale - translate[2];
		}
	}
	static void scaleUpBBox(float  bbox [], const float scale, const vector<float> & translate){
		for(int i=0; i<2; ++i){
			for(int j=0; j<3; ++j){
				bbox[i*3+j] = (bbox[i*3+j] + translate[j])*scale;
			}
		}
	}
	static void scaleUpParas(vector<float > & paras, const float scale, const vector<float> & translate){
		for(int i=0; i<paras.size()/4; ++i){
			paras[i*4+3] = (paras[i*4+3] 
			- paras[i*4]*translate[0]
			- paras[i*4+1]*translate[1]
			- paras[i*4+2]*translate[2])*scale;
		}
	}

	static inline void findTransVecAndScaleRatioOfBBox(float bbox [6], vector<float> & transVec, float & scaleRatio){
		float t_bboxLowXYZ [3];
		float t_bboxXYZScale [3];

		for(int i=0; i<3; ++i){
			t_bboxLowXYZ[i] = bbox[i];
			t_bboxXYZScale[i] = bbox[3+i]-bbox[i];
		}

		transVec.resize(3);

		float maxBBoxScale = -1;
		for(int i=0; i<3; ++i){
			maxBBoxScale = t_bboxXYZScale[i] > maxBBoxScale ? t_bboxXYZScale[i] : maxBBoxScale;
			transVec[i] = - t_bboxLowXYZ[i] - t_bboxXYZScale[i]/2;
		}
		scaleRatio = STDBBOXMAXSCALE/maxBBoxScale;
	}

	static void saveTetgenio(const char * filename, const tetgenio &tetio){
		ofstream ofs(filename, ofstream::out);
		ofs << "{";

		// points
		vector<vector<double> > Vs;
		vector<int> Vmarks;
		vector<double> oneV(3);
		Vs.reserve(tetio.numberofpoints);
		Vmarks.reserve(tetio.numberofpoints);
		for(int i=0; i<tetio.numberofpoints; ++i){
			for(int j=0; j<3; ++j){
				oneV[j] = tetio.pointlist[i*3+j];
			}
			Vs.push_back(oneV);
			Vmarks.push_back(tetio.pointmarkerlist[i]);
		}
		ofs<<"{";
		ofs<<"{"; writeVector(Vs,ofs); ofs<<"},";
		ofs<<"{"; writeVector(Vmarks,ofs); ofs<<"}";
		ofs<<"},";

		// polygons on each facets
		ofs<<"{";
		for(int fi=0; fi<tetio.numberoffacets; ++fi){
			vector<vector<int> > Polys;

			int npoly = tetio.facetlist[fi].numberofpolygons;
			for(int pi=0; pi<npoly; ++pi){
				vector<int> onePoly;
				int nPv = tetio.facetlist[fi].polygonlist[pi].numberofvertices;
				for(int vi=0; vi<nPv; ++vi){
					onePoly.push_back(tetio.facetlist[fi].polygonlist[pi].vertexlist[vi]);
				}
				Polys.push_back(onePoly);
			}
			if(fi>0){
				ofs<<",";
			}
			ofs<<"{"; writeVector(Polys,ofs); ofs<<"}";
		}
		ofs<<"}";



		ofs << "}";
		ofs.close();
	}
};

template <class T>
void Utility::writeVector(vector<vector<T>>& vec, ofstream& ofs){
	for (unsigned int i = 0; i < vec.size(); ++i){
		if (i == 0)	ofs << "{";
		else        ofs << ",{";
		for (unsigned int j = 0; j < vec[i].size(); ++j){
			if (j == 0)	ofs << vec[i][j];
			else        ofs << "," << vec[i][j];
		}
		ofs << "}";
	}
}

template <class T>
void Utility::writeAnArray(const char* filename, const T * arr, const int num){
	ofstream ofs(filename, ofstream::out);
	for(size_t i=0; i<num; ++i){
		ofs << arr[i] << "\n";
	}
	ofs.close();
}

template <class T>
void Utility::writeVector(vector<T>& vec, ofstream& ofs){
	for (unsigned int j = 0; j < vec.size(); ++j){
		if (j == 0)	ofs << vec[j];
		else        ofs << "," << vec[j];
	}
}

template <class T>
void Utility::writeVector(T & vec, const char* filename){
	ofstream ofs(filename, ofstream::out);
	ofs << "{";
	writeVector(vec, ofs);
	ofs << "}";
	ofs.close();
}

template <class T>
void Utility::writeVectorForm(vector<vector<T>>& vec, ofstream& ofs){
	ofs << vec.size() << "\n";
	for (unsigned int i = 0; i < vec.size(); ++i){
		writeVectorForm(vec[i], ofs);
	}
}

template <class T>
void Utility::writeVectorForm(vector<T>& vec, ofstream& ofs){
	for (unsigned int j = 0; j < vec.size(); ++j){
		ofs << vec[j] << " ";
	}
	ofs << "\n";
}

template <class T>
void Utility::writeVectorForm(T & vec, const char* filename){
	ofstream ofs(filename, ofstream::out);
	writeVectorForm(vec, ofs);
	ofs.close();
}

template <class T>
void Utility::printVector(vector<vector<T> >& in){
	for (int i = 0; i < in.size(); ++i){
		cout << "        <" << i << "> : ";
		for (T ele : in[i]){
			cout << ele << " ";
		}
		cout << endl;
	}
}

template <class T>
void Utility::printVector(vector<T>& in){
	for (T ele : in){
		cout << ele << " ";
	}
	cout << endl;
}

template <class T> 
string Utility::toStr(T& in){ 
	ostringstream os;
	os << in;
	return os.str();
}



template <class T>
void Utility::writeOutOneVector(const char* filename, T &vec){
	ofstream ofs(filename, ofstream::out);
	ofs << "{"; Utility::writeVector(vec, ofs); ofs << "}";
	ofs.close();
}

template <class T>
void Utility::flatten2DVectors(const vector<vector<T>>& vec, vector<T> &res) {
	if(vec.size()==0) return;
	res.reserve(vec.size()*vec[0].size());
	for (const auto& sub : vec)
		res.insert(res.end(), sub.begin(), sub.end());
}

template <class T>
void Utility::flatten2DVectors(const vector<vector<T>>& vec, T * res) {
	if(vec.size()==0) return;
	int index = 0;
	for(int i=0; i<vec.size(); ++i){
		for(int j=0; j<vec[i].size(); ++j, index++){
			res[index] = vec[i][j];
		}
	}
}

#endif //_UTILITY_H_
