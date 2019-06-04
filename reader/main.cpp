//
//  main.cpp
//  reader
//
//  Created by Zhiyang Huang on 11/13/17.
//  Copyright (c) 2017 Zhiyang Huang. All rights reserved.
//

#include "reader.h"
#include <assert.h>

//#define WRITEMATLAB

#ifdef WRITEMATLAB
	#include "/Applications/MATLAB_R2014a.app/extern/include/mat.h"
#endif

class Ctr{
	
public:
	vector<double>vertices;
	vector<unsigned int>edges2vertices;
	
	Ctr(){}
	
	bool readCtr(string filename){
		
		return readObjFile_Line(filename,vertices,edges2vertices);	
	}
	
};



class Suf{

public:
	vector<double>vertices;
	vector<unsigned int>faces2vertices;
	vector<int>facesMat;
	vector<int>CtrEdges;


	Suf(){}
	
	bool readSuf(string filename){
	
		return readSufFile(filename,vertices,faces2vertices,facesMat,CtrEdges);
	
	}
#ifdef WRITEMATLAB
	void insertSufStructMat(mxArray *target,int ifield){
		int dims[2] = {1,1};
		const char *field_names[] = {"vertices","faces2vertices","facesMat","CtrEdges"};
		mxArray *sarr  = mxCreateStructArray_700(1, dims, 4, field_names);
		/***********************************/
		mxArray *pa1 = mxCreateDoubleMatrix(3,vertices.size()/3,mxREAL);
		memcpy((void *)(mxGetPr(pa1)), (void *)vertices.data(), vertices.size()*sizeof(double));
		
		int name_field = mxGetFieldNumber(sarr,"vertices");
		mxSetFieldByNumber_700(sarr, 0, name_field, pa1);
		//mxDestroyArray(pa1);
		/***********************************/
		vector<double>doubleFaces;
		for(auto a:faces2vertices)doubleFaces.push_back(a+1);
		pa1 = mxCreateDoubleMatrix(3,doubleFaces.size()/3,mxREAL);
		memcpy((void *)(mxGetPr(pa1)), (void *)doubleFaces.data(), doubleFaces.size()*sizeof(double));
		
		name_field = mxGetFieldNumber(sarr,"faces2vertices");
		mxSetFieldByNumber_700(sarr, 0, name_field, pa1);
		//mxDestroyArray(pa1);
		/***********************************/
		vector<double>doublefacesMat;
		for(auto a:facesMat)doublefacesMat.push_back(a+1);
		pa1 = mxCreateDoubleMatrix(2,doublefacesMat.size()/2,mxREAL);
		memcpy((void *)(mxGetPr(pa1)), (void *)doublefacesMat.data(), doublefacesMat.size()*sizeof(double));
		
		name_field = mxGetFieldNumber(sarr,"facesMat");
		mxSetFieldByNumber_700(sarr, 0, name_field, pa1);
		//mxDestroyArray(pa1);
		/***********************************/
		
		vector<double>doubleCtrEdges;
		for(auto a:CtrEdges)doubleCtrEdges.push_back(a+1);
		pa1 = mxCreateDoubleMatrix(2,doubleCtrEdges.size()/2,mxREAL);
		memcpy((void *)(mxGetPr(pa1)), (void *)doubleCtrEdges.data(), doubleCtrEdges.size()*sizeof(double));
		
		name_field = mxGetFieldNumber(sarr,"CtrEdges");
		mxSetFieldByNumber_700(sarr, 0, name_field, pa1);
		//mxDestroyArray(pa1);
		/***********************************/
		
		mxSetFieldByNumber_700(target, 0, ifield, sarr);
		

	}
#endif
};


class CellTopology {
	
public:
	
	//check the folder "topoInfo" for the definition of the topology

	double _score;
	
	//within each connected surface, which loops are contained.
	vector<vector<int>> _comps; //int->loop index
	
	CellTopology(){}
	
	
	
};


int main(int argc, const char * argv[]) {
	
	
	//path to the meta folder	
	string prepath("../data/mug/meta/");
	
#ifdef WRITEMATLAB
	string file = prepath + string("matlabfile")+string(".mat");
#endif


	bool issuf = false;
	
	//structure for the contour network
	Ctr ctr;
	
	
	//Represent the non-manifold contour network into an abstract curve network which turns the sequentical edges into hyper segments. An isolated loop will form a hyper segment, a set of sequential edges between two non-manifold points will also make up a hyper segment
	//the following structure stores the edges that makes up each hyper segments.
	vector<vector<int>>seg2edges;
	

	//within each cell, each loop are made up by 1 or several hyper segments.
	//to get the vertices for each loop, u can first get all edges involved by mapping the seg2edges, and then find the vertices involved in those edges.
	vector<vector<vector<int>>>cellloop2segs;
	vector<vector<vector<int>>>cellloop2edges;
	
	//within each cell the label of each loop
	vector<vector<int>>loop2label;
	
	//within each cell, the score of each topology.
	vector<vector<double>>topo2score;
	
	//within each cell, a structure that stores the infomation of each topology
	vector<vector<CellTopology>>cell2topology;
	
	//within each cell, how many loops
	vector<int>cell2nloop;
	
	//within each cell, how many topologies
	vector<int>cell2ntopo;
	
	
	vector<vector<Suf>>cell2topoSuf;
	
	vector<int>pickTopoInd;
	vector<double>scoreAndtime;
	vector<int>comp2genus;
	vector<int>comp2labels;
	vector<int>label2genus;
	vector<double>stage1time;
	
	bool isDP = true;

	ctr.readCtr(prepath+"ctr.obj");
	
	readVVecFile(prepath+string("seg2edges_vvec.txt"), seg2edges);
	
	readVVecFile(prepath+string("loop2label_vvec.txt"), loop2label);
	
	readVVecFile(prepath+string("topoInfoVVec/info_score_vvec.txt"), topo2score);
	
	if(isDP){
		readVecFile(prepath+string("pickTopo_vec.txt"), pickTopoInd);
	
		readVecFile(prepath+string("scoreAndtime_vec.txt"), scoreAndtime);
	
		readVecFile(prepath+string("isoComp2genus_vec.txt"), comp2genus);
	
		readVecFile(prepath+string("isoComps2label_vec.txt"), comp2labels);
	
		readVecFile(prepath+string("protocol_vec.txt"),label2genus);
		
		readVecFile(prepath+string("timeStage1_vec.txt"),stage1time);
	
		for(auto &a:label2genus)a-=1;
		//label2genus.resize(comp2genus.size());
		//for(int i=0;i<label2genus.size();++i)label2genus[comp2labels[i]] = comp2genus[i]-1;
	}
	
	int n_cell = loop2label.size();
	assert(n_cell == topo2score.size());
	
	cell2nloop.resize(n_cell);
	for(int i = 0;i<n_cell;++i)cell2nloop[i] = loop2label[i].size();
	
	cell2ntopo.resize(n_cell);
	for(int i = 0;i<n_cell;++i)cell2ntopo[i] = topo2score[i].size();
	
	
	
	cellloop2segs.resize(n_cell);	
	for(int i = 0;i<n_cell;++i){
		string filename = prepath+"loop2Segs/cell"+to_string(i)+"_vvec.txt";
		readVVecFile(filename, cellloop2segs[i]);
		assert(cellloop2segs[i].size()==cell2nloop[i]);
	}
	
	cellloop2edges.resize(n_cell);
	for(int i = 0;i<n_cell;++i){
		cellloop2edges[i].resize(cellloop2segs[i].size());
		for(int j=0;j<cellloop2edges[i].size();++j){
			int countedges = 0;
			vector<bool>pickededges(ctr.edges2vertices.size()/2,false);
			for(int k=0;k<cellloop2segs[i][j].size();++k){
				for(auto b:seg2edges[cellloop2segs[i][j][k]])pickededges[b] = true;
			}
			cellloop2edges[i][j].clear();
			for(int b = 0;b<pickededges.size();++b)if(pickededges[b]){cellloop2edges[i][j].push_back(b); countedges++;}

		}

	}
	
	vector<int>testseg2edge;
	for(auto a:seg2edges[0]){
		testseg2edge.push_back(ctr.edges2vertices[a*2]);
		testseg2edge.push_back(ctr.edges2vertices[a*2+1]);
	}
	
	writeObjfile_line(prepath+string("test"),ctr.vertices,testseg2edge);
	
		
	cell2topology.resize(n_cell);
	for(int i = 0;i<n_cell;++i){
		cell2topology[i].resize(cell2ntopo[i]);
		string prefilename = prepath+"topoInfoVVec/info_"+to_string(i)+"_";
		for(int j=0;j<cell2ntopo[i];++j){
			string filename = prefilename+to_string(j)+"_vvec.txt";
			readVVecFile(filename, cell2topology[i][j]._comps);
			cell2topology[i][j]._score = topo2score[i][j];
		
		}
		
	}
	
	if(issuf){
		cell2topoSuf.resize(n_cell);
		for(int i = 0;i<n_cell;++i){
			cell2topoSuf[i].resize(cell2ntopo[i]);
			//string prefilename = prepath+"topoSuf/suf_"+to_string(i)+"_";
			string prefilename = prepath+"../suf/suf_"+to_string(i)+"_";
			for(int j=0;j<cell2ntopo[i];++j){
				string filename = prefilename+to_string(j)+".suf";
				cell2topoSuf[i][j].readSuf(filename);
			
			}
		
		}
	}

	
	
	cout<<"read finish"<<endl;

	// ***********************************************//

		
	//   do your calculation here  //
	
	// ***********************************************//
	
#ifdef WRITEMATLAB
	
	MATFile *pmat;

	printf("Creating file %s...\n\n", file.data());
	pmat = matOpen(file.data(), "w");
	if (pmat == NULL) {
		printf("Error creating file %s\n", file.data());
		printf("(Do you have write permission in this directory?)\n");
		return(EXIT_FAILURE);
	}
	
	auto writeMatrixInfoMatfile = [](string varName, MATFile *pmatfile, vector<double>&datas, int dim1, int dim2){
	
		
		mxArray *pa1;
		pa1 = mxCreateDoubleMatrix(dim1,dim2,mxREAL);
		if (pa1 == NULL) {
			printf("%s : Out of memory on line %d\n", __FILE__, __LINE__); 
			printf("Unable to create mxArray.\n");
			return(EXIT_FAILURE);
		}
		
		memcpy((void *)(mxGetPr(pa1)), (void *)datas.data(), datas.size()*sizeof(double));
		
		int status = matPutVariable(pmatfile, varName.data(), pa1);
		if (status != 0) {
			printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
			return(EXIT_FAILURE);
		}
  
		
		mxDestroyArray(pa1);
		return(EXIT_SUCCESS);

	
	
	};
	
	auto writeMatrixInfoMatfileUint = [writeMatrixInfoMatfile](string varName, MATFile *pmatfile, vector<unsigned int>&datas, int dim1, int dim2){
		
		vector<double>doubledatas(datas.size());
		for(int i=0;i<datas.size();++i)doubledatas[i] = datas[i]+1;
		return writeMatrixInfoMatfile(varName,pmatfile,doubledatas,dim1,dim2);		
	};
	
	auto writeMatrixInfoMatfileInt = [writeMatrixInfoMatfile](string varName, MATFile *pmatfile, vector<int>&datas, int dim1, int dim2){
		
		vector<double>doubledatas(datas.size());
		for(int i=0;i<datas.size();++i)doubledatas[i] = datas[i]+1;
		return writeMatrixInfoMatfile(varName,pmatfile,doubledatas,dim1,dim2);				
		
		
	};
	
	auto writeVVVecInfoMatfileInt = [writeMatrixInfoMatfile](string varName, MATFile *pmatfile, vector<vector<vector<int>>>&datas){
		
		vector<int>dims(2,1);
		dims[0] = datas.size();
		mxArray* globalCell = mxCreateCellArray_700(2, (const int*)(dims.data()));
		for(int icell = 0;icell<datas.size();++icell){
			dims[0] = datas[icell].size();
			mxArray* localCell = mxCreateCellArray_700(2, (const int*)(dims.data()));
			for(int iloop = 0;iloop<datas[icell].size();++iloop){
				vector<double>tmpconvertD;
				for(auto a:datas[icell][iloop])tmpconvertD.push_back(a+1);
				mxArray *pa1;
				pa1 = mxCreateDoubleMatrix(datas[icell][iloop].size(),1,mxREAL);
				if (pa1 == NULL) {
					printf("%s : Out of memory on line %d\n", __FILE__, __LINE__); 
					printf("Unable to create mxArray.\n");
					return(EXIT_FAILURE);
				}
				
				memcpy((void *)(mxGetPr(pa1)), (void *)tmpconvertD.data(), tmpconvertD.size()*sizeof(double));
				mxSetFieldByNumber_700(localCell, 0, iloop, pa1);
				
			}
			

			mxSetFieldByNumber_700(globalCell, 0, icell, localCell);
		
		
		}
		int status = matPutVariable(pmatfile,varName.data(), globalCell);
		
		return(EXIT_SUCCESS);
		
	};
	
	auto writeLoopConnectionsMatfileInt = [writeMatrixInfoMatfile](string varName, MATFile *pmatfile, vector<vector<CellTopology>>&datas){
		
		vector<int>dims(2,1);
		dims[0] = datas.size();
		mxArray* globalCell = mxCreateCellArray_700(2, (const int*)(dims.data()));
		for(int icell = 0;icell<datas.size();++icell){
			dims[0] = datas[icell].size();
			mxArray* localCell = mxCreateCellArray_700(2, (const int*)(dims.data()));
			for(int iTopo = 0;iTopo<datas[icell].size();++iTopo){
				dims[0] = datas[icell][iTopo]._comps.size();
				mxArray* localCell_2 = mxCreateCellArray_700(2, (const int*)(dims.data()));
				for(int iComp = 0;iComp<dims[0] ;++iComp){
					vector<double>tmpconvertD;
					for(auto a:datas[icell][iTopo]._comps[iComp])tmpconvertD.push_back(a+1);
					mxArray *pa1;
					pa1 = mxCreateDoubleMatrix(datas[icell][iTopo]._comps[iComp].size(),1,mxREAL);
					if (pa1 == NULL) {
						printf("%s : Out of memory on line %d\n", __FILE__, __LINE__); 
						printf("Unable to create mxArray.\n");
						return(EXIT_FAILURE);
					}
					memcpy((void *)(mxGetPr(pa1)), (void *)tmpconvertD.data(), tmpconvertD.size()*sizeof(double));
					mxSetFieldByNumber_700(localCell_2, 0, iComp, pa1);
				}
				
				
				mxSetFieldByNumber_700(localCell, 0, iTopo, localCell_2);
				
			}
			
			
			mxSetFieldByNumber_700(globalCell, 0, icell, localCell);
			
			
		}
		int status = matPutVariable(pmatfile,varName.data(), globalCell);
		
		return(EXIT_SUCCESS);
		
	};
	
	auto writeVVecInfoMatfileInt = [writeMatrixInfoMatfile](string varName, MATFile *pmatfile, vector<vector<int>>&datas){
		
		vector<int>dims(2,1);
		dims[0] = datas.size();
		mxArray* globalCell = mxCreateCellArray_700(2, (const int*)(dims.data()));
		for(int icell = 0;icell<datas.size();++icell){
			dims[0] = datas[icell].size();
			vector<double>tmpconvertD;
			for(auto a:datas[icell])tmpconvertD.push_back(a+1);
			mxArray *pa1;
			pa1 = mxCreateDoubleMatrix(datas[icell].size(),1,mxREAL);
			if (pa1 == NULL) {
					printf("%s : Out of memory on line %d\n", __FILE__, __LINE__); 
					printf("Unable to create mxArray.\n");
					return(EXIT_FAILURE);
			}
				
			memcpy((void *)(mxGetPr(pa1)), (void *)tmpconvertD.data(), tmpconvertD.size()*sizeof(double));							
			mxSetFieldByNumber_700(globalCell, 0, icell, pa1);
					
		}
		int status = matPutVariable(pmatfile,varName.data(), globalCell);
		
		return(EXIT_SUCCESS);
		
	};
	auto writeVVecInfoMatfile = [writeMatrixInfoMatfile](string varName, MATFile *pmatfile, vector<vector<double>>&datas){
		
		vector<int>dims(2,1);
		dims[0] = datas.size();
		mxArray* globalCell = mxCreateCellArray_700(2, (const int*)(dims.data()));
		for(int icell = 0;icell<datas.size();++icell){
			dims[0] = datas[icell].size();
			vector<double>tmpconvertD;
			for(auto a:datas[icell])tmpconvertD.push_back(a);
			mxArray *pa1;
			pa1 = mxCreateDoubleMatrix(datas[icell].size(),1,mxREAL);
			if (pa1 == NULL) {
				printf("%s : Out of memory on line %d\n", __FILE__, __LINE__); 
				printf("Unable to create mxArray.\n");
				return(EXIT_FAILURE);
			}
			
			memcpy((void *)(mxGetPr(pa1)), (void *)tmpconvertD.data(), tmpconvertD.size()*sizeof(double));							
			mxSetFieldByNumber_700(globalCell, 0, icell, pa1);
			
		}
		int status = matPutVariable(pmatfile,varName.data(), globalCell);
		
		return(EXIT_SUCCESS);
		
	};

	
	auto writeVVecInfoMatfileSuf = [writeMatrixInfoMatfile](string varName, MATFile *pmatfile, vector<vector<Suf>>&datas){
		
		vector<int>dims(2,1);
		dims[0] = datas.size();
		mxArray* globalCell = mxCreateCellArray_700(2, (const int*)(dims.data()));
		for(int icell = 0;icell<datas.size();++icell){
			dims[0] = datas[icell].size();
			mxArray* localCell = mxCreateCellArray_700(2, (const int*)(dims.data()));
			for(int iTopo = 0;iTopo<datas[icell].size();++iTopo){
						
				datas[icell][iTopo].insertSufStructMat(localCell,iTopo);
				
			}
			
			
			mxSetFieldByNumber_700(globalCell, 0, icell, localCell);
			
			
		}
		int status = matPutVariable(pmatfile,varName.data(), globalCell);
		
		return(EXIT_SUCCESS);
		
	};

	mxArray *mxCreateStructArray(mwSize ndim, const mwSize *dims,
								 int nfields, const char **fieldnames);
	writeMatrixInfoMatfile(string("ctrVertices"),pmat,ctr.vertices,3,ctr.vertices.size()/3);
	writeMatrixInfoMatfileUint(string("ctrEdges"),pmat,ctr.edges2vertices,2,ctr.edges2vertices.size()/2);
	
	writeMatrixInfoMatfileInt(string("OptimalCandidates"),pmat,pickTopoInd,1,pickTopoInd.size());
	writeMatrixInfoMatfile(string("OptimalScoreAndTime"),pmat,scoreAndtime,1,scoreAndtime.size());
	writeMatrixInfoMatfileInt(string("TargetGenusPerlabel"),pmat,label2genus,1,label2genus.size());
	writeMatrixInfoMatfile(string("TimeStage1"),pmat,stage1time,1,stage1time.size());
	
		
	writeVVVecInfoMatfileInt(string("cellloop2edges"),pmat,cellloop2edges);
	
	writeVVecInfoMatfileInt(string("cellloop2label"),pmat,loop2label);
	
	writeLoopConnectionsMatfileInt(string("loopconnection"),pmat,cell2topology);
	
	
	writeVVecInfoMatfile(string("celltopo2score"),pmat,topo2score);
	
	
	writeVVecInfoMatfileSuf(string("cell2topoSuf"),pmat,cell2topoSuf);
	
	
	if (matClose(pmat) != 0) {
		printf("Error closing file %s\n",file.data());
		return(EXIT_FAILURE);
	}
	
	cout<<"Write mat finish!"<<endl;
	
#endif
	return 0;
}
