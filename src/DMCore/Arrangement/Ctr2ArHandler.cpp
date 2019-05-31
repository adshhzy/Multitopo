/*
Ctr2ArHandler.cpp
Ming Zou @ wustl 2014

This handler reads a .contour file and outputs a space arrangement, which is 
the basic data structure used later by the cycle grouping algorithm.
File reading and space partition are done by using the same functions in Lu Liu's
Ctr2Suf project (http://www.cs.wustl.edu/~taoju/lliu/paper/ctr2suf/program.html)

Note: add _CRT_SECURE_NO_WARNINGS to the project setting
(C++ preprocessor definitions) to run Lu's code
*/

#include "Ctr2ArHandler.h"

Ctr2ArHandler::Ctr2ArHandler(){}
Ctr2ArHandler::~Ctr2ArHandler(){}

void Ctr2ArHandler::readContour(const char* filename, const char* volfilename, const char* bboxfilename, Arrangement &ar){
	
	/* ------------------- phase 1: read in .contour -------------------- */
	floatvector param;				//the parameters of the planes
	vector<floatvector> ctrvers;	//vertices of the contours
	vector<intvector> ctredges;		//edges of the contours
	float bbox[6];					//bounding box of the contours
	bool bboxset;					//if bounding box is set or not
	char** filenames = new char*[1];

	filenames[0] = new char[100];
	strcpy(filenames[0], filename);
	_global_IsCtrGraph = false;
	bboxset = false;

	// !! start reading !!
	ContourHandler::readContour(1, filenames, param, ctrvers, ctredges, bbox, bboxset);
#if _CG_LOG
	if(ar.nCell!=0){
		cout << "[+] Contour file loaded... " << endl;
	}else{
		cout << "[-] Contour file is not loaded correctly. " << endl;
	}
#endif
	// ar.bbox_org is always consistent with bbox after following processing
	ar.loadBBox(bboxfilename);
	if(ar.bbox_org.size()==6){
		for(int i=0; i<6; ++i){
			bbox[i] = ar.bbox_org[i];
		}
	}else{
		ar.bbox_org.resize(6);
		for(int i=0; i<6; ++i){
			ar.bbox_org[i] = bbox[i];
		}
	}

	// scale up the model
	Utility::findTransVecAndScaleRatioOfBBox(bbox, ar.unifyTransVec, ar.unifyScaleRatio);

	Utility::scaleUpVers(ctrvers, ar.unifyScaleRatio, ar.unifyTransVec);
	Utility::scaleUpBBox(bbox, ar.unifyScaleRatio, ar.unifyTransVec);
	Utility::scaleUpParas(param, ar.unifyScaleRatio, ar.unifyTransVec);

	/* ------------------- phase 2: preprocess contour -------------------- */
	/*pctr: plane contour*/
	/*read from .contour, get the V and E and Para of each contour on each plane*/
	int planenum;			//number of the planes
	float* pparam;			//parameters of the plane after process
	float** pctrvers;		//contour vertices after process
	int* pctrvernum;		//contour vertices number on each plane
	int** pctredges;		//contour edges after process
	int* pctredgenum;		//number of contour edges on each plane
	int** ver2planelistpos;
	vector<intvector> ver2planelist; 
	planenum = (int)ctrvers.size();
	pctrvers = new float*[planenum];
	pctrvernum = new int[planenum]; 
	pctredges = new int*[planenum];
	pctredgenum = new int[planenum];
	pparam = new float[planenum * 4 + 24];

	for (int i = 0; i < planenum; i++){
		int len = (int)ctrvers[i].size();
		pctrvernum[i] = len / 3;
		//ver
		pctrvers[i] = new float[len];
		for (int j = 0; j < len; j++){
			pctrvers[i][j] = ctrvers[i][j];
		}
		//edge
		len = (int)ctredges[i].size();
		pctredgenum[i] = len / 4;
		pctredges[i] = new int[len];
		for (int j = 0; j < len; j++){
			pctredges[i][j] = ctredges[i][j];
		}
		//param
		for (int j = 0; j < 3; j++){
			pparam[4 * i + j] = param[4 * i + j];
		}
		pparam[4 * i + 3] = MyMath::dotProduct(pparam + 4 * i, pctrvers[i]);

		float paramlen = MyMath::vectorlen(pparam + 4 * i);
		for (int j = 0; j< 4; j++)
			pparam[4 * i + j] = pparam[4 * i + j] / paramlen;
	}

	// !! start pre processing !!
	ContourHandler::preProcDataSingleMat(planenum, pparam, pctrvers, pctrvernum, pctredges, pctredgenum,
		ver2planelistpos, ver2planelist);

    cout<<"------------------- phase 3: partiting the space --------------------"<<endl;
	/* ------------------- phase 3: partiting the space -------------------- */
	/*tss, ss: (tmp) subspace*/
	/*divide the 3D bounding box into subspaces, each ss has several oriented faces, each face is part of some plane*/
	/*here V and E are the ones from the bbox, but faces are later shared with curves*/
	//int ssvernum;				//# of all V
	//float* ssver;				//coord of all V
	//int ssedgenum;			//# of all E
	//int* ssedge;				//two ends of all E
	//int ssfacenum;			//# of all F
	//int* ssfaceedgenum;		//# of Es in each F
	//int** ssface;				//E ids of each face
	//int* ssface_planeindex;	//plane ID cooresponding each face
	//int ssspacenum;			//# of spaces
	//int* ssspacefacenum;		//# of faces in each space
	//int** ssspace;			//F ids of each space
	//int** ssspace_planeside;	//which side of the oriented plane is each face in each space

	float pbbox[6];			//bouding box of the processed contour
	float enlargeratio[3];
	//temporary data
	floatvector tssver;			// tmp subspace ver
	intvector tssedge;
	vector<intvector> tssface;
	intvector tssface_planeindex;
	vector<intvector> tssspace;
	vector<intvector> tssspace_planeside;

	float scale [3];
	float minScale = INFINITY;
	float maxScale = -1;
	for (int i = 0; i < 3; i++){
		pbbox[i] = bbox[i];
		pbbox[i + 3] = bbox[i + 3];

		if (pbbox[i] == pbbox[i + 3]){
			scale [i] = 1;
		} else {
			scale [i] = abs(pbbox[i + 3]-pbbox[i]);
			if (scale[i] < minScale){
				minScale = scale[i];
			}
			if (scale[i] > maxScale){
				maxScale = scale[i];
			}
		}
	}

	for (int i = 0; i < 3; ++i){
		enlargeratio[i] = POLYMENDERENLARGERATIO;
	}

	// !! start partition !!
	SpacePartitioner partitioner;
	partitioner.partition(planenum, pparam, pbbox, enlargeratio,
		tssver, tssedge, tssface, tssface_planeindex, tssspace, tssspace_planeside, 0);

	// copy subspace info out 
	int ssvernum;
	float* ssver;
	int ssedgenum;
	int* ssedge;
	int ssfacenum;
	int* ssfaceedgenum;
	int** ssface;
	int* ssface_planeindex;
	int ssspacenum;
	int* ssspacefacenum;
	int** ssspace;
	int** ssspace_planeside;
	int versize = (int)tssver.size();
	ssvernum = versize / 3;
	ssver = new float[versize];
	for (int i = 0; i < versize; i++){
		ssver[i] = tssver[i];
	}
	int edgesize = (int)tssedge.size();
	ssedgenum = edgesize / 2;
	ssedge = new int[edgesize];
	for (int i = 0; i < edgesize; i++)
		ssedge[i] = tssedge[i];

	ssfacenum = (int)tssface.size();
	ssface = new int*[ssfacenum];
	ssfaceedgenum = new int[ssfacenum];
	ssface_planeindex = new int[ssfacenum];
	for (int i = 0; i < ssfacenum; i++){
		ssface_planeindex[i] = tssface_planeindex[i];
		ssfaceedgenum[i] = (int)tssface[i].size();
		ssface[i] = new int[ssfaceedgenum[i]];
		for (int j = 0; j < ssfaceedgenum[i]; j++){
			ssface[i][j] = tssface[i][j];
		}
	}

	ssspacenum = (int)tssspace.size();
	ssspace = new int*[ssspacenum];
	ssspace_planeside = new int*[ssspacenum];
	ssspacefacenum = new int[ssspacenum];
	for (int i = 0; i < ssspacenum; i++){
		ssspacefacenum[i] = (int)tssspace[i].size();
		ssspace[i] = new int[ssspacefacenum[i]];
		ssspace_planeside[i] = new int[ssspacefacenum[i]];
		for (int j = 0; j < ssspacefacenum[i]; j++){
			ssspace[i][j] = tssspace[i][j];
			ssspace_planeside[i][j] = tssspace_planeside[i][j];
		}
	}

//    partitioner.wfilePartition((ar.t_savedir+"partition.txt").c_str(), ssvernum, ssver, ssedgenum, ssedge, ssfacenum, ssfaceedgenum,
//            ssface, ssface_planeindex, ssspacenum, ssspacefacenum, ssspace, ssspace_planeside);
#if _DEBUG_LU
	//partitioner.wfilePartition("mmdebug/partition.txt", ssvernum, ssver, ssedgenum, ssedge, ssfacenum, ssfaceedgenum,
		//ssface, ssface_planeindex, ssspacenum, ssspacefacenum, ssspace, ssspace_planeside);
	/*partitioner.wfilePartition((ar.t_savedir+"partition.txt").c_str(), ssvernum, ssver, ssedgenum, ssedge, ssfacenum, ssfaceedgenum,
		ssface, ssface_planeindex, ssspacenum, ssspacefacenum, ssspace, ssspace_planeside);*/
#endif

    cout<<"-------------- phase 4: partition contours into space faces ---------------"<<endl;
	/* -------------- phase 4: partition contours into space faces --------------- */
	/*ctrf: contour face*/
	/*put conours into space faces constructed in last phase*/
	//int* ctrfvernum;			//# of V in each face
	//float** ctrfverpos;		//coord of V in each face
	//int** ctrfvertype;		//type of V in each face (completely in, or on the shared edge)
	//int** ctrfverval;			//?
	//int* ctrfedgenum;			//# of E in each face
	//int** ctrfedge;			//two ends of each E in each face, and left/right materials
	//int** ctrfedgetype;		//type of E in each face
	//int** ctrfedgeval;		//?
	//int** ctrfedgeancestor;	//?

	int* ctrfvernum;
	float** ctrfverpos;
	int** ctrfvertype;
	int** ctrfverval;
	int* ctrfedgenum;
	int** ctrfedge;
	int** ctrfedgetype;
	int** ctrfedgeval;
	int** ctrfedgeancestor;
	ctrfvernum = new int[ssfacenum];
	ctrfverpos = new float*[ssfacenum];
	ctrfvertype = new int*[ssfacenum];
	ctrfverval = new int*[ssfacenum];

	ctrfedgenum = new int[ssfacenum];
	ctrfedge = new int*[ssfacenum];
	ctrfedgetype = new int*[ssfacenum];
	ctrfedgeval = new int*[ssfacenum];
	ctrfedgeancestor = new int*[ssfacenum];
	for (int i = 0; i < ssfacenum; i++){
		ctrfvernum[i] = ctrfedgenum[i] = 0;
		ctrfverpos[i] = NULL;
		ctrfvertype[i] = NULL;
		ctrfverval[i] = NULL;

		ctrfedge[i] = NULL;
		ctrfedgetype[i] = NULL;
		ctrfedgeval[i] = NULL;
		ctrfedgeancestor[i] = NULL;
	}

	// !! start partition contours !!
	ContourHandler::putContourIntoFace(planenum, ssvernum,
		ssver, ssedgenum, ssedge, ssfacenum, ssfaceedgenum, ssface, ssface_planeindex, ssspacenum,
		ssspacefacenum, ssspace, pctrvernum, pctrvers, pctredgenum, pctredges, pparam, ver2planelistpos, ver2planelist,
		ctrfvernum, ctrfverpos, ctrfvertype, ctrfverval, ctrfedgenum, ctrfedge, ctrfedgetype, ctrfedgeval, ctrfedgeancestor,
		0);

#if _DEBUG_LU
	ContourHandler::writeContourInFace(ssfacenum, ctrfvernum, ctrfverpos, ctrfvertype,
		ctrfverval, ctrfedgenum, ctrfedge, ctrfedgetype, ctrfedgeval);
#endif

	float enlargedCellFrame [6];
	float tmin, tmax;
	for( int i = 0; i < 3; i ++)
	{
		tmin = pbbox[ i ];
		tmax = pbbox[ i + 3 ];
		if( tmin == tmax)
		{
			tmin = tmin  - 10;
			tmax = tmax + 10;
		}
		enlargedCellFrame[ i ] = (( 1 + enlargeratio[ i ])*tmin + ( 1 - enlargeratio[ i])*tmax)/2;
		enlargedCellFrame[ i + 3 ] = (( 1 - enlargeratio[ i ])*tmin + ( 1 + enlargeratio[ i])*tmax)/2;
	}

	// load 3D volume
	ar.loadVol(volfilename, ar.unifyTransVec, ar.unifyScaleRatio);

#if _CG_LOG
	if(ar.inputBBox){
		cout << "[+] User given BBox loaded... " << endl;
	}else{
		cout << "[-] No user-specified BBox. " << endl;
	}
	if(ar.inputVol){
		cout << "[+] User given Volume loaded... " << endl;
	}else{
		cout << "[-] No Volume info. " << endl;
	}
	cout << endl;
#endif

    cout<<"------------------ MING: construct Arrangement --------------------"<<endl;
	///* ------------------- MING: construct Arrangement -------------------- */
	ar.constructArrangement(
		ssvernum, ssedgenum, ssfacenum, ssspacenum,
		ssver, ssedge, ssface, ssspace,
		ssfaceedgenum, ssspacefacenum, ssspace_planeside, ssface_planeindex, 
		
		planenum, pparam, enlargedCellFrame,

		ctrfvernum, ctrfedgenum,
		ctrfverpos, ctrfvertype, ctrfverval,
		ctrfedge, ctrfedgetype, ctrfedgeval
		);
	///* -------------------------------------------------------------------- */


    cout<<"------------- release all memory used here by Lu's code ---------------"<<endl;
	/* -------------- release all memory used here by Lu's code --------------- */
	delete[] filenames[0];
	delete[] filenames;
	
	for (int i = 0; i < planenum; ++i){
		delete[] pctrvers[i];
		delete[] pctredges[i];
		delete[] ver2planelistpos[i];
	}
	delete[] pparam;
	delete[] pctrvernum;
	delete[] pctredgenum;
	delete[] pctrvers;
	delete[] pctredges;
	delete[] ver2planelistpos;

	for (int i = 0; i < ssfacenum; ++i){
		delete[] ssface[i];
	}
	for (int i = 0; i < ssspacenum; ++i){
		delete[] ssspace[i];
		delete[] ssspace_planeside[i];
	}
	delete[] ssver;
	delete[] ssedge;
	delete[] ssfaceedgenum;
	delete[] ssface_planeindex;
	delete[] ssspacefacenum;
	delete[] ssface;
	delete[] ssspace;
	delete[] ssspace_planeside;

	for (int i = 0; i < ssfacenum; ++i){
		delete[] ctrfverpos[i];
		delete[] ctrfvertype[i];
		delete[] ctrfverval[i];
		delete[] ctrfedge[i];
		delete[] ctrfedgetype[i];
		delete[] ctrfedgeval[i];
		delete[] ctrfedgeancestor[i];
	}
	delete[] ctrfvernum;
	delete[] ctrfverpos;
	delete[] ctrfvertype;
	delete[] ctrfverval;
	delete[] ctrfedgenum;
	delete[] ctrfedge;
	delete[] ctrfedgetype;
	delete[] ctrfedgeval;
	delete[] ctrfedgeancestor;
}
