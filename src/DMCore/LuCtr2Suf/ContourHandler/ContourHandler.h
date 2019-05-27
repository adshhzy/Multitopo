#ifndef _CONTOURHANDLER_H
#define _CONTOURHANDLER_H

#include "../config.h"
#include "../Util/SetOperation.h"
#include "../Util/sort.h"
#include "../Util/VerEdgePlaneOp.h"
#include <cmath>

class ContourHandler{
public:
	ContourHandler();
	~ContourHandler();

	//read contours
	void static readOneContour(char* filename, 
		floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
		float bbox[ 6 ], bool& bboxset);
	void static readContour(const int filenum, char** filenames, 
			floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
			float bbox[ 6 ], bool& bboxset);
	void static writeOneContourMM(
		FILE* fout, float* pctrvers, int pctrvernum,
		int* pctredges,	int pctredgenum, float* param);
	void static writeContourMM(
		int planenum,
		float** pctrvers ,
		int* pctrvernum,
		int** pctredges,
		int* pctredgenum,
		float* pparam,
		char* fname
	);

	//read contours, with intersection point informatoin. 
	void static readOneContour2(char* filename, 
		floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
		float bbox[ 6 ], bool& bboxset,
		vector<intvector>& 	ver2planelistpos,
		vector<intvector>& ver2planelist
		);
	void static readContour2(const int filenum, char** filenames, 
		floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
		float bbox[ 6 ], bool& bboxset,
		vector<intvector>& 	ver2planelistpos,
		vector<intvector>& ver2planelist);

	void static readOneContour2_v2(char* filename, 
		floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
		float bbox[ 6 ], bool& bboxset,
		vector<intvector>& 	ver2planelistpos,
		vector<intvector>& ver2planelist,
		int & radialPlaneNum
		);
	void static readContour2_v2(const int filenum, char** filenames, 
		floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
		float bbox[ 6 ], bool& bboxset,
		vector<intvector>& 	ver2planelistpos,
		vector<intvector>& ver2planelist, int& radialPlaneNum);


	void static readOneContour3_v3(char* filename, float bbox[ 6 ], bool& bboxset);
	void static readContour3_v3(const int filenum, char** filenames, 
		floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
		float bbox[ 6 ], bool& bboxset,
		vector<intvector>& 	ver2planelistpos,
		vector<intvector>& ver2planelist);

	//preprocess contour
	static void computeInterPtList(floatvector& param,vector<floatvector>& ctrvers,
		vector<intvector>& ctredges,int becut,	int tocut,
		floatvector& interptpos,intvector& interptind);
	static bool makeConsistent_InterPt_2Planes(	floatvector& param,vector<floatvector>& ctrvers,
		vector<intvector>& ctredges,int referplane,int adjustplane,	vector<intvector>& pairlist	);
    static bool makeConsistent_InterPt_2Planes2(	floatvector& param,vector<floatvector>& ctrvers,
        vector<intvector>& ctredges,int referplane,int adjustplane,	vector<intvector>& pairlist	);
	static void makeConsistent_InterPt(
		int& planenum, float*& pparam,	int*& pctrvernum, float**& pctrvers,
		int*& pctredgenum, int**& pctredges,int**& ver2planelistpos, vector<intvector>& ver2planelist);
	static void preProcDataSingleMat(
		int& planenum,float*& pparam,float**& pctrvers,int*& pctrvernum,int**& pctredges,int*& pctredgenum,
		int**& ver2planelistpos, vector<intvector>& ver2planelist
		);
	
	//make the contours on the planes sharing one line consistent
	static void interCommonLineWithOnePlane(
		float param[ 4 ],					//parameter of the plane
		floatvector& ctrvers, intvector& ctrvermark,		//the contour vertices on this cutting plane, and their marks relative to the common edge
		intvector& ctredges,	//contour edges on this cutting plane
		float comndir[ 3 ], float comnpt[ 3 ],		//common line parameters
		intvector& interptPos			//remeber the position of the vertices which lie on the common line
		);
	static void reserverPlanes
		(
		int tplanenum,
		vector<intvector> interptpos,
		int*& reservePlaneList, int& reserveNum
		);
	static void clearTheInputData(
		int& planenum, float*& pparam,
		int*& pctrvernum, float**& pctrvers,
		int*& pctredgenum, int**& pctredges);
	static bool adjustInterPt(
		int* reservelist, int reservenum,
		vector<intvector>& interptpos,
		vector<floatvector>& ctrvers,
		float comndir[ 3 ], float comnpt[ 3 ] );
	static int flipNormalToHalfSpace
		( int* reserverlist, int reservenum,
		float*& pparam, 
		vector<intvector>& ctredges,
		float comndir[ 3 ] );
	static void sortSelectMarkList
		( int planenum, int markval,	bool inc, 
		int*& markreservelist,	int*& reservelist, float*& dotProduct,
		int*& sortedplanelist,	int startind	);
	static void sortReservedPlanes
		(int& planenum, float*& pparam,
		int*& pctrvernum, float**& pctrvers,int**& pctrvermark,
		int*& pctredgenum, int**& pctredges,vector<intvector>& vermark,
		vector<floatvector>& ctrvers, vector<intvector>& ctredges,	//
		int* reservelist, int reservenum, int pickind,
		float comndir[ 3 ] );
	static  int** makeConsistentSingleMatCommonLine(
		int& planenum, float*& pparam,
		int*& pctrvernum, float**& pctrvers, 
		int*& pctredgenum, int**& pctredges,
		float comndir[ 3 ], float comnpt[ 3 ]);

	//after space partition, put contours into their faces of all the subspaces
	static inline void getElement2Planes(
		int planenum,
		int ssvernum, int ssedgenum, int* ssedge,
		int ssfacenum, int* ssfaceedgenum, int** ssface,int* ssface_planeindex,
		int ssspacenum, int* ssspacefacenum, int** ssspace,
		vector<intset>& ver2planes,	vector<intset>& edge2planes,vector<intset>& plane2planes,
		vector<intvector>& planefaces,vector<intvector>& planeedges,vector<intvector>& planevers
		,int radialPlaneNum);
	static void inline putContourIntoFace_OP_Mark(	//one plane mark
		int**& markprop,
		int planei, 
		int* ver2planesnum,	int** ver2planes,	//sorted planes that the vertex is on
		int plane2planesnum, int* plane2planes,	//sorted planes that intersects with planei
		int planefacesnum, int* planefaces,
		int planeedgesnum, int* planeedges,
		int planevernum, int* planevers,
		int planenum,	int pctrvernum, float* pctrver,
		int pctredgenum, int* pctredges,float* pparam,float* ssver,int* ssedge, int* ssfaceedgenum, 
		int** ssface);
	static int markCompare(int* mark1,int* mark2,int columnnum	);
//	static bool isSmaller(	int* mark1, int* mark2, int columnnum);
	static void inline sortMark( int**& markprop, int rownum, int columnum);
	static int inline findPos(int** markprop, int rownum, int columnnum, int* vermark);
	static void inline markVer(
		int**& markprop, int rownum,
		int planei,
		int planevernum, float*& planever,
		int**& ver2planelistpos, vector<intvector>& ver2planelist,
		int plane2planesnum, int* plane2planes,
		float* pparam,
		//result
		int*& vertype,
		int*& verval
		);
	static void inline gatherTopology(
		int* planefaces, int planefacesnum,
		int* planeedges, int planeedgesnum,
		int* planevers,int planeversnum,
		int** ssface, int* ssfaceedgenum,
		int* ssedge,

		int**& ver2edges, int*& ver2edgesnum,
		int**& ver2faces, int*& ver2facesnum,
		int**& edge2faces, int*& edge2facesnum
		);
	static void inline markEdge(
		int pctrvernum,float* pctrver,
		int pctredgenum, int* pctredges,
		int* vertype, int* verval, 

		int* planevers,	int planeversnum,		//in the increasing order of vertices indices
		int* planeedges,	int planeedgesnum,	//in the increasing order of edges indices
		int* planefaces,	int planefacesnum,	//in the increasing order of faces indices

		int* ssedge,
		int** ssface, int* ssfaceedgenum,

		vector<floatvector>& ctrfverposvec, vector<intvector>& ctrfvertypevec,vector<intvector>&ctrfvervalvec, 
		vector<intvector>&ctrfedgevec, vector<intvector>& ctrfedgetypevec, vector<intvector>& ctrfedgevalvec,
		vector<intvector>& ctrfedgeancestorvec
		);
	static void inline putContourIntoFace_OnePlane(
		int planei,
		int* ver2planesnum,	int** ver2planes,	//sorted planes that the vertex is on
		int* edge2planesnum, int** edge2planes,	//sorted planes that the edge is on
		int plane2planesnum, int* plane2planes,	//sorted planes that intersects with planei
		int planefacesnum, int* planefaces,
		int planeedgesnum, int* planeedges,
		int planevernum, int* planevers,

		int planenum,
		int pctrvernum, float* pctrver,
		int pctredgenum, int* pctredges,
		float* pparam,
		int**& ver2planelistpos, vector<intvector>& ver2planelist,

		int ssvernum,float* ssver,
		int ssedgenum, int* ssedge, int ssfacenum, int* ssfaceedgenum, 
		int** ssface,int* ssface_planeindex,

		//result
		vector<floatvector>& ctrfverposvec, vector<intvector>& ctrfvertypevec,vector<intvector>&ctrfvervalvec, 
		vector<intvector>&ctrfedgevec, vector<intvector>& ctrfedgetypevec, vector<intvector>& ctrfedgevalvec,
		vector<intvector>& ctrfedgeancestorvec
		/*int*& ctrfvernum,float**& ctrfverpos,int**& ctrfvertype,int**& ctrfverval,
		int*& ctrfedgenum, int**& ctrfedge, int**& ctrfedgetype, int**& ctrfedgeval*/
		);
   static void putContourIntoFace
		(
		int planenum,
		int ssvernum, float* ssver,int ssedgenum, int* ssedge,
		int ssfacenum, int* ssfaceedgenum, int** ssface,int* ssface_planeindex,
		int ssspacenum, int* ssspacefacenum, int** ssspace,

		int* pctrvernum, float** pctrverpos, int* pctredgenum, int** pctredge, float* pparam, 
		int**& ver2planelistpos, vector<intvector>& ver2planelist,

		//result
		int*& ctrfvernum,float**& ctrfverpos,int**& ctrfvertype,int**& ctrfverval,
		int*& ctrfedgenum, int**& ctrfedge, int**& ctrfedgetype, int**& ctrfedgeval,
		int**& ctrfedgeancestor, int radialPlaneNum
		);

	//for the cutting plane sharing common line case
   static void gatherFacesOnPlaneCMNL(int planenum,  int*& planeFaces, int ssfacenum, int* ssface_planeindex );
   static  void markCtrVerCMNL
	   (
	   int planefaces[ 2 ],		//the two faces on the plane in processing
	   float* ssver, int* ssedge, int* ssfaceedgenum, int** ssface, 
	   float comndir[ 3 ], float comnpt[ 3 ], //the common line for these cutting planes
	   int comnedgei, int comnveri[ 2 ],		//the two vertices of the common edge in ssver
	   int pctrvernum, float* pctrverpos,  
	   int* ctrvermark,		//the marked vertices

	   int*& ctrvertype, int *& ctrverval,
	   vector<floatvector>& fctrverpos,
	   vector<intvector> &fctrvertype,
	   vector<intvector>& fctrverval,
	   int* vermap[ 2 ]
	   );
    static void putControuEdgeIntoFaceCMNL
	   (
	   int planefaces[ 2 ],	
	   int* vertype, int* verval,
	   int pctredgenum, int* pctredges,	//four is a group

	   vector<intvector>& fctredges,
	   vector<intvector>& fctreedgetype,
	   vector<intvector>& fctredgeval,
	   vector<intvector>& fctredgeancestor,

	   int* vermap[ 2 ]
	   );
    static void putContourIntoFaceCMNL
	   (
	   int planenum,
	   int ssvernum, float* ssver,int ssedgenum, int* ssedge,
	   int ssfacenum, int* ssfaceedgenum, int** ssface,int* ssface_planeindex,
	   int ssspacenum, int* ssspacefacenum, int** ssspace,

	   int* pctrvernum, float** pctrverpos, int* pctredgenum, int** pctredge, float* pparam, 
	   int** pctrvermark,

	   float comndir[ 3 ], float comnpt[ 3 ], int comnedgei, int comnveri[ 2 ],
	   //result
	   int*& ctrfvernum,float**& ctrfverpos,int**& ctrfvertype,int**& ctrfverval,
	   int*& ctrfedgenum, int**& ctrfedge, int**& ctrfedgetype, int**& ctrfedgeval,
	   int**& ctrfedgeancestor
	   );

	static void writeOneContourOut
		(floatvector & ctrvers,	intvector& ctredges,
		int* ctrvermark,		intvector interptind,		
		int becut,		int tocut,		floatvector& param		);

   static void writeElement2Planes_db(
	   int ssvernum, int ssedgenum, int ssfacenum, int planenum,
	   int* ver2planesnum,int** ver2planes,int* edge2planesnum,
	   int** edge2planes,int* plane2planesnum,int** plane2planes,
	   int* planefacesnum,int** planefaces,int* planeedgesnum,
	   int** planeedges,int* planeversnum,int** planevers
	   );
  static void writeContourInFace(
	   int facenum,
	   int* ctrfvernum,float** ctrfverpos,int** ctrfvertype,int** ctrfverval,
	   int* ctrfedgenum, int** ctrfedge, int** ctrfedgetype, int** ctrfedgeval
	   );

  static void writeInterPtsOnCommonLine(
	  int* reservelist, int reservenum,  vector<intvector>& interptpos,	  vector<floatvector>& ctrvers  );

  static void writeContoursOut
	  (
	  float* param,
	  vector<floatvector>& ctrvers,
	  vector<intvector>& ctredges,
	  vector<intvector> interptpos,
      float comndir[ 3 ], float comnpt[ 3 ]
	  );

 static void writeParamOut(int* reservelist, int reservenum, float* param  , int pickind);

 static void writeFaceCtr(int facei, 
	int fctrvernum,	float* fctrverpos,	int* fvertype, int* fverval,
	int fctredgenum,	int* fedge, int* fedgetype, int* fedgeval	 );

 static void writePlaneCtr
	 (int planei, 
	 floatvector& ctrvers,
	 intvector& vermark,
	 intvector& ctredges
	 );
};

#endif
