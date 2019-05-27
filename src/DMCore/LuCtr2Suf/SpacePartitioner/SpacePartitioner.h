#ifndef _SPACEPARTITIONER_H
#define _SPACEPARTITIONER_H

#include "../config.h"
#include "../Util/VerEdgePlaneOp.h"

class SpacePartitioner
{
public:
	SpacePartitioner();
	~SpacePartitioner();

    void inline getVerMark( float planeparam[ 4 ],
		floatvector& ssver, int*& vermark);
	void inline getEdgeMark( float planeparam[ 4 ],
		floatvector& ssver, int*& vermark,
		intvector& ssedge, int*& edgenewver, int*& edgenewedge);
	void inline getFaceMark(int*& vermark,int*& edgenewver, int*& edgenewedge,intvector& ssedge,
		vector<intvector>&ssface, intvector& ssface_planeindex,
		int*& facenewedge, int*&facenewface );
	void processSubspace(int planei,intvector& ssedge,vector<intvector>& ssface, intvector& ssface_planeindex,vector<intvector>&ssspace, vector<intvector>& ssspace_planeside,
		int*& vermark,		int*& facenewedge, int*& facenewface	);
	void insertOnePlane(float planeparam[ 4 ], int planei,
		floatvector& ssver, intvector& ssedge,
		vector<intvector>&ssface, intvector& ssface_planeindex, 
		vector<intvector>& ssspace, vector<intvector>& ssspace_planeside,
		int cmnVer1 = -1, int cmnVer2 = -1);
	void partition(const int planenum, float*& param,const float boundingbox[ 6 ], const float enlarge[ 3 ],
		floatvector& ssver, intvector& ssedge, vector<intvector>&ssface, intvector& ssface_planeindex, 
		vector<intvector>& ssspace, vector<intvector>& ssspace_planeside,
		int radialPlanesNum);	

	//for contour graph, mark those vertices on common line as 0 in vermark
	void markVerOnCmnLZero( int plane0,
		floatvector& ssver, intvector& ssedge,
		vector<intvector>&ssface, intvector& ssface_planeindex, 
		vector<intvector>& ssspace, vector<intvector>& ssspace_planeside,
		int*& verMark
		);
	void insertOnePlane_CtrGraph(float planeparam[ 4 ], int planei,
		floatvector& ssver, intvector& ssedge,
		vector<intvector>&ssface, intvector& ssface_planeindex, 
		vector<intvector>& ssspace, vector<intvector>& ssspace_planeside									  
		);
	void partition_CtrGraph(const int planenum, float*& param,const float boundingbox[ 6 ], const float enlarge[ 3 ],floatvector& ssver, intvector& ssedge, vector<intvector>&ssface,
		intvector& ssface_planeindex, 
		vector<intvector>& ssspace, vector<intvector>& ssspace_planeside
		);

	//for cutting planes sharing with common line case
	 inline void reorganizeCMNL
		(
		int& comnedge,
		vector<intvector>&ssface, intvector& ssface_planeindex, 
		vector<intvector>& ssspace, vector<intvector>& ssspace_planeside
		);

	inline void mapOneSubspaceCMNL
		(  int spacei, intvector& ssedge,int& comnedge,
		vector<intvector>&ssface, vector<intvector>& ssspace,
		//reset the edges in the face, index starts from 0 and condensed
		vector<intvector>& nssface,  intvector& oedgelist,
		intvector& nssedge, intvector&overlist, //reset vertex in edgelist
		int& ncomnedge
		);
	inline void markVerCMNL
		(
		float param[ 4 ],	
		floatvector& ssver, intvector& overlist,
		intvector& vermark);

	inline void splitEdgeCMNL
		(
		float param[ 4 ],		//parameter of the plane

		floatvector& ssver,		//vertices positions
		intvector& ssedge,			//the original edge list
		intvector& oedgelist,		//the edges in current subspace
		intvector& nssedge,		//the edges composed of new vertices index
		intvector& overlist,		//vertices that are in current subspace
		intvector& vermark,			//vertex mark of the vertices

		intvector& edgetype,		//mark the edge type of the old edge list
		intvector& edgeval			//mark the edge value of the old edge list
		);

	inline void splitFaceCMNL
		(
		//int planei,
		int spacei, vector<intvector>& ssspace, 
		vector<intvector>& ssface,	//all the faces
		intvector& ssface_planeindex,
		vector<intvector>& nssface,
		intvector& ssedge,
		intvector& oedgelist, 
		intvector& edgemark, intvector& edgeval,

		intvector& facemark, intvector& faceval
		);

	void splitSubspaceCMNL
		(
		int planei,
		int spacei,
		bool above,		//true - push the half one above the plane at the end, false- the one below at the end
		vector<intvector>& ssspace,  
		vector<intvector>& ssspace_planeside,
		vector<intvector>& ssface,
		intvector& ssface_planeindex,
		intvector& facemark, intvector& faceval
		);
	void intersectPlaneSubspaceCMNL
		(
		int planei,
		float param[ 4 ],		//the parameter of the plane
		bool above,				//true intersect with the subspace above previous cutting plane, false - below
		floatvector& ssver, intvector& ssedge,  int& comnedge,
		vector<intvector>&ssface, intvector& ssface_planeindex, 
		vector<intvector>& ssspace, vector<intvector>& ssspace_planeside
		);


	void partition_ComnLine
		(
		const int planenum, float*& param,
		const float boundingbox[ 6 ], const float enlarge[ 3 ],
		floatvector& ssver, intvector& ssedge,  int& comnedge,
		vector<intvector>&ssface, intvector& ssface_planeindex, 
		vector<intvector>& ssspace, vector<intvector>& ssspace_planeside,
		float comndir[ 3 ], float comnpt[ 3 ]  //common line of the cutting planes
		);

	void wfilePartition(const char* fname,int ssvernum,float* ssver,int ssedgenum,int* ssedge,int ssfacenum,
	int* ssfaceedgenum,	int** ssface,int* ssface_planeindex,int ssspacenum,int* ssspacefacenum,	int** ssspace,	int** ssspace_planeside);

	void writemapOneSubspaceCMNL(int spacei, 
		intvector& ssedge,  int comnedge, 
		vector<intvector>& ssface,  vector<intvector>& ssspace,
		vector<intvector>& nssface, 		intvector& oedgelist,	
		intvector& nssedge, intvector& overlist,	
		int ncomnedge);	
	void writeOVerOEdgeList(intvector& overlist, intvector& oedgeslist);
	void writeOneIntvector(intvector& vec,  FILE* fout );
};

void inline SpacePartitioner::getVerMark( float planeparam[ 4 ],
										 floatvector& ssver, int*& vermark)
{
	//go through all the vertices and mark them all
	int vernum = (int)ssver.size()/ 3;
	float* vers = new float[ vernum * 3 ];
	for( int i = 0; i < 3 * vernum ;i ++)
	{
		vers[ i ] = ssver[ i ];
	}
	for( int i = 0; i < vernum ;i ++)
	{
		float val = MyMath::dotProduct( planeparam, vers + 3*i );
		if( MyMath::isEqualInToler(val, planeparam[ 3 ] , POINT_ON_PLANE_TOLERANCE)) 
		{
			vermark[ i ] = 0;
		}
		else if( val < planeparam[ 3 ] )
			vermark[ i ] = -1;
		else
			vermark[ i ] = 1;		
	}

	delete []vers;
};

void inline SpacePartitioner::getEdgeMark( float planeparam[ 4 ], 
										  floatvector& ssver, int*& vermark,
											intvector& ssedge, int*& edgenewver, int*& edgenewedge)
{
	float ver1[ 3 ];
	float ver2[ 3 ];
	float interpt[ 3 ];
	int vernum = ssver.size()/3;
	int oldedgenum = ssedge.size()/2;
	int edgenum = oldedgenum;

	//go through each edge, split it when needed, and add the intersection point and the new edge into edgenew**
	int vind[ 2 ];
	int j = 0;	//the index of edgenewver, to avoid multiplication.
	for( int i = 0 ;i < oldedgenum ; i ++)
	{
		j += 2;	//index of the edgenewver and the vertices of the edge in ssedge

		vind[ 0 ] = ssedge[ j - 2 ]; vind[  1 ]= ssedge[ j - 1 ];
		int mul = vermark[ vind[ 0 ]] * vermark[ vind[ 1 ]];

		//no intersection point
		if( mul == 1)
			continue;

		//intersection point exists
		if( mul == -1 )
		{
			for( int k = 0; k < 3; k ++ )
			{
				ver1[ k ] = ssver[ 3 * vind[ 0 ] + k];
				ver2[ k ] = ssver[ 3 * vind[ 1 ] + k];
			}
			//compute intersection point
			interPt_PlaneEdge( planeparam, ver1, ver2, interpt);

			//split the edge
			ssver.push_back( interpt[ 0 ]);
			ssver.push_back( interpt[ 1 ]);
			ssver.push_back( interpt[ 2 ]);
			ssedge.push_back( vernum );
			if( vermark[ vind[ 0 ]] == -1 )
			{
				ssedge.push_back( vind[ 0 ]);
				ssedge[ j - 2 ] = vernum;
			}
			else
			{
				ssedge.push_back( vind[ 1 ] );
				ssedge[ j - 1 ] = vernum;
			}
			//set the two mark arrays
			edgenewver[ j - 2 ] = vernum;
			edgenewedge[ i ] = edgenum;
            
			vernum ++;
			edgenum ++;
			continue;
		}
        
		//at least one of them is 0
		if( vermark[ vind [ 0 ]] == 0)
		{
			edgenewver[ j - 2 ] = vind[ 0 ];
			if( vermark[ vind[ 1 ] ] == 0 )
				edgenewver[ j - 1 ] = vind[ 1 ];
		}
		else	//vermark[vind[1]] must be 0 and vind[ 0] is not 0
		{
			edgenewver[ j - 2 ] = vind[ 1 ];
		}
	}
};

void inline SpacePartitioner::getFaceMark(int*& vermark,int*& edgenewver, int*& edgenewedge,
										  intvector& ssedge,vector<intvector>&ssface, intvector& ssface_planeindex,
								int*& facenewedge, int*&facenewface )
{
	int edgenum = ssedge.size()/2;
	int oldfacenum= ssface.size();
	int facenum = oldfacenum;
	
    intvector singletouchedge;		//ont vertex of the edge is on the plane
	intvector doubletouchedge;	//both of the vertices of the edge are on the plane
	intvector splitedge;			//the edge intersects with the plane
	
    //go through all the faces, split it when needed
	for( int i = 0; i <  oldfacenum; i ++)
	{
		//gather all the touch edges, and gather all the intersection edges
		int tedgelen = ssface[ i ].size();

		for( int j = 0;  j < tedgelen ; j ++) 
		{
			int edgei = ssface[ i ][ j ];

			if( edgenewver[ edgei * 2 ] != -1 )	//v1, -1, or (v1, v2) or (newv, -1)
			{
				if( edgenewver[  edgei * 2 + 1 ] != -1)	//(v1, v2)
					doubletouchedge.push_back( edgei);
				else			//(v1, -1), (newv, -1)
				{
					if( edgenewedge[ edgei ] != -1)	//(newv, -1)
					{
						splitedge.push_back( edgei );
					}	
					else	//(v1, -1)
						singletouchedge.push_back( edgei );
				}					
			}
			else //(-1.-1), (-1, v2)
			{
				if( edgenewver[ edgei * 2 + 1] != -1 )
					singletouchedge.push_back( edgei );
			}
		}

//#ifdef debug		//run in my debug mode
//		cout<<"single touch edge number:"<<singletouchedge.size()<<endl
//			<<"double touch edge number:"<<doubletouchedge.size()<<endl
//			<<"splitedge size:"<<splitedge.size()<<endl;
//#endif 
		//case1: there exists one double touch edge
		if( doubletouchedge.size() == 1)
		{
			//////////////////////////////////////////////////////////////////////////
			//cout<<"double touch edge index:"<<doubletouchedge[ 0 ]<<endl;
			//////////////////////////////////////////////////////////////////////////

			facenewedge[ i ] = doubletouchedge[ 0 ];
		}
		//case2: there exists one edge splitting the face
		else if( !splitedge.empty() || singletouchedge.size() == 4)
		{
			intvector newedgevers;	//new vertices of the new edge
			for (unsigned int j = 0; j < splitedge.size(); j++)
			{
				newedgevers.push_back(edgenewver[ 2*splitedge[ j ] ]);
			}
			for (unsigned int j = 0; j < singletouchedge.size(); j++)
			{
				newedgevers.push_back( edgenewver[ 2 * singletouchedge[ j ]]);
			}
			// int newedgevernum = newedgevers.size();
#ifdef debug
	//		cout<<"number of new vertices for the new edge:"<<newedgevernum<<endl;
	//		for( int j = 0; j < newedgevernum; j ++)
	//			cout<<newedgevers[ j ]<<" ";
	//		cout<<"----------------------"<<endl;
#endif
			int v1, v2;
			v1 = newedgevers[ 0 ];
			for (unsigned int j = 1; j < newedgevers.size(); j++)
			{
				v2 = newedgevers[ j ];
				if( v2 != v1)
					break;
			}
			
			//add new edge and split the old face
			ssedge.push_back( v1 );
			ssedge.push_back( v2 );
			intvector oldfaceedges;
			intvector newfaceedges;
			oldfaceedges.push_back( edgenum );
			newfaceedges.push_back( edgenum );
			facenewedge[ i ] = edgenum;
			edgenum++;
			
			for( int j = 0; j < tedgelen; j ++ )
			{
				int edgei = ssface [ i ][ j ];
				//if split edge
				if( edgenewedge[ edgei ] != -1 )
				{
					oldfaceedges.push_back( edgei );
					newfaceedges.push_back( edgenewedge[ edgei ]);
					continue;
				}
                				
				int edgevers[ 2 ] = {ssedge[ edgei * 2 ], ssedge[ edgei * 2 + 1 ]};
				int tsum = vermark[ edgevers[ 0 ]] + vermark[edgevers[ 1 ]];
				if( (tsum == 2) || (tsum == 1) )	//1+1 0+1
				{
					oldfaceedges.push_back( edgei );
				}
				else if ( (tsum == -2) || (tsum == -1))	//-1 + -1 0 + -1
				{
					newfaceedges.push_back( edgei );
				}
			}
			ssface_planeindex.push_back( ssface_planeindex[ i ] );
			ssface.push_back( newfaceedges );
			ssface[ i ].clear();
			ssface[ i ] = oldfaceedges;

			//set face mark
			facenewface[ i ] = facenum;
			facenum++;

			//clear the temporary vectors
			newedgevers.clear();
			oldfaceedges.clear();
			newfaceedges.clear();			
		}
		singletouchedge.clear();
		doubletouchedge.clear();
		splitedge.clear();
	}
};




#endif
