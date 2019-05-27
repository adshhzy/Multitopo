#include "SpacePartitioner.h"

void SpacePartitioner::markVerOnCmnLZero( int plane0,
										floatvector& ssver, intvector& ssedge,
										 vector<intvector>&ssface, intvector& ssface_planeindex, 
										 vector<intvector>& ssspace, vector<intvector>& ssspace_planeside,
										 int*& verMark
										 )
{
	//if plane i, plane j, and plane 0 belong to the same radial set, and planei < plane0 planej < plane0.
	//find all the ssver that's on both planei and planej, mark them as 0 for verMark
	int radialSetNum = _gloabl_radialPlaneSet.size();
	for( int i = 0;  i < radialSetNum; i++ )
	{
		//check if plane0 is in the set, if yes, check if it's order in the set is >= 2.
		//which means the planes before it must have already been inserted, based on the
		//fact that set is ordered in the increasing order.
		set<int>::iterator iter = _gloabl_radialPlaneSet[ i ].find( plane0 );
		if( iter == _gloabl_radialPlaneSet[ i ].end() )
			continue;

		//the plane exists in the set
		iter = _gloabl_radialPlaneSet[ i ].begin();
		bool needClear = false;
		int planes[ 2 ] = { -1, -1 };
		if( *iter != plane0 )
		{
			planes[ 0 ] = *iter;
			iter ++;
			if( *iter != plane0 )
			{
				needClear = true;
				planes[ 1 ] = *iter;
			}
		}

		if( !needClear )	//no two planes in the set have already been inserted!
			continue;

		//yes, two planes have been inserted. find all the ssver that's on the two planes		
		int ssfaceNum = ssface_planeindex.size();		
		int* verMarkPlanes[ 2 ] ;
		verMarkPlanes[ 0 ] = new int[ ssver.size()/3 ];
		verMarkPlanes[ 1 ] = new int[ ssver.size()/3 ];
		for (unsigned int j = 0; j < ssver.size() / 3; j++)
		{
			verMarkPlanes[ 0 ][ j ] = verMarkPlanes[ 1 ][ j ]  = 0;

		}
		
		for( int j = 0; j < ssfaceNum; j ++ )
		{
			int planeI = ssface_planeindex[ j ];
			if( planeI < 0 )
				continue;

			if( planeI == planes[ 0 ] )
			{
				int edgeNum = ssface[ j ].size();
				for( int k = 0; k < edgeNum; k ++ )
				{
					int edgei = ssface[ j ][ k ];
					verMarkPlanes[ 0 ][ ssedge[ 2*edgei ] ] = 1;	//mark the vertex is on plane i
					verMarkPlanes[ 0 ][ ssedge[ 2*edgei + 1]] = 1;
				}
			}
			else if( planeI == planes[ 1  ])
			{
				int edgeNum = ssface[ j ].size();
				for( int k = 0; k < edgeNum; k ++ )
				{
					int edgei = ssface[ j ][ k ];
					verMarkPlanes[ 1 ][ ssedge[ 2*edgei ] ] = 1;	//mark the vertex is on plane i
					verMarkPlanes[ 1 ][ ssedge[ 2*edgei + 1]] = 1;
				}
			}
		}

		for (unsigned int i = 0; i < ssver.size() / 3; i++)
		{
			//only both the marks are active, set the vermark to 0
			if( verMarkPlanes[ 0 ][ i ] && verMarkPlanes[ 1 ][ i ] )
				verMark[ i ] = 0;	//the vertex should be on plane0
		}
		delete []verMarkPlanes[ 0 ];
		delete []verMarkPlanes[ 1 ];
	}
}

/**
* Function that take the parameter of the planes and cut the whole space into subspaces
* @param param the parameters of the planes, ax+by+cz=d, abcd represents one plane
* @param planenum the number of all the planes
* @param ver the resulting vertices of the subspaces, x, y , z, three number is one point
* @param vernum the number of the vertices
* @param edge the resulting edges of the subspaces, v1, v2, index of the two vertices in the ver
* @param edgenum the number of the edges in the subspace
* @param face the resulting faces, edge1,.... edgen, index of the edges in edge
* @param faceedgenum, how many edges are in the corresponding face
* @param facenum the faces number in all the subspaces
* @param faceside, which side of the plane this subspace is in
* @subspacenum the number of the subspaces
*
*/
void SpacePartitioner::insertOnePlane_CtrGraph(float planeparam[ 4 ], int planei,
									  floatvector& ssver, intvector& ssedge,
									  vector<intvector>&ssface, intvector& ssface_planeindex, 
									  vector<intvector>& ssspace, vector<intvector>& ssspace_planeside									  
									  )
{
	int vernum = ssver.size()/3;
	int edgenum = ssedge.size()/2;
	int facenum = ssface.size();
	// int spacenum = ssspace.size();

	//step1. mark all the vertices
	int* vermark = new int[ vernum ];
	getVerMark( planeparam, ssver, vermark);

	//this is the part the insertion of planei special for .ctrgraph.
	//to make sure that radial set of planes go through the same 
	//common line, detect all the subspace vertices the plane should pass
	//through based on the global variable: _gloabl_radialPlaneSet
	markVerOnCmnLZero( planei, ssver, ssedge, ssface, ssface_planeindex, ssspace, ssspace_planeside, vermark );

	//step2. mark all the edges, split it when needed
	//EVENT: NEW VERTEX, NEW EDGE
	//ssver, ssedge
	//mark of them
	//intersection point index for each edge
	//case 1: one vertex + one vertex -, then intersection point index, -1
	//case 2: both of them are + or -, then -1 -1
	//case 3: one of them is on plane, then touch vertex index, -1
	//case 4: both of them are on the plane, then touch vertex indices
	int* edgenewver = new int[ 2*edgenum ];
	//new edge index if the edge is split into two
	//case1: no intersection point or at least one vertex of the edge is on the plane, -1
	//case 2: intersection point exists, the negative new edge index is in it
	int* edgenewedge = new int[ edgenum ];
	int j = 0;
	for( int i = 0 ;i < edgenum ; i ++)
	{
		edgenewver[ j++ ] = -1;
		edgenewver[ j++ ] = -1;
		edgenewedge [ i ] = -1;
	}
	getEdgeMark( planeparam, ssver, vermark,ssedge,edgenewver,edgenewedge);

	//step3. mark all the faces, split it when needed
	//EVENT: NEW EDGE, NEW FACE
	//ssedge, ssface ssface_planeindex
	//mark of them
	int* facenewedge = new int[ facenum ];
	int* facenewface = new int[ facenum ];
	for( int i = 0; i < facenum ; i++)
	{
		facenewedge[ i ] = -1;
		facenewface[ i ] = -1;			
	}
	getFaceMark(vermark,edgenewver, edgenewedge,ssedge,ssface,ssface_planeindex,facenewedge, facenewface);

	//step4. go through all the subspaces, split it when needed
	//EVENT: NEW FACE, NEW SUBSPACE
	//ssface, ssspace, ssspace_planeside
	processSubspace(planei, ssedge, ssface, ssface_planeindex,ssspace, ssspace_planeside,vermark,facenewedge, facenewface );

	delete []vermark;
	delete []edgenewedge;
	delete []edgenewver;
	delete []facenewface;
	delete []facenewedge;
}


/**
* Function partition the space with the parameters of the planes
* @param param: parameters of all the planes  ax+by+cz=d, abcd represents one plane
* @param planenum: total planenumber
* @param boundingbox the bounding box of the space,to cut
* @param enlarge, the times to enlarge the bounding box
*/
void SpacePartitioner::partition_CtrGraph(const int planenum, float*& param,const float boundingbox[ 6 ], const float enlarge[ 3 ],floatvector& ssver, intvector& ssedge, vector<intvector>&ssface,
								 intvector& ssface_planeindex, 
								 vector<intvector>& ssspace, vector<intvector>& ssspace_planeside
								)
{
	//int oldsize = planenum * 4; 
	//move them to the resize( ) part
	//step0 normalize all the planes parameters
	int oldsize = planenum * 4; 

	/*for( int i = 0; i < planenum; i ++)
	{
	float veclen = MyMath::vectorLen( param[ i * 4 ], param[ i*4 + 1 ], param[ i*4+ 2 ]);
	for(int j = 0; j < 4; j ++)
	param[ 4*i + j] = param[ 4*i + j]/veclen;
	}*/

	//step1, compute the boundingbox
	float nbx[ 6 ];
	float tmin, tmax;
	for( int i = 0; i < 3; i ++)
	{
		tmin = boundingbox[ i ];
		tmax = boundingbox[ i + 3 ];
		if( tmin == tmax)
		{
			tmin = tmin  - 10;
			tmax = tmax + 10;
		}
		nbx[ i ] = (( 1 + enlarge[ i ])*tmin + ( 1 - enlarge[ i])*tmax)/2;
		nbx[ i + 3 ] = (( 1 - enlarge[ i ])*tmin + ( 1 + enlarge[ i])*tmax)/2;
	}

	//step2. add 6 planes into the param list
	for( int i = 0; i < 24;i ++)
		param[ oldsize + i ] = 0;
	param[ oldsize ] = -1;
	param[ oldsize + 3 ] = -nbx[ 0 ];

	param[ oldsize + 5 ] = 1;
	param[ oldsize + 7 ] = nbx[ 4 ];

	param[ oldsize + 8 ] = 1;
	param[ oldsize + 11 ] = nbx[ 3 ];

	param[ oldsize + 13] = -1;
	param[ oldsize + 15 ] = -nbx[ 1 ];

	param[ oldsize + 18] = 1;
	param[ oldsize + 19] = nbx[ 5 ];

	param[ oldsize + 22] = -1;
	param[ oldsize + 23 ] = -nbx[ 2 ];

	//step3. initialize the first space
	//	floatvector tssver;
	//	intvector tssedge;
	//	vector<intvector>tssface;
	//  intvector tssface_planeindex;
	//vector<intvector> tssspace;
	//vector<intvector> tssspace_planeside;

	//ver
	ssver.resize( 24 );
	ssver[ 0 ] = nbx[ 0 ]; ssver[ 1 ] = nbx[ 1 ]; ssver[ 2 ] = nbx[ 2 ];
	ssver[ 3 ] = nbx[ 3 ]; ssver[ 4 ] = nbx[ 1 ]; ssver[ 5 ] = nbx[ 2 ];
	ssver[ 6 ] = nbx[ 3 ]; ssver[ 7 ] = nbx[ 4 ]; ssver[ 8 ] = nbx[ 2 ];
	ssver[ 9 ] = nbx[ 0 ]; ssver[ 10 ] = nbx[ 4 ]; ssver[ 11 ] = nbx[ 2 ];
	ssver[ 12 ] = nbx[ 0 ]; ssver[ 13 ] = nbx[ 1 ]; ssver[ 14 ] = nbx[ 5 ];
	ssver[ 15 ] = nbx[ 3 ]; ssver[ 16 ] = nbx[ 1 ]; ssver[ 17 ] = nbx[ 5 ];
	ssver[ 18 ] = nbx[ 3 ]; ssver[ 19 ] = nbx[ 4 ]; ssver[ 20 ] = nbx[ 5 ];
	ssver[ 21 ] = nbx[ 0 ]; ssver[ 22 ] = nbx[ 4 ]; ssver[ 23 ] = nbx[ 5 ];

	//edge
	ssedge.resize( 24 );
	ssedge[ 0 ] = 0; ssedge[ 1 ] = 1;
	ssedge[ 2 ] = 1; ssedge[ 3 ] = 2;
	ssedge[ 4 ] = 2; ssedge[ 5 ] = 3;
	ssedge[ 6 ] = 3; ssedge[ 7 ] = 0;
	ssedge[ 8 ] = 4; ssedge[ 9 ] = 5;
	ssedge[ 10 ] = 5; ssedge[ 11 ] = 6;
	ssedge[ 12 ] = 6; ssedge[ 13 ] = 7;
	ssedge[ 14 ] = 7; ssedge[ 15 ] = 4;
	ssedge[ 16 ] = 0; ssedge[ 17 ] = 4;
	ssedge[ 18 ] = 1; ssedge[ 19 ] = 5;
	ssedge[ 20 ] = 2; ssedge[ 21 ] = 6;
	ssedge[ 22 ] = 3; ssedge[ 23 ] = 7;

	//face
	ssface.resize( 6 );
	ssface_planeindex.resize( 6 );
	for( int i = 0;i  < 6; i ++)
	{
		ssface_planeindex[ i ] = -(i+1);	//which plane the face is on
		ssface[ i ].resize( 4 );
	}
	ssface[ 0 ][ 0 ] = 0;	ssface[ 0 ][ 1 ] = 1;	ssface[ 0 ][ 2 ] = 2;	ssface[ 0 ][ 3 ] = 3;
	ssface[ 1 ][ 0 ] = 4;	ssface[ 1 ][ 1 ] = 5;	ssface[ 1 ][ 2 ] = 6;	ssface[ 1 ][ 3 ] = 7;
	ssface[ 2 ][ 0 ] = 0;	ssface[ 2 ][ 1 ] = 4;	ssface[ 2 ][ 2 ] = 8;	ssface[ 2 ][ 3 ] = 9;
	ssface[ 3 ][ 0 ] = 1;	ssface[ 3 ][ 1 ] = 5;	ssface[ 3 ][ 2 ] = 9;	ssface[ 3 ][ 3 ] = 10;
	ssface[ 4 ][ 0 ] = 2;	ssface[ 4 ][ 1 ] = 6;	ssface[ 4 ][ 2 ] = 10;	ssface[ 4 ][ 3 ] = 11;
	ssface[ 5 ][ 0 ] = 3;	ssface[ 5 ][ 1 ] = 7;	ssface[ 5 ][ 2 ] = 8;	ssface[ 5 ][ 3 ] = 11;

	//subpace
	ssspace.resize( 1 );
	ssspace[ 0 ].resize( 6 );
	ssspace_planeside.resize( 1 );
	ssspace_planeside[ 0 ].resize( 6 );
	for( int i = 0; i < 6; i ++)
	{
		ssspace[ 0 ][ i ] = i;
		ssspace_planeside[ 0 ][ i ] = 0;	//negative side of the planes
	}

	for( int i = 0; i < planenum; i ++)		
	{
		insertOnePlane_CtrGraph(param + 4*i, i, ssver, ssedge, ssface, ssface_planeindex, ssspace, ssspace_planeside);
	}
}



