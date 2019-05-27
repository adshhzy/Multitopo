#include "../SpacePartitioner/SpacePartitioner.h"

/*
1.detect the common edge between all the planes
2.reorganize the subspace, so that 
first one lies above plane0, below plane1
second one above plane0  above plane1
third one - below 0 above 1
fourth one - below 0 below 1
*/
void SpacePartitioner::reorganizeCMNL
(
 int& comnedge,
 vector<intvector>&ssface, intvector& ssface_planeindex, 
 vector<intvector>& ssspace, vector<intvector>& ssspace_planeside
 )
{
	//find the faces on the planes
	int facenum = ssface.size();

	int faceind[ 4 ];
	int ind = 0;
	for( int i= 0; i< facenum; i ++ )
	{
			if( ssface_planeindex[ i ] >= 0 )
			{
				faceind[ ind ++ ] = i;
			}
	}
	
	//the same edge on the four faces
	int faceedgenum[ 2 ] = 	{(int)(ssface[ faceind[ 0 ]].size()), (int)(ssface[ faceind[ 1 ]].size())};
	bool found = false;
    for( int i = 0; i < faceedgenum[ 0 ] ; i ++ )
	{
		int edgei = ssface[ faceind[ 0  ]][  i ];

		for( int j = 0; j < faceedgenum[ 1 ]; j ++ )
		{
			if( edgei == ssface[ faceind[ 1 ] ][ j ])
			{
				found = true;
				comnedge = edgei;
				break;
			}
		}

		if(found )
			break;
	}

	//mark the subspace
	int subspacesides[ 4 ][ 2 ]; //0 - below 1 - above
	
	for( int i = 0; i < 4; i ++ )
	{
		int facesize = ssspace_planeside[  i ].size();
		for( int j = 0; j < facesize; j ++ )
		{
			int facei = ssspace[ i ][ j ];
		
			int planei = ssface_planeindex[ facei ];
			if( planei < 0 )continue;		
			
			subspacesides[ i ][ planei ] = ssspace_planeside[ i ][ j ];
		}
	}

	//sort the subspace
	int spacelist[ 4 ];
	//0 - 10 1 - 01 2 - 11 3 - 00
	for( int i = 0; i < 4; i++ )
	{
		if ( (subspacesides[ i ][ 0 ] == 1) && ( subspacesides[ i ][ 1 ] == 0 ))
		{
			spacelist[ 0 ] = i;
		}
		else if ( (subspacesides[ i ][ 0 ] == 0) && ( subspacesides[ i ][ 1 ] == 1 ))
		{
			spacelist[ 1 ] = i;
		}
		else if ( (subspacesides[ i ][ 0 ] == 1) && ( subspacesides[ i ][ 1 ] == 1 ))
		{
			spacelist[ 2 ] = i;
		}
		else if ( (subspacesides[ i ][ 0 ] == 0 ) && ( subspacesides[ i ][ 1 ] == 0 ))
		{
			spacelist[ 3 ] = i;
		}
	}
	
	vector<intvector> tssspace;
	vector<intvector> tssspace_planeside;
	tssspace.resize( 4 );
	tssspace_planeside.resize( 4 );
	for( int i = 0; i < 4; i ++ )
	{
		int size = ssspace[ i ].size() ;
		tssspace[ i ].resize( size );
		memcpy( &tssspace[ i ][ 0 ], &ssspace[ i ][ 0 ], sizeof( int ) * size);

		size = ssspace_planeside[ i ].size();
		tssspace_planeside[ i ].resize( size );
		memcpy(&tssspace_planeside[ i ][ 0 ], &ssspace_planeside[ i ][ 0 ], sizeof( int ) * size);		

		ssspace[ i ].clear();
		ssspace_planeside[ i ].clear();
	}

	for( int i = 0; i < 4; i ++)
	{
		int curspacei = spacelist[ i ];
		int size = tssspace[ curspacei ].size();
		ssspace[ i ].resize(size );
		memcpy(&ssspace[ i ][ 0 ], &tssspace[ curspacei ][ 0 ], sizeof( int ) * size);
		
		size = tssspace_planeside[ curspacei ].size();
		ssspace_planeside[ i ].resize( size );
		memcpy( &ssspace_planeside[ i ][ 0], &tssspace_planeside[ curspacei ][ 0 ], sizeof(int)*size);

		tssspace[ curspacei ].clear();
		tssspace_planeside[ curspacei ].clear();
	}	
	tssspace.clear();
	tssspace_planeside.clear();
}

void SpacePartitioner::mapOneSubspaceCMNL
(  int spacei, intvector& ssedge,int& comnedge,
 vector<intvector>&ssface, vector<intvector>& ssspace,
 //reset the edges in the face, index starts from 0 and condensed
 vector<intvector>& nssface,  intvector& oedgelist,
 intvector& nssedge, intvector&overlist, //reset vertex in edgelist
	int& ncomnedge
 )
{
	//go through the face in the subspace
	int facenum = ssspace[ spacei ].size();
	nssface.resize( facenum );
	for( int i = 0; i < facenum; i ++ )
	{
		int facei = ssspace[ spacei ][ i ];
		int fedgenum = ssface[ facei ].size();
		nssface[ i ].resize( fedgenum );
		
		//set nssface
		for( int j = 0; j < fedgenum; j ++ )
		{
			int edgei = ssface[ facei ][ j ];
			
			//trying find in the oedgelist
			int oedgelistsize = oedgelist.size();
			bool findind = false;
			for(int k = 0; k < oedgelistsize; k ++ )
			{
				if( oedgelist[ k ] == edgei )
				{
					findind = true;
					nssface[ i ][ j ] = k;	//has already been found before
					break;
				}
			}
			if( !findind )	//an old edge hasn't been seen before
			{
				oedgelist.push_back( edgei );
				nssface[ i ][ j ] = oedgelistsize;
			}	
		}
	}

	//go through the edges in this subspace
	int oedgelistsize = oedgelist.size();
	nssedge.resize( oedgelistsize * 2 );

	for( int i = 0; i < oedgelistsize; i ++ )
	{
		int edgei = oedgelist[ i ];
		int vers[ 2 ] = {ssedge[ 2*edgei], ssedge[ 2*edgei + 1]};

		//go through the overlist to see if they are old vertice not 
		//seen before or not
		int overlistsize = overlist.size();
		for( int j = 0; j < 2; j ++ )	//for these two vertices
		{
			int veri = vers[ j ];
			bool found = false;
			for( int k = 0; k < overlistsize; k ++ )
			{
				if( overlist[ k ] == veri )
				{
					found = true;
					nssedge[ 2*i + j ] = k;
					break;
				}			
			}
			if( !found)
			{
				overlist.push_back( veri );
				nssedge[ 2*i + j ] = overlistsize;
				overlistsize ++;
			}
		}
	}

	//set the ncomnedge
	oedgelistsize = oedgelist.size();
	for( int i = 0; i < oedgelistsize; i ++ )
	{
		if (  oedgelist[ i ] == comnedge )
		{
			ncomnedge = i;
			break;
		}
	}
}

void SpacePartitioner::markVerCMNL( float param[ 4 ],	 floatvector& ssver, intvector& overlist, intvector& vermark)
{	
	int vernum = overlist.size();

	for( int i = 0; i < vernum; i ++ )
	{
		int veri = overlist[ i ];
		
		if( vermark[ i ] == 0 )	//already marked on, it must be on common edge
		{			
			continue;
		}

		float val = MyMath::dotProduct( param, &ssver[ veri * 3 ]);
		val = val - param[ 3 ] ;

		if( MyMath::isEqualInToler( val, 0,  TOLERANCE_THREE ) )
		{
			vermark[ i ] = 0;
		}
		else if( val > 0 )
			vermark[ i ] = 1;
		else
			vermark[ i  ] = -1;		
	}
}

//edge type
const int EDGE_NOINTERSECT = 0;	// no intersection between the edge and the plane
const int EDGE_TOUCHONEVER_ABOVE = 1;	//the plane crosses one vertex, and the edge is above the plane
const int EDGE_TOUCHONEVER_BELOW = 4;	//similar to the one above, except that the edge is below
//const int EDGE_TOUCHONEVER = 1;	//the plane crosses one vertex of the edge
const int EDGE_INTERSECT = 2;	//the plane intersects with the edge with one intersectionp point
const int EDGE_TOUCHEDGE = 3;	//the edge is on the plane
//value corresponding to the edge type
//type = 0 --- 1 - above the plane -1 - below the plane
//type = 1 --- 0 - the first vertex 1 - the second vertex
//type = 2 --- #new edge position, the edge above the plane is at the old place, new edge pushed back,
			//and the new vertex is the second vertex of the edge at the old posiiton,and first vertex
			//of the new edge
//type = 3 --- whatever, it doesn't matter

//split the edges
void SpacePartitioner::splitEdgeCMNL
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
)
{	
	int ssvernum = ssver.size()/3;
	int ssedgenum = ssedge.size()/2;

	// int oldvernum = overlist.size();	//before splitting, vertex number
	int oldedgenum = oedgelist.size();	//before splitting, edge number
	
	edgetype.resize( oldedgenum );
	edgeval.resize( oldedgenum );

	int nveris[ 2 ];
	int overis[ 2 ];
	int vermark1, vermark2;
	float interpt[ 3 ];
	for(int i = 0; i < oldedgenum; i ++ )
	{
		int edgei = oedgelist[ i ];

		nveris[ 0 ] = nssedge[ 2*i ];
		nveris[ 1 ] = nssedge[ 2*i + 1 ];

		//above the plane	1 & 1		
		//below the plane	-1 & -1
		//touch one vertex  0 & whatever or whatever & 0
		//edge on the plane	0 & 0
		//edge intersects with the plane 1 & -1
		vermark1 = vermark[ nveris[ 0 ]];
		vermark2 = vermark[ nveris[ 1 ]];

		if( vermark1 == 0 )
		{
			if( vermark2 == 0 )
			{
				edgetype[ i ] = EDGE_TOUCHEDGE;
				edgeval[ i ] = 0;				
			}
			else
			{
				if( vermark2 == 1 )
					edgetype[ i ] = EDGE_TOUCHONEVER_ABOVE;
				else
					 edgetype[ i ] = EDGE_TOUCHONEVER_BELOW;
				edgeval[ i ] = 0;
			}
			continue;
		}
		
		//the first vertex is not on the plane
		if( vermark2 == 0 )	//the second one is on, but the first one is not!
		{
			if( vermark1 == 1 )
				edgetype[ i ] = EDGE_TOUCHONEVER_ABOVE;
			else
				edgetype[ i ] = EDGE_TOUCHONEVER_BELOW;
			edgeval[ i ] = 1;			
			continue;
		}

		//both of them are not 0!
		if( vermark1 * vermark2 == 1 )	//the are of the same sign
		{
			edgetype[ i ] = EDGE_NOINTERSECT;	//no intersection
			edgeval[ i ] = vermark1;
			continue;
		}

		//edge and plane intersect
		//ssver, ssedge, edgetype, edgeval
		//if vermark1 = -1, change the sequence of the vertices in the edge
		if( vermark1 == - 1)
		{
			//change the vertices of the edge
			int tval = ssedge[ 2 * edgei  ];
			ssedge[ 2 *edgei ] = ssedge[ 2*edgei + 1 ];
			ssedge[ 2*edgei + 1 ] = tval;		
		}

		//compute the intersection point
		edgetype[ i ] = EDGE_INTERSECT;
		edgeval[ i ] = ssedgenum;

		overis[ 0 ] = overlist[ nveris [ 0 ]] ;
		overis[ 1 ] = overlist[ nveris[ 1 ]];
		float tvals[ 2 ] ;
		for( int j = 0; j < 2; j ++ )
		{
			tvals[ j ] = MyMath::dotProduct( param, &ssver[ 3*overis[ j ]]);
			tvals[ j ] = tvals[ j ] - param[ 3 ];
		}
		MyMath::getPtOnSeg(&ssver[ 3*overis[ 0 ]], &ssver[ 3*overis[ 1 ]], tvals[ 0 ]/(tvals[ 0 ] - tvals[1]),
			interpt);

		ssver.push_back( interpt[ 0 ]);
		ssver.push_back( interpt[ 1 ]);
		ssver.push_back( interpt[ 2 ]);

		ssedge.resize( 2*ssedgenum + 2);
        memcpy( &ssedge[ 2*ssedgenum ], &ssedge[ 2*edgei ], sizeof( int ) * 2);

		ssedge[ 2*edgei + 1 ] = ssvernum;
		ssedge[2*ssedgenum ] = ssvernum;

		ssvernum ++;
		ssedgenum ++;
	}
}

const int FACE_NOINTER = 0;	//no intersection
const int FACE_INTER = 1;
const int FACE_TOUCH_ABOVE = 2;	//plane crosses one edge on it, face is above the plane
const int FACE_TOUCH_BELOW = 3;	//........ face is below the plane
//FACE_NOINTER, val: 1 / -1, above or below
//FACE_INTER, val: #face index, the place for the new half below face 
//FACE_TOUCH, val: #one edge on the face is on the plane
//split the faces of the subspace

void SpacePartitioner::splitFaceCMNL
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
 )
{
	int ssfacenum = ssspace[ spacei ].size();
	int facenum = ssface.size();
	int ssedgenum = ssedge.size()/2;

	facemark.resize(ssfacenum );
	faceval.resize( ssfacenum );

	//go through the faces, split it when necessary
	for( int i = 0; i < ssfacenum; i ++ )
	{
		int facei = ssspace[ spacei ][ i ];
		
		int oldfaceedgenum = nssface[ i ].size();

		intvector aboveedges;
		intvector belowedges;
		intvector interpts;

		//if not -1, it is the index of the subspace edge on this face
		//which the plane crosses
		int touchface = -1;	
		int above = 0;	//1 - above, -1 - below, 0 - not decided
		//divide the edges of the face into two groups.
		for( int j = 0; j < oldfaceedgenum; j ++ )
		{
			int edgei = nssface[ i ][ j ];	//the new edge index, in oedgelist
			int oldedgei = ssface[ facei ][ j ];	//the old edge index in ssedge
			switch( edgemark[ edgei ])
			{
			case EDGE_INTERSECT:
				{
					aboveedges.push_back( oldedgei  );
					belowedges.push_back( edgeval[ edgei ]);

					interpts.push_back( ssedge[ 2*oldedgei + 1 ] );
				}
				break;
			case EDGE_NOINTERSECT:
				{
					if( edgeval[ edgei ] == 1 )
					{
						aboveedges.push_back( oldedgei );
						above = 1;
					}
					else
					{
						belowedges.push_back( oldedgei );
						above = -1;
					}
				}
				break;
			case EDGE_TOUCHEDGE:
				{
					touchface = oldedgei;	//one edge on the face is on the plane!
					break;
				}
				break;
			case EDGE_TOUCHONEVER_ABOVE:
				{
					interpts.push_back( ssedge[ 2*oldedgei + edgeval[ edgei ]]);
					aboveedges.push_back( oldedgei );
					above = 1;
				}
				break;
			case EDGE_TOUCHONEVER_BELOW:
				{
					interpts.push_back( ssedge[ 2*oldedgei + edgeval[ edgei ]]);
					belowedges.push_back( oldedgei );
					above = -1;
				}
				break;
			}
		}

		if( touchface != -1 )	//the plane crosses one edge on the face
		{
			if( above == 0)	//not decided yet, could only happen when the first 
				//edge is the one lying on the plane
			{
				int tedgei = nssface[ i ][ 1 ];
				if( edgemark[ tedgei ] == EDGE_TOUCHONEVER_ABOVE )
				{
					above = 1;
				}
				else if ( edgemark[ tedgei ] == EDGE_TOUCHONEVER_BELOW )
				{
					above = -1;
				}
				else
				{
					cout<<"EDGE_TOUCHONEVER_.... is expected!! "<<endl;
					above = 1;
				}
			}
			if( above == 1 )
				facemark[ i ] = FACE_TOUCH_ABOVE;
			else
				facemark[ i ] = FACE_TOUCH_BELOW;
			faceval[ i ] = touchface;	//the subspace edge it crosses
		}
		else if( interpts.size() <= 1 )	//the face doesn't intersect with the plane
		{
			facemark[ i ] = FACE_NOINTER;
			if( aboveedges.size() == 0 )	//below
			{
				faceval[ i ] = -1;	//below
			}
			else
				faceval[ i ] = 1;	//above
		}
		else	//intersect!
		{
			facemark[ i ] = FACE_INTER;
			faceval[ i ] = facenum;	//the position of the second half of the face

			//find the two intersection points
			int twointerpts[ 2 ];
			twointerpts[ 0 ] = interpts[ 0 ];
			if( interpts.size() >= 2)
			{
				for (unsigned int k = 1; k < interpts.size(); k++)
				{
					if( interpts[ k ] != twointerpts[ 0 ])
					{
						twointerpts[ 1 ] = interpts[ k ];
						break;
					}
				}
			}

			//split the face
			//ssedge, ssface
			ssedge.push_back( twointerpts[ 0 ]);
			ssedge.push_back( twointerpts[ 1 ]);
			aboveedges.push_back( ssedgenum );
			belowedges.push_back ( ssedgenum );
			
			ssface.push_back( belowedges );
			ssface[ facei ].clear();
			ssface[ facei ].resize( aboveedges.size());
			memcpy( &ssface[ facei ][ 0 ], &aboveedges[ 0 ], sizeof( int ) * aboveedges.size());			

			ssface_planeindex.push_back( ssface_planeindex[ facei ] );

			ssedgenum ++;            
			facenum ++;
		}

		aboveedges.clear();
		belowedges.clear();
		interpts.clear();
	}

}
void SpacePartitioner::splitSubspaceCMNL
(
 int planei,
 int spacei,
 bool above,		//true - push the half one above the plane at the end, false- the one below at the end
 vector<intvector>& ssspace,  
 vector<intvector>& ssspace_planeside,
 vector<intvector>& ssface,
intvector& ssface_planeindex,
 intvector& facemark, intvector& faceval
 )
{
	//divide the faces in the subspace into two groups
	intvector abovefaces;	//the faces above the plane
	intvector abovefacesides;
	intvector belowfaces;	//the faces below the plane
	intvector belowfacesides;
	intvector interedges;	//intersection edge list

	int facenum = ssspace[ spacei ].size();
	for(int i = 0; i < facenum; i ++ )
	{
		int oldfacei = ssspace[ spacei ][ i ];
		int planeside = ssspace_planeside[ spacei ][ i ];
		switch( facemark[ i ] )
		{
		case FACE_INTER:	//faceval saves the below half face
			{
				abovefaces.push_back( oldfacei );
				abovefacesides.push_back( planeside );
				belowfaces.push_back( faceval[ i ]);
				belowfacesides.push_back( planeside );
                
				int faceedgenum = ssface[ oldfacei ].size();
				interedges.push_back( ssface[ oldfacei ][ faceedgenum - 1] );
			}
			break;
		case FACE_NOINTER:	//faceval saves 1/-1 above or below
			{
				if( faceval[ i ] == 1 )
				{
					abovefaces.push_back( oldfacei );
					abovefacesides.push_back( planeside );
				}
				else
				{
					belowfaces.push_back( oldfacei );
					belowfacesides.push_back( planeside );
				}
			}
			break;
		case FACE_TOUCH_ABOVE:	//faeeval saves the touch edge index
			{
				abovefaces.push_back( oldfacei );
				abovefacesides.push_back( planeside );
				interedges.push_back( faceval[ i ]);
			}			
			break;
		case FACE_TOUCH_BELOW:	
			{
				belowfaces.push_back( oldfacei);
				belowfacesides.push_back( planeside );

				//if this exists, there must exist another face, which is the case above
				//no need to do the following operation.
			//	interedges.push_back( faceval[ i ]);
			}
			break;
		}
	}

	//////////////////////////////////////////////////////////////////////////	
	if( interedges.size() == 0 )	//this should not have happened
	{
		cout<<"Wasn't able to find the intersection edges!"<<endl;
		intvector tvec;
		ssspace.push_back( tvec );
		return;
	}
	//////////////////////////////////////////////////////////////////////////

	//split current subspace into two	
	facenum = ssface.size();
	ssface.push_back( interedges );
	ssface_planeindex.push_back( planei );
	abovefaces.push_back( facenum );
	belowfaces.push_back( facenum );
	abovefacesides.push_back( 1 );
	belowfacesides.push_back ( 0 );

	if( above )
	{
		ssspace.push_back( abovefaces );
		//ssspace[ spacei ].clear();
		int belowfacenum = belowfaces.size();
		ssspace[ spacei ].resize( belowfacenum );		
        memcpy(&ssspace[ spacei ][ 0 ], &belowfaces[ 0 ], sizeof( int ) * belowfacenum );		

		ssspace_planeside.push_back( abovefacesides );
		//ssspace_planeside[ spacei ].clear();
		ssspace_planeside[ spacei ].resize( belowfacenum );
		memcpy(&ssspace_planeside[ spacei ][ 0 ], 
			&belowfacesides[ 0 ], sizeof(int) * belowfacenum );

	}
	else
	{
		ssspace.push_back( belowfaces );
		//ssspace[ spacei ].clear();
		int abovefacenum = abovefaces.size();
		ssspace[ spacei ].resize( abovefacenum );
		memcpy(&ssspace[ spacei ][ 0 ], &abovefaces[ 0 ], sizeof( int ) * abovefacenum );

		ssspace_planeside.push_back( belowfacesides );
		//ssspace_planeside[ spacei ].clear();
		ssspace_planeside[ spacei ].resize( abovefacenum );
		memcpy(&ssspace_planeside[ spacei ][ 0 ], 
			&abovefacesides[ 0 ], sizeof(int) * abovefacenum );
	}

	abovefaces.clear();
	belowfaces.clear();
	interedges.clear();
	abovefacesides.clear();
	belowfacesides.clear();
}
//intersect the cutting plane with the last subspace or the second to the last subspace
void SpacePartitioner::intersectPlaneSubspaceCMNL
(
 int planei,
 float param[ 4 ],		//the parameter of the plane
 bool above,				//true intersect with the subspace above previous cutting plane, false - below
 floatvector& ssver, intvector& ssedge,  int& comnedge,
 vector<intvector>&ssface, intvector& ssface_planeindex, 
 vector<intvector>& ssspace, vector<intvector>& ssspace_planeside
 )
{

	//map the vertices, edge, and faces in this subspace,
	vector<intvector> nssface;
	intvector oedgelist;
	intvector nssedge;
	intvector overlist;
	int ncomnedge = -1;

	int spacei = ssspace.size() - 2;	//the second to the last one
	mapOneSubspaceCMNL
		(spacei, ssedge,  comnedge, ssface,  ssspace,
		nssface, oedgelist,	nssedge, overlist,	ncomnedge);

	oedgelist.resize( oedgelist.size( )) ;
	
//	writemapOneSubspaceCMNL(spacei, ssedge,  comnedge, ssface,  ssspace,
//		nssface, oedgelist,	nssedge, overlist,	ncomnedge);

	//mark all the vertices
	//find the position for the two vertices on the common edge
	int ncomnvers[ 2 ];
	ncomnvers[ 0 ] = nssedge[ 2*ncomnedge ];
	ncomnvers[ 1 ] = nssedge[ 2*ncomnedge + 1];

	//initialize the vermark.
	intvector markVer;
	int versize = overlist.size();
	markVer.resize( versize );
	for( int  i = 0; i < versize; i ++ )
	{		
		markVer[ i ] = -1;
	}
	markVer[ ncomnvers[ 0 ]] = markVer[ ncomnvers[ 1 ] ] = 0;
	markVerCMNL(param,  ssver,  overlist, markVer);	

    //split edges 
	intvector edgemark;	//the types of the edges
	intvector edgeval;	//the values of the edges
	
	splitEdgeCMNL( param,ssver, ssedge, oedgelist, nssedge, overlist,	markVer,
		 edgemark, edgeval);

	//split the face
	intvector facemark;
	intvector faceval;

	splitFaceCMNL( spacei, ssspace, ssface, ssface_planeindex, nssface,ssedge,oedgelist,edgemark,edgeval,
		facemark,faceval);

	//////////////////////////////////////////////////////////////////////////
	/*writeOVerOEdgeList(overlist,  oedgelist	);
	FILE* fout = fopen( "mmdebug/vermark.txt", "w");
	writeOneIntvector(markVer, fout );
	fout = fopen( "mmdebug/edgetype.txt", "w");
	writeOneIntvector( edgemark, fout );
	fout = fopen( "mmdebug/edgeval.txt", "w");
	writeOneIntvector( edgeval, fout );
	fout = fopen("mmdebug/oface.txt","w");
	writeOneIntvector( ssspace[ spacei ],  fout );
	fout = fopen("mmdebug/facetype.txt","w");
	writeOneIntvector(facemark,  fout );
	fout = fopen( "mmdebug/faceval.txt", "w");
	writeOneIntvector(faceval, fout );*/
	//////////////////////////////////////////////////////////////////////////

	//split the subspace
	splitSubspaceCMNL(planei, spacei,above,	ssspace, ssspace_planeside, ssface, ssface_planeindex, facemark, faceval);

	//temp vars
	for (unsigned int i = 0; i< nssface.size(); i++)
	{
		nssface[ i ].clear();
	}
	nssface.clear();
	oedgelist.clear();
	nssedge.clear();
	overlist.clear();
	markVer.clear();
	edgemark.clear();
	edgeval.clear();
	facemark.clear();
	faceval.clear();
    
}

void SpacePartitioner::partition_ComnLine
(
 const int planenum, float*& param,
 const float boundingbox[ 6 ], const float enlarge[ 3 ],
 floatvector& ssver, intvector& ssedge,  int& comnedge,
 vector<intvector>&ssface, intvector& ssface_planeindex, 
 vector<intvector>& ssspace, vector<intvector>& ssspace_planeside,
 float comndir[ 3 ], float comnpt[ 3 ]  //common line of the cutting planes
 )
{
	//call the other general partition algorithm to partition the bounding box by
	//first two planes
	float* tparam = new float[ 32 ];
	memcpy( tparam, param, sizeof(float) * 8 );

	partition( 2, tparam,boundingbox,  enlarge, ssver, ssedge, ssface,
		ssface_planeindex,  ssspace,  ssspace_planeside, 0);
	//copy the 6 bounding planes parameters out
	memcpy( param + 4*planenum, tparam + 8, sizeof(float) * 24);
	delete []tparam;

	//detect the common edge between the first cutting planes.
	reorganizeCMNL( comnedge,ssface, ssface_planeindex, 
		ssspace, ssspace_planeside);
		
	////////////////////////////////////////////////////////////////////////////
	//for( int i = 0; i < planenum + 4; i ++ )
	//{
	//	cout<<param[ 4* i ]<<","
	//		<<param[ 4*i + 1]<<","
	//		<<param[ 4*i + 2]<<","
	//		<<param[ 4*i + 3]<<endl;
	//}
	//cout<<endl;
	//////////////////////////////////////////////////////////////////////////

	//insert the rest cutting planes
	for( int i = 2; i < planenum; i ++ )
	{
		intersectPlaneSubspaceCMNL(i,	
			param + 4 * i ,	true, ssver, ssedge, comnedge, ssface, 
			ssface_planeindex, 	 ssspace,  ssspace_planeside);		

		intersectPlaneSubspaceCMNL(	i,	
			param + 4 * i ,	false, ssver, ssedge, comnedge, ssface, 
			ssface_planeindex, 	 ssspace,  ssspace_planeside);		
	}
}
