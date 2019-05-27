#include "../SpacePartitioner/SpacePartitioner.h"

//////////////////////////////////////////////////////////////////////////
bool dbOutput = false;
//////////////////////////////////////////////////////////////////////////
void SpacePartitioner::processSubspace(
									   int planei,
									   intvector& ssedge, vector<intvector>& ssface, intvector& ssface_planeindex,
									     vector<intvector>&ssspace, vector<intvector>& ssspace_planeside,
										 int*& vermark,
									   int*& facenewedge, int*& facenewface  )
{	
	int facenum = ssface.size();
	int spacenum = ssspace.size();

	//go through each subspace, split it when needed
	set<int> newfaceedgesSet;
	intvector newfaceedges;

	for( int subspacei = 0; subspacei < spacenum;  subspacei ++)
	{
		//////////////////////////////////////////////////////////////////////////
		/*if( planei == 4 && subspacei == 7 )
			dbOutput = true;
		else
			dbOutput = false;*/
		//////////////////////////////////////////////////////////////////////////

		int tspacefacenum = ssspace[ subspacei ].size();

		//////////////////////////////////////////////////////////////////////////
	/*	if( dbOutput )
			cout<<"BEFORE new face size:"<<newfaceedgesSet.size()<<endl;*/
		//////////////////////////////////////////////////////////////////////////

		for( int i = 0; i < tspacefacenum; i ++)
		{	
			int tcurfacenewedge = facenewedge[ ssspace[ subspacei ] [ i ]] ;
			if(tcurfacenewedge == -1)	
				continue;
			//newfaceedges.push_back( tcurfacenewedge );
			newfaceedgesSet.insert( tcurfacenewedge );

			//////////////////////////////////////////////////////////////////////////
			//if ( dbOutput )
			//{
			//	//output the face configuration first
			//	int faceI = ssspace[ subspacei ][ i ];
			//	cout<<"face vertices:";
			//	for( int dbi = 0; dbi < ssface[ faceI ].size(); dbi++ )
			//	{
			//		int edgeI = ssface[ faceI ][ dbi ];
			//		cout<<ssedge[ 2 * edgeI ]<<" "<<ssedge[ 2 * edgeI + 1]<<"; ";
			//	}
			//	cout<<endl;

			//	cout<<"add new edge: ind"<<tcurfacenewedge<<" composition:"
			//	<<ssedge[ 2 * tcurfacenewedge ] <<","
			//	<<ssedge[ 2 * tcurfacenewedge +  1] <<"\n";
			//}
			//////////////////////////////////////////////////////////////////////////

		}
		//////////////////////////////////////////////////////////////////////////
		/*if( dbOutput )
			cout<<"new face size:"<<newfaceedgesSet.size()<<endl;*/
		//////////////////////////////////////////////////////////////////////////
		
		if( newfaceedgesSet.size() <= 2 ) //no intersection edge or passing through one of the edge
		{
			newfaceedgesSet.clear();
			continue;
		}

		set<int>::iterator iter = newfaceedgesSet.begin();

		while( iter != newfaceedgesSet.end() )
		{
			newfaceedges.push_back( *iter );
			iter ++;
		}
		newfaceedgesSet.clear();
		//////////////////////////////////////////////////////////////////////////
		/*if( dbOutput )
		{
			cout<<"newfacedges size:"<<newfaceedges.size()<<endl;
			for( int dbI = 0; dbI < newfaceedges.size(); dbI ++ )
			{
				int edgeI = newfaceedges[ dbI ];
				cout<<ssedge[ 2 * edgeI ]<<","<<ssedge[ 2 * edgeI + 1 ]<<"; ";
			}
		}
		cout<<endl;*/
		//////////////////////////////////////////////////////////////////////////

		//new face
		ssface.push_back( newfaceedges);
		ssface_planeindex.push_back( planei );
		newfaceedges.clear();

		//split the old subspace
		intvector oldsspacefaces;
		intvector newsspacefaces;
		intvector oldspacefacesides;
		intvector newspacefacesides;
		oldsspacefaces.push_back( facenum );		//index starts from 0, that's why before adding facenum
		newsspacefaces.push_back( facenum );
		facenum ++;											//one new face is added
		oldspacefacesides.push_back( 1 );		//the part of current subspace that is above current plane
		newspacefacesides.push_back( 0 );		//below the plane
		for( int i = 0; i < tspacefacenum; i ++)
		{
			int tfacei = ssspace[ subspacei ][ i ];
			int tnfacei = facenewface[ tfacei ];
			int tfaceside = ssspace_planeside[ subspacei ][ i ];
			if( tnfacei != -1 ) // the face is split into two
			{
				oldsspacefaces.push_back( tfacei );
				newsspacefaces.push_back( tnfacei );
				oldspacefacesides.push_back( tfaceside );
				newspacefacesides.push_back( tfaceside );
			}
			else	//the face is not split into two: case1, maybe touch some vertices, case2, no vertex is on the plane
			{
				//go through all the vertex on the face, stop until some vermark is 1 or -1
				int tedgenum = ssface[ tfacei ].size();
				int tvers[ 2 ];
#ifdef debug
		//		bool isset = false;
#endif
				for( int j = 0; j < tedgenum ;j ++)
				{
					int tedgei = ssface[ tfacei ][ j ];
					tvers[0] = ssedge[ 2*tedgei ];
					tvers[ 1 ] = ssedge[ 2*tedgei + 1];
					int tsum = vermark[ tvers[ 0 ]] +vermark[ tvers[ 1 ] ] ;
					if( tsum < 0 )	//the face should be in the newspacefaces
					{
						newsspacefaces.push_back( tfacei );
						newspacefacesides.push_back( tfaceside );
#ifdef debug
		//				isset = true;
#endif
						break;
					}
					else if( tsum > 0 )	//the face should be in the oldfspacefaces
					{
						oldsspacefaces.push_back( tfacei );
						oldspacefacesides.push_back( tfaceside );
#ifdef debug
		//				isset = true;
#endif
						break;
					}
				}
#ifdef debug
			//	if( !isset)	cout<<"face "<<tfacei <<" was not put into either old subspace or new subspace!"<<endl;
#endif
			}
		}
		ssspace.push_back( newsspacefaces );
		ssspace_planeside.push_back( newspacefacesides );
		ssspace[ subspacei ].clear();
		ssspace[ subspacei ] = oldsspacefaces;
		ssspace_planeside[ subspacei ].clear();
		ssspace_planeside[ subspacei ] = oldspacefacesides;

		newfaceedges.clear();
		oldspacefacesides.clear();
		oldsspacefaces.clear();
		newspacefacesides.clear();
		newsspacefaces.clear();
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
void SpacePartitioner::insertOnePlane(float planeparam[ 4 ], int planei,
									  floatvector& ssver, intvector& ssedge,
									  vector<intvector>&ssface, intvector& ssface_planeindex, 
									  vector<intvector>& ssspace, vector<intvector>& ssspace_planeside,
									  int cmnVer1, int cmnVer2
									  )
{
	int vernum = ssver.size()/3;
	int edgenum = ssedge.size()/2;
	int facenum = ssface.size();
	// int spacenum = ssspace.size();
	
	//step1. mark all the vertices
	int* vermark = new int[ vernum ];
	getVerMark( planeparam, ssver, vermark);

	if ( cmnVer1 != -1)
	{
		vermark[ cmnVer1 ] = 0;
		vermark[ cmnVer2 ] = 0;
	}

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
void SpacePartitioner::partition(const int planenum, float*& param,const float boundingbox[ 6 ], const float enlarge[ 3 ],floatvector& ssver, intvector& ssedge, vector<intvector>&ssface,
					  intvector& ssface_planeindex, 
					  vector<intvector>& ssspace, vector<intvector>& ssspace_planeside,
					  int radialPlanesNum)
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

	//MING: LUBBOX

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

	//insert one plane by one
	if( radialPlanesNum == -1 || radialPlanesNum == 0 )
	{
		for( int i = 0; i < planenum; i ++)
	//	for( int i = 0 ; i< 1; i ++)
		{
			insertOnePlane(param + 4*i, i, ssver, ssedge, ssface, ssface_planeindex, ssspace, ssspace_planeside);
		}
	}
	else
	{
		for( int i = 0 ;i < 2; i ++ )
		{
			insertOnePlane(param + 4*i, i, ssver, ssedge, ssface, ssface_planeindex, ssspace, ssspace_planeside);

			/////////////
			//////////////

			//return;
		}

		//find the commmon intersection points for the two inserted planes
		int cmnVer1, cmnVer2;
		cmnVer1 = cmnVer2 = -1;
		int ssfaceNum = ssface_planeindex.size();
		int* verCounter = new int[ ssver.size()/3 ];
		for (unsigned int j = 0; j < ssver.size() / 3; j++)
		{
			verCounter[ j ] = 0;
		}
		for( int j = 0; j < ssfaceNum; j ++ )
		{
			if( ssface_planeindex[ j ] >= 0 )
			{
				int edgeNum = ssface[ j ].size();
				for( int k = 0; k < edgeNum; k ++ )
				{
					int edgei = ssface[ j ][ k ];
					verCounter[ ssedge[ 2*edgei ] ] ++;
					verCounter[ ssedge[ 2*edgei + 1]]++;
				}
			}
		}


		for (unsigned int i = 0; i < ssver.size() / 3; i++)
		{
			if( verCounter[ i ] != 8 )
				continue;

            if( cmnVer1 == -1 )
				cmnVer1 = i;
			else
				cmnVer2 = i;
		}

		//////////////////////////////////////////////////////////////////////////
		/*for(int i = 0; i < ssver.size()/3; i ++ )
		{
			cout<<verCounter[ i ]<<",";
		}
		cout<<endl;
		cout<<"cmnver1:"<<cmnVer1<<"cmnver2:"<<cmnVer2<<endl;
		return;*/
		//////////////////////////////////////////////////////////////////////////

		for( int i = 2; i < radialPlanesNum; i ++ )
		{
			insertOnePlane(param + 4*i, i, ssver, ssedge, ssface, ssface_planeindex, ssspace, ssspace_planeside,
				cmnVer1, cmnVer2);

			
		}	

		for( int i = radialPlanesNum; i < planenum; i ++ )
		{
			insertOnePlane(param + 4*i, i, ssver, ssedge, ssface,
				ssface_planeindex, ssspace, ssspace_planeside);

		//	return;
		}

	}
//	for( int i = 0; i < planenum * 4 ; i++)
//		cout<<param[ i ]<<" ";
	cout<<endl;
}




void SpacePartitioner::wfilePartition(const char* fname,int ssvernum,float* ssver,int ssedgenum,int* ssedge,int ssfacenum,
									  int* ssfaceedgenum,	int** ssface,int* ssface_planeindex,int ssspacenum,int* ssspacefacenum,	int** ssspace,	int** ssspace_planeside)
{
	//{Length[resultVer], Length[resultEdge], Length[resultFace], Length[resultSubspace], resultVer,resultEdge, resultFace, resultSubspace} >> "partitionResult.txt";
	FILE*fout = fopen( fname, "w");
	if( fout == NULL )
	{
		cout<<"Unable to open the file : "<<fname<<"to write!"<<endl;
	}
	
	fprintf( fout, "{");
	//vertex number, edgenumber, facenumber, subspacenumber
	fprintf( fout,"%d,%d,%d,%d,{", ssvernum,ssedgenum, ssfacenum,  ssspacenum );
	//ver
	fprintf( fout, "{%f,%f,%f}", ssver[ 0 ], ssver[ 1 ], ssver[ 2 ]);
	for( int i = 1; i < ssvernum; i ++)
	{
		fprintf( fout, ",{%f,%f,%f}", ssver[3*i ], ssver[ 3*i + 1 ], ssver[ 3*i + 2 ]);
	}
	//edge
	fprintf( fout, "},{{%d,%d}", ssedge[ 0 ] + 1, ssedge[ 1 ] + 1 );
	for( int i = 1; i < ssedgenum; i ++)
	{
		fprintf( fout, ",{%d,%d}", ssedge[ 2*i ] +1, ssedge[ i * 2 + 1 ] + 1);
	}
	//face
	fprintf( fout, "},{{{%d", ssface[ 0 ][ 0 ] + 1);
	for( int i = 1; i < ssfaceedgenum[ 0 ];  i ++)
	{
		fprintf(fout, ",%d", ssface[ 0 ][ i ] + 1);
	}
	int faceside = ssface_planeindex[ 0 ];
	if( faceside >= 0 )
		faceside += 1;
	fprintf( fout, "},{%d}}", faceside);
	for( int i = 1; i < ssfacenum; i ++)
	{
		fprintf( fout, ",{{%d", ssface[ i ][ 0 ] + 1);
		for( int j = 1; j < ssfaceedgenum[ i ]; j ++)
		{
			fprintf( fout, ",%d", ssface[ i ][ j ] + 1);
		}
		faceside = ssface_planeindex[ i ];
		if( faceside >= 0 )
			faceside += 1;
		fprintf( fout,  "},{%d}}", faceside );
	}
	//subspace
	fprintf( fout, "},{{{%d", ssspace[ 0 ][ 0 ]+ 1);
	for( int i = 1; i < ssspacefacenum[ 0 ]; i ++)
	{
		fprintf( fout, ",%d", ssspace[0][ i ] + 1);
	}
	fprintf( fout, "},{%d", ssspace_planeside[ 0 ][ 0 ]);
	for( int i = 1; i < ssspacefacenum[ 0 ]; i ++)
	{
		fprintf( fout, ",%d", ssspace_planeside[ 0 ][ i ]);
	}
	fprintf( fout, "}}");
	for( int i = 1; i < ssspacenum; i ++)
	{
		fprintf( fout, ",{{%d", ssspace[ i ][ 0 ] + 1 );
		for( int j = 1; j < ssspacefacenum[ i ]; j ++)
			fprintf( fout, ",%d", ssspace[ i ][ j ] + 1);
		fprintf( fout, "},{%d", ssspace_planeside[ i ][ 0 ]);
		for( int j = 1; j < ssspacefacenum[ i ]; j ++)
			fprintf(fout, ",%d", ssspace_planeside[ i ][ j ]);
		fprintf(fout, "}}");
	}
	fprintf( fout, "}}");
	fclose( fout );
}

SpacePartitioner::SpacePartitioner()
{

}

SpacePartitioner::~SpacePartitioner()
{

}

