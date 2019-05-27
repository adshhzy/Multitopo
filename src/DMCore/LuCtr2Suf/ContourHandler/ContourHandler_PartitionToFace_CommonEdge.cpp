#include "../ContourHandler/ContourHandler.h"

//for each plane ,find the two faces on it
void ContourHandler::gatherFacesOnPlaneCMNL(int planenum,  int*& planeFaces, int ssfacenum, int* ssface_planeindex )
{
	for( int i = 0; i < planenum * 2; i ++ )
	{
		planeFaces[ i ] = -1;
	}
	for( int i = 0 ;i < ssfacenum; i ++)
	{
		int planei = ssface_planeindex[ i ];
		if ( planei < 0)		
			continue;
		if( planeFaces[2*planei] == -1)
		{
			planeFaces[ 2*planei ] = i;	
		}
		else
		{
			planeFaces[ 2*planei + 1 ] = i;			
		}
	}	
}

//divide the contour vertices on each plane,
//mark the vertices type and values on it
//set the map from the original index of the vertex to the index of the vertex on the face
void ContourHandler::markCtrVerCMNL
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
 )
{
	ctrvertype = new int [ pctrvernum ];
	ctrverval = new int[ pctrvernum ];

	//decide vermark = 1 for on the first face, or the second face of the plane
	bool isPos = true;	// find a contour vertex with mark 1
	int verind = -1;
	for( int i = 0; i < pctrvernum; i ++ )
	{
		if( ctrvermark[ i ] == 0)
			continue;
		if( ctrvermark [ i ] == 1 )
		{
			isPos = true;
			verind = i;
			break;
		}

		isPos = false;
		verind = i;
	}

	//direction, crossproduct( comnpt_ctrver, comndiir )
	//direction, crossproduct( comnpt_one ssver on the first face, comndir)
	//find one vertex on the first face
	int face1 = planefaces[ 0 ];
	int ssverind = -1;
	for( int i = 0; i < ssfaceedgenum[ face1 ]; i ++ )
	{
		int edgei = ssface[ face1 ][ i ];
		if( (ssedge[ 2*edgei ] != comnveri[ 0 ] ) && (ssedge[ 2*edgei ] != comnveri[ 1 ]) )
		{
			ssverind = ssedge[ 2 * edgei  ];
			break;
		}
		if( (ssedge[ 2*edgei +1 ] != comnveri[ 0 ] ) && (ssedge[ 2*edgei + 1] != comnveri[ 1 ]) )
		{
			ssverind = ssedge[ 2 * edgei  + 1];
			break;
		}
	}
	//////////////////////////////////////////////////////////////////////////
	/*if( ssverind == -1 )
	{
		cout<<"Unable to find a vertex of the processing face which doesn't lie on the common edge!"<<endl;
	}*/
	//////////////////////////////////////////////////////////////////////////
	float dir1[ 3 ] , dir2[ 3 ];
	float vec [ 3 ];
	for( int i = 0; i < 3; i ++ )
	{
		vec[ i ] = pctrverpos[ 3*verind + i  ] - comnpt[ i ];
	}
	MyMath::crossProduct(  vec , comndir , dir1);
	for( int i = 0; i < 3; i ++ )
	{
		vec[ i  ] = ssver[ 3*ssverind + i  ] - comnpt[ i ];
	}
	MyMath::crossProduct(  vec, comndir , dir2);

	bool POS_FACE1 = true;	//mark = 1 is on the first face
	float val = MyMath::dotProduct( dir1, dir2 );
	if( (val > 0 && isPos)  || ( val < 0 && !isPos ))
	{
		POS_FACE1 = true;
	}
	else
		POS_FACE1 = false;

	//go through each vertex and mark it
	if( !POS_FACE1 )
	{
		planefaces[ 0 ] = planefaces[ 1 ];
		planefaces[ 1 ] = face1;
	}
	//////////////////////////////////////////////////////////////////////////
	//cout<<"before looping, check the mark values:"<<endl;
	//for( int i = 0;i < pctrvernum;i ++ )
	//{
	//	cout<<ctrvermark[ i ]<<" ";
	//}
	//cout<<endl;
	//////////////////////////////////////////////////////////////////////////
	int faceversize[ 2 ] = {0 , 0 };

	for( int i = 0; i < pctrvernum; i ++ )
	{		
		switch ( ctrvermark[ i ])
		{
		case 1:	//pos
			{
				ctrvertype[ i ] = CTRTYPE_SUBFACE;
				ctrverval[ i  ]  = planefaces[ 0 ];

				for( int j = 0; j < 3; j ++   )
					fctrverpos[ planefaces[ 0 ]].push_back( pctrverpos[ 3*i + j ]);
				fctrvertype[ planefaces[ 0 ]].push_back( CTRTYPE_SUBFACE );
				fctrverval[ planefaces[ 0 ]].push_back( planefaces[ 0 ]);				

				vermap[ 0 ][ i ] = faceversize[ 0 ];
				faceversize[ 0 ] ++;
			}
			break;
		case -1:
			{
				ctrvertype[ i ] = CTRTYPE_SUBFACE;
				ctrverval[ i  ]  =  planefaces[ 1 ];

				for( int j = 0; j < 3; j ++   )
					fctrverpos[ planefaces[ 1 ]].push_back( pctrverpos[ 3*i + j ]);
				fctrvertype[ planefaces[ 1 ]].push_back( CTRTYPE_SUBFACE );
				fctrverval[ planefaces[ 1 ]].push_back( planefaces[ 1  ]);

				vermap[ 1 ][ i ] = faceversize[ 1 ];
				faceversize[ 1 ] ++;
			}
			break;
		case 0:
			{
				ctrvertype[ i ] = CTRTYPE_SUBEDGE;
				ctrverval[ i  ]  =  comnedgei;

				for( int j = 0; j < 3; j ++   )
					fctrverpos[ planefaces[ 0 ]].push_back( pctrverpos[ 3*i + j ]);
				fctrvertype[ planefaces[ 0 ]].push_back( CTRTYPE_SUBEDGE);
				fctrverval[ planefaces[ 0 ]].push_back ( comnedgei );

				for( int j = 0; j < 3; j ++   )
					fctrverpos[ planefaces[ 1 ]].push_back( pctrverpos[ 3*i + j ]);
				fctrvertype[ planefaces[ 1 ]].push_back( CTRTYPE_SUBEDGE );
				fctrverval[ planefaces[ 1 ]].push_back( comnedgei );

				vermap[ 0 ][ i ] = faceversize[ 0 ];
				faceversize[ 0 ] ++;

				vermap[ 1 ][ i ] = faceversize[ 1 ];
				faceversize[ 1 ] ++;
			}
			break;
		}
	}
}

void ContourHandler::putControuEdgeIntoFaceCMNL
(
	int planefaces[ 2 ],	
	int* vertype, int* verval,
	int pctredgenum, int* pctredges,	//four is a group

    vector<intvector>& fctredges,
	vector<intvector>& fctreedgetype,
	vector<intvector>& fctredgeval,
	vector<intvector>& fctredgeancestor,

	int* vermap[ 2 ]
 )
{	
	//go through each edge, 
	//check the types of the two verties, and put it into the right face
	for( int i = 0; i < pctredgenum; i ++ )
	{
		int vers[ 2 ];
		memcpy( vers, pctredges + 4* i , sizeof( int ) * 2);	

		int twovertype[ 2 ];
		int twoverval[ 2 ];
		for( int j = 0; j < 2; j ++ )
		{
			twovertype[ j  ] = vertype[ vers[ j ]];
			twoverval[ j ] = verval[ vers [  j] ];
		}

        int whichface;
		if( twovertype[ 0 ] == CTRTYPE_SUBEDGE )
		{
			whichface = twoverval[ 1 ];
		}
		else
			whichface = twoverval[ 0 ];

		//put the edge into the right face
		//int fctredgessize = fctredges[ whichface ].size();

        //fctredges[ whichface ].resize( fctredgessize  + 4 );
		//int ind = 0;
		//if( whichface == planefaces[ 1 ])
		//	ind = 1;
		//fctredges[ whichface ][ fctredgessize ] = vermap[ ind ][ pctredges[ 4*i ] ];
		//fctredges[ whichface ][ fctredgessize + 1] = vermap[ ind ][ pctredges[4*i + 1]];
		//memcpy( &fctredges[ whichface][ fctredgessize + 2],  &pctredges[ 4*i + 2] , sizeof( int ) * 4 );

		int ind = 0;
		if( whichface == planefaces[ 1 ])
			ind = 1;
		fctredges[ whichface ].push_back( vermap[ ind ][ pctredges[ 4*i ] ] );
		fctredges[ whichface ] .push_back( vermap[ ind ][ pctredges[4*i + 1]] );
		fctredges[ whichface].push_back( pctredges[ 4*i + 2] );
		fctredges[ whichface].push_back( pctredges[ 4*i +  3] );


		fctreedgetype[ whichface ].push_back( CTRTYPE_SUBFACE );
		fctredgeval[ whichface ].push_back ( whichface );
		fctredgeancestor[ whichface ].push_back( i );
	}		
}

//partition the contours into the faces
//set the contour vertex type
//set the contour edge type
void ContourHandler::putContourIntoFaceCMNL
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
 )
{	
	vector<floatvector> fctrver;
	vector<intvector> fctrvertype;
	vector<intvector> fctrverval;
	vector<intvector> fctredge;
	vector<intvector> fctredgetype;
	vector<intvector> fctredgeval;
	vector<intvector> fctredgeancestor;
	
	fctrver.resize( ssfacenum );
	fctrvertype.resize( ssfacenum) ;
	fctrverval.resize( ssfacenum );
	fctredge.resize( ssfacenum );
	fctredgetype.resize( ssfacenum );
	fctredgeval.resize( ssfacenum );
	fctredgeancestor.resize( ssfacenum );

	int* planefaces = new int[ 2* planenum ];
	//find the two faces for each plane
	gatherFacesOnPlaneCMNL(planenum,  planefaces,  ssfacenum,  ssface_planeindex );	

	//mark the contour vertices
	int** ctrvertype;
	int** ctrverval;
	ctrvertype = new int*[ planenum ];
	ctrverval = new int*[ planenum ];

	for( int i = 0 ; i < planenum; i++ )
	{
		int* vermap[ 2 ];

		vermap[ 0 ] = new int[ pctrvernum[ i ] ];
		vermap[ 1 ] = new int[ pctrvernum[ i ] ];

		markCtrVerCMNL(
		 planefaces + 2*i, ssver,  ssedge, ssfaceedgenum,  ssface, 
		comndir, comnpt, comnedgei, comnveri,		
		pctrvernum[ i ], pctrverpos[ i ], pctrvermark[ i ], ctrvertype[ i ],ctrverval[ i ], fctrver,
		fctrvertype,		fctrverval, vermap);		

		//////////////////////////////////////////////////////////////////////////
	/*	for( int j = 0; j < pctrvernum[ i ]; j ++ )
		{
			cout<<vermap[ 0 ][ j ]<<" ";
		}
		cout<<endl;
		for( int j = 0; j < pctrvernum[ i ]; j ++ )
		{
			cout<<vermap[ 1 ][ j ]<<" ";
		}
		cout<<endl;*/
		//////////////////////////////////////////////////////////////////////////

		//mark the contour edges
		putControuEdgeIntoFaceCMNL( planefaces + 2*i,
			ctrvertype[ i  ],  ctrverval[ i ], pctredgenum[ i ],  pctredge[ i  ],	
			fctredge,fctredgetype, fctredgeval, fctredgeancestor, vermap );

		delete []vermap[ 0 ];
		delete []vermap[ 1 ];
	}		

		
	//copy the result out
	for( int i = 0; i < ssfacenum; i ++ )
	{
		int fversize = fctrver[ i ].size();
		ctrfvernum[ i ] = fversize/3; 
		ctrfverpos[ i ] = new float[ fversize ];		
		memcpy(ctrfverpos[ i ], &fctrver[ i ][ 0 ], sizeof( float ) * fversize );

		fversize /= 3;
        ctrfvertype[ i ] = new int[ fversize ];
		memcpy( ctrfvertype[ i ], &fctrvertype[ i ][ 0 ] , sizeof( int )* fversize );

		ctrfverval[ i  ]= new int[ fversize ];
		memcpy( ctrfverval[ i ], &fctrverval[ i ][0] , sizeof( int )* fversize );

		int fedgesize = fctredge[ i ].size();
		ctrfedgenum[ i ] = fedgesize/4;
		ctrfedge[ i ] = new int[ fedgesize ];
		memcpy( ctrfedge[ i ], &fctredge[ i ][ 0 ], sizeof(int)*fedgesize );

		fedgesize /= 4;
		ctrfedgetype[ i ] = new int[ fedgesize ];
		memcpy( ctrfedgetype[ i ] , &fctredgetype[ i ][ 0 ], sizeof(int ) * fedgesize );

		ctrfedgeval [ i ] = new int[ fedgesize ];	
		memcpy( ctrfedgeval[ i ], &fctredgeval[ i ][ 0 ], sizeof( int )*fedgesize );

		ctrfedgeancestor[ i ] = new int[ fedgesize ];
		memcpy( ctrfedgeancestor[ i ], &fctredgeancestor[ i ][ 0], sizeof(int) * fedgesize);
	}
	//////////////////////////////////////////////////////////////////////////
	//for( int i = 0; i < planenum * 2; i ++)
	//{
	//	cout<<planefaces[ i ]<<"  ";
	//	if( i % 2 == 0 )
	//		cout<<", ";
	//}
	//cout<<endl;
	//cout<<"vernum, edgenum for fctredge, fctredgetype, fctreedgeval"<<endl;
	//for( int i = 0; i < ssfacenum; i ++ )
	//{
	//		cout<<fctrver[ i ].size()/3<< "  "<<fctredge[ i ].size()/4<<"  "
	//			<<fctredgetype[ i ].size()<<"  "<<fctredgeval[ i ].size()<<endl;
	//}
	//cout<<endl;
	//////////////////////////////////////////////////////////////////////////

	//release the temp vars
	for( int i = 0 ; i < planenum ;i ++ )
	{
		fctrver[ i ].clear();
		fctrvertype[ i ].clear();
		fctrverval[ i ].clear();

		fctredge[ i ].clear();
		fctredgetype[ i ].clear();
		fctredgeval[ i ].clear();
		fctredgeancestor[ i ].clear();
	}
	fctrver.clear();
	fctrvertype.clear();
	fctrverval.clear();
	fctredge.clear();
	fctredgetype.clear();
	fctredgeval.clear();
	fctredgeancestor.clear();

    delete	[]planefaces;
	for( int i = 0; i < planenum; i ++ )
	{
		delete []ctrverval[ i ];
		delete []ctrvertype[ i ];
	}
	delete []ctrverval;
	delete []ctrvertype;

	//////////////////////////////////////////////////////////////////////////
	/*for( int i = 0; i < ssfacenum; i ++ )
	{
		if( ctrfvernum[ i ] != 0)
		writeFaceCtr(i, ctrfvernum[ i ], ctrfverpos[ i ],ctrfvertype[ i ], ctrfverval[ i ],
		ctrfedgenum[i ], ctrfedge[ i ] , ctrfedgetype[ i ], ctrfedgeval[ i ]);	
	}*/
	//////////////////////////////////////////////////////////////////////////
}
