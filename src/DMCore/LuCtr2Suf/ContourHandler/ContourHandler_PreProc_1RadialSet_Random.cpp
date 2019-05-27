/************************************************************************/
/* 

//.contour2 file format, this special case is only for single material object. 
<planenum(i.e. contour_number)><contour>*
<contour> := <#vertex><vertex>*<#edges><edge>*<index of the points on plane_i>*

<#vertex> : number of vertices
<vertex>*:  positions of the vertices, one line for one vertex, 3 numbers per line, for x, y and z coordiantes
<#edge>:	number of edges
<edge>*:	indices of the two ending vertices of the edge in the vertex list, index starts at 0.
			indices of the materials on left and right side of the edge, 4 numberes per line.
<index of the points on plane_i>*:	i goes from 0 to planenum-1
			2 numbers per line, indices of the two vertices in the contour that are on planei- only for 
				the case that two contours share two intersection points.
			use -1 if there are no intersections. if i equals to the index of the contour itself, both are -1.

read in the contours and sort them.

Sort the planes:
1 set of radial contours first,
followed by a set of random oriented planes
*/
/************************************************************************/

#include "../ContourHandler/ContourHandler.h"
#include "../math/mymath.h"

//return if the contour share common points with multiple planes
//if so, the result is in the commonPlaneList
void ContourHandler::readOneContour2_v2(char* filename, 
						floatvector& paramOut,	
						vector<floatvector>& ctrversOut,vector<intvector>& ctredgesOut,
						float bbox[ 6 ], bool& bboxset,
						vector<intvector>& 	ver2planelistposOut,
						vector<intvector>& ver2planelist,
						int& radialPlaneNum
											   )
{
	FILE* fin = fopen( filename, "r");

	if( fin == NULL )
	{
		cout<<"Unable to open file:" <<filename<<endl;
		return;
	}

	floatvector param;
	vector<floatvector> ctrvers;
	vector<intvector> ctredges;
		
	vector<intvector> 	ver2planelistpos;
	//vector<intvector> ver2planelist;

	//contour number
	int ctrnum = 0;
	fscanf( fin, "%d", &ctrnum );

	int posInver2planelistpos = ver2planelistpos.size();
	ver2planelistpos.resize( posInver2planelistpos + ctrnum );
	int ver2planelistSize = ver2planelist.size();

	//////////////////////////////////////////////////////////////////////////
	//	cout<<"ctrnum:"<<ctrnum<<endl;
	//////////////////////////////////////////////////////////////////////////

	//read contour one by one
	float tparam[ 4 ];
	float tver[ 3 ];
	int tedge[ 4 ];
	int vernum, edgenum;
	floatvector tvers;
	intvector tedges;

	//this saves all the intersection points index between planei and planej
	//vector<int> cmnPtPlaneList;
	set<int> cmnPtPlaneList;
	int** cmnPtInd;
	cmnPtInd = new int*[ ctrnum ];
	for( int j = 0; j < ctrnum; j ++ )
	{
		cmnPtInd[ j ] = new int[ 3 ];		
	}
	bool found = false;


	for( int i = 0; i < ctrnum; i ++)
	{
		//plane parameter
		fscanf( fin, "%f %f %f %f\r\n", tparam, tparam + 1, tparam+2, tparam+3);
		//////////////////////////////////////////////////////////////////////////
		//cout<<"parma:"<<tparam[ 0 ]<<" "<<tparam[ 1 ]<<" "<<tparam[ 2 ]<<" "<<tparam[ 3 ]<<endl;
		//////////////////////////////////////////////////////////////////////////
		param.push_back( tparam[ 0 ]);
		param.push_back( tparam[ 1 ]);
		param.push_back( tparam[ 2 ]);
		param.push_back( tparam[ 3 ]);

		//vertex number, edge number
		fscanf( fin, "%d %d\r\n", &vernum,&edgenum);

		//////////////////////////////////////////////////////////////////////////
		//	cout<<"vernum"<<vernum<<"edgenum"<<edgenum<<endl;
		//////////////////////////////////////////////////////////////////////////
		//vertices
		for( int j = 0; j < vernum; j ++)
		{
			fscanf( fin, "%f %f %f\r\n", tver, tver+1, tver+2);
			tvers.push_back( tver[ 0 ]);
			tvers.push_back( tver[ 1 ]);
			tvers.push_back( tver[ 2 ]);

			if( !bboxset )
			{
				bbox[ 0 ] = bbox[ 3 ] = tver[ 0 ];
				bbox[ 1 ] = bbox[ 4 ] = tver[ 1 ];
				bbox[ 2 ] = bbox[ 5 ] = tver[ 2 ];
				bboxset = true;
			}
			else
			{
				if( bbox[ 0 ] > tver[ 0 ])
				{
					bbox[ 0 ] = tver[ 0 ];
				}
				else if( bbox[ 3 ] < tver[ 0 ])
				{
					bbox[ 3 ] = tver[ 0 ];
				}
				if( bbox[ 1 ] > tver[ 1 ])
				{
					bbox[ 1 ] = tver[ 1 ];
				}
				else if( bbox[ 4 ] < tver[ 1 ])
				{
					bbox[ 4 ] = tver[ 1 ];
				}
				if( bbox[ 2 ] > tver[ 2 ])
				{
					bbox[ 2 ] = tver[ 2 ];
				}
				else if( bbox[ 5 ] < tver[ 2 ])
				{
					bbox[ 5 ] = tver[ 2 ];
				}				
			}
			//////////////////////////////////////////////////////////////////////////
			//	cout<<tver[ 0 ]<<" "<<tver[ 1 ]<<" "<<tver[ 2 ]<<endl;
			//////////////////////////////////////////////////////////////////////////
		}
		ctrvers.push_back( tvers );

		//edges
		for( int j = 0; j < edgenum; j ++ )
		{
			fscanf( fin, "%d %d %d %d\r\n", tedge, tedge + 1, tedge + 2, tedge + 3 );
			tedges.push_back( tedge[ 0 ]);
			tedges.push_back( tedge[ 1 ]);
			tedges.push_back( tedge[ 2 ]);
			tedges.push_back( tedge[ 3 ]);

			//////////////////////////////////////////////////////////////////////////
			//	cout<<tedge[ 0 ]<<" "<<tedge[ 1]<<" "<<tedge[ 2 ] <<" "<<tedge[3]<<endl;
			//////////////////////////////////////////////////////////////////////////
		}
		ctredges.push_back( tedges );

		//ver pair between current planei, and plane_0-> plane_{n-1}
		int twovers[ 2 ];
		ver2planelistpos[ posInver2planelistpos ].resize( vernum, -1);

		//initialize the intersection point array
		if( !found )
		{
			for( int j = 0; j < ctrnum; j ++ )
			{
				for( int k = 0; k < 2; k ++ )
					cmnPtInd[ j ][ k ]  = -1; //will save the two intersection points indices
				cmnPtInd[ j ][ 2 ] = j;	//the plane index
			}
		}

		for( int j = 0; j < ctrnum; j ++)
		{	
			fscanf( fin, "%d %d\r\n", twovers, twovers + 1);

			//first one saves the smaller index
			if( !found )
			{
				if( twovers[ 0  ] < twovers[ 1 ] )
				{
					cmnPtInd[ j ][ 0 ] = twovers[ 0 ];
					cmnPtInd[ j ][ 1 ] = twovers[ 1 ];			
				}
				else
				{
					cmnPtInd[ j ][ 0 ] = twovers[ 1 ];
					cmnPtInd[ j ][ 1 ] = twovers[ 0 ];			
				}
			}
			

			if( j == i || twovers[ 0 ] == -1 )
				continue;


			for( int k = 0; k < 2; k ++ )
			{
				int posInVer2PlaneList =  ver2planelistpos[ posInver2planelistpos ][ twovers[ k ]] ;
				if( posInVer2PlaneList == -1 )//hasn't been set for vertex twovers[ k ]
				{
					ver2planelist.resize( ver2planelistSize + 1);
					posInVer2PlaneList = ver2planelistSize;
					ver2planelistpos[ posInver2planelistpos ][ twovers[ k ]] = posInVer2PlaneList;

					ver2planelistSize ++;
				}

				ver2planelist[ posInVer2PlaneList ].push_back( j );	//this vertex has corresponding plane, j plane
			}

		}

		//sort the intersection points based on the first index, if equal, second index
		if ( !found )
		{
			//////////////////////////////////////////////////////////////////////////
			/*cout<<"before:"<<endl;
			for( int i= 0; i < ctrnum; i ++ )
			{
				cout<<cmnPtInd[ i ][ 0 ]<<","<<cmnPtInd[ i ][ 1 ]<<","<<cmnPtInd[ i ][ 2]<<endl;
			}*/
			//////////////////////////////////////////////////////////////////////////

			for( int j = 0; j < ctrnum; j ++ )
			{
				int mini = j;		

				for( int k = j + 1; k < ctrnum; k ++ )
				{
					if( cmnPtInd[ k ][ 0 ] < cmnPtInd[ mini ][ 0 ])
					{
						mini = k;
					}
					else if( cmnPtInd[ k ][ 0 ] == cmnPtInd[ mini ][ 0 ])
					{
						if( cmnPtInd[ k ][ 1 ] < cmnPtInd[ mini ][ 1 ])
						{
							mini = k;
						}
					}
				}

				if( mini == j )	//no need to switch position
					continue;

				//switch the min with current position
				int* tptr = cmnPtInd[ j ];
				cmnPtInd[ j ] = cmnPtInd[ mini ];
				cmnPtInd[ mini ] = tptr;
			}
			//////////////////////////////////////////////////////////////////////////
			/*cout<<"after:"<<endl;
			for( int i= 0; i < ctrnum; i ++ )
			{
				cout<<cmnPtInd[ i ][ 0 ]<<","<<cmnPtInd[ i ][ 1 ]<<","<<cmnPtInd[ i ][ 2]<<endl;
			}*/
			//////////////////////////////////////////////////////////////////////////\

			//check if there are a pair of intersection points, that are shared by more than 3 planes
			int pos = 0;
			while(  pos < ctrnum && (cmnPtInd[ pos ][ 0 ] == -1))
			{
				pos ++;
			}

			while( pos < ctrnum - 1 && ( cmnPtInd[ pos ][ 0 ] != cmnPtInd[ pos + 1][ 0 ]))
				pos ++;
			
			if( pos != ctrnum - 1)	//before searching to the end, finding two pairs of interpt the same
			{
				found = true;

//				cmnPtPlaneList.push_back( i ); //the ith plane
				cmnPtPlaneList.insert( i );
				
				while( pos < ctrnum - 1&& ( cmnPtInd[ pos][ 0 ] == cmnPtInd[ pos + 1][ 0 ]) &&
					(cmnPtInd[ pos][ 1 ] == cmnPtInd[ pos ][ 1 ]))
				{
					cmnPtPlaneList.insert( cmnPtInd[ pos ][ 2]);
					pos ++;
				}
				
				if( pos == ctrnum - 1 ) //the second to the last one
				{
					cmnPtPlaneList.insert( cmnPtInd[ pos ][ 2 ]);
				}			
			}
			
		}		


		posInver2planelistpos++;

		tvers.clear();
		tedges.clear();
	}

	//delete the temp variables
	for( int i = 0; i < ctrnum; i ++ )
		delete []cmnPtInd[ i];
	delete []cmnPtInd;

	//////////////////////////////////////////////////////////////////////////
	//cout the sorted information
	//cout<<"the planes are:";
	//set<int>::iterator iter;
	//for( iter = cmnPtPlaneList.begin(); iter != cmnPtPlaneList.end(); iter++)
	////for( int i = 0; i < cmnPtPlaneList.size(); i ++ )
	//{
	//	cout<<*iter<<",";
	//}
	//cout<<endl;
	//////////////////////////////////////////////////////////////////////////

	radialPlaneNum = cmnPtPlaneList.size();

	//reorder the contours
	paramOut.resize( param.size());
	ctrversOut.resize( ctrvers.size() );
	ctredgesOut.resize( ctredges.size() );
	ver2planelistposOut.resize( ver2planelistpos.size() );
	
	int* planeMark = new int[ ctrnum ];
	for( int i = 0; i < ctrnum; i ++ )
		planeMark[ i ] = 0;

	int* planeIndList = new int[ ctrnum ];

	set<int>::iterator Iter;
	Iter = cmnPtPlaneList.begin();
	for( int i = 0; i < radialPlaneNum; i ++ )
	{
		planeMark[ *Iter ] = 1;
		planeIndList[ i ] = *Iter;

		Iter ++;
	}

	int pos = radialPlaneNum;
	int pos1 = 0;
	while( planeMark[ pos1 ] )
		pos1 ++;

	while( pos < ctrnum )
	{
		planeIndList[ pos ] = pos1;
		pos1 ++;
		while ( planeMark[ pos1 ] )
		{
			pos1 ++;
		}
		pos++;
	}

	//////////////////////////////////////////////////////////////////////////
	/*for( int i  = 0; i < ctrnum; i ++ )
	{
		cout<<planeIndList[ i ]<<",";
	}
	cout<<endl;*/
	//////////////////////////////////////////////////////////////////////////

	for( int i = 0; i < ctrnum; i ++ )
	{
		int planeInd = planeIndList[ i ];

		for( int j = 0; j < 4; j ++ )
			paramOut[ 4* i + j ] = param[ 4* planeInd + j ] ;
		
		int size = ctrvers[ planeInd ].size();
		ctrversOut[ i ].resize( size );
		for( int j = 0; j < size; j ++ )
		{
			ctrversOut[ i ][ j ] = ctrvers[ planeInd ][ j ];
		}
		
		size = ctredges[ planeInd ].size();
		ctredgesOut[ i ].resize( size );
		for( int j = 0; j < size; j++ )
		{
			ctredgesOut[ i ][ j ] = ctredges[ planeInd ][ j ];
		}

		size = ver2planelistpos[ planeInd ].size();
		ver2planelistposOut[ i ].resize( size );
		for( int j = 0; j < size; j ++ )
		{
			ver2planelistposOut[ i ][ j ] = ver2planelistpos[ planeInd ][ j ];
		}
	}	

	cmnPtPlaneList.clear();
	delete []planeMark;
	delete []planeIndList;

	

	//////////////////////////////////////////////////////////////////////////
	/*for( int i= 0; i < ver2planelist.size(); i++ )
	{
	for( int j = 0; j < ver2planelist[ i ].size(); j ++ )
	cout<<ver2planelist[ i ][ j ]<<",";
	cout<<endl;
	}*/
	//////////////////////////////////////////////////////////////////////////

	fclose( fin );
}
void ContourHandler::readContour2_v2(const int filenum, char** filenames, 
							floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
							float bbox[ 6 ], bool& bboxset,
							vector<intvector>& 	ver2planelistpos,
							vector<intvector>& ver2planelist, 
							int& radialPlaneNum)
{


	//read file one by one
	for( int i = 0; i < filenum; i ++)
	{
		readOneContour2_v2(filenames[ i ], param,ctrvers,ctredges,bbox,
			bboxset,ver2planelistpos, ver2planelist, radialPlaneNum	);

		////////////////
		/*for( int  j = 0; j < ver2planelistpos.size(); j ++ )
		{
			for( int k = 0; k < ver2planelistpos[ j ].size(); k ++ )
			{
				if( ver2planelistpos[ j ][ k ] != -1)
				{
					int pos= ver2planelistpos[  j ][ k ];
					for( int k2 = 0; k2 < ver2planelist[ pos ].size(); k2 ++ )
					{
						cout<<ver2planelist[ pos ][ k2 ]<<",";
					}
					cout<<endl;
				}
			}
		}*/
		/////////////////
	}

	//normalize the normals of the planes
	int planenum = param.size()/4;
	for( int i = 0; i < planenum; i ++)
	{
		float veclen = MyMath::vectorLen( param[ i * 4 ], param[ i*4 + 1 ], param[ i*4+ 2 ]);
		for(int j = 0; j < 4; j ++)
			param[ 4*i + j] = param[ 4*i + j]/veclen;
	}
}
