#include "../ContourHandler/ContourHandler.h"
#include "../math/mymath.h"

/**
* Function read the data from .contour file
* @param filename: the path of the file to read in
* @param param: the parameters of the planes
* @param ctrvers: all the vertices of the contours
* @param ctredges: all the edges of the contours
*/
void ContourHandler::readOneContour(char* filename, 
						   floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
						   float bbox[ 6 ], bool& bboxset)
{
	FILE* fin = fopen( filename, "r");

	if( fin == NULL )
	{
		cout<<"Unable to open file:" <<filename<<endl;
		return;
	}

	//contour number
	int ctrnum = 0;
	fscanf( fin, "%d", &ctrnum );

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

		tvers.clear();
		tedges.clear();
	}

	fclose( fin );
}
/**
* Function read the data from .contour file
* @param filenum: number of the files to read contour in
* @param filenames: the paths of the files to read in
* @param param: the parameters of the planes
* @param ctrvers: all the vertices of the contours
* @param ctredges: all the edges of the contours
*/
void ContourHandler::readContour(const int filenum, char** filenames, 
			floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
			float bbox[ 6 ], bool& bboxset)
{
	//read file one by one
	for( int i = 0; i < filenum; i ++)
	{
		readOneContour(filenames[ i ], param,ctrvers,ctredges,bbox,bboxset);
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

void ContourHandler::writeOneContourMM
(  FILE* fout, float* pctrvers, int pctrvernum,  int* pctredges,	int pctredgenum, float* param)
{
	fprintf( fout, "{");
	
	fprintf( fout, "{%f, %f, %f, %f},{", param[ 0 ], param [ 1 ], param[ 2 ], param[ 3 ]);
	
	for( int i = 0; i < pctrvernum; i ++ )
	{
		fprintf( fout, "{%f,%f,%f}", pctrvers[ 3 * i ], pctrvers[ 3*i + 1], pctrvers[ 3*i + 2 ]);
		if( i != pctrvernum - 1 )
			fprintf( fout, ",");
	}
    
	fprintf(fout, "},{")	;

	for( int i = 0; i < pctredgenum; i++)
	{
		fprintf( fout, "{{%d,%d},{%d,%d}}", pctredges[ 4 * i ], pctredges[ 4 * i + 1],
			pctredges[ 4 * i + 2 ], pctredges[ 4 * i + 3 ]);

		if ( i  != pctredgenum - 1 )
			fprintf( fout, ",");
	}

	fprintf( fout, "}}");	
}

void ContourHandler::writeContourMM(
						   int planenum,
						   float** pctrvers ,
						   int* pctrvernum,
						   int** pctredges,
						   int* pctredgenum,
						   float* pparam,
						   char* fname 
						   )
{
	FILE* fout = fopen( fname, "wb");
	if( fout == NULL )
	{
		cout<<"unable to open file "<<fname<<" to write!"<<endl;
		return;
	}

	fprintf( fout, "{");
	for( int i = 0; i < planenum; i ++ )
	{
		writeOneContourMM(fout, pctrvers[ i ], pctrvernum[ i ], pctredges[ i ],
			pctredgenum[ i ], pparam + 4 * i );

		if( i != planenum - 1)
			fprintf( fout, ",");
		else
			fprintf( fout, "}");
	}
	
	fclose( fout );
}

void ContourHandler::readOneContour2(char* filename, 
									floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
									float bbox[ 6 ], bool& bboxset,
									vector<intvector>& 	ver2planelistpos,
									vector<intvector>& ver2planelist
									)
{
	FILE* fin = fopen( filename, "r");

	if( fin == NULL )
	{
		cout<<"Unable to open file:" <<filename<<endl;
		return;
	}

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
		
		for( int j = 0; j < ctrnum; j ++)
		{	
			fscanf( fin, "%d %d\r\n", twovers, twovers + 1);

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

		posInver2planelistpos++;

		tvers.clear();
		tedges.clear();
	}

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

void ContourHandler::readContour2(const int filenum, char** filenames, 
						 floatvector& param,	vector<floatvector>& ctrvers,vector<intvector>& ctredges,
						 float bbox[ 6 ], bool& bboxset,
						 vector<intvector>& 	ver2planelistpos,
						 vector<intvector>& ver2planelist
						 )
{
	//read file one by one
	for( int i = 0; i < filenum; i ++)
	{
		readOneContour2(filenames[ i ], param,ctrvers,ctredges,bbox,bboxset,ver2planelistpos, ver2planelist);
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
