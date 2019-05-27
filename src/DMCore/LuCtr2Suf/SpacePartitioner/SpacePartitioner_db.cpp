#include "../SpacePartitioner/SpacePartitioner.h"

void SpacePartitioner::writemapOneSubspaceCMNL(int spacei, 
											   intvector& ssedge,  int comnedge, 
											   vector<intvector>& ssface,  vector<intvector>& ssspace,
											   vector<intvector>& nssface, 		intvector& oedgelist,	
											   intvector& nssedge, intvector& overlist,	
											   int ncomnedge)
{
	FILE* fout = fopen("mmdebug/mapOneSubspaceCMNL.txt" ,"w");

	fprintf(fout, "{");

	//print the old faces
	int facenum = ssspace[ spacei ].size();
	fprintf(fout, "{");
	for( int i = 0 ;i < facenum; i ++ )
	{
		int facei = ssspace[ spacei ][ i ];
		fprintf(fout, "{");
		int faceedgenum = ssface[ facei  ] .size();
		for( int j = 0; j < faceedgenum; j ++ )
		{
			int edgei = ssface[ facei ][ j ];
			fprintf(fout, "%d", edgei);
			if( j != faceedgenum - 1 )
				fprintf(fout, ",");
			else
				fprintf(fout, "}");
		}
		if( i != facenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},");	
	}

	//print the new faces
	fprintf(fout, "{");
	for( int i = 0 ;i < facenum; i ++ )
	{
		int faceedgenum = nssface[ i ].size();
		fprintf(fout, "{");
		for( int j = 0 ; j < faceedgenum; j ++ )
		{
			fprintf(fout, "%d", nssface[ i ][ j ]);
			if( j != faceedgenum - 1 )
				fprintf(fout, ",");
			else
				fprintf(fout, "}");
		}
		if( i != facenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},");
	}

	//print the old edge list
	fprintf(fout, "{");
	int oedgenum = oedgelist.size();
	for( int i = 0; i < oedgenum; i ++ )
	{
		fprintf(fout, "%d", oedgelist[ i ]);
		if( i != oedgenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},");
	}

	//print the old edges
	fprintf(fout, "{");
	for( int i = 0; i < oedgenum; i ++ )
	{
		int edgei = oedgelist[ i ];
		fprintf(fout, "{%d,%d}", ssedge[ 2* edgei ], ssedge[ 2*edgei + 1 ]);
		if( i != oedgenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},");
	}

	//print the new edges
	fprintf(fout, "{");
	for( int i = 0; i <oedgenum; i ++ )
	{
		fprintf(fout, "{%d,%d}", nssedge[ i * 2 ], nssedge[ 2* i + 1]);
		if( i != oedgenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},");
	}

	//print the vertices
	fprintf(fout, "{");
	int oversize = overlist.size();
	for(int i = 0; i < oversize; i ++ )
	{
		fprintf(fout, "%d", overlist[  i ]);
		if( i != oversize - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "}}");
	}

	fclose( fout );
}

void SpacePartitioner::writeOVerOEdgeList(intvector& overlist, intvector& oedgeslist)
{
	FILE* fout = fopen("mmdebug/Overedge.txt", "w");

	fprintf(fout, "{");
	fprintf(fout, "{");
	for(unsigned int i = 0 ;i < overlist.size(); i ++ )
	{
		fprintf(fout, "%d", overlist[ i ]);
		if( i != overlist.size() -1)
			fprintf(fout, ",");
		else
			fprintf(fout, "},{");
	}

	for (unsigned int i = 0; i < oedgeslist.size(); i++)
	{
		fprintf( fout, "%d", oedgeslist[ i ]);
		if( i != oedgeslist.size() - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "}");
	}
	fprintf(fout, "}");

	fclose( fout );
}

void SpacePartitioner::writeOneIntvector(intvector& vec,  FILE* fout )
{
	fprintf(fout, "{");
	for (unsigned int i = 0; i < vec.size(); i++)
	{
		fprintf(fout, "%d", vec[ i ]);
		if( i != vec.size() - 1)
			fprintf(fout, ",");
		else
			fprintf(fout, "}");
	}
	fclose( fout );
}