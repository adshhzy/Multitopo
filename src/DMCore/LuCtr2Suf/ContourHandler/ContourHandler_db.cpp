#include "../ContourHandler/ContourHandler.h"

void ContourHandler::writeOneContourOut
(floatvector & ctrvers,
 intvector& ctredges,
 int* ctrvermark,
 intvector interptind,
 int becut,
 int tocut,
 floatvector& param
 )
{
	char fname[ 1024 ];
	strcpy( fname, "mmdebug/becut_");
	char numstr[ 10 ];
	/*itoa( becut, numstr, 10 );*/
	sprintf(numstr, "%d", becut);
	strcat(fname, numstr);
	strcat( fname, "_tocut_");
	/*itoa( tocut, numstr, 10 );*/
	sprintf(numstr, "%d", tocut);
	strcat( fname, numstr);
	strcat( fname, ".txt");

    FILE* fout= fopen( fname, "w")	;
   	
	//write out the vertices
	int vernum = ctrvers.size()/3;
	fprintf(fout, "{{");
	for( int i = 0; i < vernum; i ++ )
	{
		fprintf(fout, "{%f,%f,%f}", ctrvers[ 3*i ], ctrvers[ 3* i + 1], ctrvers[ 3*i + 2]);
		if( i !=vernum - 1)
			fprintf(fout, ",");
		else
			fprintf(fout, "},{");
	}

	//write out the edges
	int edgenum = ctredges.size() /4;
	for( int i = 0; i < edgenum; i ++ )
	{
		fprintf(fout, "{%d,%d}", ctredges[ 4*i ], ctredges[ 4*i + 1]) ;
		if( i != edgenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},{");
	}

	//write out the marks
	for( int i = 0; i < vernum; i ++ )
	{
		fprintf(fout, "%d", ctrvermark[ i ]);
		if( i != vernum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},{");
	}

	//write out the intersection index
	int interptsize = interptind.size();
	for( int i = 0 ;i < interptsize; i++ )
	{
		fprintf(fout, "%d", interptind[ i ]);
		if( i != interptsize - 1 )
			fprintf(fout,",");
		else
			fprintf(fout ,"},{");
	}

	int planenum = param.size()/4;
	for( int i = 0; i < planenum; i ++ )
	{
		fprintf(fout, "{%f,%f,%f,%f}", param[ 4*i ], param[ 4*i+ 1],
			param[ 4*i + 2], param[ 4*i + 3]);
		if( i != planenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "}}");
	}
	fclose( fout  );
}

void ContourHandler::writeFaceCtr
(int facei, 
 int fctrvernum, float* fctrverpos,	int* fvertype, int* fverval,
 int fctredgenum,	int* fedge, int* fedgetype, int* fedgeval	 )
{
	char num[ 10 ];
	/*itoa(facei + 1, num, 10 );*/
	sprintf(num, "%d", facei + 1);

	char fname[ 50 ];
	strcpy(fname, "mmdebug/facectr/facectr");
	strcat(fname, num );
	strcat( fname, ".txt");

	FILE* fout =fopen( fname, "w");

	fprintf( fout, "{");

	//write the verticesout
	fprintf(fout, "{");
	for( int i = 0; i < fctrvernum; i ++ )
	{
		fprintf(fout, "{{%f,%f,%f},{%d,%d}}", 
			fctrverpos[ 3*i ] , fctrverpos[ 3*i + 1 ], fctrverpos[ 3*i + 2],
			fvertype[ i ], fverval[ i ]	);
		if( i != fctrvernum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},{");
	}

	//write the edge out
	for( int i = 0; i < fctredgenum; i ++ )
	{
		fprintf(fout, "{{%d,%d},{%d,%d},{%d,%d}}", 
			fedge[ 4*i ], fedge[ 4 * i +1  ], fedge[ 4*i + 2], fedge[ 4*i + 3],
			fedgetype[ i ] , fedgeval[ i ]	);
		if( i != fctredgenum - 1)
			fprintf(fout, ",");
		else
			fprintf(fout, "}");
	}
	fprintf(fout, "}");
	fclose( fout );
}

void ContourHandler::writePlaneCtr
(int planei, floatvector& ctrvers,intvector& vermark,intvector& ctredges   )
{
	char numstr[ 10 ] ;
	/*itoa( planei + 1, numstr, 10 );*/
	sprintf(numstr, "%d", planei + 1);
	char fname[ 80 ] ;
	strcpy( fname,  "mmdebug/plane");
	strcat( fname, numstr );
	strcat( fname, ".txt");

	FILE* fout = fopen( fname, "w");

	fprintf(fout, "{");

	//vertex
	fprintf(fout, "{");
	int vernum = ctrvers.size() / 3;
	for( int i = 0; i < vernum; i ++ )
	{
		fprintf(fout, "{%f,%f,%f}", ctrvers[ i*3 ], 
			ctrvers[ 3*i + 1], ctrvers[ 3* i + 2 ]);
		if( i != vernum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},{");
	}

	//edge
	int edgenum = ctredges.size()/4;
	for( int i = 0; i < edgenum; i ++ )
	{
		fprintf(fout, "{{%d,%d},{%d,%d}}",
			ctredges[ 4*i] + 1, ctredges[ 4*i + 1 ] + 1,
			ctredges[ 4*i + 2 ], ctredges[ 4*i + 3 ]);
		if( i != edgenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},{");
	}
	
	//vertex mark
	for( int i = 0; i < vernum; i ++)
	{	
		fprintf(fout, "%d", vermark[ i ]);
		if( i != vernum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "}}");
	}

	fclose( fout );
}