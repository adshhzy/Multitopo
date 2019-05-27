#include "../ContourHandler/ContourHandler.h"
#include "../config.h"

//p must not be in vec, but in st. others in st must exist in vec.
//only preserve one in st in vec
void preserveOneFromSetNotP(vector<int>& vec, intset& st, int p)
{
	//////////////////////////////////////////////////////////////////////////
	cout<<"entering func to remove set:"<<endl;
	cout<<"vector:";
	for (unsigned int i = 0; i < vec.size(); i++)
		cout<<vec[i ]<<" ";
	cout<<endl;
	intset::iterator titer = st.begin();
	while( titer != st.end()	)
	{
		cout<<*titer<<" ";
		titer ++;
	}
	cout<<endl;
	//////////////////////////////////////////////////////////////////////////

	intset::iterator iter  = st.begin();
	if( *iter == p)	//preserve the second one, ignore the first one, get rid of things in st from third on from vec.
	{
        iter ++;
		iter ++;
	}
	else	//preserve the first one, get rid of things in st from the second on from vec.
		iter ++;
	intset::iterator end = st.end(); 

	//get rid of all the ones in st from vec	
	vector<int> temp( vec.size() );
	for (unsigned int i = 0; i < temp.size(); i++)
		temp[ i ] = vec[ i ];    
	vec.clear();

	int iter1 = 0;
	int end1 = temp.size();	
	while( iter1 < end1 && iter != end )
	{
		if( temp[iter1 ] < *iter )	//move the first pointer on, add the things that are different into vec
			while( temp[ iter1 ] < *iter && iter1 < end1)
			{
				vec.push_back( temp[iter1] );
				iter1 ++;
			}
		else if( temp[ iter1 ] > *iter )	//move the second pointer on, ignore all the things that are not in vec
		{
			while( temp[ iter1 ] > *iter && iter != end )
				iter ++;
		}

		if (iter1 < end1 && iter != end && (temp[ iter1 ] == *iter) )	//if the two are the same, move both pointers on, until they are different
		{
            while(iter1 < end1 && iter != end && (temp[ iter1 ] == *iter) )
			{
				iter1 ++;
				iter ++;
			}
		}
	}
	//if iter1 hasn't reached the end, but iter has, push back all the ones after iter1
	while( iter1 < end1 )
	{
		vec.push_back( temp[ iter1 ]);
		iter1 ++;
	}

	//////////////////////////////////////////////////////////////////////////
	cout<<"after: ";
	cout<<"vector:";
	for (unsigned int i = 0; i < vec.size(); i++)
		cout<<vec[i ]<<" ";
	cout<<endl;
	//////////////////////////////////////////////////////////////////////////
}

//find out the planes that have relationship with these elements
//if the plane index is negative, then change it to be positive!
inline void ContourHandler::getElement2Planes(
	int planenum,
	int ssvernum, int ssedgenum, int* ssedge,
	int ssfacenum, int* ssfaceedgenum, int** ssface,int* ssface_planeindex,
	int ssspacenum, int* ssspacefacenum, int** ssspace,
	vector<intset>& ver2planes,	vector<intset>& edge2planes,vector<intset>& plane2planes,
	vector<intvector>& planefaces,vector<intvector>& planeedges,vector<intvector>& planevers,	
	int radialPlaneNum)
{
	ver2planes.resize( ssvernum );
	edge2planes.resize( ssedgenum );
	plane2planes.resize( planenum );
	planefaces.resize( planenum );	//not include the added 6 bounding planes
	planevers.resize( planenum );
	planeedges.resize( planenum );

	//1. go through each face, set the planes for all the edges in it
	//and set the planefaces
	for( int facei = 0; facei < ssfacenum; facei ++)
	{
		int tedgenum = ssfaceedgenum[ facei ];
		for( int edgej = 0; edgej < tedgenum; edgej++ )
		{
			int planei = ssface_planeindex[ facei ];
			if( planei < 0 )
				planei = planenum + 6 + planei;
			edge2planes[ ssface[ facei ][ edgej] ].insert( planei );
		}
		if( ssface_planeindex[ facei ] < 0 )continue;
		planefaces[ ssface_planeindex[ facei ]].push_back( facei );
	}

	//2. go through each edge, set the planes for both vertices of it
	//for the planes the edge is in, set planeedges
	int ti = 0;
	for( int edgei = 0; edgei < ssedgenum; edgei++)
	{
		for( int veri = 0; veri < 2; veri++)
		{
			int tver = ssedge[ ti ++ ];
			merge( ver2planes[ tver ], edge2planes[ edgei ]);
		}
		
		intset::iterator iter = edge2planes[ edgei ].begin();
		while( iter!= edge2planes[ edgei ].end() )
		{
			if( *iter >= planenum )
			{
				iter++;
				continue;
			}
			planeedges[ *iter ].push_back( edgei );
			iter ++;
		}
	}

	//3. go through each vertex, set planevers
	for( int veri = 0; veri < ssvernum; veri++)
	{
		intset::iterator iter = ver2planes[ veri ].begin();
		while( iter != ver2planes[ veri ].end() )
		{
			if( *iter >= planenum) 
			{
				iter++;
				continue;
			}
			planevers[ *iter ].push_back( veri );
			iter ++;
		}	
	}

	//4. go through each plane, set planeplanes from the ver2planes of the vertices on it
	for( int planei = 0; planei < planenum; planei ++)
	{
		int vernum = planevers[ planei ].size();
		for( int i = 0; i < vernum; i ++)
		{
			int veri = planevers[ planei ][ i ];
			intset::iterator iter = ver2planes[ veri ].begin();
			while ( iter != ver2planes[ veri ].end() )
			{
				if(*iter == planei)
				{
					iter ++;
					continue;
				}
				plane2planes[ planei ].insert( *iter );
				iter++;
			}
		}
	}

	//for the radial set, keep one of them
	if( !_global_IsCtrGraph )	//old format, in 1 radial set case, get rid of some of the planes in the only 
		//radial set
	{
		for( int planei = 0; planei < radialPlaneNum; planei ++ )
		{
			//////////////////////////////////////////////////////////////////////////
			/*	intset::iterator iter2 = plane2planes[ planei ].begin();
			cout<<"plane "<<planei<<"befor:\n";
			while( iter2 != plane2planes[ planei ].end() )
			{
			cout<<*iter2<<",";
			iter2 ++;
			}
			cout<<endl;*/
			//////////////////////////////////////////////////////////////////////////

			intset::iterator iter = plane2planes[ planei ].begin();
			int curInd = *iter;

			if( curInd < radialPlaneNum )
			{
				iter ++;
				intset::iterator iterend = iter;
				while( *iterend < radialPlaneNum && ! (iterend == plane2planes[ planei ].end()))
				{
					iterend ++;
				}
				plane2planes[ planei ].erase( iter, iterend );			
			}

			//////////////////////////////////////////////////////////////////////////
			/*cout<<"after:\n";
			iter2 = plane2planes[ planei ].begin();
			while( iter2 != plane2planes[ planei ].end() )
			{
			cout<<*iter2<<",";
			iter2 ++;
			}
			cout<<endl;*/
			//////////////////////////////////////////////////////////////////////////
		}
	}
	else	//.CtrGraph case!!
	{		
		//////////////////////////////////////////////////////////////////////////
		cout<<"before removing radial sets!"<<endl;
		for (unsigned int i = 0; i < plane2planes.size(); i++)
		{
			intset::iterator iter = plane2planes[ i ].begin();
			while( iter != plane2planes[ i ].end() )
			{
				cout<<*iter<<" ";
				iter ++;
			}
			cout<<endl;
		}
		//////////////////////////////////////////////////////////////////////////

		//if currentplane is in a radial set, only one plane from the radial set should be saved in plane2planes.
		//all others should be remvoed
		//otherwise if currentplane is not in the radial set, there is no need to remove planes from the radial set
		for (unsigned int i = 0; i < plane2planes.size(); i++)
		{		
			vector<int> oldplanes( plane2planes[ i ].size() );			

			intset::iterator iter = plane2planes[ i ].begin();
			int pos = 0;
			while( iter != plane2planes[ i ].end())
			{
				oldplanes[ pos ] = *iter;				
				iter ++; 
				pos++;
			}
							
			//go through all the radialset
			for (unsigned int j = 0; j < _gloabl_radialPlaneSet.size(); j++)
			{
				if( _gloabl_radialPlaneSet[ j ].find( i ) == _gloabl_radialPlaneSet[ j ].end() )
				{
					continue;
				}

				//for the radial set, only preserve one from it other than the plane i itself
				preserveOneFromSetNotP( oldplanes, _gloabl_radialPlaneSet[ j ], i );
			}

			//replace the old set with the new one
			plane2planes[ i ].clear();
			for (unsigned int j = 0; j < oldplanes.size(); j++)
				plane2planes[ i ].insert( oldplanes[ j ]);
		}
	
	}
	//////////////////////////////////////////////////////////////////////////
	//cout<<"after removing radial sets!"<<endl;
	for (unsigned int i = 0; i < plane2planes.size(); i++)
	{
		intset::iterator iter = plane2planes[ i ].begin();
		while( iter != plane2planes[ i ].end() )
		{
			//cout<<*iter<<" ";
			iter ++;
		}
		//cout<<endl;
	}
	//////////////////////////////////////////////////////////////////////////
	
}
void inline ContourHandler::putContourIntoFace_OP_Mark(
	int**& markprop,
	int planei, 
	int* ver2planesnum,	int** ver2planes,	//sorted planes that the vertex is on
	int plane2planesnum, int* plane2planes,	//sorted planes that intersects with planei
	int planefacesnum, int* planefaces,
	int planeedgesnum, int* planeedges,
	int planevernum, int* planevers,
	int planenum,	int pctrvernum, float* pctrver,
	int pctredgenum, int* pctredges,float* pparam,float* ssver,int* ssedge, int* ssfaceedgenum, 
	int** ssface)
{
	//set mark for all the subspace vertices on this plane
	int rowi = 0;
	for( int i = 0; i < planevernum; i ++)
	{
		markprop[ rowi ][ plane2planesnum ] = CTRTYPE_SUBVER;
		//markprop[ rowi ][ plane2planesnum + 1 ] = planevers[ i ];
		markprop[ rowi ][ plane2planesnum + 1 ] = i;	//save the vertex index in planevers
		int i1 = 0;	//ver2planes
		int i2 = 0;	//plane2planes
		// int columj = 0;
		int curver = planevers[ i ];
		//set for the planes ver[ i ] is on
		while( (i2 < plane2planesnum) && (i1 < ver2planesnum[ curver ]))
		{
			if( plane2planes[ i2 ] == ver2planes[ curver ][ i1 ])
			{
				markprop[ rowi ][ i2 ] = 0;
				i1++;
				i2++;
			}
			else if( plane2planes[ i2 ] < ver2planes[ curver ][ i1 ])
				i2 ++;
			else
				i1++;
		}
		
		//set for all other planes
		for( int j = 0; j < plane2planesnum; j++)
		{
			if( markprop[ rowi ][ j ] == 0 )
				continue;
			//compute the relative position between vertex and the plane
			int tplanei = plane2planes[ j ];
			if( MyMath::dotProduct( ssver + 3*curver, pparam + 4* tplanei) > pparam[ 4*tplanei + 3 ])
			{
				markprop[ rowi ][ j ] = 1;
			}
			else
				markprop[ rowi ][ j ] = -1;
		}
		rowi ++;

		//////////////////////////////////////////////////////////////////////////
		/*if( curver == 30 || curver == 42)
		{
			for( int j = 0; j < plane2planesnum + 2; j++ )
			{
				cout<<markprop[ rowi -1][ j ]<<" ";
			}
			cout<<endl;
		}*/
		//////////////////////////////////////////////////////////////////////////
	}
	//set mark for edges on this plane
	for( int i = 0; i < planeedgesnum; i++)
	{
		int edgei = planeedges[ i ];
		int veri[ 2 ] = {ssedge[ 2*edgei ],  ssedge[ 2*edgei +  1 ]};
		int verrowi[ 2 ] = { -1, -1 };
		for( int j = 0; j < planevernum; j ++)
		{
			if(( verrowi[ 0 ] == -1) && ( veri[ 0 ] == planevers[markprop[ j ][ plane2planesnum + 1 ]]))
				verrowi[ 0 ] = j;
			else if (( verrowi[ 1 ] == -1) && ( veri[ 1 ] == planevers[markprop[ j ][ plane2planesnum + 1 ]]))
				verrowi[ 1 ] = j;

			if( (verrowi[ 0 ] != -1) && (verrowi[ 1 ] != -1))
				break;
		}
		int colum = 0;
		for(; colum < plane2planesnum; colum++)
		{
			int tmark = markprop[ verrowi[ 0 ]][ colum ] + markprop[ verrowi[ 1 ]][ colum ] ;
			if( tmark > 0) 
				tmark = 1;
			else if( tmark < 0 )
				tmark = -1;
			markprop[ rowi ][ colum ] = tmark;
		}
		markprop[ rowi ][ colum++ ] = CTRTYPE_SUBEDGE;
		markprop[ rowi ][ colum ++ ] = i;	//the edge index in planeedges
		//markprop[ rowi ][ colum++ ] = planeedges[ i ];
		rowi ++;
	}
	//set mark for the faces on this plane
	for( int i = 0; i < planefacesnum; i++)
	{
		int facei = planefaces[ i ];
		int tfacedgenum = ssfaceedgenum[ facei ];
		int* markrowlistt = new int[ tfacedgenum ];
		for( int j = 0; j < tfacedgenum; j++)
		{
			int tedgei = ssface[ facei ][ j ];
			for( int k = planevernum; k < planevernum + planeedgesnum; k++)
			{
				if( planeedges[markprop[ k ][ plane2planesnum + 1 ]] == tedgei )
				{
					markrowlistt[ j ] = k;
					break;
				}
			}
		}
		for( int j = 0; j < plane2planesnum; j++)
		{
			int tmark = 0;
			for( int k = 0; k < tfacedgenum; k++ )
			{
				tmark += markprop[ markrowlistt[ k ]][ j ];
			}
			if( tmark > 0 ) tmark = 1;
			else if( tmark < 0 ) tmark = -1;
			markprop[ rowi ][ j ] = tmark;
		}
		markprop[ rowi ][ plane2planesnum ] = CTRTYPE_SUBFACE;
		markprop[ rowi ][ plane2planesnum+1 ] = i;
		//markprop[ rowi ][ plane2planesnum+1 ] = facei;
		delete []markrowlistt;
		rowi ++;
	}	
}
int ContourHandler::markCompare(int* mark1,int* mark2,int columnnum	)
{
	for( int i = 0; i < columnnum; i++)
	{
		if( mark1[ i ] < mark2[ i ])	//mark1 is smaller
			return -1;
		if( mark1[ i ] > mark2[ i ])	//mark2 is smaller
			return 1;
	}
	return 0;	//they are equal, or else, already return
}

void inline ContourHandler::sortMark( 
int**& markprop, int rownum, int columnum
									)
{
	for( int i = 0; i < rownum; i ++)
	{
		//find the smallest one
		int smallest = i;
		for( int j = i+1; j < rownum; j++)
		{
			if( markCompare(markprop[ j ], markprop[smallest], columnum) ==-1 )
			{
				smallest = j ;
			}
		}
		if( smallest == i) continue;
		int* t = markprop[ i ] ;
		markprop[ i ] = markprop[ smallest ];
		markprop[ smallest ] = t;
		t = NULL;
	}
}
int inline ContourHandler::findPos(int** markprop, int rownum, int columnnum, int* vermark)
{
	int begin = 0;
	int end = rownum - 1;
	int mid = (begin+end)/2;
	while( begin <= end )
	{
		int comp = markCompare( vermark, markprop[mid] , columnnum);
		if( comp == 0 )
			return mid;
		if( comp == -1 )
			end = mid - 1;
		else
			begin = mid + 1;
		mid = (begin + end)/2;
	}
	//////////////////////////////////////////////////////////////////////////
	/*for( int k = 0; k < columnnum; k++)
		cout<<vermark[ k ]<<" ";
	cout<<endl;*/
		//////////////////////////////////////////////////////////////////////////
	cout<<"ERROR! UNABLE TO LOCATE THE VERTEX MARK IN DOUBLE MARK ARRAY!"<<endl;
	return mid;
}
void inline ContourHandler::markVer(
int**& markprop, int rownum,
int planei,
int planevernum, float*& planever,
int**& ver2planelistpos, vector<intvector>& ver2planelist,
int plane2planesnum, int* plane2planes,
float* pparam,
//result
int*& vertype,
int*& verval
)
{
	vertype = new int[ planevernum ];
	verval = new int[ planevernum ];

	int* vermark = new int[ plane2planesnum ];
	for( int i = 0; i < planevernum; i++)
	{
		//////////////////////////////////////////////////////////////////////////
		//cout<<"vertex position:"<<"{"<<planever[ 3*i ]<<","<<planever[ 3*i + 1]<<","<<planever[ 3*i + 2]<<endl;
		/*for( int k1 = 0; k1 < ver2planelist.size(); k1 ++ )
		{
			for( int k2 = 0; k2< ver2planelist[ k1 ].size(); k2++ )
				cout<<ver2planelist[ k1 ][ k2 ]<<",";
			cout<<endl;
		}*/
		//////////////////////////////////////////////////////////////////////////
		for( int j = 0; j < plane2planesnum; j++)
			vermark[ j ] = -1;
		int posinver2planelist = ver2planelistpos[ planei ][ i ];
		if( posinver2planelist != NOT_EXIST )
		{
			//set mark for those planes 0
			int planelistlen = ver2planelist[ posinver2planelist ].size();
			for( int j = 0; j < planelistlen; j++ )
			{
				int planei = ver2planelist[ posinver2planelist ][ j ];
				for( int k = 0; k < plane2planesnum; k++ )
				{
					if( plane2planes[ k ] == planei )
					{
						vermark[ k ] = 0;
						break;
					}
				}
			}
		}
		//for the planes mark haven't been set
		for( int j = 0; j < plane2planesnum; j++ )
		{
			if( vermark[ j ] == 0 )continue;

			//compute the relative position betweent this vertex and the plane
			if( MyMath::dotProduct( planever + 3*i, pparam + 4*plane2planes[ j ] ) < pparam[ 4*plane2planes[ j ] + 3 ])
				vermark[ j ] = -1;
			else
				vermark[ j ] = 1;
		}

		//with the mark, find the mark in markprop, to get the property for current vertex
		int markpos = findPos( markprop, rownum, plane2planesnum, vermark );
		vertype[ i ] = markprop[ markpos ][ plane2planesnum ];
		verval[ i ] = markprop[ markpos ][ plane2planesnum + 1 ];
	}

    delete []vermark;
}

void inline ContourHandler::gatherTopology(
int* planefaces, int planefacesnum,
int* planeedges, int planeedgesnum,
int* planevers,int planeversnum,
int** ssface, int* ssfaceedgenum,
int* ssedge,

int**& ver2edges, int*& ver2edgesnum,
int**& ver2faces, int*& ver2facesnum,
int**& edge2faces, int*& edge2facesnum
	)
{
	vector<intset> edge2face;	//the result makes all the faces corresponding to one edge in the increasing order
	vector<intset> ver2edge;
	vector<intset> ver2face;
	edge2face.resize( planeedgesnum );
	ver2edge.resize( planeversnum );
	ver2face.resize( planeversnum );

	//step1. go through each face, and set the edge2face
	for( int i = 0; i < planefacesnum; i++)
	{
		int facei = planefaces[ i ];
		int edgenum = ssfaceedgenum[ facei ] ;
		for( int j = 0; j < edgenum; j ++)
		{
			int edgei = ssface[ facei ][ j ];
			int tpos = findPosInIncSortedArray( planeedges, planeedgesnum, edgei );
			if( tpos == -1 )
			{
				cout<<"ERROR! Unable to find the edge in planeedges!"<<endl;
				continue;
			}
			edge2face[ tpos ].insert( i );
		}
	}
	intset::iterator iter;
	for( int i = 0; i < planeedgesnum; i ++)
	{
		int edgei = planeedges[ i ];
		for( int j = 0; j < 2; j ++)
		{
			int veri = ssedge[ 2*edgei + j ];
			int tpos = findPosInIncSortedArray( planevers, planeversnum, veri );
			if( tpos == -1 )
			{
				cout<<"ERROR! Unable to find the vertex in the planevers"<<endl;
				continue;
			}

			ver2edge[ tpos ].insert( i );
			iter = edge2face[ i ].begin();
			while( iter != edge2face[ i ].end())
			{
				ver2face[ tpos ].insert( *iter );
				iter ++;
			}
		}
	}

	//set those arrays
	ver2facesnum = new int[ planeversnum ];
	ver2faces = new int*[ planeversnum ];
	ver2edgesnum = new int[ planeversnum];
	ver2edges = new int*[ planeversnum ];
	edge2facesnum = new int[ planeedgesnum];
	edge2faces = new int*[ planeedgesnum];
	for( int i = 0; i < planeversnum; i ++)
	{
		ver2facesnum[ i ] = ver2face[ i ].size();
		ver2faces[ i ] = new int[ ver2facesnum[ i ]];
		iter = ver2face[ i ].begin();
		for( int j = 0; j < ver2facesnum[ i ]; j ++)
		{
			ver2faces[ i ][ j ] = *iter;
			iter++;
		}
		ver2edgesnum[ i ] = ver2edge[ i ].size();
		ver2edges[ i ] = new int[ ver2edgesnum[ i ]];
		iter = ver2edge[ i ].begin();
		for( int j = 0; j < ver2edgesnum[ i ]; j ++)
		{
			ver2edges[ i ][ j ] = *iter;
			iter ++;
		}
		ver2face[ i ].clear();
		ver2edge[ i ].clear();
	}
	ver2face.clear();
	ver2edge.clear();

	for( int i = 0; i< planeedgesnum; i++)
	{
		edge2facesnum[ i ] = edge2face[ i ].size();
		edge2faces[ i ] = new int[ edge2facesnum[ i ]];
		iter = edge2face[ i ].begin();
		for( int j = 0; j < edge2facesnum[ i ]; j ++)
		{
			edge2faces[ i ][ j ] = *iter;
			iter++;
		}

		edge2face[ i ].clear();
	}
	edge2face.clear();

	//////////////////////////////////////////////////////////////////////////
	////cout those topology information
	//for(int i = 0; i < planeversnum; i ++)
	//{
	//	///int veri = planevers[ i ];
	//	cout<<"the "<<i<<"th vertex!"<<endl<<"ver2edge:"<<endl;
	//	for( int j = 0; j < ver2edgesnum[ i ]; j++)
	//		cout<<ver2edges[ i ][ j ]<<" ";
	//	cout<<endl;
	//	cout<<"ver2face"<<endl;
	//	for( int j = 0; j < ver2facesnum[ i ]; j ++)
	//	{
	//		cout<<ver2faces[ i ][ j ]<<" ";
	//	}
	//	cout<<endl;
	//}
	//for( int i = 0; i < planeedgesnum; i ++)
	//{
	//	cout<<"the "<<i<<"th edge!"<<endl<<"edge2face"<<endl;
	//	for( int j = 0; j < edge2facesnum[ i ]; j ++)
	//	{
	//		cout<<edge2faces[ i ][ j ]<<" ";
	//	}
	//	cout<<endl;
	//}
	
	//////////////////////////////////////////////////////////////////////////
}
void inline ContourHandler::markEdge(
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
)
{
	int** ver2edges;
	int* ver2edgesnum;
	int** ver2faces;
	int* ver2facesnum;
	int** edge2faces;
	int* edge2facesnum;
	//step1.find those topology information
	gatherTopology(planefaces, planefacesnum,planeedges, planeedgesnum,planevers,planeversnum,
		ssface, ssfaceedgenum,ssedge,
		ver2edges, ver2edgesnum,ver2faces, ver2facesnum,edge2faces,edge2facesnum);

	//step2. go through each edge, and from the mark of the two vertices, put the vertex and edge into the proper face
	int** mapctrverface;
	mapctrverface = new int*[ pctrvernum ];
	for( int i = 0; i < pctrvernum; i++)
	{
		mapctrverface[ i ] = new int[ planefacesnum];
		for( int j = 0; j < planefacesnum; j++)
			mapctrverface[ i ][ j] = -1;
	}
	int tverval[ 2 ];
	int tvertype[ 2 ];
	int* tvers;
	int nverpos[ 2 ];
	int comedge;	//common subspace edge two vertices are on
	int comfacenum;		//how many common faces current two vertices are on
	int comface[ 2 ];	//common subspace face two vertices are on
	//int maprow[ 2]; //MING: commented because of unreferenced local variable
	for( int i = 0; i< pctredgenum; i ++)
	{
		tvers = pctredges + 4*i;	//two for indices and two for material, so toally 4 is one edge
		for( int j = 0; j < 2; j ++)
		{
			tverval[ j ] = verval[ tvers[ j ]];
			tvertype[ j ] = vertype[ tvers[ j ]];
		}
		int mul = tvertype[ 0 ] * tvertype[ 1 ];
		if ( (mul == 1) || ( mul == 2 ) || ( mul == 4 ))
		{
			//set common edge, if not exists, -1, otherwise, the real common edge index
			if ( mul == 1 )
			{
				comedge = commonOfTwoIncSortedArray_PosArray( 
					ver2edges[ tverval[ 0 ]], ver2edgesnum[ tverval[ 0 ]],
					ver2edges[ tverval[ 1 ]], ver2edgesnum[ tverval[ 1 ]]);
			}
			else if( mul == 2)
			{
				// if ( tvertype[ 0 ] = CTRTYPE_SUBVER )
				// MING
				if ( tvertype[ 0 ] == CTRTYPE_SUBVER )
				{
					comedge = findPosInIncSortedArray( ver2edges[ tverval[ 0 ]], ver2edgesnum[ tverval[ 0 ]], tverval[ 1 ]);
					if( comedge != -1 )
						comedge = tverval[ 1 ];
				}
				else if( tvertype[ 1 ] == CTRTYPE_SUBVER )
				{
					comedge = findPosInIncSortedArray( ver2edges[ tverval[ 1 ]], ver2edgesnum[ tverval[ 1 ]], tverval[ 0 ]);
					if( comedge != -1 )
						comedge = tverval[ 0 ];
				}
			}
			else
			{
				if( tverval[ 0 ] == tverval[ 1 ])
					comedge = tverval[ 0 ];
				else
					comedge = -1;
			}

			//set common face, and set the number of common faces.
			if( comedge != -1 )	//two common faces for this edge
			{
				comfacenum = edge2facesnum[ comedge ];
				for( int j = 0; j < edge2facesnum[ comedge ]; j++)
				{
					comface[ j ] = edge2faces[ comedge ][ j ];
				}
			}
			else
			{
				comfacenum = 1;
				if( mul == 1 )
				{
					comface[ 0 ] = commonOfTwoIncSortedArray_PosArray( ver2faces[ tverval[ 0 ]], ver2facesnum[ tverval[0 ] ] ,
					ver2faces[ tverval[1 ]], ver2facesnum[ tverval[ 1 ]]);

					//////////////////////////////////////////////////////////////////////////
					/*if( comface[ 0 ] == -1)
					{
						cout<<"vertex:"<< tverval[ 0 ]<<"vertex:"<< tverval[1 ]<<endl;
						for( int ti = 0; ti < ver2facesnum[ tverval[ 0 ]]; ti++)
						{
							cout<<ver2faces[ tverval[ 0 ]][ ti ]<<" ";
						}
						cout<<endl;
                        for( int ti = 0; ti < ver2facesnum[ tverval[ 1 ]]; ti++)
							cout<<ver2faces[ tverval[ 1 ]][ ti ]<<" ";
						cout<<endl;
					}*/
					//////////////// //////////////////////////////////////////////////////////
				}
				else if ( mul == 2 )
				{
					if( tvertype[ 0 ] == CTRTYPE_SUBVER )
					{
						comface[ 0 ] = commonOfTwoIncSortedArray_PosArray( ver2faces[ tverval[ 0 ]], ver2facesnum[ tverval[0 ] ] ,
							edge2faces[ tverval[1 ]], edge2facesnum[ tverval[ 1 ]]);
						//////////////////////////////////////////////////////////////////////////
						/*if( comface[ 0 ] == -1)
						{
							cout<<"vertex:"<< tverval[ 0 ]<<"edge:"<< tverval[1 ]<<endl;
							for( int ti = 0; ti < ver2facesnum[ tverval[ 0 ]]; ti++)
							{
								cout<<ver2faces[ tverval[ 0 ]][ ti ]<<" ";
							}
							cout<<endl;
							for( int ti = 0; ti < edge2facesnum[ tverval[ 1 ]]; ti++)
								cout<<edge2faces[ tverval[ 1 ]][ ti]<<" ";
							cout<<endl;
						}*/
						//////////////// //////////////////////////////////////////////////////////
					}
					else
					{
						comface[ 0 ] = commonOfTwoIncSortedArray_PosArray( ver2faces[ tverval[ 1 ]], ver2facesnum[ tverval[1 ] ] ,
							edge2faces[ tverval[0 ]], edge2facesnum[ tverval[ 0 ]]);
						//////////////////////////////////////////////////////////////////////////
						/*if( comface[ 0 ] == -1)
						{
							cout<<"edge:"<< tverval[ 0 ]<<"vertex:"<< tverval[1 ]<<endl;
							for( int ti = 0; ti < edge2facesnum[ tverval[ 0 ]]; ti++)
								cout<<edge2faces[ tverval[ 0 ]][ti]<<" ";
							cout<<endl;
							for( int ti = 0; ti < ver2facesnum[ tverval[ 1 ]]; ti++)
							{
								cout<<ver2faces[ tverval[ 1 ]][ ti ]<<" ";
							}
							cout<<endl;
						}*/
						//////////////// //////////////////////////////////////////////////////////
					}
					
				}
				else
				{
					comface[ 0 ] = commonOfTwoIncSortedArray_PosArray( edge2faces[ tverval[ 0 ]], edge2facesnum[ tverval[0 ] ] ,
						edge2faces[ tverval[1 ]], edge2facesnum[ tverval[ 1 ]]);
					//////////////////////////////////////////////////////////////////////////
					/*if( comface[ 0 ] == -1)
					{
						cout<<"edge:"<< tverval[ 0 ]<<"edge:"<< tverval[1 ]<<endl;
						for( int ti = 0; ti < edge2facesnum[ tverval[ 0 ]]; ti++)
							cout<<edge2faces[ tverval[ 0 ]][ ti ]<<" ";
						cout<<endl;
						for( int ti = 0; ti < edge2facesnum[ tverval[ 1 ]]; ti++)
						{
							cout<<edge2faces[ tverval[ 1 ]][ ti ]<<" ";
						}
						cout<<endl;
					}*/
					//////////////// //////////////////////////////////////////////////////////
				}
				if( comface[ 0 ] == -1)
				{
					cout<<"ERROR! Can not find the common face of the two vertices!"<<endl;
					cout<<"edge: vers"<<tvers[ 0 ]<<","<<tvers[1]<<" | "<< 
						tvertype[ 0 ] << ","<< tvertype[  1]<<"|"<<
						tverval[ 0 ] <<","<<tverval[ 1 ]<<endl;

					continue;
				}
			}
		}
		else	//mul = 3, 6 or 9
		{
			comedge = -1;
			comfacenum = 1;
			if( mul == 3)
			{
				if( tvertype[ 0 ] == CTRTYPE_SUBVER )
				{
					comface[ 0 ] = findPosInIncSortedArray( ver2faces[ tverval[ 0 ]], 
						ver2facesnum[ tverval[ 0 ]], tverval[ 1 ]);
					//////////////////////////////////////////////////////////////////////////
					/*if( comface[ 0 ] == -1)
					{
						cout<<"vertex:"<< tverval[ 0 ]<<"face:"<<tverval[1 ]<<endl;
						for( int ti = 0; ti < ver2facesnum[ tverval[ 0 ]]; ti++)
							cout<<ver2faces[ tverval[ 0 ]][ ti ] <<" ";
						cout<<endl;
						cout<<tverval[ 1 ]<<endl;
					}*/
					//////////////// //////////////////////////////////////////////////////////
					if( comface[ 0 ] == -1 )
					{
						cout<<"ERROR! Can not find the common face of the two contour vertices!"<<endl;
						cout<<"edge: vers"<<tvers[ 0 ]<<","<<tvers[1]<<" | "<< 
							tvertype[ 0 ] << ","<< tvertype[  1]<<"|"<<
							tverval[ 0 ] <<","<<tverval[ 1 ]<<endl;
						continue;
					}
					comface[ 0 ] = tverval[ 1 ];
				}
				else
				{
					comface[ 0 ] = findPosInIncSortedArray( ver2faces[ tverval[ 1 ]], 
						ver2facesnum[ tverval[ 1 ]], tverval[ 0 ]);
					//////////////////////////////////////////////////////////////////////////
					/*if( comface[ 0 ] == -1)
					{
						cout<<"face:"<< tverval[ 0 ]<<"vertex:"<< tverval[1 ]<<endl;
						cout<<tverval[ 0 ]<<endl;
						for( int ti = 0; ti < ver2facesnum[ tverval[ 1 ]]; ti++)
							cout<<ver2faces[ tverval[ 1 ]][ ti ]<<" ";
						cout<<endl;
						
					}*/
					//////////////// //////////////////////////////////////////////////////////
					if( comface[ 0 ] == -1 )
					{
						cout<<"ERROR! Can not find the common face of the two contour vertices!"<<endl;
						cout<<"edge: vers"<<tvers[ 0 ]<<","<<tvers[1]<<" | "<< 
							tvertype[ 0 ] << ","<< tvertype[  1]<<"|"<<
							tverval[ 0 ] <<","<<tverval[ 1 ]<<endl;
						continue;
					}
					comface[ 0 ] = tverval[ 0 ];
				}
			}
			else if ( mul == 6 )
			{
				if( tvertype[ 0 ] == CTRTYPE_SUBEDGE)
				{
					comface[ 0 ] = findPosInIncSortedArray( edge2faces[ tverval[ 0 ]], 
						edge2facesnum[ tverval[ 0 ]], tverval[ 1 ]);
					//////////////////////////////////////////////////////////////////////////
					/*if( comface[ 0 ] == -1)
					{
						cout<<"the index of the two vertices:"<<tvers[ 0 ]<<"\t"<<tvers[1]<<endl;
						cout<<"edge:"<<tverval[ 0 ]<<"face:"<< tverval[1 ]<<endl;
						for( int ti = 0; ti < edge2facesnum[ tverval[ 0 ]]; ti++)
							cout<<edge2faces[ tverval[ 0 ]][ ti ]<<" ";
						cout<<endl;
						cout<<tverval[ 1 ]<<endl;

					}*/
					//////////////// //////////////////////////////////////////////////////////
					if( comface[ 0 ] == -1 )
					{
						cout<<"ERROR! Can not find the common face of the two contour vertices!"<<endl;
						cout<<"edge: vers"<<tvers[ 0 ]<<","<<tvers[1]<<" | "<< 
							tvertype[ 0 ] << ","<< tvertype[  1]<<"|"<<
							tverval[ 0 ] <<","<<tverval[ 1 ]<<endl;
						continue;
					}
					comface[ 0 ] = tverval[ 1 ];
				}
				else
				{
					comface[ 0 ] = findPosInIncSortedArray( edge2faces[ tverval[ 1 ]], 
						edge2facesnum[ tverval[ 1 ]], tverval[ 0 ]);
					//////////////////////////////////////////////////////////////////////////
					/*if( comface[ 0 ] == -1)
					{
						cout<<"face:"<< tverval[ 0 ]<<"edge:"<< tverval[1 ]<<endl;
						cout<<tverval[ 0 ]<<endl;
						for( int ti = 0; ti < edge2facesnum[ tverval[ 1 ]]; ti++)
							cout<<edge2faces[ tverval[ 1 ]][ ti ]<<" ";
						cout<<endl;
						

					}*/
					//////////////// //////////////////////////////////////////////////////////
					if( comface[ 0 ] == -1 )
					{
						cout<<"ERROR! Can not find the common face of the two contour vertices!"<<endl;
						cout<<"edge: vers"<<tvers[ 0 ]<<","<<tvers[1]<<" | "<< 
							tvertype[ 0 ] << ","<< tvertype[  1]<<"|"<<
							tverval[ 0 ] <<","<<tverval[ 1 ]<<endl;
						continue;
					}
					comface[ 0 ] = tverval[ 0 ];
				}
			}
			else
			{
				if( tverval[ 0 ] != tverval[ 1 ] )
				{
					//////////////////////////////////////////////////////////////////////////
					/*cout<<"face:"<< tverval[ 0 ]<<"face:"<< tverval[1 ]<<endl;
					cout<<tverval[ 0 ]<<endl;
					cout<<tverval[ 1 ]<<endl;*/
					//////////////////////////////////////////////////////////////////////////
					cout<<"ERROR! Can not find the common face of the two contour vertices!"<<endl;
					cout<<"edge: vers"<<tvers[ 0 ]<<","<<tvers[1]<<" | "<< 
						tvertype[ 0 ] << ","<< tvertype[  1]<<"|"<<
						tverval[ 0 ] <<","<<tverval[ 1 ]<<endl;
					continue;
				}
				comface[ 0 ] = tverval[ 0 ];
			}
		}

		
		//for this contour edge, put it into the proper faces.
	//	for( int k = 0; k < 2; k ++)
	//	{
	//		maprow[ k ] = findPosInIncSortedArray( planevers, planeversnum, tvers[ k ]);
	//	}
		for( int j = 0; j < comfacenum; j ++)
		{
			int facei_planeface = comface[ j ];
			int facei_subspace = planefaces[ facei_planeface ];
			//check if the two vertices have already been put into current facectrver
			//int mapcol = findPosInIncSortedArray( planefaces, planefacesnum, facei );
			for( int k = 0; k < 2; k ++ )
			{
				//nverpos[ k ] = mapctrverface[ maprow[ k ]][ mapcol ];
				nverpos[ k ] = mapctrverface[ tvers[ k ] ][ facei_planeface];
				if( nverpos[ k ] == -1 )
				{
					for( int k1 = 0; k1 < 3; k1++ )
					{
						ctrfverposvec[ facei_subspace ].push_back( pctrver[ 3 * tvers[ k ] + k1 ]);
					}
					ctrfvertypevec[ facei_subspace ].push_back( tvertype[ k ]);
					if( tvertype[ k ] == CTRTYPE_SUBVER )
					{
							ctrfvervalvec[ facei_subspace].push_back( planevers[tverval[ k ]]);
					}
					else if( tvertype[ k ] == CTRTYPE_SUBEDGE )
					{
						ctrfvervalvec[ facei_subspace ].push_back( planeedges[ tverval[ k ]]);
					}
					else
						ctrfvervalvec[ facei_subspace ].push_back( facei_subspace );
					
					mapctrverface[ tvers[ k ]][ facei_planeface ] = nverpos[ k ] = ctrfverposvec[ facei_subspace ].size()/3 - 1;
				}
			}

			//put the new edge into the face
			ctrfedgevec[ facei_subspace ].push_back( nverpos[ 0 ]);
			ctrfedgevec[ facei_subspace].push_back( nverpos[ 1 ]);
			ctrfedgevec[ facei_subspace ].push_back( tvers[ 2 ]);	//material
			ctrfedgevec[ facei_subspace ].push_back( tvers[ 3 ]);
			ctrfedgeancestorvec[ facei_subspace ].push_back( i );
			if( comedge  !=  -1 )
			{
				ctrfedgetypevec[ facei_subspace].push_back( CTRTYPE_SUBEDGE );
				ctrfedgevalvec[ facei_subspace].push_back( planeedges[comedge]);
			}
			else
			{
				ctrfedgetypevec[ facei_subspace].push_back( CTRTYPE_SUBFACE );
				ctrfedgevalvec[ facei_subspace].push_back( facei_subspace);
			}
		}
	
	}
	tvers = NULL;
	for(int i = 0; i < planeversnum; i ++)
	{
		delete []ver2edges[ i ];
		delete []ver2faces[ i ];
	}
	delete []ver2edges;
	delete []ver2faces;
	delete []ver2edgesnum;
	delete []ver2facesnum;

	for( int i = 0; i < planeedgesnum; i ++)
	{
		delete []edge2faces[ i ];
	}
	delete []edge2faces;
	delete []edge2facesnum;

	//////////////////////////////////////////////////////////////////////////
	//cout<<"---\n";
	//////////////////////////////////////////////////////////////////////////
}

void inline ContourHandler::putContourIntoFace_OnePlane(
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
	)
{
	//step1. mark - ver, edge, face from the relative position against, last two are the real property
	int elemnum = planevernum + planeedgesnum + planefacesnum;
	int** markprop = new int*[ elemnum ];
	for( int i = 0; i < elemnum; i ++)
	{
		markprop[ i ] = new int[ plane2planesnum + 2 ];
		for( int j = 0; j < plane2planesnum; j++)
		{
			markprop[ i ][ j ] = -2;
		}
	}
	putContourIntoFace_OP_Mark(markprop,planei,ver2planesnum,ver2planes,plane2planesnum,  plane2planes,
		 planefacesnum,  planefaces, planeedgesnum,  planeedges,
		 planevernum,  planevers, planenum, pctrvernum,  pctrver, pctredgenum,  pctredges, pparam,
		 ssver,ssedge,  ssfaceedgenum,ssface);

	////////////////////////////////////////////////////////////////////////////
	/*cout<<"markprop:"<<endl;
	for( int i = 0; i < elemnum; i++ )
	{
		for( int j = 0; j < plane2planesnum + 2; j ++ )
			cout<<markprop[ i ][j]<<" ";
		cout<<endl;
	}*/
	////////////////////////////////////////////////////////////////////////////

	//step2. go through each contour vertex, and mark the property of them
	////2.1. sort these properties first
	sortMark(markprop, elemnum, plane2planesnum);	//only pass the number of column joining in comparison

	////////////////////////////////////////////////////////////////////////////
	/*char fname[ 100 ];
	itoa( planei, fname, 10);
	strcat( fname, ".mark");
	FILE* fout = fopen( fname,"w" );
	if( fout != NULL )
	{
		//print out all the marks
		fprintf( fout, "{{%d", markprop[ 0 ][ 0 ]);
		for( int i = 1; i < plane2planesnum; i ++)
		{
			fprintf( fout, ",%d", markprop[ 0 ][ i ]);
		}
		fprintf( fout, ",%d,%d", markprop[ 0 ][ plane2planesnum ] , markprop[ 0 ][ plane2planesnum + 1]+ 1);
		for( int i = 1; i < elemnum; i ++)
		{
			fprintf( fout, "},{%d", markprop[ i ][ 0 ]);
			for( int j = 1; j < plane2planesnum; j ++)
			{
				fprintf(fout, ",%d", markprop[ i ][ j ]);

			}
			fprintf( fout, ",%d,%d", markprop[ i][ plane2planesnum ], markprop[ i ][ plane2planesnum + 1]+ 1);
		}
		fprintf( fout, "}}");
		fclose( fout );
	}
	cout<<"markprop:"<<endl;
	for( int i = 0; i < elemnum; i++ )
	{
		for( int j = 0; j < plane2planesnum + 2; j ++ )
		cout<<markprop[ i ][j]<<" ";
		cout<<endl;
	}*/
/*	for( int i = 0 ;i < planenum; i ++ )
	{
		for( int j = 0; j < pctrvernum; j ++ )
		{
			if( ver2planelistpos[ i ][ j ] != -1)
			{
				int pos = ver2planelistpos[ i ][ j ];
				for( int k = 0; k < ver2planelist[ pos ].size(); k ++ )
				{
					cout<<ver2planelist[ pos][ k ]<<",";
				}
				cout<<endl;
			}
		}
	}*/
	////////////////////////////////////////////////////////////////////////////

	//2.2 go through each contour vertex, and set the property
	int* vertype;
	int* verval;
	markVer(markprop, elemnum,planei,pctrvernum, pctrver,ver2planelistpos, ver2planelist,
		plane2planesnum, plane2planes,pparam,vertype,verval);
	
	//////////////////////////////////////////////////////////////////////////
	//FILE* fout = fopen("ctrver.txt","w");
	//if( fout!= NULL )
	//{
	//	//contour vertex positions
	//	fprintf( fout, "{{{%f,%f,%f}",pctrver[0], pctrver[1],pctrver[2]);
	//	for(int i = 1; i < pctrvernum; i ++)
	//	{
	//		fprintf( fout, ",{%f,%f,%f}", pctrver[ 3*i ], pctrver[ 3*i + 1], pctrver[ 3*i + 2 ]);
	//	}

	//	//contour vertex properties
	//	fprintf( fout, "},{{%d,%d}", vertype[ 0 ], verval[ 0 ] + 1);
	//	for( int i = 1; i < pctrvernum ; i++)
	//		fprintf(fout, ",{%d,%d}", vertype[ i ], verval[ i ] + 1 );

	//	//planevers
	//	fprintf( fout, "},{%d", planevers[ 0 ] + 1);
	//	for( int i = 1; i < planevernum; i ++)
	//		fprintf( fout, ",%d", planevers[ i ] + 1);

	//	//planeedges
	//	fprintf( fout, "},{%d", planeedges[ 0 ] + 1);
	//	for( int i = 1; i < planeedgesnum; i ++)
	//		fprintf( fout, ",%d", planeedges[ i ] + 1 );

	//	//planefaces
	//	fprintf( fout ,"},{%d", planefaces[ 0 ] +1 );
	//	for( int i = 1; i < planefacesnum; i ++)
	//		fprintf( fout, ",%d", planefaces[ i ] + 1);

	//	//ctredge
	//	//fprintf( fout, "},{{%d,%d}", pctredges[ 0 ] + 1, pctredges[ 1 ] + 1 );
	//	
	//	fprintf( fout, "}}");
	//	fclose( fout );
	//}
	//////////////////////////////////////////////////////////////////////////
	//step3. go through each edge, and mark the property of them, at the same time, put them into correct face
	markEdge(
		pctrvernum,pctrver,
		pctredgenum, pctredges,
		vertype,  verval, 

		 planevers,	planevernum,		//in the increasing order of vertices indices
		 planeedges,	planeedgesnum,	//in the increasing order of edges indices
		 planefaces,	 planefacesnum,	//in the increasing order of faces indices

		 ssedge,
		 ssface,  ssfaceedgenum,

		ctrfverposvec,  ctrfvertypevec,ctrfvervalvec, 
		ctrfedgevec,  ctrfedgetypevec, ctrfedgevalvec,ctrfedgeancestorvec);
}

void ContourHandler::writeElement2Planes_db(
	int ssvernum, int ssedgenum, int ssfacenum, int planenum,
int* ver2planesnum,int** ver2planes,int* edge2planesnum,
int** edge2planes,int* plane2planesnum,int** plane2planes,
int* planefacesnum,int** planefaces,int* planeedgesnum,
int** planeedges,int* planeversnum,int** planevers
	)
{
	FILE* fout = fopen ("mmdebug/planeelem.txt", "w");
	if( fout == NULL )

	{
		cout<<"Unable to open file planeelme.txt to write!"<<endl;
		return;
	}

	//write out map from ver, edge and face 2 planes
	//ver2planes
	fprintf( fout , "{{{%d", ver2planes[ 0 ][ 0 ] + 1);
	for( int i = 1; i < ver2planesnum[ 0 ]; i++ )
	{
		fprintf( fout, ",%d", ver2planes[ 0 ][ i ] + 1);
	}
	for( int i = 1; i < ssvernum; i ++)
	{
		fprintf( fout, "},{%d", ver2planes[ i ][ 0 ] + 1 );
		for(int j = 1; j < ver2planesnum[ i ]; j ++)
			fprintf( fout, ",%d", ver2planes[ i ][ j ] + 1);
	}
	
	//edge2planes
	fprintf(fout, "}},{{%d", edge2planes[ 0 ][ 0 ] + 1);
	for( int i = 1; i < edge2planesnum[ 0 ]; i ++)
	{
		fprintf( fout, ",%d", edge2planes[ 0 ][ i ] + 1);
	}
	for( int i = 1; i < ssedgenum; i ++)
	{
		fprintf( fout, "},{%d", edge2planes[ i ][ 0 ] + 1);
		for( int j = 1; j < edge2planesnum[ i ]; j ++)
			fprintf( fout, ",%d", edge2planes[ i ][ j ] + 1) ;
	}

	//write out map from planes to vertex, edge and faces
	//planevers
	fprintf( fout, "}},{{%d",  planevers[ 0 ][ 0 ] + 1);
	for( int i = 1; i < planeversnum[ 0 ]; i ++)
	{
		fprintf( fout, ",%d", planevers[ 0 ][ i ] + 1);
	}

	for( int i = 1; i < planenum; i++)
	{
		fprintf( fout, "},{%d", planevers[ i ][ 0 ] + 1);
		for( int j = 1; j < planeversnum[ i ]; j ++)
		{
			fprintf( fout, ",%d", planevers[ i ][ j ] + 1 );
		}
	}

	//planeedges
	fprintf( fout, "}},{{%d", planeedges[ 0 ][ 0 ]+ 1);
	for( int i = 1; i < planeedgesnum[ 0 ]; i ++)
	{
		fprintf(fout, ",%d", planeedges[ 0 ][ i ] + 1) ;
	}
	for( int i = 1; i < planenum; i ++)
	{
		fprintf( fout, "},{%d", planeedges[ i ][ 0 ] + 1);
		for( int j = 1; j < planeedgesnum[ i ]; j ++)
			fprintf( fout, ",%d", planeedges[ i ][ j ] + 1);
	}
    
	//planeface
	fprintf( fout, "}},{{%d", planefaces[ 0 ][ 0 ] + 1);
	for( int i = 1; i < planefacesnum[0 ]; i++)
		fprintf( fout, ",%d", planefaces[ 0 ][ i ] + 1);
	for( int i =1 ;  i < planenum; i ++)
	{
		fprintf( fout, "},{%d", planefaces[ i ][ 0 ] + 1);
		for( int j = 1; j < planefacesnum[ i ]; j ++)
			fprintf( fout, ",%d", planefaces[ i ][ j ] + 1);
	}

	//plane2planes
	fprintf( fout, "}},{{%d", plane2planes[ 0 ][ 0 ] + 1);
	for( int i = 1; i < planefacesnum[ 0 ]; i ++)
	{
		fprintf(fout, ",%d", plane2planes[ 0 ][ i ] + 1);
	}
	for (int i = 1; i < planenum; i ++)
	{
		fprintf( fout, "},{%d", plane2planes[ i ][ 0 ] + 1);
		for( int j = 1; j < plane2planesnum[ i ]; j++ )
		{
			fprintf( fout, ",%d", plane2planes[ i ][ j ] + 1);
		}
	}
	fprintf( fout, "}}}");
	fclose( fout );
}

void ContourHandler::writeContourInFace(
	int facenum,
int* ctrfvernum,float** ctrfverpos,int** ctrfvertype,int** ctrfverval,
int* ctrfedgenum, int** ctrfedge, int** ctrfedgetype, int** ctrfedgeval
)
{
    FILE* fout = fopen( "/Users/Research/Geometry/MM/ConvertFolder/facectr.txt","w");
	if( fout== NULL )
	{
		cout<<"Unable to open facectr.txt to write!"<<endl;
		return;
	}

	//write face by face
	for( int facei = 0; facei < facenum; facei ++)
	{
		if( ctrfverpos[ facei ] == NULL )
		{
			if( facei == 0 )
				fprintf( fout, "{{{},{},{},{");
			else
				fprintf( fout, "}},{{},{},{},{");
			continue;
		}
		//write contour vertex
		if( facei == 0)
		{
			fprintf( fout, "{{{{%f,%f,%f}", ctrfverpos[ 0 ][ 0 ], ctrfverpos[ 0 ][1], ctrfverpos[ 0 ][ 2 ]);
		}
		else
			fprintf( fout, "}},{{{%f,%f,%f}", ctrfverpos[ facei ][ 0 ], ctrfverpos[ facei][1], ctrfverpos[ facei ][ 2 ]);

		for( int i = 1; i < ctrfvernum[ facei ]; i ++)
		{
			fprintf( fout, ",{%f,%f,%f}", ctrfverpos[ facei ][ 3*i ],ctrfverpos[ facei ][ 3*i + 1],ctrfverpos[ facei ][ 3*i +2 ]);
		}

		//write vertex property, { vertype, verval }
		fprintf( fout, "},{{%d,%d}", ctrfvertype[ facei ][ 0 ], ctrfverval[ facei ][ 0 ] + 1);
		for( int i = 1; i < ctrfvernum[ facei ]; i ++)
		{
			fprintf( fout, ",{%d,%d}", ctrfvertype[ facei ][ i ], ctrfverval[ facei ][ i ] + 1);
		}

		//write contour edges {{ver1,ver2}, {mat1,mat2}}
		fprintf( fout, "},{{{%d,%d},{%d,%d}}", ctrfedge[ facei ][ 0 ] + 1, ctrfedge[ facei ][ 1 ] + 1, 
			ctrfedge[facei][ 2 ] + 1, ctrfedge[ facei ][3 ] + 1	);
		for( int i = 1; i < ctrfedgenum[ facei ]; i ++)
			fprintf( fout, ",{{%d,%d},{%d,%d}}", ctrfedge[ facei ][ 4*i ] + 1, ctrfedge[ facei ][ 4*i + 1 ] + 1, 
			ctrfedge[facei][ 4*i + 2 ] + 1, ctrfedge[ facei ][ 4*i + 3 ] + 1);

		//write edge property {edgetype, edgeval}
		fprintf( fout, "},{{%d,%d}", ctrfedgetype[ facei ][ 0 ] , ctrfedgeval[ facei ][ 1 ] + 1);
		for( int i = 1; i < ctrfedgenum[ facei ]; i ++)
			fprintf( fout, ",{%d,%d}", ctrfedgetype[ facei ][ i ], ctrfedgeval[ facei ][ i ] + 1);

	}
	fprintf ( fout, "}}}");
	fclose( fout );
}
void ContourHandler::putContourIntoFace
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
int**& ctrfedgeancestor,

int radialPlaneNum
)
{
	//////////////////////////////////////////////////////////////////////////
	/*for( int i = 0; i < ver2planelist.size(); i ++ )
	{
		for( int j = 0; j < ver2planelist[ i ].size(); j ++)
		{
			cout<<ver2planelist[ i ][ j ]<<"j"<<",";
			if( j == (ver2planelist[ i ].size() - 1) )
				cout<<"   |||   ";
		}
	}
	cout<<endl;*/
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	/*FILE* fout = fopen( "param.txt","w");
	if( fout != NULL )
	{
		fprintf( fout, "{{%f", pparam[ 0 ]);
		for( int i = 1; i < 4; i ++ )
			fprintf( fout, ",%f", pparam[ i ]);
		for( int i = 1; i < planenum + 6; i ++ )
		{
			fprintf( fout, "},{%f", pparam[ 4*i ]);
			for( int j = 1; j < 4; j ++)
				fprintf( fout, ",%f", pparam[ 4*i + j]);
		}
		fprintf( fout, "}}");
			fclose( fout );
	}*/
	//////////////////////////////////////////////////////////////////////////
	//step1. for each vertex, each edge,get the planes that they are on
	//for each plane, get the planes that they cross
	vector<intset> ver2planes;
	vector<intset> edge2planes;
	vector<intset> plane2planes;
    vector<intvector> planefaces;
	vector<intvector> planeedges;
	vector<intvector> planevers;
	getElement2Planes(planenum,ssvernum, ssedgenum, ssedge,ssfacenum, ssfaceedgenum,
		ssface,ssface_planeindex,ssspacenum,ssspacefacenum,ssspace,		
		ver2planes,	edge2planes, plane2planes, planefaces, planeedges, planevers,
		radialPlaneNum);
	int* ver2planesnum;
    int** daver2planes;
	int* edge2planesnum;
	int** daedge2planes;
	int* plane2planesnum;
	int** daplane2planes;
	int* planefacesnum;
	int** daplanefaces;
	int* planeedgesnum;
	int** daplaneedges;
	int* planeversnum;
	int** daplanevers;
	ver2planesnum = new int[ ssvernum ];
	daver2planes = new int*[ ssvernum ];
	edge2planesnum = new int[ ssedgenum ];
	daedge2planes	= new int*[ ssedgenum];		
	plane2planesnum	= new int[ planenum];
	daplane2planes = new int*[ planenum];
	planefacesnum	= new int[ planenum ];
	daplanefaces	= new int*[ planenum ];
	planeedgesnum= new int[ planenum ];
	daplaneedges= new int*[ planenum ];
	planeversnum= new int[ planenum ];
	daplanevers= new int*[ planenum ];
	intset::iterator iter;
	for( int i = 0; i < ssvernum; i ++)
	{
		ver2planesnum[ i ] = ver2planes[ i ].size();
		daver2planes[ i ] = new int[ ver2planesnum[ i ]];
		int j = 0;
		for( iter = ver2planes[ i ].begin(); iter!=ver2planes[ i ].end(); iter++)
		{
			daver2planes[ i ][ j++ ] = *iter;
		}
		ver2planes[ i ].clear();
	}
	ver2planes.clear();

	for( int i = 0; i < ssedgenum; i++)
	{
		edge2planesnum[ i ] = edge2planes[i].size();
		daedge2planes[ i ] = new int[ edge2planesnum[ i ]];
		int j = 0;
		for( iter = edge2planes[ i ].begin(); iter!= edge2planes[ i ].end(); iter++)
		{
			daedge2planes[ i ][ j++ ] = *iter;
		}
		edge2planes[ i ].clear();
	}
	edge2planes.clear();

	for( int i = 0; i < planenum; i ++)
	{
		plane2planesnum[ i ] = plane2planes[ i ].size();
		daplane2planes[ i ] = new int[ plane2planesnum[ i ]];

		planefacesnum[ i ] = planefaces[ i ].size();
		daplanefaces[ i ] = new int[ planefacesnum[ i ]];

		planeedgesnum[ i ] =  planeedges[ i ].size();
		daplaneedges[ i ] = new int[ planeedgesnum[ i ]];

		planeversnum[ i ] = planevers[ i ].size();
		daplanevers[ i ] = new int[ planeversnum[ i ]];

		iter = plane2planes[ i ].begin();
		for( int j = 0; j < plane2planesnum[ i ]; j++)
		{
			daplane2planes[ i ][ j ] = *iter;
			iter++;
		}
			
		for( int j = 0; j < planefacesnum[ i ]; j++)
		{
			daplanefaces[ i ][ j ] = planefaces[ i ][ j ];
		}

		for( int j = 0; j < planeedgesnum[ i ]; j ++)
		{
			daplaneedges[ i ][ j ] = planeedges[ i ][ j ];
		}

		for( int j = 0; j <planeversnum[ i ]; j ++)
		{
			daplanevers[ i ][ j ] = planevers[ i ][ j ];
		}

		plane2planes[ i ].clear();
		planefaces[ i ].clear();
		planeedges[ i ].clear();
		planevers[ i ].clear();
	}
	plane2planes.clear();
	planefaces.clear();
	planeedges.clear();
	planevers.clear();
	for( int i = 0; i < planenum; i ++)
	{
		selectSortInc( daplanevers[ i ], planeversnum[ i  ]);
		selectSortInc( daplaneedges[ i ], planeedgesnum[ i ]);
		selectSortInc( daplanefaces[ i ], planefacesnum[ i ]);
	}
	
	//////////////////////////////////////////////////////////////////////////
	/*writeElement2Planes_db(
		ssvernum, ssedgenum, ssfacenum, planenum,
		ver2planesnum, daver2planes,edge2planesnum,
		daedge2planes,plane2planesnum, daplane2planes,
		planefacesnum, daplanefaces, planeedgesnum,
		 daplaneedges, planeversnum, daplanevers		);*/
	//////////////////////////////////////////////////////////////////////////

	//step2. process each plane and put the contour into each face
	vector<floatvector> ctrfverposvec;
	vector<intvector> ctrfvertypevec;
	vector<intvector>ctrfvervalvec;
	vector<intvector>ctrfedgevec;
	vector<intvector> ctrfedgetypevec;
	vector<intvector> ctrfedgevalvec;
	vector<intvector> ctrfedgeancestorvec;
	ctrfverposvec.resize( ssfacenum );
	ctrfvertypevec.resize( ssfacenum );
	ctrfvervalvec.resize( ssfacenum );
	ctrfedgevec.resize( ssfacenum );
	ctrfedgetypevec.resize( ssfacenum );
	ctrfedgevalvec.resize( ssfacenum );
	ctrfedgeancestorvec.resize( ssfacenum);
	//////////////////////////////////////////////////////////////////////////
	//for( int i = 0; i < planenum; i ++)
		//////////////////////////////////////////////////////////////////////////
	for( int i = 0; i < planenum; i ++)
	{
		putContourIntoFace_OnePlane(i, 
		ver2planesnum, daver2planes, edge2planesnum, daedge2planes, 
		plane2planesnum[ i ],daplane2planes[ i ], planefacesnum[ i ],daplanefaces[ i ],
		planeedgesnum[ i ], daplaneedges[ i ], planeversnum[ i ], daplanevers[ i ],
		planenum, pctrvernum[ i ], pctrverpos[ i ], pctredgenum[ i ], pctredge[i], pparam,
		ver2planelistpos,ver2planelist, ssvernum, ssver, ssedgenum, ssedge, ssfacenum, ssfaceedgenum,
		ssface,ssface_planeindex, ctrfverposvec, ctrfvertypevec, ctrfvervalvec, ctrfedgevec, ctrfedgetypevec,
		ctrfedgevalvec, ctrfedgeancestorvec);
	}
	
	//set the ** array ,and clear those vectors
	/*int*& ctrfvernum,float**& ctrfverpos,int**& ctrfvertype,int**& ctrfverval,
		int*& ctrfedgenum, int**& ctrfedge, int**& ctrfedgetype, int**& ctrfedgeval*/
	for( int i = 0; i < ssfacenum; i ++)
	{
		ctrfvernum[ i ] = ctrfverposvec[ i ].size()/3;
		if( ctrfvernum[ i ] == 0)
		{
			ctrfedgenum[ i ] = 0;
			ctrfverpos[ i ] = NULL;
			ctrfvertype[ i ] = NULL;
			ctrfverval[ i ] = NULL;
			ctrfedge[ i ] = NULL;
			ctrfedgetype[ i ] = NULL;
			ctrfedgeval[ i ] = NULL;
			ctrfedgeancestor[ i ] = NULL;
			continue;
		}
		ctrfverpos[ i ] = new float[ctrfverposvec[ i ].size()];
		ctrfvertype[ i ] = new int[ ctrfvertypevec[ i ].size()];
		ctrfverval[ i ] = new int[ ctrfvervalvec[ i ].size()];
		for (unsigned int j = 0; j < ctrfverposvec[i].size(); j++)
		{
			ctrfverpos[ i ][ j ] = ctrfverposvec[ i][j ];
		}
		for (unsigned int j = 0; j < ctrfvertypevec[i].size(); j++)
			ctrfvertype[ i ][ j ] = ctrfvertypevec[ i ][ j ];
		for (unsigned int j = 0; j < ctrfvervalvec[i].size(); j++)
		{
			ctrfverval[ i ][ j ] = ctrfvervalvec[ i ][j ];
		}
		ctrfedgenum[ i ] = ctrfedgevec[ i ].size()/4;
		ctrfedge[ i ] = new int[ctrfedgevec[ i ].size()];
		for (unsigned int j = 0; j < ctrfedgevec[i].size(); j++)
		{
			ctrfedge[ i ][ j ]= ctrfedgevec[ i ][ j ];
		}
		ctrfedgetype[ i ] = new int[ ctrfedgetypevec[ i ].size()];
		for (unsigned int j = 0; j < ctrfedgetypevec[i].size(); j++)
		{
			ctrfedgetype[ i ][ j ] = ctrfedgetypevec[ i ][ j ];
		}
		ctrfedgeval[ i ] = new int[ ctrfedgevalvec[ i ].size()];
		for (unsigned int j = 0; j < ctrfedgevalvec[i].size(); j++)
		{
			ctrfedgeval[ i ][ j ] = ctrfedgevalvec[ i ][ j ];
		}
		ctrfedgeancestor[ i ] = new int[ ctrfedgeancestorvec[ i ].size()];
		for (unsigned int j = 0; j < ctrfedgeancestorvec[i].size(); j++)
			ctrfedgeancestor[ i ][ j ]= ctrfedgeancestorvec[ i][ j ];
	}

	//clear temporary arrays
	for( int i = 0 ; i < ssvernum; i ++)
	{
		delete []daver2planes[ i ];
	}
	delete []ver2planesnum;
	delete []daver2planes;
	for( int i = 0; i < ssedgenum; i++)
		delete []daedge2planes[ i ];
	delete []daedge2planes;
	delete []edge2planesnum;
	for( int i = 0 ;i < planenum; i ++)
	{
		delete []daplane2planes[ i ];
		delete []daplanefaces[ i ];
		delete []daplaneedges[ i ];
		delete []daplanevers[ i ];
	}
	delete []daplane2planes;
	delete []daplanefaces;
	delete []daplaneedges;
	delete []daplanevers;

	//////////////////////////////////////////////////////////////////////////
	/*writeContourInFace(		ssfacenum,
		 ctrfvernum,ctrfverpos, ctrfvertype, ctrfverval,
		 ctrfedgenum,  ctrfedge,  ctrfedgetype,  ctrfedgeval
		);*/
	//////////////////////////////////////////////////////////////////////////
}
