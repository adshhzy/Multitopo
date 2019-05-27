/*
 * The input contour planes share one common edge, i.e, all the cutting planes rotate around a common line
 * Make the contours on these planes consistent with each other.
 * 
 */
#include "../ContourHandler/ContourHandler.h"

void ContourHandler::interCommonLineWithOnePlane(
	float param[ 4 ],					//parameter of the plane
	floatvector& ctrvers, intvector& ctrvermark,		//the contour vertices on this cutting plane, and their marks relative to the common edge
	intvector& ctredges,	//contour edges on this cutting plane
	float comndir[ 3 ], float comnpt[ 3 ],		//common line parameters
	intvector& interptPos			//remeber the position of the vertices which lie on the common line
	)
{
	floatvector fvec_ctrvermark;

	float dir [ 3 ];	
	MyMath::crossProduct(param, comndir , dir);

	//step1, mark all the contour vertices
	int ctrvernum = ctrvers.size() / 3;
	fvec_ctrvermark.resize( ctrvernum );
	float comnpt_ctrver_vec[ 3 ];
	int k = 0;
	for( int i = 0; i < ctrvernum; i ++ )
	{
		for( int j = 0; j < 3; j ++ )
			comnpt_ctrver_vec[ j ] = ctrvers[ k + j ] - comnpt[ j ];
		k += 3;

		float val = MyMath::dotProduct(dir, comnpt_ctrver_vec);
		fvec_ctrvermark[ i ] = val ;
		if( MyMath::isEqualInToler(val , 0, TOLERANCE_FOUR))
		{
			fvec_ctrvermark[ i ] = 0;	//on the common line	
			ctrvermark[ i ] = 0;
			interptPos.push_back( i );
		}
		else if( val > 0 )
		{
			ctrvermark[ i ] = 1;
		}
		else
			ctrvermark[ i ] = -1;
	}

	//step2, go through all the edges, split it when necessary
	int ctredgenum = ctredges.size()/4;
	int ctredgesize = ctredges.size();

	k = 0;
	for( int i = 0; i < ctredgenum; i ++ )
	{
		if( ctrvermark[ ctredges[ k ] ] * ctrvermark[ ctredges[ k + 1]] != -1 )
		{
			k += 4;
			continue;
		}

		//1*-1 intersect with the common line, split the contour edge
		//computer intersection point,
		//ctrvers, ctrvermark, interptpos
		float interpt[ 3 ];
		MyMath::getPtOnSeg( &ctrvers[ ctredges[ k ] * 3],
			&ctrvers[ ctredges[ k + 1] * 3],
			fvec_ctrvermark[ ctredges[ k ] ]/(fvec_ctrvermark[ ctredges[ k ] ] - fvec_ctrvermark[ ctredges[ k  + 1] ]),
			interpt);
		for( int tk = 0; tk < 3; tk ++ )
		{
			ctrvers.push_back( interpt[ tk ]);
		}
		ctrvermark.push_back( 0 );
		interptPos.push_back( ctrvernum );

		//split the contour edge
		//ctredges
		ctredges.resize( ctredgesize + 4 );
		memcpy( &ctredges[ ctredgesize], &ctredges[ k ],  sizeof( int ) * 4 );
		//for( int tk = 0; tk < 4; tk ++ )
		//	ctredges.push_back( ctredges[ k + tk]);
		ctredges[ k + 1 ] = ctrvernum;
		ctredges[ctredgesize ] = ctrvernum;
		ctredgesize += 4;

		ctrvernum ++;

		k += 4;
	}

	fvec_ctrvermark.clear();

	//////////////////////////////////////////////////////////////////////////
	////check code
	//if( ctrvermark.size() != ctrvers.size()/3 )
	//{
	//	cout<<"Vertex size for the mark of the vertex doesn't match with the actual vertex size!"<<endl;
	//}
	//else
	//{
	//	cout<<"GOOD! Vertex size for the mark of the vertex matches with the actual vertex size!"<<endl;
	//}
	//int versize = interptPos.size();
	//for( int i = 0; i < versize; i ++ )
	//{
	//	int intervi = interptPos[ i ];
	//	if( ctrvermark[ intervi ] != 0 )
	//	{
	//		cout<<"Oops, pushed in a vertex which doesn't lie on the common edge!"<<endl;
	//	}
	//}
	//////////////////////////////////////////////////////////////////////////
}

void ContourHandler::clearTheInputData(
									   int& planenum, float*& pparam,
									   int*& pctrvernum, float**& pctrvers,
									   int*& pctredgenum, int**& pctredges)
{
	for( int i = 0; i < planenum - 6; i ++ )
	{
		delete []pctrvers[ i ];
		delete []pctredges[ i ];
	}
	delete []pctrvers;
	delete []pctredges;
	delete []pctrvers;
	delete []pctredges;
	delete []pctrvernum;
	delete []pctredgenum;
	delete []pparam;

	planenum = 0;
}

/*
reserve the planes with the same amout intersection point
and reserve those with the most hitted intersection point
*/
void ContourHandler::reserverPlanes
(
 int tplanenum,
 vector<intvector> interptpos,
 int*& reservePlaneList,
 int& reserveNum
 )
{
	int* interptcount = new int[ tplanenum ];
	intvector ptcount;
	int ptcountsize = 0;
	intvector crspplanenum ;

	for( int i = 0; i < tplanenum; i ++ )
	{
		interptcount[ i ] = interptpos[ i ].size();

		if( ptcountsize == 0  )	//the first one
		{
			ptcount.push_back( interptcount[ i ]);
			ptcountsize++;
			crspplanenum.push_back( 1 );
			continue;
		} 

		//not the first one
		int foundind = -1;
		for( int j = 0; j < ptcountsize;  j ++)
		{
			if( ptcount[ j ] == interptcount[ i ]) 
			{
				foundind = j;
				break;
			}
		}

		//the plane i has different number of intersection points than the 
		//previous planes
		if( foundind == -1  ) 
		{
			ptcount.push_back( interptcount[ i ]);
			ptcountsize ++;
			crspplanenum.push_back( 1 );
		}
		else	//the plane has the same number of intersection points as some ones
		{
			crspplanenum[ foundind ] ++;
		}
	}

	//find the most hitted number, and preserve those cutting planes
	//////////////////////////////////////////////////////////////////////////
	//for debug reason
	//if( ptcountsize != ptcount.size()) 
	//{
	//	cout<<"wrong ! ptcountsize is not equal to the real size!"<<endl;		
	//	//	clearTheInputData(
	//	//		planenum, pparam,	pctrvernum, pctrvers, pctredgenum, pctredges);

	//	delete []interptcount;
	//	ptcount.clear();
	//	crspplanenum.clear();
	//	reservePlaneList = NULL;

	//	reserveNum = 0;
	//	return;
	//}
	//////////////////////////////////////////////////////////////////////////

	if( ptcountsize == 1)	//only one kind of intersection points number
	{
		reservePlaneList = new int[ tplanenum ];
		for( int i = 0; i< tplanenum; i ++ )
			reservePlaneList[ i ] = i;

		delete []interptcount;
		ptcount.clear();
		crspplanenum.clear();
		reserveNum = tplanenum;
		return;
	}

	//some planes have to be discarded
	int max_interptcount = crspplanenum[ 0 ];
	int max_pos = 0;
	for( int i = 1; i < ptcountsize; i ++ )	//find the one with most cutting planes
	{
		if( crspplanenum[ i ] > max_interptcount )
		{
			max_interptcount = crspplanenum[ i ];
			max_pos = i;
		}
	}

	reservePlaneList = new int[ crspplanenum[ max_pos ]];
	int countshouldbe = ptcount[ max_pos ];
	int k = 0;
	for( int i = 0; i < tplanenum; i++ )
	{
		if( interptcount[ i ] == countshouldbe )
			reservePlaneList[ k ++ ] = i;
	}
	reserveNum = crspplanenum[ max_pos ];

	delete []interptcount;
	ptcount.clear();
	crspplanenum.clear();		
}

bool ContourHandler::adjustInterPt
(
	int* reservelist, int reservenum,
	vector<intvector>& interptpos,

	vector<floatvector>& ctrvers,

	float comndir[ 3 ], float comnpt[ 3 ] 
 )
{
	bool result = true;

	int interptnum = interptpos[ reservelist[ 0 ]].size();

	float* oldtOnRay = new float[ interptnum ];
	float* tOnRay = new float[ interptnum ];

	int plane0 = reservelist[ 0 ];

	for( int i = 0; i < reservenum; i ++ )
	{
		int planei = reservelist[ i ];

		//compute the ts for the intersection points
		for( int j = 0; j < interptnum; j ++ )
		{
			int verstarti = interptpos[ planei ][ j ] * 3;
			tOnRay[ j ] = getTOnRay(comnpt, comndir, &ctrvers[ planei ][ verstarti ]);
		}
		
		//sort the ts and change the sequence of interptpos as well
		//select sort
		for( int j = 0; j < interptnum - 1; j ++ )
		{
			int mini = j;
			float minval = tOnRay[ j ];
			for( int k = j + 1; k < interptnum; k ++ )
			{
				if( tOnRay[ k ] < minval )
				{
					minval = tOnRay[ k ];
					mini = k;
				}
			}

			if( mini == j )
				continue;

			//change the current one with the mininum one
			float tval = tOnRay[ mini ];
			tOnRay[ mini ] = tOnRay[ j ];
			tOnRay[ j ] = tval;

			int tpos = interptpos[ planei ][ mini ] ;
			interptpos[ planei ][ mini ] = interptpos[ planei ][ j ];
			interptpos[ planei ][ j ] = tpos;
		}

		
		if( i == 0 )
		{
			float* tpointer = oldtOnRay;
			oldtOnRay = tOnRay;
			tOnRay = tpointer;
			tpointer = NULL;
			continue;
		}

		//check if it is consistent.
		for( int j = 0; j < interptnum; j ++ )
		{
			/*if( !MyMath::isEqualInToler( oldtOnRay[ j ] , tOnRay[ j ], TOLERANCE_ONE ) )
			{
				result = false;
				break;
			}*/

			//set the intersection point to the same as the first one
			int curverstarti = interptpos[ planei ][ j ] * 3;
			int curverstarti0 = interptpos[ plane0][ j ] * 3;
			memcpy( &ctrvers[ planei ][ curverstarti ], & ctrvers[ plane0 ][ curverstarti0 ], sizeof( float ) * 3 );	
		}

		if( !result )
			break;
	}
	
	delete []oldtOnRay;
	delete []tOnRay;

	return result;
}

//return the half space, all the normals are flipped to,
//0 - flip to the x+ halfspace,
//1 - flip to the y+ half space
//2 - z+ half space
int ContourHandler::flipNormalToHalfSpace
( int* reserverlist, int reservenum,
  float*& pparam, 
  vector<intvector>& ctredges,
  float comndir[ 3 ] )
{
	// int plane0 = reserverlist[ 0 ];	
	int pickind = 0;
	
	float tmin;
	if( abs( comndir[ 0 ]) < abs( comndir[ 1 ]))
	{
		tmin = abs( comndir[ 0 ]);
		pickind = 0;
	}
	else
	{
		tmin = abs( comndir[ 1 ]);
		pickind = 1;
	}
	if( abs( comndir[ 2 ] < tmin ) )
	{
		pickind = 2;
	}

	//go through the reserved planes, and flip the normals if necessary
	for( int i = 0; i < reservenum; i ++ )
	{
		int planei = reserverlist[ i ];
		
		if( pparam[ planei * 4 + pickind ] < 0 )	//flip it
		{
			for( int j = 0 ; j < 4;j ++ )
			{
				pparam[ planei * 4 + j ] = -pparam[ planei * 4 + j ];
			}
			
			int edgesize = ctredges[ planei ].size();
			for( int j = 2; j < edgesize; j += 4 )
			{
				int tmat = ctredges[ planei ][ j + 1];
				ctredges[ planei ][ j + 1 ] = ctredges[ planei ] [ j ];
				ctredges[ planei ][ j ] = tmat;
			}
		}
	}
	return pickind;
}

void ContourHandler::sortSelectMarkList
( int planenum, int markval,	bool inc, //the number of the marked planes, and mark value 1 or -1, inc-true, in increasing order
	int*& markreservelist,			//all the mark list
	int*& reservelist, float*& dotProduct, //the reserve plane list, the dotproduct between the first normal and other normals
	int*& sortedplanelist,		//put the sorted marked plane list into the result
	int startind							//the starting ind to put into the array
 )
{
	int* tplanelist = new int[ planenum ];
	float* tdotprodlist = new float[ planenum ];

	int j = 1;
	for( int i = 0; i <planenum; i ++ )	
	{
		while( markreservelist[ j ] != markval )
		{
			j ++;
		}
		tplanelist[ i ]  = reservelist[ j ];
		tdotprodlist[ i ] = dotProduct[ j ];
		j ++;
	}

	//////////////////////////////////////////////////////////////////////////
	/*for( int i = 0 ; i < planenum;  i ++ )
	{
		cout<<tplanelist[ i ]<<" ";
	}
	cout<<endl;*/
	//////////////////////////////////////////////////////////////////////////

	//sort in the decreasing order of the dotproduct for all the mark = 1
	for( int i = 0; i < planenum - 1; i ++ )
	{
		int mini = i;				//i don't want to change the name, maybe not min, can be max too, depends on inc
		float minval = tdotprodlist[ i ];

		for( int j = i + 1; j < planenum; j ++ )
		{
			if(( inc && tdotprodlist[ j ] < minval ) ||
				( !inc && tdotprodlist[ j ] > minval ))
			{
				minval = tdotprodlist[ j ];
				mini = j;
			}
		}

		if( mini == i )
			continue;

		int ti = tplanelist[ mini ];
		tplanelist[ mini ] = tplanelist[ i ];
		tplanelist[ i ] = ti;

		float tprod = tdotprodlist[ mini ];
		tdotprodlist[ mini ] = tdotprodlist[ i ];
		tdotprodlist[ i ] = tprod;
	}

	//copy the sorted neglist out
	for( int i = 0; i < planenum; i ++ )
	{
		sortedplanelist[ startind + i ] = tplanelist[ i ];
	}

	delete []tplanelist;
	delete []tdotprodlist;
}

void removeTooCloseCtrVeres
(
 intvector& vermark,
 floatvector& ctrvers, 
 intvector& ctredges	
 )
{
	int oldvernum = vermark.size();
	int oldedgenum = ctredges.size()/4;

	intvector interpos;
	int* vermap = new int[ oldvernum ];

	//step1, find all the points on the intersection line
	for( int i = 0; i < oldvernum; i ++ )
	{
		if ( vermark[ i ] == 0 )
			interpos.push_back( i );

		vermap[ i ] = -1;
	}

	//////////////////////////////////////////////////////////////////////////
	/*cout<<"interpt pos indx:";
	for( int i = 0; i < interpos.size(); i++ )
	{
		cout<<interpos[ i ]<<" ";
	}
	cout<<endl;*/
	//////////////////////////////////////////////////////////////////////////

	//step2, go through each intersection point
	int interptnum = interpos.size();

	intvector getridind;

	for( int i= 0; i < interptnum; i ++ )
	{
		int interpti = interpos[ i ];

		//step2.1   find at most the left and right 3 secutive intersection points
		//step2.2   push these indices of the neighbors into a same vector
					//set the map information of these vertices disappearVer2
		int lastedges[ 2 ] = {-1, -1};		
		int curvers[ 2 ] = {0, 0};

		//find the first two which touches the two vertices
		for( int j = 0; j < oldedgenum; j ++ )
		{
			int veris[ 2 ] = { ctredges[ 4*j ], ctredges[ 4*j + 1]};
			for( int k = 0; k < 2; k ++)
			{
				if( veris[ k ] == interpti )
				{
					if( lastedges[ 0 ] == -1 )	//first on hasn't been found
					{
						lastedges[ 0 ] = j;
						curvers[ 0 ] = veris[ 1 - k ];
						break;
					}

					lastedges[ 1 ] = j;
					curvers[ 1 ] = veris[ 1 - k ];
					break;
				}		
			}				
		}

		//go through all the edges, to find the edge, which includes the curvex, and are not lastedges
		bool shouldcontinue[ 2 ] ={true, true };
		int count[ 2 ] = {0, 0};
		float vec[ 3 ];
		
		while( true )
		{
			count[ 0 ] ++;
			count[ 1 ] ++;

			//check if should continue
			for( int j = 0; j < 2; j ++ )
			{
				if( !shouldcontinue[ j ])
					continue;

				if( count[ j ] >= 5 )
				{
					shouldcontinue[ j ] = false;
					continue;
				}

				MyMath::getVec( &ctrvers[ 3*interpti ], &ctrvers[ 3*curvers[ j ]], vec);			

				if( MyMath::vectorlen( vec ) >= EXPECTLEN )
					shouldcontinue[ j ] = false;
				else
				{					
					getridind.push_back( curvers[ j ] );
					vermap[ curvers[ j ]] = interpti;
				}
			}

			if( !shouldcontinue[ 0 ] && !shouldcontinue[ 1 ])
				break;

			//go through the edges to find the next vertex
			bool twosidesset[ 2 ] = {!shouldcontinue[ 0 ], !shouldcontinue[ 1 ]};
			
			for( int j = 0; j < oldedgenum; j ++ )
			{
				int veris[ 2 ] = { ctredges[ 4*j ], ctredges[ 4* j + 1 ]};
				
				//check the two curvers
				for( int k = 0; k < 2; k ++ )	//k for the two sides vertices
				{
					if( twosidesset[ k ] )					
						continue;
					

					for( int k1 = 0; k1 < 2; k1++ )	//k1 for the two vertices of the edge
					{
						if( lastedges[ k ] == j )	//this edge has already been pushed in
							continue;

						if( veris[ k1 ] == curvers[ k ] )
						{
							curvers[ k ] = veris[ 1-k1 ];
							lastedges[ k ]  = j;
							twosidesset[ k ] = true;
							break;	//found!
						}
					}
				}

				if( twosidesset[ 0 ] && twosidesset[ 1 ])	//both have been set!
					break;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	/*cout<<"to remove index:"<<endl;
	for( int i = 0; i < getridind.size(); i++)
	{
		cout<<getridind[ i ]<<" ";
	}
	cout<<endl;
	cout<<"vermap:"<<endl;
	for( int  i = 0; i < oldvernum; i++ )
	{
		cout<<vermap[ i ]<<" ";
	}
	cout<<endl;*/
	//////////////////////////////////////////////////////////////////////////

	//step 3. sort these intersection points,
				//set the map for all the points, points other than intersection point will be got rid of, 
				//all the indices of the vertices( = points)  are remapped.
	int getridindsize = getridind.size();

	if( getridindsize == 0 )
		return;
	
	//select sort
	for( int i = 0; i < getridindsize - 1; i ++ )
	{
		int minval = getridind[ i ];
		int mini = i;
		for( int j = i + 1; j < getridindsize; j ++  )
		{
			if( getridind[ j ] < minval )
			{
				minval = getridind[ j ];
				mini = j;
			}
		}
		
		if( mini == i )
			continue;

		getridind[ mini ] = getridind[ i ];
		getridind[ i ] = minval;
	}

    //set the map
	//set the first segment before the first getting rid index
	for( int i = 0; i < getridind[ 0 ] ; i ++ )
	{
		vermap[ i ] = i;
	}
	int reduceamount = 1;

	for( int i = 0; i < getridindsize - 1; i ++ )
	{
		int start = getridind[ i ] + 1;
		int end = getridind[ i + 1 ] - 1;
		for( int j = start; j <= end; j ++ )
		{
			vermap[ j ] = j - reduceamount;
		}

		reduceamount ++;
	}
	for( int i = getridind[ getridindsize - 1] + 1; i < oldvernum; i ++ )
	{
		vermap[ i ] = i - reduceamount;
	}

	//step 4,	after remapping, go through all the getting rid of indices, and reset the map, they should be mapped
		//to the new index of the intersection point
	for( int i = 0; i < getridindsize; i ++ )
	{
		int indi = getridind[ i ];			//the index of get rid of index

		int oldinterpti = vermap[ indi ];	//the old position of the intersection point
		
		vermap[ indi ] = vermap[ oldinterpti ];	//set to the new postion of the intersection point
	}

	//step 5. get rid of those redundant points
	floatvector tvec ;
	tvec.resize( ctrvers.size());
	memcpy( &tvec[ 0 ] , &ctrvers[ 0 ] , oldvernum * 3 * sizeof( float ));

	int cursize = 0;
	ctrvers.clear();
	ctrvers.resize( oldvernum * 3 - getridindsize * 3 );
	int start = 0;
	int end;

	//vermark
	intvector tvec2;
	tvec2.resize( oldvernum );
	memcpy( &tvec2[ 0 ], &vermark[ 0 ], sizeof(int) * oldvernum );
	int cursize2 = 0;
	vermark.clear();
	vermark.resize( oldvernum - getridindsize );

    for( int i = 0; i < getridindsize; i ++ )
	{
		end = getridind[ i ] - 1;
		if ( start > end )
		{
			start = end + 2;
			continue;
		}

		memcpy( &ctrvers[ cursize ], &tvec[ 3 * start ], sizeof(float ) * ( end - start + 1) * 3);
		memcpy( &vermark[ cursize2 ], &tvec2[  start ], sizeof(int ) * ( end - start + 1 ));

		cursize = cursize + 3 * ( end - start + 1 );
		cursize2 = cursize2 + end - start + 1;

		start = end + 2;
	}
	//the last segment
	end = oldvernum - 1;
	if( end >= start)
	{
		memcpy( &ctrvers[ cursize ], &tvec[ 3 * start ], sizeof(float ) * ( end - start + 1) * 3);
		memcpy( &vermark[ cursize2 ], &tvec2[  start ], sizeof(int ) * ( end - start + 1 ));
	}

	tvec.clear();
	tvec2.clear();

	//step 6. for all the edges, reset the map information, if v1 and v2 are same, get rid of it, otherwise, 
	//reset the vertices of the indices
	tvec2.resize( oldedgenum * 4 );
	memcpy( &tvec2[ 0 ] , &ctredges[ 0 ] , sizeof( int ) * oldedgenum * 4);

	ctredges.clear();
	for( int i = 0; i < oldedgenum; i ++ )
	{
		int veris[ 2] = { vermap[ tvec2[ 4*i ] ] , vermap [ tvec2[ 4* i + 1] ] };

		if( veris[ 0 ] ==  veris[ 1 ])
			continue;

		ctredges.push_back( veris[ 0 ] );
		ctredges.push_back( veris[ 1 ] );
		ctredges.push_back( tvec2[ 4*i + 2 ]);
		ctredges.push_back( tvec2[ 4*i + 3 ]);
	}
	tvec2.clear();

	//////////////////////////////////////////////////////////////////////////
	/*for( int i = 0; i< oldvernum; i++)
	{
		cout<<vermap[ i ]<<" ";
	}
	cout<<endl;
	
	cout<<"vermark.size:"<<vermark.size()<<endl;
	for( int i =0 ; i < vermark.size(); i++ )
		cout<<vermark[ i ]<<" ";
	cout<<endl;*/
	//////////////////////////////////////////////////////////////////////////
	
	interpos.clear();
	delete []vermap;
	getridind.clear();	
}

//sort the reserved planes, and save them into the result
void ContourHandler::sortReservedPlanes
( int& planenum, float*& pparam,
 int*& pctrvernum, float**& pctrvers, int**& pctrvermark,
 int*& pctredgenum, int**& pctredges,vector<intvector>& vermark,
 vector<floatvector>& ctrvers, vector<intvector>& ctredges,	//
 int* reservelist, int reservenum, int pickind,
 float comndir[ 3 ] )
{
	//divide the planes into two types, on left side of the first normal, and on the right side of the normal
	int plane0 = reservelist[ 0 ];
	float normal0[ 3 ];
	memcpy( normal0, pparam + 4* plane0 , sizeof( float ) * 3 );

	int* markreservelist = new int[ reservenum ];
	float* dotProduct = new float[ reservenum ];
	int negnum , posnum;
	negnum = posnum = 0;
	for( int i = 1; i < reservenum; i ++ )
	{
		//compute the crossproduct with normal0 and current normal
		float dir[ 3 ];
		MyMath::crossProductNotNorm(normal0, pparam + 4*reservelist[ i ], dir );

		//compute dot product between the cross product and the comndir
		float dotprod = MyMath::dotProduct( dir, comndir );

		//set the mark
		if( dotprod > 0 )
		{
			markreservelist[ i ] = 1;
			posnum ++;
		}
		else
		{
			markreservelist[ i ] = -1;
			negnum ++;
		}

		//compute the dotproduct
		dotProduct[ i ] = MyMath::dotProduct(normal0, pparam + 4*reservelist[ i ]);
	}

	int* sortedplanelist = new int[ reservenum ];

	//////////////////////////////////////////////////////////////////////////
	/*for( int i = 0; i < reservenum ; i ++ )
		cout<<markreservelist[ i ]<<" ";
	cout<<endl;
	cout<<"posnum"<<posnum<<endl;*/
	//////////////////////////////////////////////////////////////////////////
	//sort in the increasing order of the dotproduct for all the mark = -1
	if( negnum != 0 )
	{
		sortSelectMarkList( negnum, -1, true,  markreservelist, reservelist, dotProduct,sortedplanelist,	0);
	}	
	
	sortedplanelist[ negnum ] = plane0;

	//sort in the decaresing order for those dotproduct with mark = 1
	if( posnum != 0 )
	{
		sortSelectMarkList( posnum, 1, false,  markreservelist, reservelist, dotProduct,sortedplanelist,	negnum + 1);
	}

	//////////////////////////////////////////////////////////////////////////
	/*cout<<"The sorted list is :";
	for(int  i = 0; i < reservenum; i ++ )
	{
		cout<<sortedplanelist[ i ]<<" ";
	}
	cout<<endl;*/
	//////////////////////////////////////////////////////////////////////////
	
	//remove some contour vertices
	for( int i = 0 ;i < reservenum; i ++ )
	{
		int planei = sortedplanelist[ i ];
		removeTooCloseCtrVeres(	vermark[ sortedplanelist[ planei ] ], ctrvers[ sortedplanelist[ planei ]] , 
			ctredges[ sortedplanelist[ planei ]]);

		/*writePlaneCtr( i, ctrvers[ sortedplanelist[ planei ]] , vermark[ sortedplanelist[ planei ] ],
			ctredges[ sortedplanelist[ planei ]]);*/
	}

	//planenum, pparam
	int oldplanenum = planenum;
	planenum = reservenum ;
	float* tpointer = pparam;
	pparam = NULL;
	pparam = new float[ planenum * 4 + 24] ;
	for( int i = 0; i  < reservenum; i ++ )
	{
		int tplanei = sortedplanelist[ i ];
		memcpy(pparam + 4*i, tpointer + 4* tplanei, sizeof(float) * 4 );
	}
	memcpy( pparam + planenum * 4,  tpointer + 4 * oldplanenum, sizeof ( float ) * 4 );
	delete []tpointer;

	//int*& pctrvernum, float**& pctrvers; 
	//int*& pctredgenum, int**& pctredges;
	delete []pctrvernum;
	for( int i = 0; i < oldplanenum; i ++ )
		delete []pctrvers[ i ];
	delete []pctrvers;
	
	delete []pctredgenum;
	for( int i = 0; i < oldplanenum; i ++ )
	{
		delete []pctredges[ i ];
	}
	delete []pctredges;

	pctrvernum = new int [ planenum ];
	pctredgenum = new int[ planenum ];
	pctrvermark = new int*[planenum];
	pctrvers = new float*[ planenum] ;
	pctredges = new int*[ planenum ];
	for( int i = 0; i < planenum; i ++ )
	{
		int tplanei = sortedplanelist[ i ];
		int versize = ctrvers[ tplanei ].size();
		pctrvernum[ i ] = versize/3;
		pctrvers[ i ] = new float[ versize ];
		for( int j = 0; j < versize; j ++ )
		{
			pctrvers[ i ][ j ] = ctrvers[ tplanei ][ j ];
		}
		versize /= 3;
		pctrvermark[ i  ] = new int[ versize ];
		for( int j = 0;  j < versize; j++ )
		{
			pctrvermark[ i ][ j ] = vermark[ tplanei ][ j ];
		}

		int edgesize = ctredges[ tplanei ].size();
		pctredgenum[ i ] = edgesize/4;
		pctredges[ i ] = new int[ edgesize ];
		for( int j = 0; j < edgesize; j ++ )
		{
			pctredges[ i ][ j ] = ctredges[ tplanei ][ j ];
		}
	}
	
	//////////////////////////////////////////////////////////////////////////
	/*for( int i = 0; i < planenum + 4; i ++ )
	{
		cout<<pparam[ 4* i ]<<","
			<<pparam[ 4*i + 1]<<","
			<<pparam[ 4*i + 2]<<","
			<<pparam[ 4*i + 3]<<endl;
	}
	cout<<endl;*/
	//////////////////////////////////////////////////////////////////////////

	delete []markreservelist;
	delete []dotProduct;
	delete []sortedplanelist;
}

/*---------------------------------
make the contours on the planes sharing one line consistent
input and also output
planenum: the real number of the cutting plane, not including the added six bounding planes
pparam: the parameters of the planes, including the six bounding planes
pctrvernum: contour vertices numbers on the cuttong planes
pctrvers: the contour vertices positions, three is a group, for one vertex, x y and z
pctredgenum: the contour edges on the planes, four is a group, first vertex index, second vertex index, left material, right material

---------------------------------*/
int** ContourHandler::makeConsistentSingleMatCommonLine(
	int& planenum, float*& pparam,
	int*& pctrvernum, float**& pctrvers,
	int*& pctredgenum, int**& pctredges,
	
	float comndir[ 3 ], float comnpt[ 3 ])
{
	//initialization
	int** pctrvermark;
	int tplanenum;	
	vector<floatvector> tctrvers;
	vector<intvector> tctredges;
	tctrvers.resize( planenum );
	tctredges.resize( planenum );
	vector<intvector> tctrvermark;	//-1 on one side of the common line, 0 on the common line 1 on the other side
	tctrvermark.resize( planenum );

	tplanenum = planenum ;
	for( int i = 0; i < tplanenum; i ++ )
	{
		tctrvers[ i ].resize( pctrvernum[ i ] * 3);
		memcpy( &tctrvers[ i ][ 0 ], pctrvers[ i ], pctrvernum[ i ] * 3 * sizeof(float));

		tctredges[ i ].resize( pctredgenum[ i ] * 4 );
		memcpy( &tctredges[ i ][ 0 ], pctredges[ i ], pctredgenum[ i ] * 4 * sizeof(int));

		tctrvermark[ i ].resize( pctrvernum[ i ]);
	}

	//step1, split the contour edges by the common line, compute intersection points on the common line
	vector<intvector> interptpos;
	interptpos.resize( tplanenum );
	for( int i = 0; i < tplanenum; i ++ )
	{
		interCommonLineWithOnePlane(pparam + 4 * i, tctrvers[ i ], tctrvermark[ i ],tctredges[ i ], comndir, comnpt, interptpos[ i ]);
	}
	//////////////////////////////////////////////////////////////////////////
	//writeContoursOut		(		pparam,		tctrvers,		tctredges,		interptpos,		comndir, comnpt		);
	//////////////////////////////////////////////////////////////////////////

	//step2, reserve the ones that have the same number of intersection points, and 
	//adjust the position of those intersection points, making them lying on the common position
	int* reservelist;
	int reservenum = 0;
	reserverPlanes( tplanenum, interptpos,  reservelist, reservenum );
	if ( reservelist == NULL ) //all the planes are discarded due to some reason
	{
		cout<<"Not able to make the contours on the planes consistent! Fail!"<<endl;
		clearTheInputData(planenum, pparam, pctrvernum, pctrvers,
			pctredgenum, pctredges);
		//clear temporary data
		for( int i = 0; i< tplanenum; i++)
		{
			tctrvers[ i].clear();
			tctredges[ i ].clear();
			tctrvermark[ i ].clear();
			interptpos[ i ].clear();
		}
		tctrvers.clear();
		tctredges.clear();
		tctrvermark.clear();
		interptpos.clear();
		return NULL;
	}
	//data is partially good
	//ouptut the discard ones
	for( int i = 0; i < reservenum; i ++ )
	{
		int end ;		
		if( i != reservenum -1 )
			end = reservelist[ i + 1]  - 1;
		else
			end = tplanenum - 1;
		for( int j = reservelist[ i ] + 1; j <= end; j ++)
			cout<<"Plane "<< j <<" is discarded due to the intersection points on the common edge is not right! "<<endl;
	}

	//////////////////////////////////////////////////////////////////////////
	/*writeInterPtsOnCommonLine(
		 reservelist,  reservenum,   interptpos,	 tctrvers  );*/
	//////////////////////////////////////////////////////////////////////////

	//adjust the position of the intersection points
	if( !adjustInterPt(reservelist, reservenum, interptpos, tctrvers,comndir, comnpt))
	{
		cout<<"FAIL! Unable to make the intersection points consistent with each other!"<<endl;
		clearTheInputData(planenum, pparam, pctrvernum, pctrvers,
			pctredgenum, pctredges);
		//clear temporary data
		for( int i = 0; i< tplanenum; i++)
		{
			tctrvers[ i].clear();
			tctredges[ i ].clear();
			tctrvermark[ i ].clear();
			interptpos[ i ].clear();
		}
		tctrvers.clear();
		tctredges.clear();
		tctrvermark.clear();
		interptpos.clear();
		return NULL;
	}

	//step3, set the parameters of the planes, so that the normal of them are pointing to the same half space
	int pickind = flipNormalToHalfSpace	( reservelist, reservenum,	pparam, tctredges,comndir );

	//////////////////////////////////////////////////////////////////////////
	//writeParamOut( reservelist,  reservenum, pparam  , pickind);
	//////////////////////////////////////////////////////////////////////////

	//step4, sort all the cutting planes, so that it is easier for partitioning	
	sortReservedPlanes(planenum, pparam,
		pctrvernum, pctrvers, pctrvermark, pctredgenum, pctredges, tctrvermark,
		tctrvers,tctredges,reservelist, reservenum, pickind,	comndir );

	//////////////////////////////////////////////////////////////////////////
	//writeParamOut( reservelist,  reservenum, pparam  , pickind);
	//////////////////////////////////////////////////////////////////////////

	//release the temporary memory
	for( int i = 0; i< tplanenum; i++)
	{
		tctrvers[ i].clear();
		tctredges[ i ].clear();
		tctrvermark[ i ].clear();
		interptpos[ i ].clear();
	}
	tctrvers.clear();
	tctredges.clear();
	tctrvermark.clear();
	interptpos.clear();
	delete reservelist;

	return pctrvermark;
}

//for debug reason
void ContourHandler::writeInterPtsOnCommonLine(
									 int* reservelist, int reservenum,
									  vector<intvector>& interptpos,	  vector<floatvector>& ctrvers  )
{
	FILE* fout = fopen( "mmdebug/interpt.txt","w");
	
	fprintf( fout, "{");

	for( int i = 0 ;  i < reservenum; i ++ )
	{
		int planei = reservelist[ i ];

		fprintf( fout,"{");
		int interptpossize = interptpos[ planei ].size();
		for( int j = 0; j < interptpossize ;  j++ )
		{
			int veri = interptpos[ planei ][ j ];
			fprintf( fout, "{%f,%f,%f}", ctrvers[ planei ][ 3*veri  ], ctrvers[ planei ][ 3*veri + 1 ], ctrvers[ planei ][ 3*veri + 2]);

			if( j != interptpossize - 1 )
				fprintf( fout, ",");
		}
		if( i != reservenum - 1 )
			fprintf( fout ,"},");
		else
			fprintf(fout, "}");
	}
	fprintf(fout ,"}");
	fclose( fout );
}

void ContourHandler::writeContoursOut
(
 float* param,
 vector<floatvector>& ctrvers,
 vector<intvector>& ctredges,
 vector<intvector> interptpos,
 float comndir[ 3 ], float comnpt[ 3 ]
 )
 {
	 //write out the common line parameters
	 //write out the param
	 //write out the vertices
	 //write out the edges
	 //write out the intersection points

	 FILE* fout= fopen( "mmdebug/ctr.txt", "w");

	 fprintf( fout ,"{");
	
	 int planenum = ctrvers.size();

	 //write out the common line parameters
	 fprintf( fout, "{%f,%f,%f},", comndir[ 0 ], comndir[ 1 ] , comndir[ 2 ]);
	 fprintf( fout, "{%f,%f,%f},", comnpt[ 0 ], comnpt[ 1 ], comnpt[ 2 ]);

	 //write out the param
	 fprintf(fout, "{");
	for( int i = 0; i < planenum ; i++ )
	{
		fprintf(fout, "{%f,%f,%f,%f}", param[ 4*i ], param[ 4*i+ 1], param[ 4*i +2], param[ 4*i + 3]);
		if( i != planenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},");
	}

	 //write out the vertices
	fprintf(fout, "{");
	for( int i = 0; i < planenum; i ++ )
	{
		fprintf(fout, "{");
		int vernum  = ctrvers[ i ].size()/3;
		for ( int j= 0; j < vernum ; j ++ )
		{
			fprintf(fout, "{%f,%f,%f}", ctrvers[ i ][ 3*j], ctrvers[ i ][ 3*j + 1 ] , ctrvers[ i ][ 3*j + 2 ]);
			if ( j != vernum - 1)
				fprintf(fout, ",");
			else
				fprintf(fout, "}");
		}
		if( i != planenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "},");
	}

	 //write out the edges
	fprintf(fout, "{");
	for( int i =0; i <planenum; i ++ )
	{
		fprintf(fout, "{");
		int edgenum = ctredges[ i ].size() / 4; 
		for( int j = 0; j < edgenum; j ++ )
		{
			fprintf(fout, "{{%d,%d},{%d,%d}}", ctredges[ i ][ 4*j ] + 1, ctredges[ i ][ 4*j +1] + 1, 
				ctredges[ i ][ 4*j + 2 ] + 1, ctredges[ i ][ 4*j + 3 ] + 1);
			if( j != edgenum - 1 )
				fprintf(fout, ",");
			else
				fprintf(fout, "}");		
		}
		if( i != planenum - 1)
			fprintf(fout, ",");
		else
			fprintf(fout, "},");
	}

	 //write out the intersection points
	fprintf(fout, "{");
	for( int i =0 ;i < planenum; i ++ )
	{
		fprintf(fout, "{");
		int interptnum = interptpos[ i ].size();
		for( int j= 0; j < interptnum; j ++ )
		{
			fprintf(fout, "%d", interptpos[ i  ][ j ] + 1);
			if( j != interptnum - 1 )
				fprintf(fout, ",");
			else
				fprintf(fout, "}");
		}
		if( i != planenum - 1 )
			fprintf(fout, ",");
		else
			fprintf(fout, "}");
	}

	 fprintf( fout, "}" );
	 fclose( fout  );
 }

 void ContourHandler::writeParamOut(int* reservelist, int reservenum, float* param  , int pickind)
 {
	 FILE* fout = fopen( "mmdebug/planeparam.txt", "w");
	 fprintf(fout, "{");
	 fprintf(fout, "%d,", pickind );
	 fprintf(fout, "{");

	 for( int i  = 0; i < reservenum  ; i ++ )
	 {
		 int planei = reservelist[ i ];
		 fprintf(fout, "{%f,%f,%f,%f}", param[ 4*planei ], param[ 4*planei + 1 ], param[ 4*planei +2],
			 param[ planei * 4 + 3 ]);
		 if( i != reservenum - 1 )
			 fprintf(fout, ",");
		 else
			 fprintf(fout, "}");
	 }
	 fprintf(fout, "}");
	 fclose( fout  );
  }
