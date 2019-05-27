#include "../ContourHandler/ContourHandler.h"

//void ContourHandler::computeInterPtList(
//floatvector param,
//vector<floatvector>& ctrvers,
//vector<intvector>& ctredges,
//int becut,	//the index of the plane to be cut by the other plane
//int tocut,	//the index of the plane that is going to cut another plane
//floatvector& interptpos,	//result intersection points list
//intvector& interptind
//)	
//{	
//	/////////////////////////////////////////////////////////////////////////
//	//int tbecut = becut;
//	//int ttocut = tocut;
//	//////////////////////////////////////////////////////////////////////////
//	//step1. mark all the vertices
//	int tvernum = ctrvers[ becut ].size()/3;
//	int* tvermark = new int[ tvernum ];
//	float* tctrvers = new float[ tvernum * 3 ];
//	float tparam[ 4 ];
//	for( int i = 0; i < 4; i ++)
//		tparam[ i ] = param[ 4*tocut + i ];
//	for( int i = 0; i < tvernum*3 ; i++ )
//	{
//		tctrvers[ i ] = ctrvers[ becut ][ i ];
//	}
//	//////////////////////////////////////////////////////////////////////////
//	//cout<<"tocut plane "<<ttocut<<": param:"<<tparam[ 0 ]<<" " <<tparam[ 1 ]<<" "<<tparam[ 2 ]<<tparam[ 3 ]<<endl;
//	//cout<<"the fist vertex:"<< tctrvers[ 0 ]<<" "<<tctrvers[ 1 ]<<" "<<tctrvers[2]<<endl;
//	////////////////////////////////////////////////////////////////////////////
//
//	for( int i = 0 ;i < tvernum; i ++)
//	{
//		//because it is resized to 1000, and the plane normal is unit normal,
//		//if distance is 0.05, it is though as lying on the plane
//		float dist = MyMath::dotProduct( tctrvers + 3*i, tparam) - tparam[ 3 ];
//        if (MyMath::isEqualInToler(dist, 0, 0.05))
//		{
//			tvermark[ i ] = 0;
//			interptpos.push_back( tctrvers[ 3*i ]);
//			interptpos.push_back( tctrvers[ 3*i + 1 ] );
//			interptpos.push_back( tctrvers[ 3*i + 2 ]);
//			interptind.push_back( i );
//		}
//		else if( dist > 0 )
//			tvermark[ i ] = 1;
//		else
//			tvermark[ i ] = -1;
//	}
//
//	//step2. go through each edge, split it when needed
//	int tedgenum = ctredges[ becut ].size()/4;
//	int tvers[ 2 ];
//	float interpt[ 3 ];
//	for( int i = 0; i< tedgenum; i ++ )
//	{
//		tvers[ 0 ] = ctredges[ becut ][ 4*i ];
//		tvers[ 1 ] = ctredges[ becut ][ 4*i + 1 ];
//
//		if( tvermark[tvers[ 0 ]] * tvermark[tvers[ 1 ]] == -1 )	//intersection point exists
//		{
//			interPt_PlaneEdge( tparam, tctrvers+3*tvers[ 0 ], tctrvers + 3*tvers[ 1 ], interpt);
//			//////////////////////////////////////////////////////////////////////////
//		/*	printf("{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}", 
//				tctrvers[ 3*tvers[ 0 ]],
//				tctrvers[ 3*tvers[ 0 ] + 1 ],
//				tctrvers[ 3*tvers[ 0 ] + 2],
//				interpt[ 0 ],
//				interpt[ 1 ],
//				interpt[ 2 ],
//				tctrvers[ 3*tvers[ 1 ]],
//				tctrvers[ 3*tvers[ 1 ] + 1 ],
//				tctrvers[ 3*tvers[ 1 ] + 2]
//				);*/
//			//////////////////////////////////////////////////////////////////////////
//			interptind.push_back( tvernum );
//			ctrvers[ becut].push_back( interpt[ 0 ]);
//			ctrvers[ becut ].push_back( interpt[ 1 ]);
//			ctrvers[ becut ].push_back( interpt[ 2 ]);
//			ctredges[ becut ].push_back( tvernum );
//			ctredges[ becut ].push_back( tvers[ 1 ]);
//			ctredges[ becut ].push_back( ctredges[ becut ][ 4*i + 2 ]);
//			ctredges[ becut ].push_back( ctredges[ becut ][ 4*i + 3 ]);
//			ctredges[ becut ][ 4*i + 1 ] = tvernum;
//			tvernum ++;
//			interptpos.push_back( interpt[ 0 ]);
//			interptpos.push_back( interpt[  1 ]);
//			interptpos.push_back( interpt[ 2 ]);			
//		}
//	}
//
//	delete []tvermark;
//	delete []tctrvers;
//}
//bool ContourHandler::makeConsistent_InterPt_2Planes(
//	floatvector param,
//	vector<floatvector>& ctrvers,
//	vector<intvector>& ctredges,
//	//vector<intvector>& v2plistpos,
//	//vector<intvector>& v2plist,
//	int referplane,
//	int adjustplane,
//	vector<intvector>& pairlist
//	)
//{
//	//////////////////////////////////////////////////////////////////////////
//	/*int trerferplane = referplane;
//	int tadjustplane = adjustplane;*/
//	//////////////////////////////////////////////////////////////////////////
//	intvector interptind1;
//	intvector interptind2;
//	floatvector interptpos1;
//	floatvector interptpos2;
//	int len;
//
//	//step1.process the reference plane: get the intersection points
//	computeInterPtList( param, ctrvers, ctredges, referplane, adjustplane, interptpos1, interptind1);
//	len = interptind1.size();
//
//	//step2.process the new plane, compute the intersection points
//	computeInterPtList( param, ctrvers, ctredges, adjustplane, referplane, interptpos2, interptind2);
//	if( interptind2.size() != len )
//	{
//		interptind1.clear();
//		interptind2.clear();
//		interptpos1.clear();
//		interptpos2.clear();
//	//	pairlist.resize( pairlist.size() + 1 );
//		return false;
//	}
//
//	if( len == 0 )	//no intersection points at all
//	{
//		pairlist.resize( pairlist.size() + 1 );
//		return true;
//	}
//    //step3. compute the parameters of these intersection points, compute the parameters for them and sort them out
//
//	//////////////////////////////////////////////////////////////////////////
//	/*cout<<"========================"<<endl;
//	for( int i= 0 ; i < interptpos1.size()/3; i ++)
//	{
//		cout<<interptind1[ i ]<<"\t"<<interptpos1[ 3*i ]<<"\t"<<interptpos1[ 3*i + 1 ]<<"\t"<<interptpos1[ 3*i +2 ]<<endl;
//	}
//	cout<<"----------"<<endl;
//	for( int i= 0 ; i < interptpos2.size()/3; i ++)
//	{
//		cout<<interptind2[ i ]<<"\t"<<interptpos2[ 3*i ]<<"\t"<<interptpos2[ 3*i + 1 ]<<"\t"<<interptpos2[ 3*i +2 ]<<endl;
//	}
//	cout<<"========================"<<endl;*/
//	//////////////////////////////////////////////////////////////////////////
//	float dir[ 3 ];
//	float param1[ 4 ], param2[ 4 ];
//	for( int i = 0; i < 4; i++)
//	{
//		param1[ i ] = param[ 4*referplane + i ];
//		param2[ i ] = param[ 4*adjustplane + i ];
//	}
//	MyMath::crossProduct( param1, param2, dir);
//
//	float startpt[ 3 ];
//	float curpt[ 3 ];
//	startpt[ 0 ] = interptpos1[ 0 ];
//	startpt[ 1 ] = interptpos1[ 1 ];
//	startpt[ 2 ] = interptpos1[ 2 ];
//	float* interptparam1 = new float[ len ];
//	float* interptparam2 = new float[ len  ];
//	for( int i = 0; i < len; i ++ )
//	{
//		for( int j = 0 ; j < 3; j ++)
//			curpt[ j ] = interptpos1[  i*3 + j ];
//		interptparam1[ i ] = getTOnRay( startpt, dir, curpt);
//		for( int j = 0; j < 3; j ++)
//		{
//			curpt[ j ] = interptpos2[ i *3 + j];
//		}
//		interptparam2[ i ] = getTOnRay( startpt, dir, curpt );
//	}
//	
//	//sort the points accroding to their parameters
//	for( int i = 0 ; i < len ; i ++)
//	{
//		int minpos = i;
//		for( int j = i + 1; j < len ;j ++)
//		{
//			if( interptparam1[ j ] < interptparam1[ minpos ])
//			{
//				minpos = j;
//			}
//		}
//		if( minpos == i )
//			continue;
//		//switch the two values there
//		//param
//		float t = interptparam1[ minpos ];
//		interptparam1[ minpos ] = interptparam1[ i ];
//		interptparam1[ i ] = t;
//		//position
//		for( int k = 0; k < 3; k ++)
//		{
//			t = interptpos1[ 3*minpos + k ];
//			interptpos1[ 3*minpos + k ] = interptpos1[ 3*i + k ];
//			interptpos1[ 3*i + k ] = t;
//		}
//		//ind
//		int it = interptind1[ minpos ];
//		interptind1[ minpos ] = interptind1[ i ];
//		interptind1[ i ] = it;
//	}
//
//	for( int i = 0 ; i < len ; i ++)
//	{
//		int minpos = i;
//		for( int j = i + 1; j < len ;j ++)
//		{
//			if( interptparam2[ j ] < interptparam2[ minpos ])
//			{
//				minpos = j;
//			}
//		}
//		if( minpos == i )
//			continue;
//		//switch the two values there
//		//param
//		float t = interptparam2[ minpos ];
//		interptparam2[ minpos ] = interptparam2[ i ];
//		interptparam2[ i ] = t;
//		//position
//		for( int k = 0; k < 3; k ++)
//		{
//			t = interptpos2[ 3*minpos + k ];
//			interptpos2[ 3*minpos + k ] = interptpos2[ 3*i + k ];
//			interptpos2[ 3*i + k ] = t;
//		}
//		//ind
//		int it = interptind2[ minpos ];
//		interptind2[ minpos ] = interptind2[ i ];
//		interptind2[ i ] = it;
//	}
//
//	//check if the parameters are close enough
//	for( int i = 0; i < len ; i ++)
//	{
//		if(MyMath::isEqualInToler( interptparam2[ i ] , interptparam1[ i ], 1 ) )	//the distance between them is 1
//		{
//			//adjust the adjustcon
//			for( int k = 0; k < 3; k ++)
//				ctrvers[ adjustplane ][ 3*interptind2[ i ] + k ] = interptpos1[ i*3 + k ];
//		}
//		else	//not correct return false
//		{
//			interptind1.clear();
//			interptind2.clear();
//			interptpos1.clear();
//			interptpos2.clear();
//			delete []interptparam1;
//			delete []interptparam2;
//		//	pairlist.resize( pairlist.size() + 1 );
//			return false;
//		}
//	}
//
//	//step3.set the pairlist
//	int pairlisti = pairlist.size();
//	pairlist.resize( pairlisti + 1 );
//	for( int i = 0; i < len ;i ++)
//	{
//		pairlist[ pairlisti ].push_back( interptind1[ i ]);
//		pairlist[ pairlisti ].push_back( interptind2[ i ]);
//	}
//
//	interptind1.clear();
//	interptind2.clear();
//	interptpos1.clear();
//	interptpos2.clear();
//	delete []interptparam1;
//	delete []interptparam2;
//	return true;
//}
/*--------------------------------
make the contours on the cutting plane consistent
the following functions assume the contours are only of single material
and only make the intersection point on the common line at the same position
--------------------------------*/
void ContourHandler::computeInterPtList(
        floatvector& param,
        vector<floatvector>& ctrvers,
        vector<intvector>& ctredges,
        int becut,	//the index of the plane to be cut by the other plane
        int tocut,	//the index of the plane that is going to cut another plane
        floatvector& interptpos,	//result intersection points list
        intvector& interptind
        )
{	
    //////////////////////////////////////////////////////////////////////////
    /*cout<<"In computeInterPtList:"<<endl;
    cout<<"param:"<<endl;
    for( int i = 0; i < param.size(); i ++ )
    {
        cout<<param[ i ]<<",";
        if( i +1 % 4 == 0)
            cout<<"\t";
    }
    cout<<endl;*/
    //////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //int tbecut = becut;
    //int ttocut = tocut;
    //////////////////////////////////////////////////////////////////////////
    //step1. mark all the vertices
    int tvernum = ctrvers[ becut ].size()/3;
    float* tverdis = new float[ tvernum ];
    int* tvermark = new int[ tvernum ];
    float* tctrvers = new float[ tvernum * 3 ];
    float tparam[ 4 ];
    for( int i = 0; i < 4; i ++)
        tparam[ i ] = param[ 4*tocut + i ];
    for( int i = 0; i < tvernum*3 ; i++ )
    {
        tctrvers[ i ] = ctrvers[ becut ][ i ];
    }
    //////////////////////////////////////////////////////////////////////////
    //cout<<"tocut plane "<<ttocut<<": param:"<<tparam[ 0 ]<<" " <<tparam[ 1 ]<<" "<<tparam[ 2 ]<<tparam[ 3 ]<<endl;
    //cout<<"the fist vertex:"<< tctrvers[ 0 ]<<" "<<tctrvers[ 1 ]<<" "<<tctrvers[2]<<endl;
    ////////////////////////////////////////////////////////////////////////////

    for( int i = 0 ;i < tvernum; i ++)
    {
        //because it is resized to 1000, and the plane normal is unit normal,
        //if distance is 0.05, it is though as lying on the plane
        float dist = MyMath::dotProduct( tctrvers + 3*i, tparam) - tparam[ 3 ];
        /*if (MyMath::isEqualInToler(dist, 0, 0.05))
        {
            tvermark[ i ] = 0;
            interptpos.push_back( tctrvers[ 3*i ]);
            interptpos.push_back( tctrvers[ 3*i + 1 ] );
            interptpos.push_back( tctrvers[ 3*i + 2 ]);
            interptind.push_back( i );
        }
        else */
        if( dist == 0 )
        {
            tvermark[ i ] = 0;

            interptpos.push_back( tctrvers[ 3*i ]);
            interptpos.push_back( tctrvers[ 3*i + 1 ] );
            interptpos.push_back( tctrvers[ 3*i + 2 ]);
            interptind.push_back( i );
        }
        else if( dist > 0 )
            tvermark[ i ] = 1;
        else
            tvermark[ i ] = -1;
        tverdis[ i ] = dist;
    }

    //////////////////////////////////////////////////////////////////////////
    /*cout<<"tvermark:"<<endl;
    for( int i = 0; i < tvernum; i ++ )
    {
        cout<<tvermark[ i ]<<" ";
    }*/
    //////////////////////////////////////////////////////////////////////////

    //step2. go through each edge, split it when needed
    int tedgenum = ctredges[ becut ].size()/4;
    int tvers[ 2 ];
    //float interpt[ 3 ]; //MING: commented because of unreferenced local variable
    for( int i = 0; i< tedgenum; i ++ )
    {
        tvers[ 0 ] = ctredges[ becut ][ 4*i ];
        tvers[ 1 ] = ctredges[ becut ][ 4*i + 1 ];

        if( tvermark[tvers[ 0 ]] * tvermark[tvers[ 1 ]] == -1 )	//intersection point exists
        {
            //////////////////////////////////////////////////////////////////////////
            //cout<<"two vertices are:"<<tvers[ 0 ]<<","<<tvers[ 1 ]<<endl;
            //////////////////////////////////////////////////////////////////////////
            if( abs(tverdis[ tvers[ 0 ]]) < abs(tverdis[ tvers[ 1 ]]))
            {
                interptpos.push_back( tctrvers[ 3*tvers[ 0 ] ]);
                interptpos.push_back( tctrvers[ 3*tvers[ 0 ] + 1 ] );
                interptpos.push_back( tctrvers[ 3*tvers[ 0 ] + 2 ]);
                interptind.push_back( tvers[ 0 ] );
            }
            else
            {
                interptpos.push_back( tctrvers[ 3*tvers[ 1 ] ]);
                interptpos.push_back( tctrvers[ 3*tvers[ 1 ] + 1 ] );
                interptpos.push_back( tctrvers[ 3*tvers[ 1 ] + 2 ]);
                interptind.push_back( tvers[ 1 ] );
            }
            //interPt_PlaneEdge( tparam, tctrvers+3*tvers[ 0 ], tctrvers + 3*tvers[ 1 ], interpt);
            ////////////////////////////////////////////////////////////////////////////
            ///*	printf("{{%f,%f,%f},{%f,%f,%f},{%f,%f,%f}}",
            //tctrvers[ 3*tvers[ 0 ]],
            //tctrvers[ 3*tvers[ 0 ] + 1 ],
            //tctrvers[ 3*tvers[ 0 ] + 2],
            //interpt[ 0 ],
            //interpt[ 1 ],
            //interpt[ 2 ],
            //tctrvers[ 3*tvers[ 1 ]],
            //tctrvers[ 3*tvers[ 1 ] + 1 ],
            //tctrvers[ 3*tvers[ 1 ] + 2]
            //);*/
            ////////////////////////////////////////////////////////////////////////////
            //interptind.push_back( tvernum );
            //ctrvers[ becut].push_back( interpt[ 0 ]);
            //ctrvers[ becut ].push_back( interpt[ 1 ]);
            //ctrvers[ becut ].push_back( interpt[ 2 ]);
            //ctredges[ becut ].push_back( tvernum );
            //ctredges[ becut ].push_back( tvers[ 1 ]);
            //ctredges[ becut ].push_back( ctredges[ becut ][ 4*i + 2 ]);
            //ctredges[ becut ].push_back( ctredges[ becut ][ 4*i + 3 ]);
            //ctredges[ becut ][ 4*i + 1 ] = tvernum;
            //tvernum ++;
            //interptpos.push_back( interpt[ 0 ]);
            //interptpos.push_back( interpt[  1 ]);
            //interptpos.push_back( interpt[ 2 ]);
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    //writeOneContourOut
    //	( ctrvers[becut], ctredges[becut],	 tvermark, interptind,
    //	 becut,		 tocut,		param		);
    ////////////////////////////////////////////////////////////////////////////

    delete []tvermark;
    delete []tctrvers;
}
//bool ContourHandler::makeConsistent_InterPt_2Planes(
//	floatvector& param,
//	vector<floatvector>& ctrvers,
//	vector<intvector>& ctredges,
//	//vector<intvector>& v2plistpos,
//	//vector<intvector>& v2plist,
//	int referplane,
//	int adjustplane,
//	vector<intvector>& pairlist
//	)
//{
//	//////////////////////////////////////////////////////////////////////////
//	/*int trerferplane = referplane;
//	int tadjustplane = adjustplane;*/
//	//////////////////////////////////////////////////////////////////////////
//	intvector interptind1;
//	intvector interptind2;
//	floatvector interptpos1;
//	floatvector interptpos2;
//	int len;

//	//step1.process the reference plane: get the intersection points
//	computeInterPtList( param, ctrvers, ctredges, referplane, adjustplane, interptpos1, interptind1);
//	len = interptind1.size();

//	//////////////////////////////////////////////////////////////////////////
//	/*cout<<"interptpos1:";
//	for( int i = 0; i < len; i ++ )
//	{
//		cout<<interptind1[i]<<",";
//	}
//	cout<<endl;
//	cout<<"len is ---"<<len;*/
//	//////////////////////////////////////////////////////////////////////////

//	//step2.process the new plane, compute the intersection points
//	computeInterPtList( param, ctrvers, ctredges, adjustplane, referplane, interptpos2, interptind2);

//	//////////////////////////////////////////////////////////////////////////
//	//cout<<"len is ====="<<len;
//	//////////////////////////////////////////////////////////////////////////

//	if( interptind2.size() != len )
//	{
//		//////////////////////////////////////////////////////////////////////////
//		/*cout<<"referplane:"<<referplane<<"adjustplane"<<adjustplane
//			<<"intersection points len1:"<<len<<"intersection points len2:"<<interptind2.size()
//			<<endl;*/
//		//////////////////////////////////////////////////////////////////////////
//		interptind1.clear();
//		interptind2.clear();
//		interptpos1.clear();
//		interptpos2.clear();
//		//	pairlist.resize( pairlist.size() + 1 );
//		return false;
//	}

//	if( len == 0 )	//no intersection points at all
//	{
//		pairlist.resize( pairlist.size() + 1 );
//		return true;
//	}
//	//step3. compute the parameters of these intersection points, compute the parameters for them and sort them out

//	//////////////////////////////////////////////////////////////////////////
//	/*cout<<"========================"<<endl;
//	for( int i= 0 ; i < interptpos1.size()/3; i ++)
//	{
//	cout<<interptind1[ i ]<<"\t"<<interptpos1[ 3*i ]<<"\t"<<interptpos1[ 3*i + 1 ]<<"\t"<<interptpos1[ 3*i +2 ]<<endl;
//	}
//	cout<<"----------"<<endl;
//	for( int i= 0 ; i < interptpos2.size()/3; i ++)
//	{
//	cout<<interptind2[ i ]<<"\t"<<interptpos2[ 3*i ]<<"\t"<<interptpos2[ 3*i + 1 ]<<"\t"<<interptpos2[ 3*i +2 ]<<endl;
//	}
//	cout<<"========================"<<endl;*/
//	//////////////////////////////////////////////////////////////////////////
//	float dir[ 3 ];
//	float param1[ 4 ], param2[ 4 ];
//	for( int i = 0; i < 4; i++)
//	{
//		param1[ i ] = param[ 4*referplane + i ];
//		param2[ i ] = param[ 4*adjustplane + i ];
//	}
//	MyMath::crossProduct( param1, param2, dir);

//	float startpt[ 3 ];
//	float curpt[ 3 ];
//	startpt[ 0 ] = interptpos1[ 0 ];
//	startpt[ 1 ] = interptpos1[ 1 ];
//	startpt[ 2 ] = interptpos1[ 2 ];
//	float* interptparam1 = new float[ len ];
//	float* interptparam2 = new float[ len  ];
//	for( int i = 0; i < len; i ++ )
//	{
//		for( int j = 0 ; j < 3; j ++)
//			curpt[ j ] = interptpos1[  i*3 + j ];
//		interptparam1[ i ] = getTOnRay( startpt, dir, curpt);
//		for( int j = 0; j < 3; j ++)
//		{
//			curpt[ j ] = interptpos2[ i *3 + j];
//		}
//		interptparam2[ i ] = getTOnRay( startpt, dir, curpt );
//	}

//	//sort the points accroding to their parameters
//	for( int i = 0 ; i < len ; i ++)
//	{
//		int minpos = i;
//		for( int j = i + 1; j < len ;j ++)
//		{
//			if( interptparam1[ j ] < interptparam1[ minpos ])
//			{
//				minpos = j;
//			}
//		}
//		if( minpos == i )
//			continue;
//		//switch the two values there
//		//param
//		float t = interptparam1[ minpos ];
//		interptparam1[ minpos ] = interptparam1[ i ];
//		interptparam1[ i ] = t;
//		//position
//		for( int k = 0; k < 3; k ++)
//		{
//			t = interptpos1[ 3*minpos + k ];
//			interptpos1[ 3*minpos + k ] = interptpos1[ 3*i + k ];
//			interptpos1[ 3*i + k ] = t;
//		}
//		//ind
//		int it = interptind1[ minpos ];
//		interptind1[ minpos ] = interptind1[ i ];
//		interptind1[ i ] = it;
//	}

//	for( int i = 0 ; i < len ; i ++)
//	{
//		int minpos = i;
//		for( int j = i + 1; j < len ;j ++)
//		{
//			if( interptparam2[ j ] < interptparam2[ minpos ])
//			{
//				minpos = j;
//			}
//		}
//		if( minpos == i )
//			continue;
//		//switch the two values there
//		//param
//		float t = interptparam2[ minpos ];
//		interptparam2[ minpos ] = interptparam2[ i ];
//		interptparam2[ i ] = t;
//		//position
//		for( int k = 0; k < 3; k ++)
//		{
//			t = interptpos2[ 3*minpos + k ];
//			interptpos2[ 3*minpos + k ] = interptpos2[ 3*i + k ];
//			interptpos2[ 3*i + k ] = t;
//		}
//		//ind
//		int it = interptind2[ minpos ];
//		interptind2[ minpos ] = interptind2[ i ];
//		interptind2[ i ] = it;
//	}

//	//check if the parameters are close enough
//	for( int i = 0; i < len ; i ++)
//	{
//        //if(MyMath::isEqualInToler( interptparam2[ i ] , interptparam1[ i ], 1 ) )	//the distance between them is 1
//        if(true){
//			//adjust the adjustcon
//			for( int k = 0; k < 3; k ++)
//				ctrvers[ adjustplane ][ 3*interptind2[ i ] + k ] = interptpos1[ i*3 + k ];
//		}
//		else	//not correct return false
//		{
//			interptind1.clear();
//			interptind2.clear();
//			interptpos1.clear();
//			interptpos2.clear();
//			delete []interptparam1;
//			delete []interptparam2;
//			//	pairlist.resize( pairlist.size() + 1 );
//			return false;
//		}
//	}

//	//step3.set the pairlist
//	int pairlisti = pairlist.size();
//	pairlist.resize( pairlisti + 1 );
//	for( int i = 0; i < len ;i ++)
//	{
//		pairlist[ pairlisti ].push_back( interptind1[ i ]);
//		pairlist[ pairlisti ].push_back( interptind2[ i ]);
//	}

//	interptind1.clear();
//	interptind2.clear();
//	interptpos1.clear();
//	interptpos2.clear();
//	delete []interptparam1;
//	delete []interptparam2;
//	return true;
//}
bool ContourHandler::makeConsistent_InterPt_2Planes(
        floatvector& param,
        vector<floatvector>& ctrvers,
        vector<intvector>& ctredges,
        //vector<intvector>& v2plistpos,
        //vector<intvector>& v2plist,
        int referplane,
        int adjustplane,
        vector<intvector>& pairlist
        )
{
    //////////////////////////////////////////////////////////////////////////
    /*int trerferplane = referplane;
    int tadjustplane = adjustplane;*/
    //////////////////////////////////////////////////////////////////////////
    intvector interptind1;
    intvector interptind2;
    floatvector interptpos1;
    floatvector interptpos2;
    int len;

    //step1.process the reference plane: get the intersection points
    computeInterPtList( param, ctrvers, ctredges, referplane, adjustplane, interptpos1, interptind1);
    len = interptind1.size();

    //////////////////////////////////////////////////////////////////////////
    /*cout<<"interptpos1:";
    for( int i = 0; i < len; i ++ )
    {
        cout<<interptind1[i]<<",";
    }
    cout<<endl;
    cout<<"len is ---"<<len;*/
    //////////////////////////////////////////////////////////////////////////

    //step2.process the new plane, compute the intersection points
    computeInterPtList( param, ctrvers, ctredges, adjustplane, referplane, interptpos2, interptind2);

    //////////////////////////////////////////////////////////////////////////
    //cout<<"len is ====="<<len;
    //////////////////////////////////////////////////////////////////////////

    if( interptind2.size() != len )
    {
        //////////////////////////////////////////////////////////////////////////
        /*cout<<"referplane:"<<referplane<<"adjustplane"<<adjustplane
            <<"intersection points len1:"<<len<<"intersection points len2:"<<interptind2.size()
            <<endl;*/
        //////////////////////////////////////////////////////////////////////////
        interptind1.clear();
        interptind2.clear();
        interptpos1.clear();
        interptpos2.clear();
        //pairlist.resize( pairlist.size() + 1 );
        return false;

    }

    if( len == 0 )	//no intersection points at all
    {
        pairlist.resize( pairlist.size() + 1 );
        return true;
    }
    //step3. compute the parameters of these intersection points, compute the parameters for them and sort them out

    //////////////////////////////////////////////////////////////////////////
    /*cout<<"========================"<<endl;
    for( int i= 0 ; i < interptpos1.size()/3; i ++)
    {
    cout<<interptind1[ i ]<<"\t"<<interptpos1[ 3*i ]<<"\t"<<interptpos1[ 3*i + 1 ]<<"\t"<<interptpos1[ 3*i +2 ]<<endl;
    }
    cout<<"----------"<<endl;
    for( int i= 0 ; i < interptpos2.size()/3; i ++)
    {
    cout<<interptind2[ i ]<<"\t"<<interptpos2[ 3*i ]<<"\t"<<interptpos2[ 3*i + 1 ]<<"\t"<<interptpos2[ 3*i +2 ]<<endl;
    }
    cout<<"========================"<<endl;*/
    //////////////////////////////////////////////////////////////////////////
    float dir[ 3 ];
    float param1[ 4 ], param2[ 4 ];
    for( int i = 0; i < 4; i++)
    {
        param1[ i ] = param[ 4*referplane + i ];
        param2[ i ] = param[ 4*adjustplane + i ];
    }
    MyMath::crossProduct( param1, param2, dir);

    float startpt[ 3 ];
    float curpt[ 3 ];
    startpt[ 0 ] = interptpos1[ 0 ];
    startpt[ 1 ] = interptpos1[ 1 ];
    startpt[ 2 ] = interptpos1[ 2 ];
    float* interptparam1 = new float[ len ];
    float* interptparam2 = new float[ len  ];
    for( int i = 0; i < len; i ++ )
    {
        for( int j = 0 ; j < 3; j ++)
            curpt[ j ] = interptpos1[  i*3 + j ];
        interptparam1[ i ] = getTOnRay( startpt, dir, curpt);
        for( int j = 0; j < 3; j ++)
        {
            curpt[ j ] = interptpos2[ i *3 + j];
        }
        interptparam2[ i ] = getTOnRay( startpt, dir, curpt );
    }

    //sort the points accroding to their parameters
    for( int i = 0 ; i < len ; i ++)
    {
        int minpos = i;
        for( int j = i + 1; j < len ;j ++)
        {
            if( interptparam1[ j ] < interptparam1[ minpos ])
            {
                minpos = j;
            }
        }
        if( minpos == i )
            continue;
        //switch the two values there
        //param
        float t = interptparam1[ minpos ];
        interptparam1[ minpos ] = interptparam1[ i ];
        interptparam1[ i ] = t;
        //position
        for( int k = 0; k < 3; k ++)
        {
            t = interptpos1[ 3*minpos + k ];
            interptpos1[ 3*minpos + k ] = interptpos1[ 3*i + k ];
            interptpos1[ 3*i + k ] = t;
        }
        //ind
        int it = interptind1[ minpos ];
        interptind1[ minpos ] = interptind1[ i ];
        interptind1[ i ] = it;
    }

    for( int i = 0 ; i < len ; i ++)
    {
        int minpos = i;
        for( int j = i + 1; j < len ;j ++)
        {
            if( interptparam2[ j ] < interptparam2[ minpos ])
            {
                minpos = j;
            }
        }
        if( minpos == i )
            continue;
        //switch the two values there
        //param
        float t = interptparam2[ minpos ];
        interptparam2[ minpos ] = interptparam2[ i ];
        interptparam2[ i ] = t;
        //position
        for( int k = 0; k < 3; k ++)
        {
            t = interptpos2[ 3*minpos + k ];
            interptpos2[ 3*minpos + k ] = interptpos2[ 3*i + k ];
            interptpos2[ 3*i + k ] = t;
        }
        //ind
        int it = interptind2[ minpos ];
        interptind2[ minpos ] = interptind2[ i ];
        interptind2[ i ] = it;
    }

    //check if the parameters are close enough
    for( int i = 0; i < len ; i ++)
    {
        if(MyMath::isEqualInToler( interptparam2[ i ] , interptparam1[ i ], 1 ) )	//the distance between them is 1
            if(true){
                //adjust the adjustcon
                for( int k = 0; k < 3; k ++)
                    ctrvers[ adjustplane ][ 3*interptind2[ i ] + k ] = interptpos1[ i*3 + k ];
            }
            else	//not correct return false
            {
                interptind1.clear();
                interptind2.clear();
                interptpos1.clear();
                interptpos2.clear();
                delete []interptparam1;
                delete []interptparam2;
                //	pairlist.resize( pairlist.size() + 1 );
                return false;
            }
    }

    //step3.set the pairlist
    int pairlisti = pairlist.size();
    pairlist.resize( pairlisti + 1 );
    for( int i = 0; i < len ;i ++)
    {
        pairlist[ pairlisti ].push_back( interptind1[ i ]);
        pairlist[ pairlisti ].push_back( interptind2[ i ]);
    }

    interptind1.clear();
    interptind2.clear();
    interptpos1.clear();
    interptpos2.clear();
    delete []interptparam1;
    delete []interptparam2;
    return true;
}





static void projectAllOnEdgeVertices(float *vpos,float *pparam1,float *pparam2){


    auto projectErr = [](float *pp,float* ppara){
        return abs(pp[0]*ppara[0]+pp[1]*ppara[1]+pp[2]*ppara[2]-ppara[3]);
    };
    auto project2p = [](float *pp,float* ppara){
        float err = (pp[0]*ppara[0]+pp[1]*ppara[1]+pp[2]*ppara[2]-ppara[3]);
        for(int i=0;i<3;++i)pp[i]-=ppara[i]*err;
    };

    float error = 100;
    while(error>1e-2){
        project2p(vpos,pparam1);
        project2p(vpos,pparam2);
        error = 0;
        error+=projectErr(vpos,pparam1);
        error+=projectErr(vpos,pparam2);
    }

}
static void projectAllOnEdgeVertices(vector<float>&vpos,vector<vector<int>>&v2ps,float *pparam){


    auto projectErr = [pparam](float *pp,int pind){
        auto ppara = pparam+pind*4;
        return abs(pp[0]*ppara[0]+pp[1]*ppara[1]+pp[2]*ppara[2]-ppara[3]);
    };
    auto project2p = [pparam](float *pp,int pind){
        auto ppara = pparam+pind*4;
        float err = (pp[0]*ppara[0]+pp[1]*ppara[1]+pp[2]*ppara[2]-ppara[3]);
        for(int i=0;i<3;++i)pp[i]+=ppara[i]*err;
    };

    int npoints = v2ps.size();
    for(int i=0;i<npoints;++i){
        auto &pl = v2ps[i];
        float *pp = vpos.data()+i*3;
        float error = 100;
        while(error>1e-2){
            for(auto a:pl)project2p(pp,a);
            error = 0;
            for(auto a:pl)error+=projectErr(pp,a);
        }
    }


}
bool ContourHandler::makeConsistent_InterPt_2Planes2(
        floatvector& param,
        vector<floatvector>& ctrvers,
        vector<intvector>& ctredges,
        //vector<intvector>& v2plistpos,
        //vector<intvector>& v2plist,
        int referplane,
        int adjustplane,
        vector<intvector>& pairlist
        )
{
    //////////////////////////////////////////////////////////////////////////
    /*int trerferplane = referplane;
    int tadjustplane = adjustplane;*/
    //////////////////////////////////////////////////////////////////////////
    intvector interptind1;
    intvector interptind2;
    floatvector interptpos1;
    floatvector interptpos2;
    int len;

    //step1.process the reference plane: get the intersection points
    computeInterPtList( param, ctrvers, ctredges, referplane, adjustplane, interptpos1, interptind1);
    len = interptind1.size();


    //////////////////////////////////////////////////////////////////////////
    /*cout<<"interptpos1:";
    for( int i = 0; i < len; i ++ )
    {
        cout<<interptind1[i]<<",";
    }
    cout<<endl;
    cout<<"len is ---"<<len;*/
    //////////////////////////////////////////////////////////////////////////

    //step2.process the new plane, compute the intersection points
    computeInterPtList( param, ctrvers, ctredges, adjustplane, referplane, interptpos2, interptind2);

    //////////////////////////////////////////////////////////////////////////
    //cout<<"len is ====="<<len;


    if( interptind1.size() == 0 && interptind2.size() ==0 )	//no intersection points at all
    {
        pairlist.resize( pairlist.size() + 1 );
        return true;
    }

    if(0)if( true)	//no intersection points at all
    {
        pairlist.resize( pairlist.size() + 1 );
        return true;
    }
    //////////////////////////////////////////////////////////////////////////

    for(int i=0;i< interptind1.size();++i)projectAllOnEdgeVertices(interptpos1.data()+i*3,param.data()+adjustplane*4,param.data()+referplane*4);

    for(int i=0;i< interptind2.size();++i)projectAllOnEdgeVertices(interptpos2.data()+i*3,param.data()+adjustplane*4,param.data()+referplane*4);

    for(int i=0;i< interptind1.size();++i)for(int j=0;j<3;++j)ctrvers[ referplane ][ interptind1[i]*3+j ] = interptpos1[i*3+j];
    for(int i=0;i< interptind2.size();++i)for(int j=0;j<3;++j)ctrvers[ adjustplane ][ interptind2[i]*3+j ] = interptpos2[i*3+j];
    //////////////////////////////////////////////////////////////////////////
    //step3. compute the parameters of these intersection points, compute the parameters for them and sort them out

    //////////////////////////////////////////////////////////////////////////
    /*cout<<"========================"<<endl;
    for( int i= 0 ; i < interptpos1.size()/3; i ++)
    {
    cout<<interptind1[ i ]<<"\t"<<interptpos1[ 3*i ]<<"\t"<<interptpos1[ 3*i + 1 ]<<"\t"<<interptpos1[ 3*i +2 ]<<endl;
    }
    cout<<"----------"<<endl;
    for( int i= 0 ; i < interptpos2.size()/3; i ++)
    {
    cout<<interptind2[ i ]<<"\t"<<interptpos2[ 3*i ]<<"\t"<<interptpos2[ 3*i + 1 ]<<"\t"<<interptpos2[ 3*i +2 ]<<endl;
    }
    cout<<"========================"<<endl;*/
    //////////////////////////////////////////////////////////////////////////
    float dir[ 3 ];
    float param1[ 4 ], param2[ 4 ];
    for( int i = 0; i < 4; i++)
    {
        param1[ i ] = param[ 4*referplane + i ];
        param2[ i ] = param[ 4*adjustplane + i ];
    }
    MyMath::crossProduct( param1, param2, dir);

    float startpt[ 3 ];
    float curpt[ 3 ];
    if(interptind1.size() != 0){
        startpt[ 0 ] = interptpos1[ 0 ];
        startpt[ 1 ] = interptpos1[ 1 ];
        startpt[ 2 ] = interptpos1[ 2 ];
    }else{
        startpt[ 0 ] = interptpos2[ 0 ];
        startpt[ 1 ] = interptpos2[ 1 ];
        startpt[ 2 ] = interptpos2[ 2 ];
    }
    //    float* interptparam1 = new float[ interptind1.size() ];
    //    float* interptparam2 = new float[ interptind2.size()  ];
    vector<float>interptparam1(interptind1.size());
    vector<float>interptparam2(interptind2.size());
    for( int i = 0; i < interptind1.size(); i ++ )
    {
        for( int j = 0 ; j < 3; j ++)
            curpt[ j ] = interptpos1[  i*3 + j ];
        interptparam1[ i ] = getTOnRay( startpt, dir, curpt);
    }
    for( int i = 0; i < interptind2.size(); i ++ ){
        for( int j = 0; j < 3; j ++)
        {
            curpt[ j ] = interptpos2[ i *3 + j];
        }
        interptparam2[ i ] = getTOnRay( startpt, dir, curpt );
    }

    //sort the points accroding to their parameters
    for( int i = 0 ; i < interptind1.size() ; i ++)
    {
        int minpos = i;
        for( int j = i + 1; j < interptind1.size() ;j ++)
        {
            if( interptparam1[ j ] < interptparam1[ minpos ])
            {
                minpos = j;
            }
        }
        if( minpos == i )
            continue;
        //switch the two values there
        //param
        float t = interptparam1[ minpos ];
        interptparam1[ minpos ] = interptparam1[ i ];
        interptparam1[ i ] = t;
        //position
        for( int k = 0; k < 3; k ++)
        {
            t = interptpos1[ 3*minpos + k ];
            interptpos1[ 3*minpos + k ] = interptpos1[ 3*i + k ];
            interptpos1[ 3*i + k ] = t;
        }
        //ind
        int it = interptind1[ minpos ];
        interptind1[ minpos ] = interptind1[ i ];
        interptind1[ i ] = it;
    }

    for( int i = 0 ; i < interptind2.size() ; i ++)
    {
        int minpos = i;
        for( int j = i + 1; j < interptind2.size() ;j ++)
        {
            if( interptparam2[ j ] < interptparam2[ minpos ])
            {
                minpos = j;
            }
        }
        if( minpos == i )
            continue;
        //switch the two values there
        //param
        float t = interptparam2[ minpos ];
        interptparam2[ minpos ] = interptparam2[ i ];
        interptparam2[ i ] = t;
        //position
        for( int k = 0; k < 3; k ++)
        {
            t = interptpos2[ 3*minpos + k ];
            interptpos2[ 3*minpos + k ] = interptpos2[ 3*i + k ];
            interptpos2[ 3*i + k ] = t;
        }
        //ind
        int it = interptind2[ minpos ];
        interptind2[ minpos ] = interptind2[ i ];
        interptind2[ i ] = it;
    }

    //check if the parameters are close enough
    vector<float>paradist(interptind2.size() * interptind1.size());
    int len1 = interptind1.size(),len2 = interptind2.size();
    for(int i=0;i<len1;++i){
        for(int j=0;j<len2;++j){
            paradist[i*len2+j] = fabs(interptparam2[j] - interptparam1[i]);
        }
    }

    len = len1>len2?len1:len2;

    vector<int>pick1(len1,-1);
    vector<int>pick2(len2,-1);

    if(paradist.size()!=0)while(true){
        auto pmin = min_element(paradist.begin(),paradist.end());

        int minind = pmin - paradist.begin();

        if(*pmin>2.0)break;

        *pmin = 999999.9;

        int m2 = minind%len2, m1 = (minind-m2)/len2;

        assert(pick1[m1]==-1);
        assert(pick2[m2]==-1);
        pick1[m1] = m2;
        pick2[m2] = m1;

        for(int i=0;i<len2;++i)paradist[m1*len2+i] = 999999.9;
        for(int i=0;i<len1;++i)paradist[i*len2+m2] = 999999.9;

    }

    int pairlisti = pairlist.size();
    pairlist.resize( pairlisti + 1 );
    auto &ppairlist = pairlist[ pairlisti ];

    for(int i=0;i<len1;++i){
        if(pick1[i]!=-1){
            for( int k = 0; k < 3; k ++)ctrvers[ adjustplane ][ 3*interptind2[ pick1[i] ] + k ] = interptpos1[ i*3 + k ];
            ppairlist.push_back( interptind1[ i ]);
            ppairlist.push_back( interptind2[ pick1[i] ]);
        }else{
           ppairlist.push_back( interptind1[ i ]);
           ppairlist.push_back( ctrvers[ adjustplane ].size()/3);
           for( int k = 0; k < 3; k ++)ctrvers[ adjustplane ].push_back(interptpos1[ i*3 + k ]);
        }
    }

    for(int i=0;i<len2;++i){
        if(pick2[i]==-1){
          ppairlist.push_back( ctrvers[ referplane ].size()/3);
          ppairlist.push_back( interptind2[ i ]);
          for( int k = 0; k < 3; k ++)ctrvers[ referplane ].push_back(interptpos2[ i*3 + k ]);
        }
    }



    //    for( int i = 0; i < len ; i ++)
    //    {
    //        if(MyMath::isEqualInToler( interptparam2[ i ] , interptparam1[ i ], 1 ) )	//the distance between them is 1
    //        if(true){
    //            //adjust the adjustcon
    //            for( int k = 0; k < 3; k ++)
    //                ctrvers[ adjustplane ][ 3*interptind2[ i ] + k ] = interptpos1[ i*3 + k ];
    //        }
    //        else	//not correct return false
    //        {
    //            interptind1.clear();
    //            interptind2.clear();
    //            interptpos1.clear();
    //            interptpos2.clear();
    //            delete []interptparam1;
    //            delete []interptparam2;
    //            //	pairlist.resize( pairlist.size() + 1 );
    //            return false;
    //        }
    //    }



    //step3.set the pairlist

    //    for( int i = 0; i < len ;i ++)
    //    {
    //        pairlist[ pairlisti ].push_back( interptind1[ i ]);
    //        pairlist[ pairlisti ].push_back( interptind2[ i ]);
    //    }

    interptind1.clear();
    interptind2.clear();
    interptpos1.clear();
    interptpos2.clear();
    //delete []interptparam1;
    //delete []interptparam2;
    return true;
}
/*
ver2planelistpos:specify the position in ver2planelist
ver2planelist: array of the planes, the contour vertex is on
because of numeric issue, when preprocessing, remember the planes the vertex should be on
so that when the contour vertices are partitioned into the face, they will give us right result, we
specify the location instead of checking if it is on or not, which might not be consistent with this step.
*/
void ContourHandler::makeConsistent_InterPt(
        int& planenum, float*& pparam,
        int*& pctrvernum, float**& pctrvers,
        int*& pctredgenum, int**& pctredges,
        int**& ver2planelistpos, vector<intvector>& ver2planelist)
{
    // int tplanenum = 0;
    floatvector tparam;
    vector<floatvector> tctrvers;
    //intvector tctrvernum;
    vector<intvector> tctredges;
    //intvector tctredgenum;
    vector<intvector> tv2plistpos;
    vector<intvector> pairlist;

    for( int i = 0; i < planenum; i ++)
    {
        int rowi = tctrvers.size();	//the row index of the new row
        floatvector tvec;
        tctrvers.push_back( tvec );
        tctrvers[ rowi ].resize( pctrvernum[ i ]*3);
        for( int j = 0; j < pctrvernum[ i ]*3; j++ )
        {
            tctrvers[ rowi ][ j ] = pctrvers[ i ][ j ];
        }
        intvector tvec2;
        tctredges.push_back( tvec2 );
        tctredges[ rowi ].resize( pctredgenum[ i ]*4 );
        for( int j = 0; j < pctredgenum[  i ]*4; j ++)
        {
            tctredges[ rowi ][ j ] = pctredges[ i ][ j ];
        }
        for( int j = 0; j < 4; j ++)
            tparam.push_back( pparam[ 4*i + j] );

        bool suc = true;
        for( int j = 0; j < rowi; j ++)
        {
            suc = makeConsistent_InterPt_2Planes2(tparam,tctrvers,tctredges,j,rowi,pairlist);
            if( !suc )
            {
                suc = false;
                break;
            }
        }

        //suc=true;

        if( suc )
        {
            //set the ver2planelist
            tv2plistpos.resize( rowi + 1 );
            for( int j = 0; j <= rowi; j ++)
            {
                int tvernum = tctrvers[ j ].size()/3;
                tv2plistpos[ j ].resize( tvernum, -1);
            }

            //for those pair list, set the vertex to plane list
            for( int j = 0; j < rowi; j ++)	//for each row in pairlist
            {
                int pairnum = pairlist[ j ].size()/2;
                int plane2[ 2 ] = { j, rowi };
                for( int k = 0; k < pairnum; k++)
                {
                    int tver2[ 2 ] = { pairlist[ j ][ 2*k ], pairlist[ j ][ k*2 + 1 ]};
                    for( int k1 = 0; k1 < 2; k1 ++)
                    {
                        int posinv2p = tv2plistpos[ plane2[ k1 ] ][ tver2[ k1 ]];
                        if(	posinv2p == -1 )	//this vertex has not been set before
                        {
                            intvector tempvec;
                            ver2planelist.push_back( tempvec);
                            ver2planelist[ ver2planelist.size() - 1].push_back( plane2[ 1 - k1] );

                            tv2plistpos[ plane2[ k1 ] ][ tver2[ k1 ]] = ver2planelist.size() - 1;
                        }
                        else
                        {
                            ver2planelist[ posinv2p ].push_back( plane2[ 1 - k1 ] );
                        }
                    }
                }
            }
        }
        else
        {

            //////////////////////////////////////////////////////////////////////////
            cout<<"plane "<<i<<" is discarded!"<<endl;
            //////////////////////////////////////////////////////////////////////////
            //pop out the already pushed in vertex and edges
            tparam.pop_back();
            tparam.pop_back();
            tparam.pop_back();
            tparam.pop_back();
            tctrvers[ rowi ].clear();
            tctredges[ rowi ].clear();
            tctrvers.resize( rowi );
            tctredges.resize( rowi );

        }
        //clear pairlist
        //////////////////////////////////////////////////////////////////////////
        //cout<<"pairlist size"<<pairlist.size()<<endl;
        //////////////////////////////////////////////////////////////////////////
        for (unsigned int j = 0; j < pairlist.size(); j++)
            pairlist[ j ].clear();
        pairlist.clear();
    }

    //reset the contour parameters
    //for( int i = 0; i < planenum; i ++)
    //{
    //	delete []pctrvers[ i ];
    //	delete []pctredges[ i ];
    //}
    //delete []pctrvers;
    //delete []pctredges;
    //delete []pctrvernum;
    //delete []pctredgenum;
    //
    ////before delete the old parameters, we should first push back all the added six planes parameters
    //int starti = planenum * 4;
    //   for(int j = 0; j < 24; j ++)
    //{
    //	tparam.push_back( pparam[ starti + j ] );
    //}
    //delete []pparam;

    planenum = tctrvers.size();
    pparam = new float[ 4*planenum + 24];
    pctrvernum = new int[ planenum ];
    pctrvers = new float*[ planenum ];
    pctredgenum = new int[ planenum ];
    pctredges = new int*[ planenum ];

    for( int i = 0; i < planenum; i ++)
    {
        pctrvernum[ i ] = tctrvers[ i ].size()/3;
        pctrvers[ i ] = new float[ tctrvers[ i].size() ];
        for (unsigned int j = 0; j < tctrvers[i].size(); j++)
        {
            pctrvers[ i ][ j ]= tctrvers[ i ][ j ];
        }
        pctredgenum[ i ] = tctredges[ i ].size()/4;
        pctredges[ i ] = new int[ tctredges[ i ].size()];
        for (unsigned int j = 0; j < tctredges[i].size(); j++)
        {
            pctredges[ i ][ j ] = tctredges[ i ][ j ];
        }

        tctrvers[ i ].clear();
        tctredges[ i ].clear();
    }
    tctrvers.clear();
    tctredges.clear();

    for (unsigned int i = 0; i < tparam.size(); i++)
        pparam[ i ] = tparam[ i ];
    tparam.clear();

    //set the ver2planelistpos, and ver2planelist
    ver2planelistpos = new int*[ planenum ];
    for( int i = 0 ;i < planenum; i ++)
    {
        ver2planelistpos[ i ] = new int[ pctrvernum[ i ]];
        for( int j = 0; j < pctrvernum[ i ]; j ++)
        {
            ver2planelistpos[ i ][ j ] = tv2plistpos [ i ][ j ];
        }
        tv2plistpos[ i ].clear();
    }
    tv2plistpos.clear();

}
void ContourHandler::preProcDataSingleMat(
        int& planenum,float*& pparam,float**& pctrvers,int*& pctrvernum,int**& pctredges,int*& pctredgenum,
        int**& ver2planelistpos, vector<intvector>& ver2planelist
        )
{
    //////////////////////////////////////////////////////////////////////////
    /*cout<<"plane 1:"<<endl;
    cout<<"ver1:"<<pctrvers[ 1 ][ 0 ]<<"  "<<pctrvers[ 1 ][ 1 ]<<"  "<<pctrvers[ 1 ][ 2 ]<<endl;
    cout<<"ver2:"<<pctrvers[ 1 ][ 3 ]<<"  "<<pctrvers[ 1 ][ 4 ]<<"  "<<pctrvers[ 1 ][ 5 ]<<endl;
    cout<<"plane 2: "<<endl;
    cout<<"ver1:"<<pctrvers[ 2 ][ 0 ]<<"  "<<pctrvers[ 2 ][ 1 ]<<"  "<<pctrvers[ 2 ][ 2 ]<<endl;
    cout<<"ver2:"<<pctrvers[ 2 ][ 3 ]<<"  "<<pctrvers[ 2 ][ 4 ]<<"  "<<pctrvers[ 2 ][ 5 ]<<endl;

    cout<<"plane1 param:"<<pparam[ 4 ]<<" "<<pparam[ 5 ] <<"  "<<pparam[ 6 ]<<" "<<pparam[ 7 ]<<endl;
    cout<<"ver1 relative to this param:"<<MyMath::dotProduct( pparam + 4, pctrvers[ 1 ])<<endl;*/
    //////////////////////////////////////////////////////////////////////////
    makeConsistent_InterPt(
                planenum, pparam,pctrvernum, pctrvers,pctredgenum, pctredges,
                ver2planelistpos, ver2planelist);
}

