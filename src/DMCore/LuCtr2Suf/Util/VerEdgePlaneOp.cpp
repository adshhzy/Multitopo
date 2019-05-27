#include "../Util/VerEdgePlaneOp.h"
#include "../Math/mymath.h"

//for debugging reason
void writeComnLineDB(float param1[ 4 ], float param2[ 4 ], float comndir[ 3 ], float comnpt[ 3 ])
{
	FILE* fout = fopen("comnparam.txt", "w");
	
	fprintf( fout, "{{%f,%f,%f,%f},", param1[ 0 ] , param1[ 1 ], param1[ 2 ], param1[ 3 ]);
	fprintf(fout, "{%f,%f,%f,%f},", param2[ 0 ], param2[ 1 ], param2[ 2 ], param2[ 3 ]);
	fprintf( fout, "{%f,%f,%f},", comndir[ 0 ], comndir[ 1 ], comndir[ 2 ]);
	fprintf(fout, "{%f,%f,%f}}", comnpt[ 0 ], comnpt[ 1 ], comnpt[ 2 ]);

	fclose( fout );
}

/*
* function: compute the common line between two intersecting planes
*	@param:param1, ax+by+cz = d	for plane 1, normal should be normalized before passing in
* @param: param2, for plane 2
* @comndir: direction of the common line	-result
* @comnpt: point on the common line			-result
*/ 
void computeComnLine( float param1[ 4 ], float param2[ 4 ], float comndir[ 3 ], float comnpt[ 3 ] )
{
	//compute the direction of the line
	MyMath::crossProduct( param1, param2, comndir );

	//check if there is one plane, b and c are both 0
	int bczero = -1;
	if( MyMath::isEqualInToler(param1[ 1 ], 0, TOLERANCE_FOUR) &&
			MyMath::isEqualInToler(param1[ 2 ], 0, TOLERANCE_FOUR))
			bczero = 0;
	else if(MyMath::isEqualInToler(param2[ 1 ], 0, TOLERANCE_FOUR) &&
		MyMath::isEqualInToler(param2[ 2 ], 0, TOLERANCE_FOUR))
		bczero = 1;
	if( bczero != -1 )	//the plane bczero, with b and c both zero
	{
		float* param;
		if( bczero == 0 )
			param = param1;
		else
			param = param2;
		//set the common x
		comnpt[ 0 ] = param[ 3 ]/param[ 0 ];

		if( param == param1 ) 
			param = param2;
		else
			param = param1;
		
		float b = param[ 1 ];
		float c = param[ 2 ];
		float val = param[ 3 ] - param[ 0 ] * comnpt[ 0 ];
		param = NULL;

		//check if b is equal to zero or not
		if( !MyMath::isEqualInToler( b, 0, TOLERANCE_FOUR ))
		{
			comnpt[ 2 ] = 0;
			comnpt[ 1 ] =  val / b;
		}
		else
		{
			comnpt[ 1 ] = 0;
			comnpt[ 2 ] = val / c;
		}

		return;
	}

	//there is no plane with b and c both are zero
	comnpt[ 0 ] = 0;
	double temp =  ((double) param1[ 3 ] * 1000 * param2[ 2 ] -(double) param2[ 3 ] * 1000 *  param1[ 2 ] )
		/((double)param2[ 2 ] * 1000 *  param1[ 1 ] - (double)param1[ 2 ] * 1000* param2[ 1 ]);
	comnpt[ 1 ] = (float)temp;
	temp =  ((double)param1[ 3 ] * 1000 * param2[ 1 ] - (double) param2[ 3 ] * 1000 *  param1[ 1 ] )
		/((double)param1[ 2 ] *  1000 * param2[ 1 ] - (double) param2[ 2 ] *  1000 * param1[ 1 ]);
	comnpt[2] = (float)temp;

	//////////////////////////////////////////////////////////////////////////
	cout<<"plane 0 "<<comnpt[ 1 ] * param1[ 1 ] + comnpt[ 2 ] * param1 [ 2 ] - param1[ 3 ] <<endl;
cout<<"plane 1 "<<comnpt[ 1 ] * param2[ 1 ] + comnpt[ 2 ] * param2 [ 2 ] - param2[ 3 ] <<endl;

	cout<<param1[ 0 ]<<","<<param1[ 1 ]<<","<<param1[ 2 ]<<","<<param1[ 3 ]<<endl;
	cout<<param2[ 0 ]<<","<<param2[ 1]<<","<<param2[ 2 ] <<","<<param2[ 3 ]<<endl;	
	cout<<"comndir:"<<comndir[0]<<","<<comndir[ 1 ]<<","<<comndir[ 2 ]<<endl;
	cout<<"comnpt:"<<comnpt[ 0 ]<<","<<comnpt[ 1 ]<<","<<comnpt[ 2 ]<<endl;
	
	//////////////////////////////////////////////////////////////////////////
}

float getAngleVecX2D(float vec[ 2 ], float xaxis[2] )
{
	//compute cos
	float cosval = vec[ 0 ]*xaxis[ 0 ] + vec[ 1 ]*xaxis[ 1 ];
	cosval = cosval/(sqrt((float)(vec[ 0 ]*vec[ 0 ] + vec[1 ]*vec[1])));

	if( (xaxis[ 0 ]*vec[1 ] - xaxis[ 1 ]*vec[ 0 ]) < 0)
		cosval += 2;
	else
		cosval = -cosval;
	return cosval;
}
//the dotproduct between vec and xaxis,
//0-2pi, to make it monotonously increaing, 0-pi, flip, pi - 2pi, translate by 2
//xaxis must be unit vector, or else, not correct answer
float getAngleVecX(float vec[ 3 ], float xaxis[ 3 ], float zaxis[ 3 ])
{
	//compute cos
	float cosval = MyMath::dotProduct( vec, xaxis );
	cosval /= MyMath::vectorlen( vec );

	//range
	float dir[ 3 ];
	MyMath::crossProductNotNorm( xaxis, vec, dir );	
	if( MyMath::dotProduct( dir, zaxis ) < 0 )	//pi - 2pi
		cosval += 2;
	else
		cosval = -cosval;
	return cosval;
}

void getDirOfTrian(float pt1[ 3 ], float pt2[ 3 ], float pt3[ 3 ],
				   float dir[ 3 ])
{
	float vec1[ 3 ], vec2[ 3 ];
	MyMath::getVec( pt1, pt2, vec1);
	MyMath::getVec( pt1, pt3, vec2);
	MyMath::crossProductNotNorm( vec1, vec2, dir);
}

//decide if the two vectors are in the same direction
//the two direction must both be NON-ZERO
//criteria: 1, angle between them are small,
//2, the distance from the end point of short edge to the long is small
//in order to make it consistent with criteria for vertex lying on
//edge, 1. is removed!
bool isSameDir( float dir1[ 3 ], float dir2[ 3 ])
{
	float len1 = MyMath::vectorlen( dir1 );
	float len2 = MyMath::vectorlen( dir2 );
	float cosv = MyMath::dotProduct( dir1, dir2 )/(len1*len2);
	if( cosv < 0 )
		return false;
//	if( cosv < (1 - TOLERANCE_ANGLE ) )
//		return false;

    float sinv = sqrt( 1 - cosv*cosv );
	if( sinv * min( len1, len2 ) < TOLERANCE_SAME_VER )
		return true;
	return false;
}

int getMajorDirec( float vec[ 3 ])
{
	int maxi = 0;
	if( abs( vec[ 1 ] ) > abs( vec[ maxi ] ))
		maxi = 1;
	if( abs( vec[ 2 ]) > abs( vec[ maxi ]))
		maxi = 2;
	return maxi;
}

//pick the axes, along which the difference is biggest
int getMajorDirec( float pt1[ 3 ], float pt2[ 3 ])
{
	int maxi = 0;
	if( abs( pt1[ 1 ] - pt2[ 1 ] ) > abs( pt1[ maxi ] - pt2[ maxi ]))
		maxi = 1;
	if( abs( pt1[ 2 ] - pt2[ 2 ] ) > abs( pt1[ maxi ] - pt2[ maxi ]))
		maxi = 2;
	return maxi;
}
bool pointEqualOriginCube( float pos1[ 3 ], float toler )
{
	for( int i = 0;i < 3; i ++)
	{
		if( !MyMath::isEqualInToler( pos1[ i ],0, toler ))
			return false;
	}
	return true;
}

bool pointEqualCube( float pos1[ 3 ], float pos2[ 3 ], float toler )
{
	for( int i = 0;i < 3; i ++)
	{
		if( !MyMath::isEqualInToler( pos1[ i ], pos2[ i ], toler ))
			return false;
	}
	return true;
}	

bool pointEqualSphere(float pos1[ 3 ], float pos2[ 3 ], float toler)
{
	float diff[ 3 ];
	MyMath::getVec( pos1, pos2, diff);
	return MyMath::isEqualInToler( MyMath::vectorlen( diff ), 0, toler );
}
//return true, if the third point lies in the range of the two translated TOLERANCE_FOUR lines
bool onSameLine( float pt1[3], float pt2[ 3 ], float pt3[ 3 ])
{
	float vec1[ 3 ];
	float vec2[ 3 ];
	MyMath::getVec( pt1, pt3, vec1);
	MyMath::getVec( pt2, pt3, vec2 );
	//MyMath::normalize( vec1 );
	//MyMath::normalize( vec2 );
	float dir[ 3 ];
	MyMath::crossProductNotNorm( vec1, vec2, dir );
	if( MyMath::vectorlen(dir)/MyMath::vectorlen( pt1, pt2 ) < TOLERANCE_FOUR )
		return true;
	return false;
	//return pointEqualOriginCube( dir, TOLERANCE_SIX);
}

//compute the parameter of the point on this segment
//make sure that the two points of the segment do not overlap
float getParamOnSeg(float pt1[ 3 ], float pt2[ 3 ], float pt[ 3 ])
{
	int maxind = 0;
	if( abs( pt1[ 1 ] - pt2[ 1 ]) > abs(pt1[ maxind] - pt2[ maxind ]) )		maxind = 1;
	if( abs( pt1[ 2 ] - pt2[ 2 ]) > abs(pt1[ maxind] - pt2[ maxind ]) )		maxind = 2;

	return ( pt[ maxind ] - pt1[ maxind ])/(pt2[ maxind ] - pt1[ maxind ]);
}
//compute the projection of a point on the plane
//planeparam: ax + by + cz = d {a,b,c,d}
void projectionPtOnPlane(float pt[ 3 ],  float planeparam[ 4 ], bool normisunit, float projectionpt[ 3 ])
{
	float t = (planeparam[ 3 ]-MyMath::dotProduct( pt, planeparam ))/MyMath::dotProduct( planeparam, planeparam );
	MyMath::getPtOnRay( pt, planeparam, t,projectionpt);
}
//no check is done for it, so before passing your parameters in
//make sure, 1: the dir is not 0 vector, and the dir is unit vector
float getTOnRay( float startpt[ 3 ], float dir[ 3 ], float pt[ 3 ])
{
	int maxind = 0;
	if( abs(dir[ 1 ]) > abs(dir[ maxind ])) maxind = 1;
	if( abs(dir[ 2 ]) > abs(dir[ maxind ])) maxind = 2;
	return (pt[ maxind ] - startpt[ maxind ])/dir[ maxind ];
}
/**
* compute the distance from the point to the plane along direction dir
* if the direction is parallel to the plane, no matter pt is on the plane or not return -1
* otherwise, return the distance from pt to the intersection point, might be negative, if the
* intersection point is in the opposite direction of dir
*/
float distPt2PlaneAlongDir( float pt[ 3 ], float dir[ 3 ],bool dirisunit, float planept[ 3 ], float planenorm[ 3 ], bool normisunit)
{
	if( !dirisunit )
		MyMath::normalize( dir );
	if( !normisunit )
		MyMath::normalize( planenorm );

	float dotval = MyMath::dotProduct( dir, planenorm);
	if( MyMath::isEqualInToler( dotval, 0, TOLERANCE_SIX ))	//dir is parallel to the plane
	{
		return -1;
	}

	float dir2[ 3 ];
	MyMath::getVec( pt, planept, dir2 );
	float dir3[ 3 ];
	MyMath::normalize( dir2, dir3 );
	dotval = MyMath::dotProduct( dir3, planenorm );
	if( MyMath::isEqualInToler( dotval, 0, TOLERANCE_SIX ))	//point is on the plane
		return 0;

	return (MyMath::dotProduct(dir2, planenorm)/MyMath::dotProduct(dir, planenorm));
}

/**
* compute the intersection point between the plane and the segment. 
* intersection must exist, otherwise, may crash, no check for it.
*/
void interPt_PlaneEdge( float param[ 4 ], float ver1[ 3 ], float ver2[ 3 ], float interpt[ 3 ])
{
	float t1 = MyMath::dotProduct( param, ver1 );
	t1 -= param[ 3 ];
	float t2 = MyMath::dotProduct( param, ver2);
	t2 -= param[ 3 ];

	for( int i = 0; i < 3; i ++)
	{
		interpt[ i ] = (ver1[ i ]* t2  - ver2[ i ] * t1)/(t2 - t1);
	}
}

/**
*	decide if the line with two end points ver1 and ver2 is parallel to the plane with parameter param
*/
bool isParalLinePlane( float ver1[ 3 ], float ver2[ 3 ], float normal[ 3 ], bool isunit )
{
	if( !isunit )
	{
		MyMath::normalize( normal );
	}

	//compute the dot product between the direction of the line and the plane
	float dir[ 3 ];
	MyMath::getVec( ver1, ver2, dir );
	MyMath::normalize( dir );
	
	float dotp = MyMath::dotProduct( dir, normal );
	return MyMath::isEqualInToler( dotp,0, TOLERANCE_SIX);
}
