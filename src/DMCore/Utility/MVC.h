#ifndef MVC_H
#define MVC_H

/************************************************************************/
/* Mean Value Coordinates                                               */
/************************************************************************/

#include <math.h>
#include <vector>
#include <stdlib.h>
#include <assert.h>

using namespace std;




class MVC
{
	// Mesh
	vector<vector<float> > & vs;
	vector<vector<int> > & ts;
	vector<bool> & flags;

	int nv, nt ;
	float ** vert;
	float * mags;
	float * weights;
	

public:
	// Constructor
	MVC( vector<vector<float> > & _verts, vector<vector<int> > & _tris, vector<bool> & _trisFlag )	: vs(_verts), ts(_tris), flags(_trisFlag)
	{	
		nv = vs.size() ;
		nt = ts.size() ;

		vert = new float* [nv] ;
		mags = new float [ nv ] ;
		weights = new float [ nv ];
		for ( int i = 0 ; i < nv ; i ++ )
		{
			vert[ i ] = new float[ 3 ] ;
		}
	}

	~MVC(){
		for(int i=0; i<nv; ++i){
			delete [] vert[i];
		}
		delete [] vert;
		delete [] mags;
		delete [] weights;
	}

	// Compute weights
	float getValue ( const float x[3])
	{
		int i, j;
		float* v [ 3 ];
		float a [ 3 ], w [ 3 ], totalW = 0, totalWF = 0, totalA, vol;
		float l [ 3 ], t [ 3 ], s, sinS, sinT [ 3 ], cs [ 3 ], sinG [ 3 ];
		int *f;
		double EPSILN = 0.00001;
		float mag;

		float numerator = 0;
		
		for ( i = 0; i < nv; i++ )
		{
			weights [ i ] = 0;
		}
		
		for ( i = 0; i < nv; i++ )
		{
			vert [ i ][ 0 ] = vs[ i ][ 0 ] - x[ 0 ];
			vert [ i ][ 1 ] = vs[ i ][ 1 ] - x[ 1 ];
			vert [ i ][ 2 ] = vs[ i ][ 2 ] - x[ 2 ];

			mag = vert[i][0] * vert[i][0] + vert[i][1] * vert[i][1] + vert[i][2] * vert[i][2] ;
			if ( mag < EPSILN * EPSILN )
			{
				weights [ i ] = 1;
				return 0;
			}
			mag = (float)sqrt ( mag );
			vert [ i ][ 0 ] /= mag;
			vert [ i ][ 1 ] /= mag;
			vert [ i ][ 2 ] /= mag;
			mags[ i ] = mag ;
		}
		
		for ( i = 0; i < nt; i++ )
		{
			v [ 0 ] = vert [ ts[ i ] [ 0 ] ];
			v [ 1 ] = vert [ ts[ i ] [ 1 ] ];
			v [ 2 ] = vert [ ts[ i ] [ 2 ] ];
			
			l[0] = l[1] = l[2] = 0 ;
			for ( j = 0 ; j < 3 ; j ++ )
			{
				l [ 0 ] += ( v [ 2 ][ j ] - v [ 1 ][ j ] ) * ( v [ 2 ][ j ] - v [ 1 ][ j ] ) ;
				l [ 1 ] += ( v [ 0 ][ j ] - v [ 2 ][ j ] ) * ( v [ 0 ][ j ] - v [ 2 ][ j ] ) ;
				l [ 2 ] += ( v [ 1 ][ j ] - v [ 0 ][ j ] ) * ( v [ 1 ][ j ] - v [ 0 ][ j ] ) ;
			}
			l[ 0 ] = sqrt( l[ 0 ] ) ;
			l[ 1 ] = sqrt( l[ 1 ] ) ;
			l[ 2 ] = sqrt( l[ 2 ] ) ;
			//assert ( _finite ( l [ 0 ] ) && _finite ( l [ 1 ] ) && _finite ( l [ 2 ] ) );
			
			t [ 0 ] = 2 * (float)asin ( l [ 0 ] / 2 );
			t [ 1 ] = 2 * (float)asin ( l [ 1 ] / 2 );
			t [ 2 ] = 2 * (float)asin ( l [ 2 ] / 2 );
			//assert ( _finite ( t [ 0 ] ) && _finite ( t [ 1 ] ) && _finite ( t [ 2 ] ) );
			
			s = 0.5f * ( t [ 0 ] + t [ 1 ] + t [ 2 ] );
			//assert ( _finite ( s ) );
			
			if ( !(3.141592654 - s < EPSILN) )
			{
				sinS = (float)sin ( s );
				//assert ( _finite ( sinS ) );
				sinT [ 0 ] = (float)sin ( t [ 0 ] );
				sinT [ 1 ] = (float)sin ( t [ 1 ] );
				sinT [ 2 ] = (float)sin ( t [ 2 ] );
				//assert ( _finite ( sinT [ 0 ] ) && _finite ( sinT [ 1 ] ) && _finite ( sinT [ 2 ] ) );
				
				cs [ 0 ] = (float)(( 2 * sinS * sin ( s - t [ 0 ] ) ) / ( sinT [ 2 ] * sinT [ 1 ] ) - 1);
				cs [ 1 ] = (float)(( 2 * sinS * sin ( s - t [ 1 ] ) ) / ( sinT [ 0 ] * sinT [ 2 ] ) - 1);
				cs [ 2 ] = (float)(( 2 * sinS * sin ( s - t [ 2 ] ) ) / ( sinT [ 1 ] * sinT [ 0 ] ) - 1);
				//assert ( _finite ( cs [ 0 ] ) && _finite ( cs [ 1 ] ) && _finite ( cs [ 2 ] ) );
				
				//			if ( fabs ( vol ) > EPSILN * EPSILN )
				if ( ( 1 - cs [ 0 ] ) > EPSILN && ( 1 - cs [ 1 ] ) > EPSILN && ( 1 - cs [ 2 ] ) > EPSILN )
				{
					float cv[ 3 ] ;
					cv[ 0 ] = v[ 1 ][ 1 ] * v[ 2 ][ 2 ] - v[ 1 ][ 2 ] * v[ 2 ][ 1 ] ;
					cv[ 1 ] = v[ 1 ][ 2 ] * v[ 2 ][ 0 ] - v[ 1 ][ 0 ] * v[ 2 ][ 2 ] ;
					cv[ 2 ] = v[ 1 ][ 0 ] * v[ 2 ][ 1 ] - v[ 1 ][ 1 ] * v[ 2 ][ 0 ] ;

					vol = v[ 0 ][ 0 ] * cv[ 0 ] + v[ 0 ][ 1 ] * cv[ 1 ] + v[ 0 ][ 2 ] * cv[ 2 ] ;
						
					if ( vol < 0 )
					{
						sinG [ 0 ] = -1 * sqrt ( 1 - cs [ 0 ] * cs [ 0 ] );
						sinG [ 1 ] = -1 * sqrt ( 1 - cs [ 1 ] * cs [ 1 ] );
						sinG [ 2 ] = -1 * sqrt ( 1 - cs [ 2 ] * cs [ 2 ] );
					}
					else
					{
						sinG [ 0 ] = sqrt ( 1 - cs [ 0 ] * cs [ 0 ] );
						sinG [ 1 ] = sqrt ( 1 - cs [ 1 ] * cs [ 1 ] );
						sinG [ 2 ] = sqrt ( 1 - cs [ 2 ] * cs [ 2 ] );
					}
					
					w [ 0 ] = (float)(( t [ 0 ] - cs [ 1 ] * t [ 2 ] - cs [ 2 ] * t [ 1 ] ) / ( sinT [ 1 ] * sinG [ 2 ] ));
					w [ 1 ] = (float)(( t [ 1 ] - cs [ 2 ] * t [ 0 ] - cs [ 0 ] * t [ 2 ] ) / ( sinT [ 2 ] * sinG [ 0 ] ));
					w [ 2 ] = (float)(( t [ 2 ] - cs [ 0 ] * t [ 1 ] - cs [ 1 ] * t [ 0 ] ) / ( sinT [ 0 ] * sinG [ 1 ] ));
					//assert ( _finite ( w [ 0 ] ) && _finite ( w [ 1 ] ) && _finite ( w [ 2 ] ) );
					
					w [ 0 ] /= mags[ ts[ i ][ 0 ] ] ;
					w [ 1 ] /= mags[ ts[ i ][ 1 ] ] ;
					w [ 2 ] /= mags[ ts[ i ][ 2 ] ] ;
					
					weights [ ts[ i ][ 0 ] ] += w [ 0 ];
					weights [ ts[ i ][ 1 ] ] += w [ 1 ];
					weights [ ts[ i ][ 2 ] ] += w [ 2 ];

					if (flags[i])
					{
						numerator += (w[0] + w[1] + w[2]);
					}
					totalW += ( w [ 0 ] + w [ 1 ] + w [ 2 ] );
				}
				else
				{
					continue;
				}
			}
			else
			{
				a [ 0 ] = (float)sin ( t [ 0 ] ) * mags[ ts[ i ][ 2 ] ] * mags[ ts[ i ][ 1 ] ];
				a [ 1 ] = (float)sin ( t [ 1 ] ) * mags[ ts[ i ][ 0 ] ] * mags[ ts[ i ][ 2 ] ];
				a [ 2 ] = (float)sin ( t [ 2 ] ) * mags[ ts[ i ][ 1 ] ] * mags[ ts[ i ][ 0 ] ];
				
				totalA = a [ 0 ] + a [ 1 ] + a [ 2 ];
				
				// we're on a face
				a [ 0 ] /= totalA;
				a [ 1 ] /= totalA;
				a [ 2 ] /= totalA;
				
				for ( j = 0; j < nv; j++ )
				{
					weights [ j ] = 0;
				}
				
				weights [ ts[ i ][ 0 ] ] = a [ 0 ];
				weights [ ts[ i ][ 1 ] ] = a [ 1 ];
				weights [ ts[ i ][ 2 ] ] = a [ 2 ];
				return 0;
			}
		}
		
		/*for ( j = 0; j < nv; j++ )
		{
			weights [ j ] /= totalW;
		}*/
		
		return numerator / totalW ;
		/*if ( totalW > 0 )
		{
			return numerator / totalW ;
		}
		else
		{
			return -1 ;
		}*/
	}

	int _finite( float x )
	{
		if ( fabs(x) < 1000 )
		{
			return 1 ;
		}
		return 0 ;
	}
};


class MVC_TAO
{
	// Mesh
	float ** vs;
	int ** ts;
	int nv, nt ;
	float ** vert;
	float * mags;
	bool * flags;

	float * weights;

public:
	// Constructor
	MVC_TAO( int numv, int numt, float ** verts, int ** tris, bool * trisFlag )	
	{	
		vs = verts ;
		ts = tris ;
		nv = numv ;
		nt = numt ;
		flags = trisFlag;

		vert = new float* [numv] ;
		mags = new float [ numv ] ;

		weights = new float [ numv ] ;
		for ( int i = 0 ; i < numv ; i ++ )
		{
			vert[ i ] = new float[ 3 ] ;
		}
	} ;

	~MVC_TAO(){
		for(int i=0; i<nv; ++i){
			delete [] vert[i];
		}
		delete [] vert;
		delete [] mags;
		delete [] weights;
	}

	// Compute weights
	float /*getWeights*/getValue ( float x[3]/*, float* weights */)
	{
		int i, j;
		float* v [ 3 ];
		float a [ 3 ], w [ 3 ], totalW = 0, totalWF = 0, totalA, vol;
		float l [ 3 ], t [ 3 ], s, sinS, sinT [ 3 ], cs [ 3 ], sinG [ 3 ];
		int *f;
		double EPSILN = 0.00001;
		float mag;

		float numerator = 0.0f;
		
		for ( i = 0; i < nv; i++ )
		{
			weights [ i ] = 0;
		}
		
		for ( i = 0; i < nv; i++ )
		{
			vert [ i ][ 0 ] = vs[ i ][ 0 ] - x[ 0 ];
			vert [ i ][ 1 ] = vs[ i ][ 1 ] - x[ 1 ];
			vert [ i ][ 2 ] = vs[ i ][ 2 ] - x[ 2 ];

			mag = vert[i][0] * vert[i][0] + vert[i][1] * vert[i][1] + vert[i][2] * vert[i][2] ;
			if ( mag < EPSILN * EPSILN )
			{
				weights [ i ] = 1;
				return 0;
			}
			mag = (float)sqrt ( mag );
			vert [ i ][ 0 ] /= mag;
			vert [ i ][ 1 ] /= mag;
			vert [ i ][ 2 ] /= mag;
			mags[ i ] = mag ;
		}
		
		for ( i = 0; i < nt; i++ )
		{
			v [ 0 ] = vert [ ts[ i ] [ 0 ] ];
			v [ 1 ] = vert [ ts[ i ] [ 1 ] ];
			v [ 2 ] = vert [ ts[ i ] [ 2 ] ];
			
			l[0] = l[1] = l[2] = 0 ;
			for ( j = 0 ; j < 3 ; j ++ )
			{
				l [ 0 ] += ( v [ 2 ][ j ] - v [ 1 ][ j ] ) * ( v [ 2 ][ j ] - v [ 1 ][ j ] ) ;
				l [ 1 ] += ( v [ 0 ][ j ] - v [ 2 ][ j ] ) * ( v [ 0 ][ j ] - v [ 2 ][ j ] ) ;
				l [ 2 ] += ( v [ 1 ][ j ] - v [ 0 ][ j ] ) * ( v [ 1 ][ j ] - v [ 0 ][ j ] ) ;
			}
			l[ 0 ] = sqrt( l[ 0 ] ) ;
			l[ 1 ] = sqrt( l[ 1 ] ) ;
			l[ 2 ] = sqrt( l[ 2 ] ) ;
			//assert ( _finite ( l [ 0 ] ) && _finite ( l [ 1 ] ) && _finite ( l [ 2 ] ) );
			
			t [ 0 ] = 2 * (float)asin ( l [ 0 ] / 2 );
			t [ 1 ] = 2 * (float)asin ( l [ 1 ] / 2 );
			t [ 2 ] = 2 * (float)asin ( l [ 2 ] / 2 );
			//assert ( _finite ( t [ 0 ] ) && _finite ( t [ 1 ] ) && _finite ( t [ 2 ] ) );
			
			s = 0.5f * ( t [ 0 ] + t [ 1 ] + t [ 2 ] );
			//assert ( _finite ( s ) );
			
			if ( !(3.141592654 - s < EPSILN) )
			{
				sinS = (float)sin ( s );
				//assert ( _finite ( sinS ) );
				sinT [ 0 ] = (float)sin ( t [ 0 ] );
				sinT [ 1 ] = (float)sin ( t [ 1 ] );
				sinT [ 2 ] = (float)sin ( t [ 2 ] );
				//assert ( _finite ( sinT [ 0 ] ) && _finite ( sinT [ 1 ] ) && _finite ( sinT [ 2 ] ) );
				
				cs [ 0 ] = (float)(( 2 * sinS * sin ( s - t [ 0 ] ) ) / ( sinT [ 2 ] * sinT [ 1 ] ) - 1);
				cs [ 1 ] = (float)(( 2 * sinS * sin ( s - t [ 1 ] ) ) / ( sinT [ 0 ] * sinT [ 2 ] ) - 1);
				cs [ 2 ] = (float)(( 2 * sinS * sin ( s - t [ 2 ] ) ) / ( sinT [ 1 ] * sinT [ 0 ] ) - 1);
				//assert ( _finite ( cs [ 0 ] ) && _finite ( cs [ 1 ] ) && _finite ( cs [ 2 ] ) );
				
				//			if ( fabs ( vol ) > EPSILN * EPSILN )
				if ( ( 1 - cs [ 0 ] ) > EPSILN && ( 1 - cs [ 1 ] ) > EPSILN && ( 1 - cs [ 2 ] ) > EPSILN )
				{
					float cv[ 3 ] ;
					cv[ 0 ] = v[ 1 ][ 1 ] * v[ 2 ][ 2 ] - v[ 1 ][ 2 ] * v[ 2 ][ 1 ] ;
					cv[ 1 ] = v[ 1 ][ 2 ] * v[ 2 ][ 0 ] - v[ 1 ][ 0 ] * v[ 2 ][ 2 ] ;
					cv[ 2 ] = v[ 1 ][ 0 ] * v[ 2 ][ 1 ] - v[ 1 ][ 1 ] * v[ 2 ][ 0 ] ;

					vol = v[ 0 ][ 0 ] * cv[ 0 ] + v[ 0 ][ 1 ] * cv[ 1 ] + v[ 0 ][ 2 ] * cv[ 2 ] ;
						
					if ( vol < 0 )
					{
						sinG [ 0 ] = -1 * sqrt ( 1 - cs [ 0 ] * cs [ 0 ] );
						sinG [ 1 ] = -1 * sqrt ( 1 - cs [ 1 ] * cs [ 1 ] );
						sinG [ 2 ] = -1 * sqrt ( 1 - cs [ 2 ] * cs [ 2 ] );
					}
					else
					{
						sinG [ 0 ] = sqrt ( 1 - cs [ 0 ] * cs [ 0 ] );
						sinG [ 1 ] = sqrt ( 1 - cs [ 1 ] * cs [ 1 ] );
						sinG [ 2 ] = sqrt ( 1 - cs [ 2 ] * cs [ 2 ] );
					}
					
					w [ 0 ] = (float)(( t [ 0 ] - cs [ 1 ] * t [ 2 ] - cs [ 2 ] * t [ 1 ] ) / ( sinT [ 1 ] * sinG [ 2 ] ));
					w [ 1 ] = (float)(( t [ 1 ] - cs [ 2 ] * t [ 0 ] - cs [ 0 ] * t [ 2 ] ) / ( sinT [ 2 ] * sinG [ 0 ] ));
					w [ 2 ] = (float)(( t [ 2 ] - cs [ 0 ] * t [ 1 ] - cs [ 1 ] * t [ 0 ] ) / ( sinT [ 0 ] * sinG [ 1 ] ));
					//assert ( _finite ( w [ 0 ] ) && _finite ( w [ 1 ] ) && _finite ( w [ 2 ] ) );
					
					w [ 0 ] /= mags[ ts[ i ][ 0 ] ] ;
					w [ 1 ] /= mags[ ts[ i ][ 1 ] ] ;
					w [ 2 ] /= mags[ ts[ i ][ 2 ] ] ;
					
					weights [ ts[ i ][ 0 ] ] += w [ 0 ];
					weights [ ts[ i ][ 1 ] ] += w [ 1 ];
					weights [ ts[ i ][ 2 ] ] += w [ 2 ];

					if (flags[i])
					{
						numerator += (w[0] + w[1] + w[2]);
					}
					totalW += ( w [ 0 ] + w [ 1 ] + w [ 2 ] );
				}
				else
				{
					continue;
				}
			}
			else
			{
				a [ 0 ] = (float)sin ( t [ 0 ] ) * mags[ ts[ i ][ 2 ] ] * mags[ ts[ i ][ 1 ] ];
				a [ 1 ] = (float)sin ( t [ 1 ] ) * mags[ ts[ i ][ 0 ] ] * mags[ ts[ i ][ 2 ] ];
				a [ 2 ] = (float)sin ( t [ 2 ] ) * mags[ ts[ i ][ 1 ] ] * mags[ ts[ i ][ 0 ] ];
				
				totalA = a [ 0 ] + a [ 1 ] + a [ 2 ];
				
				// we're on a face
				a [ 0 ] /= totalA;
				a [ 1 ] /= totalA;
				a [ 2 ] /= totalA;
				
				for ( j = 0; j < nv; j++ )
				{
					weights [ j ] = 0;
				}
				
				weights [ ts[ i ][ 0 ] ] = a [ 0 ];
				weights [ ts[ i ][ 1 ] ] = a [ 1 ];
				weights [ ts[ i ][ 2 ] ] = a [ 2 ];
				return 0;
			}
		}
		
		for ( j = 0; j < nv; j++ )
		{
			weights [ j ] /= totalW;
		}
		
		if ( totalW > 0 )
		{
			return numerator / totalW ;
		}
		else
		{
			return -1 ;
		}
	};

	int _finite( float x )
	{
		if ( fabs(x) < 1000 )
		{
			return 1 ;
		}
		return 0 ;
	}
};




#endif