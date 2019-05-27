#include "../Util/SetOperation.h"

//if there exists one point in the ptlist with mark as validmark and close to the pt, return its position
//if no such point, return -1
int findOverLapPoint( float* ptlist,int* ptmark, int ptnum, float pt[ 3 ], int validmark)
{
	for( int i = 0; i< ptnum; i ++)
	{
		if( ptmark[ i ] != validmark )
			continue;
		if(pointEqualSphere( ptlist + 3*i, pt, TOLERANCE_SAME_VER))
			return i;
	}
	return -1;
}
//find position in the array with length len. return the position of it if found, otherwise -1
int findPosInArray( int* array, int len, int val)
{
	for( int i = 0;i < len; i ++)
		if( array[ i ] == val )
			return i;
	return -1;
}
int commonOfTwoIncSortedArray_PosArray( int* array1, int len1, int* array2, int len2)
{
	int i = 0;
	int j = 0;
	while( (i < len1) && (j < len2) )
	{
		if( array1[ i ] == array2[ j ] )
			return array1[ i ];
		if( array1[ i ] < array2[ j ])
			i++;
		else
			j++;
	}
	return -1;
}
int commonOfTwoIncSortedArray_Ind1( int* array1, int len1, int* array2, int len2)
{
	int i = 0;
	int j = 0;
	while( (i < len1) && (j < len2) )
	{
		if( array1[ i ] == array2[ j ] )
			return i;
		if( array1[ i ] < array2[ j ])
			i++;
		else
			j++;
	}
	return -1;
}

int findPosInIncSortedArray( int* desarray, int len, int val )
{
	int begin = 0;
	int end = len - 1;
	int mid = (begin + end )/2;
	while ( begin <= end )		
	{
		if( desarray[ mid ] == val )
			return mid;
		if( desarray[ mid ] < val )
		{
			begin = mid + 1;
		}
		else
			end = mid - 1;
		mid = ( begin + end )/2;
	}
	return -1;
}
void merge( set<int>& set1, set<int>& set2)
{
	set<int>::iterator iter = set2.begin();
	for( ; iter != set2.end(); iter++)
	{
		set1.insert( *iter );
	}
}

bool findOneCommon(int* list1, int* list2, int listlen1, int listlen2, int& common)
{
	for( int i = 0; i < listlen1; i ++)
	{
		common = list1[ i ];
		for( int j = 0; j < listlen2; j ++)
		{
			if( common == list2[ j ])
				return true;
		}
	}
	return false;
}
void merge( vector<int>& set1, vector<int>& set2)
{
	int size1 = set1.size();
	int size2 = set2.size();
	for( int i = 0; i < size2; i ++)
	{
		bool isin = false;
		for(int j = 0; j < size1; j++)
		{
			if( set2[ i ] == set1[ j ])
			{
				isin = true;
				break;
			}
		}
		if( !isin )
		{
			set1.push_back( set2[ i ]);
		}
	}
}
void merge( vector<int>& set1, int* set2, int size2)
{
	int size1 = set1.size();
	for( int i = 0; i < size2; i ++)
	{
		bool isin = false;
		for(int j = 0; j < size1; j++)
		{
			if( set2[ i ] == set1[ j ])
			{
				isin = true;
				break;
			}
		}
		if( !isin )
		{
			set1.push_back( set2[ i ]);
		}
	}
}
bool isSetEqual( vector<int>& set, int startindex, int len, int* subset)
{
	int* vermark = new int[ len ];
	for( int i = 0; i < len; i ++)
		vermark[ i ] = 0;

	for( int i = 0; i < len; i++)
	{
		bool equalone = false;
		for( int j = 0; j < len; j++)
		{
			if( subset[ i ] == set[ j + startindex ])
			{
				vermark[ j ] = 1;
				equalone = true;
				break;
			}
		}
		if( !equalone )
		{
			delete []vermark;
			return false;
		}
	}
	for( int i = 0; i < len; i ++)
		if( vermark[ i ] == 0 )
		{
			delete []vermark;
			return false;
		}
	delete []vermark;
	return true;
}