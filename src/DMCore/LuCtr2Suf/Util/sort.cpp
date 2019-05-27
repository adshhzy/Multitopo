#include "../Util/Sort.h"

void selectSortInc( int*& sa, int size)
{
	for( int i = 0; i < size; i ++)
	{
		int mini = i;
		for( int j = i + 1; j < size; j ++)
		{
			if( sa[ mini ] > sa[ j ])
			{
				mini = j;
			}
		}
		if( mini == i) continue;
		int t = sa[ mini ];
		sa[ mini ] = sa[ i ];
		sa[ i ] = t;
	}
}