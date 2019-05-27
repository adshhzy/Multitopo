#ifndef _SETOPERATION_H_
#define _SETOPERATION_H_

#include <vector>
#include <set>
using namespace std;

#include "../util/veredgeplaneop.h"

const double TOLERANCE_POINT_SAME = 0.000001;
//const TOLERANCE_SIX_ = 0.000001;

int findOverLapPoint( float* ptlist,int* ptmark, int ptnum, float pt[ 3 ], int validmark);
int findPosInArray( int* array, int len, int val);
int commonOfTwoIncSortedArray_PosArray( int* array1, int len1, int* array2, int len2);
//return the position in array1, if exist, return the common value, otherwise, -1

int commonOfTwoIncSortedArray_Ind1( int* array1, int len1, int* array2, int len2);
//return the position in array1, if exist, the position is non-negative, otherwise, -1

int findPosInIncSortedArray( int* desarray,int len, int val );
//if exists, return the position, otherwise, return -1

bool findOneCommon(int* list1, int* list2, int listlen1, int listlen2, int& common);

void merge( vector<int>& set1, vector<int>& set2);
//merge set1 and set2, and save the result into set1

void merge( vector<int>& set1, int* set2, int size2);
void merge( set<int>& set1, set<int>& set2);
bool isSetEqual( vector<int>& set, int startindex, int len, int* subset);

#endif