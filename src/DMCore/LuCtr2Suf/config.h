#ifndef _CONFIG_H_
#define _CONFIG_H_

// this header fixes a number of platform specific problems
// make sure you include it at the top of every file you write

//#ifdef _MSC_VER
//// turn off symbol length warnings
//#pragma warning (disable: 4786)
//#pragma warning (disable: 4503)
////fix non-compliant 'for' scoping
//#define for if(false) {} else for 
//// conversion double -> float
//#pragma warning (disable: 4244)
//// truncation double -> float
//#pragma warning (disable: 4305)
//#endif

// define usefull constants

#ifndef M_PI
#define M_PI               3.1415926535897932384626433832795
#endif

#include <cassert>
#include <iostream>
//#include <FL/gl.h>
//#include <fl/math.h>
#include <iostream>
#include <fstream>
using namespace std;

//#include "./math/mymath.h"
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <algorithm>

typedef vector<int> intvector;
typedef vector<float> floatvector;
typedef set<int> intset;
const float TOLERANCE_HALF = 0.5f;
const float TOLERANCE_POS_ONE = 1.0f;
const float TOLERANCE_SIX = 0.000001f;
const float TOLERANCE_THREE = 0.001f;
const float TOLERANCE_FOUR = 0.0001f;
const float TOLERANCE_TWO = 0.01f;
const float TOLERANCE_ONE = 0.1f;
const float TOLERANCE_ANGLE = 0.01f;
const float POINT_ON_PLANE_TOLERANCE = 0.000001f;
const float PLANE_PARAM_PARALLEL_TOLERANCE = 0.000001f;

const int CTRTYPE_SUBVER = 1;
const int CTRTYPE_SUBEDGE = 2;
const int CTRTYPE_SUBFACE = 3;
const int CTRTYPE_MAJPT = 1;
const int CTRTYPE_MASEAM = 2;
const int CTRTYEP_MASHEET = 3;
const int NOT_EXIST = -1;

//expected length for contour edge when preprocessing
const float EXPECTLEN = 20;
const float TOLERANCE_SAME_VER = 0.01f;

#ifdef LINUX
using std::min;
using std::max;
#endif

#define debug
const float PROCSIZE = 1000;

extern void RedrawWindow();

//global variables
extern float DIM;
extern bool _global_IsCtrGraph;
extern vector<set<int> > _gloabl_radialPlaneSet;

#endif /* _CONFIG_H_ */
