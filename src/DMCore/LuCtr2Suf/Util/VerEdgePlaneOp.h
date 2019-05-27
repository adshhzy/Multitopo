#ifndef _VEREDGEPLANEOP_H_
#define _VEREDGEPLANEOP_H_

#include "../Math/mymath.h"
#include "../config.h"

//compute the common line between two intersecting planes
void computeComnLine( float param1[ 4 ], float param2[ 4 ], float comndir[ 3 ], float comnpt[ 3 ] );

//for 2D
float getAngleVecX2D(float vec[ 2 ], float xaxis[2] );

float getAngleVecX(float vec[ 3 ], float xaxis[ 3 ], float zaxis[ 3 ]);
void getDirOfTrian(float pt1[ 3 ], float pt2[ 3 ], float pt3[ 3 ], float dir[ 3 ]);
bool isSameDir( float dir1[ 3 ], float dir2[ 3 ]);
int getMajorDirec( float vec[ 3 ]);
int getMajorDirec( float pt1[ 3 ], float pt2[ 3 ]);
bool pointEqualOriginCube( float pos1[ 3 ], float toler );
bool pointEqualCube( float pos1[ 3 ], float pos2[ 3 ], float toler );
bool pointEqualSphere(float pos1[ 3 ], float pos2[ 3 ], float toler);
bool onSameLine( float pt1[3], float pt2[ 3 ], float pt3[ 3 ]);
float getParamOnSeg(float pt1[ 3 ], float pt2[ 3 ], float pt[ 3 ]);
void projectionPtOnPlane(float pt[ 3 ],  float planeparam[ 4 ], bool normisunit, float projectionpt[ 3 ]);
//void projectionPtOnPlane(float pt[ 3 ],  float planenorm[ 3 ], bool normisunit,float projectionpt[ 3 ]);
float getTOnRay( float startpt[ 3 ], float dir[ 3 ], float pt[ 3 ]);
float distPt2PlaneAlongDir( float pt[ 3 ], float dir[ 3 ],bool dirisunit, float planept[ 3 ], float planenorm[ 3 ], bool normisunit);
void interPt_PlaneEdge( float param[ 4 ], float ver1[ 3 ], float ver2[ 3 ], float interpt[ 3 ]);
bool isParalLinePlane( float ver1[ 3 ], float ver2[ 3 ], float normal[ 3 ], bool isunit );


//for debugging reason
void writeComnLineDB(float param1[ 4 ], float param2[ 4 ], float comndir[ 3 ], float comnpt[ 3 ]);

#endif