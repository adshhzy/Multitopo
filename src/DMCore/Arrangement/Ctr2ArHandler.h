
/*
Ctr2ArHandler.h
Ming Zou @ wustl 2014

This handler reads a .contour file and outputs a space arrangement, which is
the basic data structure used later by the cycle grouping algorithm.
File reading and space partition are done by using the same functions in Lu Liu's
Ctr2Suf project (http://www.cs.wustl.edu/~taoju/lliu/paper/ctr2suf/program.html)

Note: add _CRT_SECURE_NO_WARNINGS to the project setting
(C++ preprocessor definitions) to run Lu's original code
*/

#ifndef _CTR2ARHANDLER_H_
#define _CTR2ARHANDLER_H_

#ifndef _DEBUG_LU
#define _DEBUG_LU 0
#endif

#include <vector>
#include <string.h>
#include "Arrangement.h"
#include "../LuCtr2Suf/config.h"
#include "../LuCtr2Suf/ContourHandler/ContourHandler.h"
#include "../LuCtr2Suf/SpacePartitioner/SpacePartitioner.h"

using namespace std;

class Ctr2ArHandler{
public:
	Ctr2ArHandler();
	~Ctr2ArHandler();

	void static readContour(const char* filename, const char* volfilename, const char* bboxfilename, Arrangement &ar);
    void static processContourMM(const char* filename, Arrangement &ar);

private:
};


#endif //_Ctr2ArHandler_H_
