#ifndef _TAGMAPPING_H_
#define _TAGMAPPING_H_

class TagMapping{
public:

	enum GROUPINFO {GI_ISO = 0, GI_SCORE = 1};
	enum ETETMARK {E_INTERIOR = 0, E_BD = 1, E_FACE = 3};
	enum NXTSEGV {NXTSEGV_SHARE = -2, NXTSEGV_NAN = -1};
	enum MAT {MAT_IN=1, MAT_OUT=0, MAT_NAN=-1, MAT_SEG=2};

	static inline bool isEdgeOnSeg(int edgemark){ return edgemark < 0; }
	static inline bool isFaceOnBB(int facemark){ return facemark > 0; }
	static inline int getSegiFromEdgemark(int edgemark){ return -edgemark-10;}
	static inline int getFacetiFromFacemark(int facemark){ return facemark-100;}

	static inline int getSegiTag(int segi){return segi*2 + 1;}
	static inline int tetCellVTag(){return -1;}
	static inline int tetSegVTag(int segii){return - 10 - segii;}
	static inline int tetSegETag(int segi){return - 10 - segi;}
	static inline int tetFaceTag(int fi){return 100 + fi;}

	static inline bool isTetVsOnSeg(int tag){return tag == MAT_SEG;}
	static inline int InteriorSegVTag2Segi(int tag){return  (int)(( - tag - 11 ) / 2); }
	static inline bool isTetVsInterior(int tag){return tag == 0;}

	static inline int getOtherMat(int matmark, vector<int> & matMarks){ return matMarks[0]+matMarks[1]-matmark;}
};


#endif //_TAGMAPPING_H_
