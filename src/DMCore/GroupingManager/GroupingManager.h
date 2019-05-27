#ifndef _GROUPINGMANAGER_H_
#define _GROUPINGMANAGER_H_

#include "../Arrangement/Ctr2ArHandler.h"
#include "../Arrangement/Arrangement.h"
#include "../Utility/Utility.h"
#include "../Utility/UnionFind_o.h"
#include "../Utility/AlgoX.h"
#include <unordered_map>
#include <queue>
#include <algorithm>

using namespace std;

typedef BigInt GState;
typedef BigInt GGrouping;
typedef int GSolID;

typedef unordered_map<GState, GSolID, myHash, myEq> GStateHash;

class SufAttribute{
public:
	int sum_Bi;
	int sum_gi;
	int L;
	int k;

	SufAttribute() : sum_Bi(0), sum_gi(0), L(0), k(0){};
	inline void operator+= (const SufAttribute &att);

	SufAttribute(const SufAttribute& _sa) :
		sum_Bi(_sa.sum_Bi),
		sum_gi(_sa.sum_gi),
		L(_sa.L), k(_sa.k){};
	inline SufAttribute & operator = (const SufAttribute &_sa);
	inline int calcGenus(int B);
};

void SufAttribute::operator+=(const SufAttribute &att){
	sum_Bi += att.sum_Bi;
	sum_gi += att.sum_gi;
	L += att.L;
	k += att.k;
}
SufAttribute & SufAttribute::operator = (const SufAttribute &_sa){
	sum_Bi = _sa.sum_Bi;
	sum_gi = _sa.sum_gi;
	L = _sa.L;
	k = _sa.k;
	return *this;
};
int SufAttribute::calcGenus(int B){
	return (int)(sum_gi + 1 - k + 0.5*(sum_Bi - B + L));
}

class GroupingManager{
public:
	Arrangement ar;
	vector<int> Order;
	int goalGenus;
	int maxGenus;

	int maxGenus_log;
	vector<int> nfrontcycs_log;
	vector<int> ncellcycs_log;

	void Contour2Arrangement(const char* ctrfilename, const char* volfilename, const char* bboxfilename);
    void Contour2ArrangementMM(const char* ctrfilename, const char* volfilename, const char* bboxfilename);

	void CalcSearchOrder();
	void Grouping();
	void GenSurface();
	void SmoothSurface();
	void SaveResBefSmooth();

	GroupingManager();
	~GroupingManager();

	vector<vector<int> > solGrouping;


	void writeSearchOrder(const char* filename);
	void writeSolGrouping(const char* filename);
	void writeOutputs(const char* filename, vector<float> &data);
	inline void setGenus(int _goalGenus);

    inline void setPara(const char *_outDir, float _teta, float _randwb, int _nSmBefLoop, int _nLoop, int _nSmInLoop);
	inline bool solFound(){
		return !solGrouping.empty();
	}
	int t_celli;
	float t_teta;
	float t_randwb;
	int t_nSmBefLoop;
	int t_nLoop;
	int t_nSmInLoop;
	string outDir;
	vector<int> resGrouping;
private:

	GState initState;
	GState lastState;
	void Propagate(GState &f0State, int fronti, vector<GState> &f1States,
		GStateHash &f0StateHash, GStateHash &f1StateHash,
		vector<vector<int> > &f0SolGrouping,
		vector<vector<int> > &f1SolGrouping,
		vector<float> &f0SolGroupingCost,
		vector<float> &f1SolGroupingCost);

	void GenerateAllSingleDoubleTripleCurveSets(int nCyc, vector<vector<int> > &sets);
	void CalculateSetCost(int celli, vector<vector<int> > &sets, vector<float> &weights);
	void DetectConflicts(int celli, vector<vector<int> > &sets, vector<vector<int> > &conflicts);

	void UnPackGState(int fronti, GState &state, vector<vector<int> > &suf2cycid, vector<int> &suf2genus);
	void PackGState(int fronti, vector<int> &finalSufs, vector<int> &finalSufs_genus, vector<BigInt> &sufBoundaries, GState &nxtState);

	void FindUniqueOrderMap(int maxSi, vector<int> &Ss, vector<int> &map);

	void PackGGrouping(int celli, vector<int> &curGrouping, vector<vector<int> > &sets, GGrouping &ggroup);
	void UnPackGGrouping(int celli, GGrouping &ggroup, vector<vector<int> > &grouping);
};

inline void GroupingManager::setGenus(int _goalGenus){
	goalGenus = _goalGenus;
	ar.setPara(outDir, t_teta, t_randwb, t_nSmBefLoop, t_nLoop, t_nSmInLoop, goalGenus);
}

inline void GroupingManager::setPara(const char* _outDir, float _teta, float _randwb, int _nSmBefLoop, int _nLoop, int _nSmInLoop){
    //outDir = string(_outDir)+"/";
    outDir = string(_outDir);
	t_teta = _teta;
	t_randwb = _randwb;
	t_nSmBefLoop = _nSmBefLoop;
	t_nLoop = _nLoop;
	t_nSmInLoop = _nSmInLoop;
}

#endif //_GROUPINGMANAGER_H_
