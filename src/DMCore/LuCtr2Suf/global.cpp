#include <vector>
#include <set>
using namespace std;

float DIM = 1.5;
bool _global_IsCtrGraph;
//the set of planes are share one common line
vector<set<int> > _gloabl_radialPlaneSet;	//finally only those with length >= 3 are considered.