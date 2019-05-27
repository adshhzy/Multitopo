#include "MMLevelSet.h"
#include <algorithm>
bool operator<(const EventPoint& lh, const EventPoint& rh) {
    return lh._lamda < rh._lamda;
}
bool operator>(const EventPoint& lh, const EventPoint& rh) {
    return lh._lamda > rh._lamda;
}
bool operator==(const EventPoint& lh, const EventPoint& rh) {
    return lh._lamda == rh._lamda;
}
static bool operator==(const CellTopology& lh, const CellTopology& rh) {
    // return (lh._comps == rh._comps) && (lh._conns == rh._conns);
    //return (lh._comps == rh._comps);
    //return false;
    //return (lh._comps == rh._comps) &&(lh.comp2VEF==rh.comp2VEF) ;
    return (lh._comps == rh._comps) &&(lh.comp2VEF==rh.comp2VEF) &&(lh._junctionpoints == rh._junctionpoints);
}
static bool operator<(const CellTopology& lh, const CellTopology& rh) {
    // TODO is this correct
    // return (lh._comps < rh._comps) && (lh._conns < rh._conns);
    //return false;
    //return (lh._comps < rh._comps ||(lh._comps == rh._comps && lh._conns < rh._conns));
    //return !((lh._comps == rh._comps) &&(lh.comp2VEF==rh.comp2VEF) &&(lh._junctionpoints == rh._junctionpoints));
    //return !((lh._comps == rh._comps) &&(lh.comp2VEF==rh.comp2VEF));


    if(lh._comps != rh._comps)return lh._comps < rh._comps;
    if(lh._junctionpoints != rh._junctionpoints)return lh._junctionpoints < rh._junctionpoints;
    return lh.comp2VEF < rh.comp2VEF;

}
struct CellTopologyCompareWOJP
{
    bool operator() (const CellTopology& lh, const CellTopology& rh) const
    {
        if(lh._comps != rh._comps)return lh._comps < rh._comps;
        return lh.comp2VEF < rh.comp2VEF;
    }
//    bool operator() (const CellTopology& lh, const CellTopology& rh) const
//    {
//        return (lh._comps == rh._comps) &&(lh.comp2VEF==rh.comp2VEF) ;
//    }
//    bool operator==(const CellTopology& lh, const CellTopology& rh) const{

//        return (lh._comps == rh._comps) &&(lh.comp2VEF==rh.comp2VEF) ;
//        //return (lh._comps == rh._comps) &&(lh.comp2VEF==rh.comp2VEF) &&(lh._junctionpoints == rh._junctionpoints);
//    }
} ;

void MMLevelSet::GenerateRays(float density) {
    RayInt ray(_K, 0);
    set<RayInt> rays;
    GenerateRays_helper(density, 0, ray, rays);

    vector<RayInt> raysInt(rays.size());
    copy(rays.begin(), rays.end(), raysInt.begin());
    _rays.resize(rays.size(), Ray(_K));
    // Convert integer-based rays to float.
    for (int i = 0; i < rays.size(); ++i) {
        copy(raysInt[i].begin(), raysInt[i].end(), _rays[i].begin());
    }
    // Modify the length of each ray to have the other end (beside the origin)
    // touching the boundary of range space.
    // It is equavalent to make the max value in the ray representation to be 1.
    for (auto& oneray : _rays) {
        float maxVal = *max_element(oneray.begin(), oneray.end());
        for (auto& ele : oneray) {
            ele /= maxVal;
        }
    }
}
void MMLevelSet::GenerateRays_helper(float density, int pos, RayInt& ray,
                                     set<RayInt>& rays) {
    if (pos == _K) {
        // Filter case[0] : (0,0,..,0).
        if (*max_element(ray.begin(), ray.end()) == 0) return;
        // Filter case[1] : v1 = v0 + lamda*(1,1,..,1). e.g. (0,0,2) and (1,1,3).
        if (*min_element(ray.begin(), ray.end()) != 0) return;
        // Filter cases[2]: v1=lamda*v0. e.g. (0,0,2)and(0,0,1);
        int gcd = MingUtility::GCD(ray);
        RayInt ray_final(ray.size(), 0);
        for (int i = 0; i < ray.size(); ++i) {
            ray_final[i] = ray[i] / gcd;
        }
        rays.insert(ray_final);
        return;
    }
    for (int i = 0; i <= density; ++i) {
        ray[pos] = i;
        GenerateRays_helper(density, pos + 1, ray, rays);
    }
}
extern bool isCJP;
int MMLevelSet::ExploreTopologiesAlongRays(bool allowHighGenus, float density,
                                           float segLenLowBound) {
    // Generate the sampling rays with user specified density.
    GenerateRays(density);
    cout << "nRay = " << _rays.size() << endl;
    //if(_rays.size() >10000)return 0;
    Labeling labeling(_N, 0);
    _criticalPtsPerRay.resize(_rays.size());
    _topologyPtsPerRay.resize(_rays.size());

    map<CellTopology, pair<int, int>,CellTopologyCompareWOJP> cellTopo2RayNPt;
    map<CellTopology, int,CellTopologyCompareWOJP> cellTopo2junctionpoint;
    map<CellTopology, CellTopology,CellTopologyCompareWOJP> cellTopo2trueTopo;
    map<CellTopology, vector<vector<int>>> cellTopo2AllRayNPt;
    _badCellTopoPtID.resize(2);
    // For each ray, find the event points and further the critical segments.
    for (int ri = 0; ri < _rays.size(); ++ri) {
        if(ri>=50 && ri%50==0)cout<<ri<<endl;
        const auto& ray = _rays[ri];
        float rayLen = MingUtility::L2Norm(ray);
        vector<EventPoint> eventPts;
        // Find all the event points which are the intersection points of the ray
        // with the C(n,2) half-hyperplane of each vertex. Each point satisfies the
        // following equations/inequalities:
        // event point = S = lamda * ray; (1)
        // V2F[vi][i] + lamda * ray[i] == V2F[vi][j] + lamda * ray[j]; (2)
        // V2F[vi][i] + lamda * ray[i] >= V2F[vi][k] + lamda * ray[k]; (3:k!=i or j)
        // S <= upbound * (1,...,1); (4:all i)
        // use equation (2) to solve lamda then check the validness with the
        // inequalities in (3) and (4).
        for (int vi = 0; vi < _N; ++vi) {
            if (_bdVs[vi]) continue;
            for (int i = 0; i < _K; ++i) {
                for (int j = i + 1; j < _K; ++j) {
                    // Degenerated case, skip.
                    if (ray[i] - ray[j] == 0) continue;
                    float lamda = (_V2F[vi][j] - _V2F[vi][i]) / (ray[i] - ray[j]);
                    // Only extract event points along the direction of the ray.
                    if (lamda < 0.0) continue;
                    for (int k = 0; k < _K; ++k) {
                        if (lamda * ray[k] > _upbound) {
                            lamda = -1.0;
                            break;
                        }
                        if (k == i || k == j) continue;
                        if (lamda * (ray[i] - ray[k]) < _V2F[vi][k] - _V2F[vi][i]) {
                            lamda = -1.0;
                            break;
                        }
                    }
                    if (lamda >= 0.0) {
                        eventPts.push_back(EventPoint(lamda, vi, i, j));
                    }
                }
            }
        }
        sort(eventPts.begin(), eventPts.end());
        // TODO find critical segments
        vector<EventPoint>& criticalPts = _criticalPtsPerRay[ri];
        for (auto pt : eventPts) {
            float lamda = pt._lamda;
            int vi = pt._vi;
            const vector<int>& nbVs = _mesh->Get1RingNeighborVs(vi);
            OffsetVector offset = ray;
            for (int i = 0; i < _K; ++i) {
                offset[i] *= lamda;
            }
            OffsetVector2LabelingForVs(offset, nbVs, labeling);
            if (_mesh->isVCritical(vi, labeling, pt._mat0, pt._mat1)) {
                //                if(vi==1948){
                //                    for (int i = 0; i < _K; ++i)cout<<offset[i]<<' ';cout<<endl;
                //                }
                pt.offset = ray;
                for (int i = 0; i < _K; ++i) {
                    pt.offset[i] *= lamda;
                }
                //OffsetVector2Labeling(pt.offset, labeling);
                pt._score = OffsetVector2Score(pt.offset);
                criticalPts.push_back(pt);
                // cout << pt._score << endl;
            }
        }
        // cout << "ray[" << ri << "]: " << endl;
        // cout << "# eventPts = " << eventPts.size() << endl;
        // cout << "# critcPts = " << criticalPts.size() << endl;
        // cout << "# critcPts = " << _criticalPtsPerRay[ri].size() << endl;

        // Collect the critical segment/topologyPts along the ray.
        // The first critical segment (segment between origin and the first
        // criticalPt, represent the same topo as the one on origin) will always be
        // skipped. Instead, the center point/topology will be pushed into the first
        // ray in the beginning.
        if (ri == 0) {
            EventPoint pt(0.0, -1, -1, -1);
            pt.offset = OffsetVector(_K, 0.0);
            //OffsetVector2Labeling(pt.offset, pt.label);
            pt._score = OffsetVector2Score(pt.offset);
            _topologyPtsPerRay[ri].push_back(pt);
        }
        // Skip rays if no critical point is found.
        vector<EventPoint>& topoPts = _topologyPtsPerRay[ri];
        if (criticalPts.size() > 0) {
            // vector<EventPoint>& topoPts = _topologyPtsPerRay[ri];
            for (int pi = 0; pi < criticalPts.size(); ++pi) {
                float score = criticalPts[pi]._score;
                float lamda0 = criticalPts[pi]._lamda;
                float lamda1 =
                        pi + 1 == criticalPts.size() ? 1.0 : criticalPts[pi + 1]._lamda;
                float span = (lamda1 - lamda0) * rayLen;
                // Filter out segs shorter than the length lower bound.
                if (span < segLenLowBound) continue;
                float lamda = (lamda0 + lamda1) / 2.0;
               // float lamda = lamda1 - (lamda1 - lamda0) / 100.0;
                EventPoint pt(lamda, -1, -1, -1);
                pt._score = score;
                pt._span = span;
                pt.offset = ray;
                for (int i = 0; i < _K; ++i) {
                    pt.offset[i] *= lamda;
                }
                //OffsetVector2Labeling(pt.offset, pt.label);
                topoPts.push_back(pt);
            }
        }
        // TODO I haven't implement the boundary preprocessing for 2D mesh.
        // I may be not gonna implenent in the future, since real cases are all 3D.
        if (_mesh->_dim == 2) continue;

        Labeling label;
        for (int pti = 0; pti < topoPts.size(); ++pti) {
            const auto& pt = topoPts[pti];
            //const Labeling& label = pt.label;
            OffsetVector2Labeling(pt.offset, label);
            CellTopology cellTopo;
            cellTopo._score = pt._score;
            //int res = _mesh->ExtractLevelSets(allowHighGenus, label, cellTopo._comps,
            //   cellTopo._conns);
            int res = _mesh->ExtractLevelSets(allowHighGenus, label, cellTopo._comps,
                                              cellTopo._conns,cellTopo.comp2VEF,cellTopo._junctionpoints);
            if (res == 0) {
                // if (_mesh->ExtractLevelSets(_K, label, cellTopo._comps,
                // cellTopo._conns)) {
                // cellTopo.Unique();
                // for (auto& comp : cellTopo._comps) {
                //   sort(comp.begin(), comp.end());
                // }
                // sort(cellTopo._comps.begin(), cellTopo._comps.end());
                // cout << "found a comp:" << endl;
                // MingUtility::printVector(cellTopo._comps);
                // If a new cell topology is found, or same cell topology exist but has
                // lower score than current one, map the current cell topology to the
                // specific topologyPt pti on ray ri.
                if (cellTopo2RayNPt.find(cellTopo) == cellTopo2RayNPt.end()) {
                    cellTopo2RayNPt[cellTopo] = make_pair(ri, pti);
                    cellTopo2junctionpoint[cellTopo] = cellTopo._junctionpoints;
                    cellTopo2trueTopo[cellTopo] = cellTopo;
                } else {
                    const auto& rayNpt = cellTopo2RayNPt[cellTopo];
                    if(isCJP){
                        if(cellTopo2junctionpoint[cellTopo]>cellTopo._junctionpoints ||
                                (cellTopo2junctionpoint[cellTopo] == cellTopo._junctionpoints && _topologyPtsPerRay[rayNpt.first][rayNpt.second]._score <
                                 pt._score)){
                            cellTopo2RayNPt[cellTopo] = make_pair(ri, pti);
                            cellTopo2junctionpoint[cellTopo] = cellTopo._junctionpoints;
                            cellTopo2trueTopo[cellTopo] = cellTopo;

                        }
                    }else{

                        if (_topologyPtsPerRay[rayNpt.first][rayNpt.second]._score <
                                pt._score) {
                            cellTopo2RayNPt[cellTopo] = make_pair(ri, pti);
                            cellTopo2trueTopo[cellTopo] = cellTopo;
                        }
                    }
                }
                // For debug use.
                cellTopo2AllRayNPt[cellTopo].push_back({ri, pti});
            } else {
                if (res == -1) {
                    _badCellTopoPtID[0].push_back({ri, pti});
                } else {
                    _badCellTopoPtID[1].push_back({ri, pti});
                }
            }
        }
    }

    if(1){
        // Sort the topologies by its score.
        map<CellTopology, int, CellTopologyCompareWOJP>WOJP;
        _cellNexplored2Topogroup.clear();int nex = 0;
        typedef pair<float, map<CellTopology, pair<int, int>>::iterator> sitpair;
        vector<sitpair> scoreOrder;
        for (auto it = cellTopo2RayNPt.begin(); it != cellTopo2RayNPt.end(); ++it) {
            const auto& rayNpt = it->second;
            //float score = _topologyPtsPerRay[rayNpt.first][rayNpt.second]._score;
            float score = cellTopo2trueTopo[it->first]._score;
            //cout<<rayNpt.first<<' '<<rayNpt.second<<endl;
            scoreOrder.push_back(make_pair(score, it));
            if(WOJP.find(it->first)==WOJP.end())WOJP[it->first] = nex++;
        }

        sort(scoreOrder.begin(), scoreOrder.end(),
             [](const sitpair& lhs,
             const sitpair& rhs) { return lhs.first > rhs.first; });

        for (const auto& onepair : scoreOrder) {

            const auto& pr = *(onepair.second);
            const auto& rayNpt = pr.second;
            int index = _cellTopologies.size();
            //_cellTopologies.push_back(pr.first);
            _cellTopologies.push_back(cellTopo2trueTopo[pr.first]);
            _cellNexplored2Topogroup.push_back(WOJP[pr.first]);
            // auto& newTopo = _cellTopologies.back();
            auto& newTopo = _cellTopologies[index];
            newTopo._score = _topologyPtsPerRay[rayNpt.first][rayNpt.second]._score;
            newTopo._span = _topologyPtsPerRay[rayNpt.first][rayNpt.second]._span;
            newTopo._offset = _topologyPtsPerRay[rayNpt.first][rayNpt.second].offset;
            //newTopo._label = _topologyPtsPerRay[rayNpt.first][rayNpt.second].label;
            // cout << "collected comp:" << endl;
            // MingUtility::printVector(newTopo._comps);
            // For debug.
            _cellTopo2AllTopoPtID.push_back(cellTopo2AllRayNPt[pr.first]);
            _cellTopo2topoPtID.push_back({rayNpt.first, rayNpt.second});
        }
        // for (const auto& pr : cellTopo2RayNPt) {
        //   const auto& rayNpt = pr.second;
        //   int index = _cellTopologies.size();
        //   _cellTopologies.push_back(pr.first);
        //   // auto& newTopo = _cellTopologies.back();
        //   auto& newTopo = _cellTopologies[index];
        //   newTopo._score = _topologyPtsPerRay[rayNpt.first][rayNpt.second]._score;
        //   newTopo._span = _topologyPtsPerRay[rayNpt.first][rayNpt.second]._span;
        //   newTopo._offset = _topologyPtsPerRay[rayNpt.first][rayNpt.second].offset;
        //   newTopo._label = _topologyPtsPerRay[rayNpt.first][rayNpt.second].label;
        //   // cout << "collected comp:" << endl;
        //   // MingUtility::printVector(newTopo._comps);
        //   // For debug.
        //   _cellTopo2AllTopoPtID.push_back(cellTopo2AllRayNPt[pr.first]);
        //   _cellTopo2topoPtID.push_back({rayNpt.first, rayNpt.second});
        // }
    }else{
        _cellNexplored2Topogroup.clear();
        int nex = 0;
        for (const auto& pr : cellTopo2AllRayNPt) {


            typedef pair<float, const vector<int>*> sitpair;
            vector<sitpair> scoreOrder;
            for (const auto &rayNpt:pr.second) {
                float score = _topologyPtsPerRay[rayNpt[0]][rayNpt[1]]._score;
                scoreOrder.push_back(make_pair(score, &rayNpt));
            }

            sort(scoreOrder.begin(), scoreOrder.end(),
                 [](const sitpair& lhs,
                 const sitpair& rhs) { return lhs.first > rhs.first; });



            for(const auto &onepair:scoreOrder){
                _cellNexplored2Topogroup.push_back(nex);
                const auto &rayNpt = *(onepair.second);
                int index = _cellTopologies.size();
                _cellTopologies.push_back(pr.first);
                // auto& newTopo = _cellTopologies.back();
                auto& newTopo = _cellTopologies[index];
                newTopo._score = _topologyPtsPerRay[rayNpt[0]][rayNpt[1]]._score;
                newTopo._span = _topologyPtsPerRay[rayNpt[0]][rayNpt[1]]._span;
                newTopo._offset = _topologyPtsPerRay[rayNpt[0]][rayNpt[1]].offset;
                //newTopo._label = _topologyPtsPerRay[rayNpt[0]][rayNpt[1]].label;
                // cout << "collected comp:" << endl;
                // MingUtility::printVector(newTopo._comps);
                // For debug.
                _cellTopo2AllTopoPtID.push_back(cellTopo2AllRayNPt[pr.first]);
                _cellTopo2topoPtID.push_back(rayNpt);
            }
            ++nex;

        }
    }
    cout << "nTopo = " << _cellTopologies.size() << endl;
    return _rays.size();
}
