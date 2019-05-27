#ifndef MYMESH_H
#define MYMESH_H


#include "utility.h"
#include<math.h>
#include<iostream>
#include<vector>
using namespace std;
using namespace MyUtility;

class CellComplex{

public:
    uint n_vertices, n_edges, n_faces, n_voxels, cell_dim;
    vector< uint >dimtable;
    vector< uint >sizetable;
    vector< double >vertices;
    vector< uint >edges;//edges2vertices
    vector< uint >faces;//faces2edges
    vector< uint >voxels;
    vector< uint >faces2vertices;
    vector< vector< uint >* >cc;
    CellComplex(uint vertices_dim = 3, uint edges_dim = 2, uint faces_dim = 3, uint voxels_dim = 0){
		n_vertices = n_edges = n_faces = n_voxels = 0;
		dimtable.push_back(vertices_dim);
		dimtable.push_back(edges_dim);
		dimtable.push_back(faces_dim);
		dimtable.push_back(voxels_dim);
		cc.resize(4, NULL);
		if (voxels_dim == 0){
			cell_dim = 3;
		}
		else{
			cell_dim = 4;
		}
		cc[1] = &edges; cc[2] = &faces; cc[3] = &voxels;
	}
    virtual void reset(){
		vertices.clear(); edges.clear(); faces.clear(); voxels.clear(); faces2vertices.clear();
		n_vertices = n_edges = n_faces = n_voxels = 0;
		sizetable.clear();
	}
    virtual void setparameters(){
		sizetable.resize(4, 0);
		if (dimtable[3] != 0){ n_voxels = voxels.size() / dimtable[3]; sizetable[3] = n_voxels; }
		n_vertices = vertices.size() / dimtable[0]; n_edges = edges.size() / dimtable[1];
		n_faces = faces.size() / dimtable[2];
		sizetable[0] = n_vertices;
		for (int i = 1; i < 3; i++)sizetable[i] = cc[i]->size() / dimtable[i];
		if (n_faces == 0 && faces2vertices.size() != 0){ n_faces = sizetable[2] = faces2vertices.size() / dimtable[2]; }
        cc[1] = &edges; cc[2] = &faces; cc[3] = &voxels;
	}



    void inline vectorize(const uint first, const uint second,double *e){
        auto vdim = dimtable[0];
        auto p_v1 = &(vertices[first*vdim]);
        auto p_v2 = &(vertices[second*vdim]);
        minusVec(p_v1,p_v2,e,vdim);
    }
    void inline vectorize(const double *v1, const uint second,double *e){
        auto vdim = dimtable[0];
        auto p_v2 = &(vertices[second*vdim]);
        minusVec(v1,p_v2,e,vdim);
    }
    double inline VerticesDistance(const uint first, const uint second)const{
        double dist = 0.0f;
		auto vdim = dimtable[0];
        auto p_v1 = &(vertices[first*vdim]);
        auto p_v2 = &(vertices[second*vdim]);
        for (uint i = 0; i < vdim; ++i)
            dist += pow((p_v1[i] - p_v2[i] ), 2);
		return sqrt(dist);
	}

    double inline VerticesDistance(const double *p_v1, const uint second)const{
        double dist = 0.0f;
        auto vdim = dimtable[0];
        auto p_v2 = &(vertices[second*vdim]);
        for(int i = 0;i<vdim;++i)dist += powf((p_v1[i] - p_v2[i] ), 2);
        return sqrt(dist);
    }
    void inline VerticesMidpoint(const uint first,const  uint second, double *p_vc)const{
        auto vdim = dimtable[0];
        auto p_v1 = &(vertices[first*vdim]);
        auto p_v2 = &(vertices[second*vdim]);
        for(int i = 0;i<vdim;++i)p_vc[i] = (p_v1[i] + p_v2[i]) /2;
    }

    double inline TriangleArea(const uint first,const  uint second,const  uint third)const{
		auto v1 = &(vertices[3 * first]);
		auto v2 = &(vertices[3 * second]);
		auto v3 = &(vertices[3 * third]);
        return _TriangleArea(v1,v2,v3);
	}

    double inline TriangleLeastAngle(const uint first,const uint second,const uint third){
        const auto vdim = dimtable[0];
        const auto v1 = &(vertices[vdim * first]);
        const auto v2 = &(vertices[vdim * second]);
        const auto v3 = &(vertices[vdim * third]);
        return _TriangleLeastAngle(v1,v2,v3,vdim);
    }


};


class Mesh : public CellComplex{
public:

    bool ismanifold = true;

    int num_max_neighor;
    int step;
	int f2f_step;
	bool isbuildneighbortable;
    bool isbuildEdges;
    vector< uint >faces2faces;
    vector< uint >vertices2faces;
    vector< unsigned long >vertices2faces_accumulate;
    vector< uint >vertices2faces_inverse;


    vector< uint >edges2faces;
    vector< double >faces_normal;
    vector< double >vertices_normal;
    vector<bool>faces_flinge;
    vector<bool>vertices_flinge;
    vector<double>vertices_defect;
    vector<uchar>faces2edgesPN;

    vector< uint >vertices2edges;
    vector< uint >vertices2vertices;
    vector< unsigned long >vertices2edges_accumulate;
    vector<uchar>vertices2edgesPN;
    vector< uint >vertices2vertices_inverse;

public:

    vector<bool>edge_non_manifold;
    vector<uint>efnum_nonmanifold;

    vector<bool>vertices_non_manifold;
    vector<bool>faces_in_manifold;

    vector< uint >vertices2vertices_nonmanifold;
    vector< unsigned long >vertices2edges_accumulate_nonmanifold;

    vector< uint >faces2faces_nonmanifold;
    vector<  unsigned long >faces2faces_accumulate_nonmanifold;

    vector< uint >edges2faces_nonmanifold;
    vector<  unsigned long >edges2faces_accumulate_nonmanifold;

public:
    double lscale, pcenter[3];

public:
    Mesh(uint vertices_dim = 3, uint edges_dim = 2, uint faces_dim = 3) : CellComplex(vertices_dim, edges_dim, faces_dim){
		f2f_step = faces_dim + 1;
		isbuildneighbortable = false;
        num_max_neighor = 16;
        step = num_max_neighor + 1;
	}

	Mesh(CellComplex& cc_in) : CellComplex(cc_in.dimtable[0], cc_in.dimtable[1], cc_in.dimtable[2]){
		vertices = cc_in.vertices;
		faces2vertices = cc_in.faces2vertices;
		setparameters();
	}
    void inline setparameters(){

        sizetable.resize(4, 0);
        if (dimtable[3] != 0){ n_voxels = voxels.size() / dimtable[3]; sizetable[3] = n_voxels; }
        n_vertices = vertices.size() / dimtable[0]; n_edges = edges.size() / dimtable[1];
        sizetable[0] = n_vertices;
        sizetable[1] = n_edges;
        n_faces = sizetable[2] = faces2vertices.size() / dimtable[2];
        cc[1] = &edges; cc[2] = &faces; cc[3] = &voxels;

    }
    void inline setCellComplex(CellComplex& cc_in){
		reset();
		dimtable = cc_in.dimtable;
		vertices = cc_in.vertices;
		faces2vertices = cc_in.faces2vertices;
		setparameters();
	}
    void inline reset(){
         vertices.clear(); edges.clear(); faces.clear(); voxels.clear(); faces2vertices.clear();
         faces2faces.clear();vertices2faces.clear();faces_normal.clear();vertices_normal.clear();
         faces_in_manifold.clear();faces_flinge.clear();
         faces2edgesPN.clear();
         vertices2edgesPN.clear();
         vertices2faces_accumulate.clear();edges2faces.clear();
         vertices2edges_accumulate.clear();vertices2edges.clear();
         vertices2vertices.clear();vertices2vertices_inverse.clear();
         isbuildneighbortable = false;isbuildEdges=false;
         n_vertices = n_edges = n_faces = n_voxels = 0;
         sizetable.clear();

         faces2faces_nonmanifold.clear();
         faces2faces_accumulate_nonmanifold.clear();

         edges2faces_nonmanifold.clear();
         edges2faces_accumulate_nonmanifold.clear();

         vertices2vertices_nonmanifold.clear();
         vertices2edges_accumulate_nonmanifold.clear();
         edge_non_manifold.clear();
         efnum_nonmanifold.clear();

         vertices_non_manifold.clear();
         faces_in_manifold.clear();

    }
public:
    void BuildUp(bool ismanifold_in);
    void CheckValidation();

	void BuildNeighborTable();
    void BuildNeighborTable_nonmanifold();
    void ComputeDefectAngle();
	void DeleteInmanifoldedFaces();
    double ComputeDihedralAngle(uint e_ind);
    void GetRescaleInfo(double &lscale, double *pcenter);
    void ReScale_uniform(double lscale = 1.0, bool isMoveCent = true);
    void ReScale(double xscale,double yscale,double zscale);
    void ReScale(double lscale,double *pcenter);
    void ReScale_uniform_directratio(double ratio);
    void MoveBottom(bool xx,bool yy,bool zz);
    //void MultiplyMatrix(double *m_transform);
	void ReOrientFaces();
    void InverseFaces();
    double ComputeTotalArea();
    void ComputeEdgefromFace(bool isRecomputeNeighborhood = false);
    void ComputeEdgefromFaceWithDuplication();
    void ComputeFaceNormal(bool isComputeVerticesNormal = true,vector<bool>*p_isinver =NULL);






    bool RayHitTriangle(double* raystart, double* rayend, uint f_ind,double *hitpoint = NULL);
    bool RayHitTarget(double* raystart,double* rayend);
	void addCellComplex(CellComplex& cc_in);
    void addMesh(Mesh& mesh_in);
    void TriangleCircumcenter(uint f_ind,double *cc,double &r);
    void TriangleMidpoint(uint f_ind,double *cc);
    double point2TriSquareDistance(const double *P,uint find,double *cp);
    double EdgeCottan(uint e_ind);
    void Fairing(bool usingedges = true, bool iskeepboundary = true, bool isRecomputeNeighborhood = false, int round = 50,
                 double lamda = 0.5, double kaipa = 0.1, vector<bool> *staticv = NULL);

    void SurfaceFFFairing(double lambda,int iter,vector<bool>* inverNor,vector<bool>* staticv);
    void SurfaceLaplacianFairing(double lambda,int iter,vector<bool>* staticv);

    bool readfile(string filename);
    bool createToy(int toyind);

    bool savePLY2File(string filename);
private:
    void SpcecialNormalComputation(vector<bool>*p_isinver,vector<double>&out_vnor,vector<double>&out_mag);


public:
    bool readObjfile(string filename);
    bool readOfffile(string filename);
    bool saveObjFile(string filename);


public:


//	inline uint* vf_begin(uint v_ind){ return &(vertices2faces[v_ind * step + 1]); }
//	inline uint* vf_end(uint v_ind){ return &(vertices2faces[v_ind * step + 1]) + vertices2faces[v_ind * step]; }
    inline uint* vf_begin(uint v_ind){ return &(vertices2faces[vertices2faces_accumulate[v_ind]]); }
    inline uint* vf_end(uint v_ind){ return &(vertices2faces[vertices2faces_accumulate[v_ind+1]]); }
    inline uint vf_num(uint v_ind){ return vertices2faces_accumulate[v_ind+1]-vertices2faces_accumulate[v_ind]; }
    inline uint* vfinv_begin(uint v_ind){ return &(vertices2faces_inverse[vertices2faces_accumulate[v_ind]]); }
    inline uint* vfinv_end(uint v_ind){ return &(vertices2faces_inverse[vertices2faces_accumulate[v_ind+1]]); }


    inline uint* ev_begin(uint e_ind){ return &(edges[e_ind * 2]); }
    inline uint* ev_end(uint e_ind){ return &(edges[e_ind * 2 + 2]); }
    inline uint* ve_begin(uint v_ind){ return &(vertices2edges[vertices2edges_accumulate[v_ind]]); }
    inline uint* ve_end(uint v_ind){ return &(vertices2edges[vertices2edges_accumulate[v_ind+1]]); }
    inline uint ve_num(uint v_ind){ return vertices2edges_accumulate[v_ind+1]-vertices2edges_accumulate[v_ind]; }
    inline uchar* vepn_begin(uint v_ind){ return &(vertices2edgesPN[vertices2edges_accumulate[v_ind]]); }
    inline uchar* vepn_end(uint v_ind){ return &(vertices2edgesPN[vertices2edges_accumulate[v_ind+1]]); }
    inline uint* vv_begin(uint v_ind){ return &(vertices2vertices[vertices2edges_accumulate[v_ind]]); }
    inline uint* vv_end(uint v_ind){ return &(vertices2vertices[vertices2edges_accumulate[v_ind+1]]); }
    inline uint* vvinv_begin(uint v_ind){ return &(vertices2vertices_inverse[vertices2edges_accumulate[v_ind]]); }
    inline uint* vvinv_end(uint v_ind){ return &(vertices2vertices_inverse[vertices2edges_accumulate[v_ind+1]]); }
    inline uint vv_num(uint v_ind){ return vertices2edges_accumulate[v_ind+1]-vertices2edges_accumulate[v_ind]; }

    inline uint* fv_begin(uint f_ind){ return &(faces2vertices[f_ind * 3]); }
    inline uint* fv_end(uint f_ind){ return &(faces2vertices[f_ind * 3 + 3]); }
    inline uint* ff_begin(uint f_ind){ return &(faces2faces[f_ind * f2f_step + 1]); }
    inline uint* ff_end(uint f_ind){ return &(faces2faces[f_ind * f2f_step + 1]) + faces2faces[f_ind * f2f_step]; }
    inline uint ff_num(uint f_ind){ return faces2faces[f_ind * f2f_step]; }
    inline uint* ff_endspace(uint f_ind){ return &(faces2faces[f_ind * f2f_step + f2f_step]); }
    inline uint* ef_begin(uint e_ind){ return &(edges2faces[e_ind * 2]); }
    inline uint* ef_end(uint e_ind){ if(edges2faces[e_ind * 2 + 1]!=UINT_MAX)return &(edges2faces[e_ind * 2 + 2]);
        else return &(edges2faces[e_ind * 2 + 1]);}
    inline uint* ef_endspace(uint e_ind){ return &(edges2faces[e_ind * 2 + 2]);}
    inline uint ef_num(uint e_ind){ if(edges2faces[e_ind * 2 + 1]!=UINT_MAX)return 2;else return 1; }
    inline bool isflingeedge(uint e_ind){ return(edges2faces[e_ind * 2 + 1]==UINT_MAX); }
    inline uint* fe_begin(uint f_ind){ return &(faces[f_ind * 3]); }
    inline uint* fe_end(uint f_ind){ return &(faces[f_ind * 3 + 3]); }
    inline uchar* fepn_begin(uint f_ind){ return &(faces2edgesPN[f_ind * 3]); }
    inline uchar* fepn_end(uint f_ind){ return &(faces2edgesPN[f_ind * 3 + 3]); }

    inline double* get_face_normal(uint f_ind){ return &(faces_normal[f_ind * 3]); }
    inline double* get_vertice(uint v_ind){ return &(vertices[v_ind * dimtable[0]]); }
    inline double* get_vertice_normal(uint v_ind){ return &(vertices_normal[v_ind * 3]); }
    inline double* v_begin(uint v_ind){ return &(vertices[v_ind * dimtable[0]]); }
    inline double* v_end(uint v_ind){ return &(vertices[(v_ind+1) * dimtable[0]]); }

    inline double* fnor_begin(uint f_ind){ return &(faces_normal[f_ind * 3]); }
    inline double* fnor_end(uint f_ind){ return &(faces_normal[f_ind * 3 + 3]); }
    inline double* vnor_begin(uint v_ind){ return &(vertices_normal[v_ind * 3]); }
    inline double* vnor_end(uint v_ind){ return &(vertices_normal[v_ind * 3 + 3]); }


    inline double v_defect(uint v_ind){ return (vertices_defect[v_ind]); }
    inline bool isflingevertives(uint v_ind){ return (vertices_flinge[v_ind]); }


public:
    inline uint* vvnon_begin(uint v_ind){ return &(vertices2vertices_nonmanifold[vertices2edges_accumulate_nonmanifold[v_ind]]); }
    inline uint* vvnon_end(uint v_ind){ return &(vertices2vertices_nonmanifold[vertices2edges_accumulate_nonmanifold[v_ind+1]]); }
    inline uint vvnon_num(uint v_ind){ return vertices2edges_accumulate_nonmanifold[v_ind+1]-vertices2edges_accumulate_nonmanifold[v_ind]; }


    inline uint* ffnon_begin(uint f_ind){ return &(faces2faces_nonmanifold[faces2faces_accumulate_nonmanifold[f_ind]]); }
    inline uint* ffnon_end(uint f_ind){ return &(faces2faces_nonmanifold[faces2faces_accumulate_nonmanifold[f_ind+1]]); }
    inline uint ffnon_num(uint f_ind){ return faces2faces_accumulate_nonmanifold[f_ind+1]-faces2faces_accumulate_nonmanifold[f_ind]; }

    inline uint* efnon_begin(uint e_ind){ return &(edges2faces_nonmanifold[edges2faces_accumulate_nonmanifold[e_ind]]);}
    inline uint* efnon_end(uint e_ind){ return &(edges2faces_nonmanifold[edges2faces_accumulate_nonmanifold[e_ind+1]]);}
    inline uint efnon_num(uint e_ind){ return edges2faces_accumulate_nonmanifold[e_ind+1]-edges2faces_accumulate_nonmanifold[e_ind]; }
};

#endif
