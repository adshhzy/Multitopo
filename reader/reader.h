#ifndef _reader_h
#define _reader_h

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

template <class T>
bool readVVecFile(string filename, vector<vector<T>> &vvec){
	ifstream reader(filename.data(), ifstream::in);
	if (!reader.good()) {
		cout << "Can not open the file " << filename << endl;
		return false;
	}
	int nnum = 0;
	reader>>nnum;
	vector<int>vvnum(nnum,0);
	for(int i=0;i<nnum;++i)reader>>vvnum[i];
	
	vvec.resize(nnum);
	for(int i=0;i<nnum;++i){
		vvec[i].resize(vvnum[i]);
		for(int j=0;j<vvnum[i];++j)reader>>vvec[i][j];
	}
	reader.close();
	return true;
}

template <class T>
bool readVecFile(string filename,vector<T>&vec){
	ifstream reader(filename.data(), ifstream::in);
	if (!reader.good()) {
		cout << "Can not open the file " << filename << endl;
		return false;
	}
	int nnum = 0;
	reader>>nnum;
	vec.resize(nnum);
	for(int i=0;i<nnum;++i)reader>>vec[i];
	
	reader.close();
	return true;
	
}

bool readObjFile_Line(string filename,vector<double>&vertices,vector<unsigned int>&edges2vertices){
	
	
	ifstream fin(filename.data());
	if(fin.fail()){
		cout<<"Fail to open input file: "<<filename<<endl;
		return false;
	}
	
	vertices.clear();
	edges2vertices.clear();
	auto readVertices = [&vertices](stringstream &strs){
		double dvalue;
		for(int i=0;i<3;++i){strs>>dvalue;vertices.push_back(dvalue);}
	};
	auto readEdges = [&edges2vertices](stringstream &strs){
		string oneset,indstring;
		unsigned int ivalue,ivalue2;
		strs>>ivalue;
		while (strs>>ivalue2){
			edges2vertices.push_back(ivalue-1);
			edges2vertices.push_back(ivalue2-1);
			ivalue = ivalue2;
		}
		
	};
	
	
	string oneline;
	
	cout<<"reading: "<<filename<<endl;
	
	while( getline( fin, oneline ) ){
		stringstream strs( oneline );
		string prefix;
		
		strs >> prefix;
		
		if( prefix == "v"  ) { readVertices( strs ); continue; } // vertex
		if( prefix == "l"  ) { readEdges( strs ); continue; } // edges
		if( prefix == "vt" ) {  continue; } // texture coordinate
		if( prefix == "vn" ) {  continue; } // vertex normal
		if( prefix == "vf" ) { /*readVerticesField( ss );*/ continue; } // tangent vector
		if( prefix == "f"  ) {  continue; } // face
		if( prefix[0] == '#' ) continue; // comment
		if( prefix == "o" ) continue; // object name
		if( prefix == "g" ) continue; // group name
		if( prefix == "s" ) continue; // smoothing group
		if( prefix == "mtllib" ) continue; // material library
		if( prefix == "usemtl" ) continue; // material
		if( prefix == "k" ) continue; // field degree
		if( prefix == "fs" ) continue; // field singularity
		if( prefix == "" ) continue; // empty string
		if( prefix == "c" ) continue;
		
		cout << "Error: not a valid curf file!" << endl;
		cout << "(Offending line: " << oneline << ")" << endl;
		return false;
	}
	
	
	fin.close();
	return true;
	
	
	
	
}









bool readSufFile(string filename,  vector<double>&vertices,  vector<unsigned int>&faces2vertices,  vector<int>&facesMat,  vector<int> &CtrEdges){
	
	ifstream fin(filename.data(), ofstream::out);
	if (!fin.good()) {
		cout << "Can not open suf file " << filename << endl;
		return false;
	}
	int n_vertices;
	int n_faces;
	int n_ctredges;
	
	fin>>n_vertices;
	fin>>n_faces;
	
	vertices.clear();
	vertices.resize(n_vertices*3);
	for(int i=0;i<n_vertices;++i){
		auto p_v = vertices.data()+i*3;
		for(int j=0;j<3;++j)fin>>p_v[j];
	}
	
	faces2vertices.clear();
	faces2vertices.resize(n_faces*3);
	facesMat.clear();
	facesMat.resize(n_faces*2);
	for(int i=0;i<n_faces;++i){
		auto p_fv = faces2vertices.data()+i*3;
		auto p_fm = facesMat.data()+i*2;
		for(int j=0;j<3;++j)fin>>p_fv[j];
		for(int j=0;j<2;++j)fin>>p_fm[j];

	}
	
	fin>>n_ctredges;
	CtrEdges.resize(n_ctredges*2);
	for(int i=0;i<n_ctredges;++i){
		auto p_ctr = CtrEdges.data()+i*2;
		for(int j=0;j<2;++j)fin>>p_ctr[j];
	}
	fin.close();
	return true;
	
}












template <class T>
bool writeVVecFile(string filename, const vector< vector<T> > &vvec){
	filename = filename + "_vvec.txt";
	ofstream fout(filename.data());
	if(fout.fail()){
		cout<<"Fail to create output file: "<<filename<<endl;
		return false;
	}
	
	int numof = vvec.size();
	
	fout<<numof<<endl;
	
	if(numof==0)return true;
	for(int i=0;i<numof-1;i++){
		fout<<vvec[i].size()<<' ';
	}
	fout<<vvec[numof-1].size()<<endl;
	
	for(auto &a:vvec){
		if(a.size()==0)continue;
		for(int i=0;i<a.size()-1;i++){
			fout<<a[i]<<' ';
		}
		fout<<a[a.size()-1]<<endl;
	}
	fout.close();
	
	cout<<"Write: "<<filename<<endl;
	return true;
}

template <class T>
bool writeVecFile(string filename, const vector<T>&vec){
	filename = filename + "_vec.txt";
	ofstream fout(filename.data());
	if(fout.fail()){
		cout<<"Fail to create output file: "<<filename<<endl;
		return false;
	}
	
	int numof = vec.size();
	fout<<numof<<endl;
	for(auto a:vec)fout<<a<<endl;
	fout.close();
	
	cout<<"Write: "<<filename<<endl;
	return true;
}

bool writeSufFile(string filename, const vector<double>&vertices, const vector<unsigned int>&faces2vertices, const vector<int>&facesMat, const vector<int> &CtrEdges){
	
	filename = filename + ".suf";
	ofstream outer(filename.data(), ofstream::out);
	if (!outer.good()) {
		cout << "Can not create output suf file " << filename << endl;
		return false;
	}
	int n_vertices = vertices.size()/3;
	int n_faces = faces2vertices.size()/3;
	
	
	outer<<n_vertices<<' '<<n_faces<<endl;
	for(int i=0;i<n_vertices;++i){
		auto p_v = vertices.data()+i*3;
		outer << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
	}
	for(int i=0;i<n_faces;++i){
		auto p_fv = faces2vertices.data()+i*3;
		auto p_fm = facesMat.data()+i*2;
		outer<< p_fv[0]<< " "<< p_fv[1] << " "<< p_fv[2] << " "<< p_fm[0]<< " "<<p_fm[1] <<endl;
	}
	
	outer<<CtrEdges.size()/2<<endl;
	for(int i=0;i<CtrEdges.size()/2;++i){
		outer<<CtrEdges[i*2]<<' '<<CtrEdges[i*2+1]<<endl;
	}
	outer.close();
	return true;
	
}



bool writeCtrGraphFile(string filename,const vector<float>&vertices,const vector<vector<int>>&edge2vertices,const vector<vector<int>>&edgeMat, const vector<vector<float>>&planepara){
	filename = filename + ".CtrGraph";
	ofstream outer(filename.data(), ofstream::out);
	if (!outer.good()) {
		cout << "Can not create output CtrGraph file " << filename << endl;
		return false;
	}
	
	int n_vertices = vertices.size()/3;
	int n_plane = edge2vertices.size();
	outer<<"n "<<n_vertices<<endl;
	auto p_vd = vertices.data();
	for(int i=0;i<n_vertices;++i){
		auto p_v = p_vd+ i*3;
		outer << "v "<<p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
	}
	
	outer<<"n "<<n_plane<<endl;
	for(int i=0;i<n_plane;++i){
		outer<<"p ";
		for(int j=0;j<4;++j)outer<<planepara[i][j]<<' ';
		
		auto &p_ev = edge2vertices[i];
		auto &p_em = edgeMat[i];
		int ne = p_ev.size()/2;
		outer<<ne<<endl;
		
		for(int j=0;j<ne;++j){
			auto ind = j*2;
			outer<<"e "<<p_ev[ind]<<' '<<p_ev[ind+1]<<' '<<p_em[ind]<<' '<<p_em[ind+1]<<endl;
		}
	}
	
	outer.close();
	
	return true;
	
	
	
	
	
}
void writeObjfile_line(string filename,const vector<double>&vertices,const vector<int>&edge2vertices){
	
	filename = filename + ".obj";
	ofstream outer(filename.data(), ofstream::out);
	if (!outer.good()) {
		cout << "Can not create output Obj file " << filename << endl;
		//return false;
	}
	
	
	
	int n_vertices = vertices.size()/3;
	
	for(int i=0;i<n_vertices;++i){
		auto p_v = vertices.data()+i*3;
		outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
	}
	
	
	{
		int n_edges = edge2vertices.size()/2;
		for(int i=0;i<n_edges;++i){
			auto p_ev = edge2vertices.data()+i*2;
			outer << "l " << p_ev[0]+1<< " "<< p_ev[1]+1 << endl;
		}
	}
	
	
	outer.close();
	cout<<"saving finish: "<<filename<<endl;
	//return true;
}


#endif
