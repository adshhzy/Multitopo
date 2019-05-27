#ifndef _VOLREADER_H_
#define _VOLREADER_H_
/* 
* Function to read/write volume data
*/
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "Utility.h"
using namespace std;

class Volume{
public:
	Volume():minValue(-1), maxValue(0),meanValue(0){}
	~Volume(){}
	void inline init(int x, int y, int z) {
		data.resize(x, vector<vector<float> >(y, vector<float>(z)));
		gridSize.resize(3);
		gridSize[0] = x;
		gridSize[1] = y;
		gridSize[2] = z;
	}
	void inline clear(){
		data.clear();
	}
	bool inline empty(){
		return data.empty();
	}
	void inline setDataAt(int i, int j, int k, float d){
		data[i][j][k] = d;
	}
	void inline setMinMaxValue(const float min, const float max){
		minValue = min;
		maxValue = max;
	}
	float inline getDataAt(int i, int j, int k){
		return data[i][j][k];
	}
	float inline getUnit(int i){
		return unitXYZ[i];
	}
	float inline getLowerCorner(int i){
		return lowCornerXYZ[i];
	}
	int inline getGridSize(int i){
		return gridSize[i];
	}
	void inline setLowerCornerXYZ(float x, float y, float z){
		lowCornerXYZ.resize(3);
		lowCornerXYZ[0] = x;
		lowCornerXYZ[1] = y;
		lowCornerXYZ[2] = z;
	}
	void inline setUnitXYZ(float x, float y, float z){
		unitXYZ.resize(3);
		unitXYZ[0] = x;
		unitXYZ[1] = y;
		unitXYZ[2] = z;
	}
	void normalize(float min, float max){
		float ratio = (max-min)/(maxValue-minValue);
		float trans = min + minValue*ratio;
		for(int i=0; i<gridSize[0]; ++i){
			for(int j=0; j<gridSize[1]; ++j){
				for(int k=0; k<gridSize[2]; ++k){
					data[i][j][k] = ratio * data[i][j][k] + trans;
				}
			}
		}
	}
	void writeVol(const char* fname){
		int x = data.size();
		int y = data[0].size();
		int z = data[0][0].size();
		ofstream ofs(fname, std::ofstream::out);
		ofs << "{";
		for (int i = 0; i < x; ++i){
			ofs << "{";
			for (int j = 0; j < y; ++j){
				ofs << "{";
				for (int k = 0; k < z; ++k){
					ofs << data[i][j][k];
					if (k < z - 1){
						ofs << ",";
					}
				}
				if (j < y - 1){
					ofs << "},";
				}
				else{
					ofs << "}";
				}
			}
			if (i < x - 1){
				ofs << "},";
			}
			else{
				ofs << "}";
			}
		}
		ofs << "}";
		ofs.close();
	}



	// Function to write out to a MRC file
	void toMRCFile( const char* fname )
	{
		FILE* fout = fopen( fname, "wb" ) ;

		// Write header
		fwrite( &gridSize[0], sizeof( int ), 1, fout ) ;
		fwrite( &gridSize[1], sizeof( int ), 1, fout ) ;
		fwrite( &gridSize[2], sizeof( int ), 1, fout ) ;

		int mode = 2 ;
		fwrite( &mode, sizeof ( int ), 1, fout ) ;

		int off[3] = {0,0,0} ;
		int intv[3] = { gridSize[0] - 1, gridSize[1] - 1, gridSize[2] - 1 } ;
		fwrite( off, sizeof( int ), 3, fout ) ;
		fwrite( intv, sizeof( int ), 3, fout ) ;

		float cella[3] = {intv[0]*unitXYZ[0],intv[1]*unitXYZ[1],intv[2]*unitXYZ[2]} ;
		float cellb[3] = {90,90,90} ;
		fwrite( cella, sizeof( float ), 3, fout ) ;
		fwrite( cellb, sizeof( float ), 3, fout ) ;

		int cols[3] = {1,2,3} ;
		fwrite( cols, sizeof( int ), 3, fout ) ;

		float ds[3] = {minValue, maxValue, meanValue} ;
		fwrite( ds, sizeof( float ), 3, fout ) ;

		int zero = 0 ;
		for ( int i = 22 ; i < 256 ; i ++ )
		{
			fwrite( &zero, sizeof( int ), 1, fout ) ;
		}

		// Write contents
		for ( int z = 0 ; z < gridSize[2] ; z ++ ){
			for ( int y = 0 ; y < gridSize[1] ; y ++ ){
				for ( int x = 0 ; x < gridSize[0] ; x ++ )
				{
					fwrite( &data[x][y][z], sizeof( float ), 1, fout ) ;
				}
			}
		}

		fclose( fout ) ;
	}

private:
	vector<vector<vector<float> > > data;
	vector<float> lowCornerXYZ;
	vector<float> unitXYZ;
	vector<int> gridSize;
	float minValue;
	float maxValue;
	float meanValue;
};

class MRCReader 
{
public:
	/* Initializer */
	MRCReader():dimx(0), dimy(0), dimz(0){}
	MRCReader( const char* fname )
	{
		sprintf( mrcfile, "%s", fname ) ;

		FILE* fin = fopen( fname, "rb" ) ;
		if (fin == NULL){
			return;
		}

		// Parse header
		fread( &totx, sizeof( int ), 1, fin ) ;
		fread( &toty, sizeof( int ), 1, fin ) ;
		fread( &totz, sizeof( int ), 1, fin ) ;

		fread( &mode, sizeof( int ), 1, fin ) ;

		fread( &offx, sizeof( int ), 1, fin ) ;
		fread( &offy, sizeof( int ), 1, fin ) ;
		fread( &offz, sizeof( int ), 1, fin ) ;

		fread( &dimx, sizeof( int ), 1, fin ) ;
		fread( &dimy, sizeof( int ), 1, fin ) ;
		fread( &dimz, sizeof( int ), 1, fin ) ;
		dimx ++ ;
		dimy ++ ;
		dimz ++ ;

		fread( &angsx, sizeof( float ), 1, fin ) ;
		fread( &angsy, sizeof( float ), 1, fin ) ;
		fread( &angsz, sizeof( float ), 1, fin ) ;

		fread( &anglex, sizeof( float ), 1, fin ) ;
		fread( &angley, sizeof( float ), 1, fin ) ;
		fread( &anglez, sizeof( float ), 1, fin ) ;

		fseek( fin, 4 * 3, SEEK_CUR ) ;

		fread( &dmin, sizeof( float ), 1, fin ) ;
		fread( &dmax, sizeof( float ), 1, fin ) ;
		fread( &dmean, sizeof( float ), 1, fin ) ;

		fseek( fin, 4 * 27, SEEK_CUR ) ;

		fread( &orgx, sizeof( float ), 1, fin ) ;
		fread( &orgy, sizeof( float ), 1, fin ) ;
		fread( &orgz, sizeof( float ), 1, fin ) ;

		fseek( fin, 4 * 2, SEEK_CUR ) ;

		fread( &drms, sizeof( float ), 1, fin ) ;
		fclose( fin ) ;

		dimx = totx ;
		dimy = toty ;
		dimz = totz ;

		if ( mode > 2 )
		{
			printf("Complex mode not supported.\n") ;
			exit(0) ;
		}
	}

	/* Read volume */
	void getVolume(Volume & vol )
	{
		FILE* fin = fopen( mrcfile, "rb" ) ;
		if (fin == NULL){
			return;
		}
		fseek( fin, 1024, SEEK_SET ) ;

		char chard ;
		short shortd ;
		float floatd ;
		float d ;

		float minValue = INFINITY;
		float maxValue = -INFINITY;
		vol.init( dimx, dimy, dimz ) ;
		for ( int i = 0 ; i < dimz ; i ++ ){
			for ( int j = 0 ; j < dimy ; j ++ ){
				for ( int k = 0 ; k < dimx ; k ++ )
				{
					switch ( mode )
					{
					case 0: 
						fread( &chard, sizeof( char ), 1, fin ) ;
						d = (float) chard ;
						break ;
					case 1:
						fread( &shortd, sizeof( short ), 1, fin ) ;
						d = (float) shortd ;
						break ;
					case 2:
						fread( &floatd, sizeof( float ), 1, fin ) ;
						d = (float) floatd ;
						break ;
					}
					minValue = d < minValue ? d : minValue;
					maxValue = d > maxValue ? d : maxValue;
					vol.setDataAt( k, j, i, d ) ;
				}
			}
		}
		vol.setMinMaxValue(minValue,maxValue);
		fclose( fin ) ;
	}

	/* Get resolution */
	void getSpacing( float& ax, float& ay, float& az )
	{
		ax = angsx / (dimx - 1);
		ay = angsy / (dimy - 1) ;
		az = angsz / (dimz - 1) ;
	}


private:

	int totx, toty, totz ;
	int offx, offy, offz ;
	int dimx, dimy, dimz ;

	float angsx, angsy, angsz ;
	float anglex, angley, anglez ;
	float dmin, dmax, dmean, drms ;
	float orgx, orgy, orgz;

	int mode ;

	char mrcfile[1024] ;
};


class VolReader 
{
public:
	/* Initializer */
	VolReader():dimx(0), dimy(0), dimz(0){}
	VolReader( const char* hdrfname, const char* datfname )
	{
		sprintf( datfile, "%s", datfname ) ;

		ifstream fin ( hdrfname );
		if(!fin.good()){
			cout << "Can not read HRD file " << hdrfname << endl;
			return;
		}

		lowCorner.resize(3);
		unitXYZ.resize(3);
		vecX.resize(3);
		vecY.resize(3);
		vecZ.resize(3);

		fin >> lowCorner[0] >> lowCorner[1] >> lowCorner[2];
		fin >> vecX[0] >> vecX[1] >> vecX[2];
		fin >> vecY[0] >> vecY[1] >> vecY[2];
		fin >> vecZ[0] >> vecZ[1] >> vecZ[2];
		fin >> dimx >> dimy >> dimz;
		fin >> unitXYZ[0] >> unitXYZ[1] >> unitXYZ[2];

		fin.close();

		int n = 0;
		while (dimx > pow(2.0,n)) {
			n++;
		}
		dimAllX = pow(2.0,n);
		n = 0;
		while (dimy > pow(2.0,n)) {
			n++;
		}
		dimAllY = pow(2.0,n);
		n = 0;
		while (dimz > pow(2.0,n)) {
			n++;
		}
		dimAllZ = pow(2.0,n);

#if _CG_LOG
		cout << "lowCorner: " << lowCorner[0] << ", "<< lowCorner[1] << ", "<< lowCorner[2] <<endl;
		cout << "unitXYZ: " << unitXYZ[0] << ", "<< unitXYZ[1] << ", "<< unitXYZ[2] <<endl;
		cout << "dimXYZ: " << dimx << ", "<< dimy << ", "<< dimz <<endl;
		cout << "dimAllXYZ: " << dimAllX << ", "<< dimAllY << ", "<< dimAllZ <<endl;
#endif
	}

	/* Read volume */
	void getVolume(Volume & vol )
	{
		ifstream fin(datfile, ios::in | ios::binary);
		if(!fin.good()){
			cout << "Can not read DAT file " << datfile << endl;
			return;
		}
		short int * data_in = new short int[ dimAllX * dimAllY * dimAllZ ];
		fin.read((char*) data_in, dimAllX * dimAllY * dimAllZ * sizeof(short int));
		fin.close();
		float minValue = INFINITY;
		float maxValue = -INFINITY;
		float curValue;
		vol.init( dimAllX, dimAllY, dimAllZ ) ;
		int index = 0;

		for ( int i = 0 ; i < dimAllZ ; i ++ ){
			for ( int j = dimAllY-1 ; j >=0 ; j -- ){
				for ( int k = 0 ; k <dimAllX ; k++, index++ ){		
					float curValue = (float) data_in[index];
					minValue = curValue < minValue ? curValue : minValue;
					maxValue = curValue > maxValue ? curValue : maxValue;
					vol.setDataAt( k, j, i, curValue ) ;
				}
			}
		}

		delete [] data_in;
		vol.setMinMaxValue(minValue,maxValue);
	}
	float inline getLowCorner(int i){
		return lowCorner[i];
	}
	float inline getUnit(int i){
		return unitXYZ[i];
	}
private:
	vector<float> lowCorner;
	vector<float> vecX;
	vector<float> vecY;
	vector<float> vecZ;
	int dimx, dimy, dimz ;
	int dimAllX, dimAllY, dimAllZ;
	vector<float> unitXYZ;

	char datfile[1024] ;
};

#endif