#include <iostream>
#include <unistd.h>
#include <vector>
#include <string>
#include<dirent.h>
#include <sys/stat.h>
using namespace std;


bool isCJP = true;

int az2(string incontour, string outdir, float _teta, float _randwb, int _nSmBefLoop, int _nLoop, int _nSmInLoop, int _genus, int MMprop_f2vmethod,
        int endstep, vector<int>&numofTopo, vector<int>&cellSeq, int &nMat, int computeSaveLoad, string iofile);
vector<int> runDynamicProgramming(string protocalname, vector<int>&TopoConstraintInd);
void outputCombineSurface(string outfolder, vector<int>&Cell2NTopo, vector<int>&pickTopoInd);


void makemethodfolder(string &foldername){
    struct stat sb;
    if(!(stat(foldername.data(), &sb) == 0 && S_ISDIR(sb.st_mode)))mkdir(foldername.data(),0777);
    return;
}

int main(int argc, char** argv)
{
    //cout << argc << endl;

    int  MMprop_f2vmethod = 1;
    //string ctrmodelname = "outTestmapping";       MMprop_f2vmethod = 1;
    //string ctrmodelname = "ringc";                   MMprop_f2vmethod = 0;
    string ctrmodelname = "mug";                     MMprop_f2vmethod = 0;

    string mediumfolder("/Users/Research/Geometry/MM/MultiTopo/ConvertFolder/");

    string inputctr = mediumfolder+ctrmodelname+string(".contour");
    string outputss = mediumfolder+ctrmodelname+string("_/");
    string ioTopofilename = string("../saveTopos/")+ctrmodelname+string(".topos");

    string s_outdir = outputss;
    string s_modelname = ctrmodelname;


    vector<double>points;vector<uint>edges;
    vector<int>CellSeq;
    int n_m,nMat;

    makemethodfolder(s_outdir);

    //original curve 0; curve net: 1;  cssegs: 2; cssegs New Ctr: 3; DP and Cell :4 ;

    isCJP = false;
    int methodid = 4;

    //compute: 0; compute and save: 1;  load: 2;
    int computeSaveLoad = 0;

    int endstep = 1;
    if(methodid < 4 ) endstep =0;


    cout<<"adsa"<<endl;
    vector<int>numofTopo;
    if(1)if(methodid)n_m = az2(inputctr,
                               outputss,
                               10,50,0,5,10,1,
                               MMprop_f2vmethod,endstep,numofTopo,CellSeq,nMat,computeSaveLoad,ioTopofilename);

    vector<int>TopoConstraintInd(numofTopo.size(),-1);

    vector<int>pickTopoInd = runDynamicProgramming(string("../protocal/")+ctrmodelname+string(".txt"),TopoConstraintInd);




    outputCombineSurface(outputss,numofTopo,pickTopoInd);

    return 0;
}





