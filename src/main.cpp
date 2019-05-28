#include <iostream>
#include <unistd.h>
#include <vector>
#include <string>
#include<dirent.h>
#include <sys/stat.h>
using namespace std;


bool isCJP = true;

int az2(string incontour, string outdir, int MMprop_f2vmethod,
        int endstep, int density, vector<int>&numofTopo, vector<int>&cellSeq, int &nMat, int computeSaveLoad, string iofile);
vector<int> runDynamicProgramming(string protocalname, vector<int>&TopoConstraintInd);
void outputCombineSurface(string outfolder, vector<int>&Cell2NTopo, vector<int>&pickTopoInd);


void makemethodfolder(string &foldername){
    struct stat sb;
    if(!(stat(foldername.data(), &sb) == 0 && S_ISDIR(sb.st_mode)))mkdir(foldername.data(),0777);
    return;
}

int main(int argc, char** argv)
{
    cout << argc << endl;


    int  MMprop_f2vmethod = 1;
    //string ctrmodelname = "outTestmapping";       MMprop_f2vmethod = 1;
    //string ctrmodelname = "ringc";                   MMprop_f2vmethod = 0;
    string ctrmodelname = "mug";                     MMprop_f2vmethod = 0;

    string mediumfolder("/Users/Research/Geometry/MM/MultiTopo/ConvertFolder/");

    string inputctr = mediumfolder+ctrmodelname+string(".contour");
    string outputss = mediumfolder+ctrmodelname+string("_/");
    string ioTopofilename = string("../saveTopos/")+ctrmodelname+string(".topos");

    string protocol_name = string("../protocal/")+ctrmodelname+string(".txt");



    //original curve 0; curve net: 1;  cssegs: 2; cssegs New Ctr: 3; DP and Cell :4 ;
    int methodid = 4;

    //compute: 0; compute and save: 1;  load: 2;
    int computeSaveLoad = 0;

    int ray_density = 2;
    isCJP = false;
    int endstep = 1;

    int c;
    optind=1;
    while ((c = getopt(argc, argv, "i:o:p:m:d:s:j")) != -1) {
        switch (c) {
        case 'i':
            inputctr = string(optarg);
            break;
        case 'o':
            outputss = string(optarg);
            break;
        case 'p':
            protocol_name = string(optarg);
            break;
        case 'm':
            MMprop_f2vmethod = atoi(optarg);
            break;
        case 'd':
            ray_density = atoi(optarg);
            break;
        case 's':
            computeSaveLoad = atoi(optarg);
            if(computeSaveLoad<0 || computeSaveLoad > 2)computeSaveLoad = 0;
            break;
        case 'j':
            isCJP = true;
            break;
        case '?':
            cout << "Bad argument setting!" << endl;
            break;
        }
    }


    cout<<"Input file: "<<inputctr<<endl;
    cout<<"Output path: "<<outputss<<endl;

    cout<<"Protocol name: "<<protocol_name<<endl;
    cout<<"Propagation method: "<<MMprop_f2vmethod<<endl;

    cout<<"Ray density: "<<ray_density<<endl;

    cout<<"Compute Save & Load "<<computeSaveLoad<<endl;



    vector<int>CellSeq;
    int n_m,nMat;

    makemethodfolder(outputss);


    if(methodid < 4 ) endstep =0;



    vector<int>numofTopo;
    if(methodid)n_m = az2(inputctr,
                               outputss,
                               MMprop_f2vmethod,endstep,ray_density,numofTopo,CellSeq,nMat,computeSaveLoad,ioTopofilename);

    vector<int>TopoConstraintInd(numofTopo.size(),-1);

    vector<int>pickTopoInd = runDynamicProgramming(protocol_name,TopoConstraintInd);



    outputCombineSurface(outputss,numofTopo,pickTopoInd);

    return 0;
}





