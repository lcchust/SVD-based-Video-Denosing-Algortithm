#include "TypeDef.h"
#include "SANF.h"
#include <time.h>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

//function declaration
bool parseCfg(string cfgFilePath);

//define global variables
string origFile, recFile, filterFile;
Int WIDTH, HEIGHT, FRAME_NUM, FRAME_SIZE;
Double THRESHOLD;
Int TEMPLATEWINDOWSIZE;
Int SEARCHWINDOWSIZE;
Int NUMBERPATCHES;
Int SLIDINGDIS;

int main(int argc, char* argv[])
{
	//parse config file
	string cfgName;
	if(argc > 1)
		cfgName = argv[1];
	parseCfg(cfgName);
	
	//open file and acquire file pointer
	ifstream pOrigFile, pRecFile;
	pOrigFile.open(origFile, ios::in | ios::binary);
	pRecFile.open(recFile, ios::in | ios::binary);
	if(pOrigFile.fail() || pRecFile.fail()){
		cout<<"OrigFile Or RecFile Error!"<<endl;
		return -1;
	}

	ofstream pFiltered;
	pFiltered.open(filterFile, ios::out | ios::binary);
	if(pFiltered.fail()){
		cout<<"FilterFile Error!"<<endl;
		return -1;
	}

	//starting time
	Double dResult;
	clock_t lBefore = clock();

	SANF sanf;
 	sanf.setParameter(COMPONENT_Y, WIDTH, HEIGHT, THRESHOLD, TEMPLATEWINDOWSIZE, SEARCHWINDOWSIZE, NUMBERPATCHES, SLIDINGDIS, FRAME_NUM);
	sanf.setParameter(COMPONENT_Cb, WIDTH/2, HEIGHT/2, THRESHOLD, TEMPLATEWINDOWSIZE, SEARCHWINDOWSIZE, NUMBERPATCHES, SLIDINGDIS,FRAME_NUM);
	sanf.setParameter(COMPONENT_Cr,WIDTH/2, HEIGHT/2, THRESHOLD, TEMPLATEWINDOWSIZE, SEARCHWINDOWSIZE, NUMBERPATCHES, SLIDINGDIS,FRAME_NUM);
 	sanf.init();
 	
	cout<<left<<setw(20)<<"Before/After";
	cout<<left<<setw(20)<<"Y";
	cout<<left<<setw(20)<<"U";
	cout<<left<<setw(20)<<"V"<<endl;
 	for(int cur_frame = 0; cur_frame < FRAME_NUM; cur_frame++){
		cout<<left<<setw(12)<<cur_frame;
 		sanf.fillBuffer(pOrigFile, pRecFile, cur_frame);
		sanf.StructureFilter(pFiltered, cur_frame);
 	}
	dResult = (Double)(clock()-lBefore) / CLOCKS_PER_SEC;
	printf("\nTotal Time: %12.3f sec.\n", dResult);
 	sanf.destroy();
	pOrigFile.close();
	pRecFile.close();
	pFiltered.close();
}

//parse input config file
bool parseCfg(string cfgFilePath){

	// read config file
	ifstream cfgFile(cfgFilePath);

	if(! cfgFile){
		cout<<"can't open cfg file"<<endl;
	}
	string value[cfgParamNum];
	string tmp;

	for(int i = 0; i < cfgParamNum; i++){
		getline(cfgFile, tmp);
		size_t pos = tmp.find('=');
		if(pos == string::npos) return false;
		value[i] = tmp.substr(pos + 2);	
	}
	cfgFile.close();

	//Assign parameters
	origFile = value[0];
	recFile = value[1];
	filterFile = value[2];

	stringstream transfer;
	// WIDTH
	transfer<<value[3];
	transfer>>WIDTH;
	transfer.clear();
	//HEIGHT
	transfer<<value[4];
	transfer>>HEIGHT;
	transfer.clear();
	//FRAME_NUM
	transfer<<value[5];
	transfer>>FRAME_NUM;
	transfer.clear();
	//Threshold
	transfer<<value[6];
	transfer>>THRESHOLD;
	transfer.clear();
	//TemplateWindowSize
	transfer<<value[7];
	transfer>>TEMPLATEWINDOWSIZE;
	transfer.clear();
	//SearchWindowSize
	transfer<<value[8];
	transfer>>SEARCHWINDOWSIZE;
	transfer.clear();
 	//NumberPatches
 	transfer<<value[9];
 	transfer>>NUMBERPATCHES;
 	transfer.clear();
	//SlidingDis
	transfer<<value[10];
	transfer>>SLIDINGDIS;
	transfer.clear();

	FRAME_SIZE = WIDTH * HEIGHT * 3 / 2;

	return true;
}
