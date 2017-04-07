#ifndef __NLMEANS__
#define __NLMEANS__
#include "TypeDef.h"
#include <iostream>
#include <vector>
using namespace std;
/************************************************************************/
/* How to use this Class ?
	1. create an instance of the class
	2. setParameter(...params): set param(Y Cb Cr)
	3. init(): initialize YUV buffer and patchSet, curArray
	4. before every frame's process, call fillBuffer(...params) to store YUV inform
	5. StructureFilter(...param) to filter rec image(Y,Cb,Cr)
	6. destroy(): destroy YUV buffer and patchSet, curArray
/************************************************************************/
class SANF{
protected:
//Parameters
	Parameter param[MAX_NUM_COMPONENT];
	Int FrameSize;
	Double PSNR_before[MAX_NUM_COMPONENT];
	Double PSNR_after[MAX_NUM_COMPONENT];
//Buffers
	Char* Rec_buffer[MAX_NUM_COMPONENT]; // Y:1 Cb:2 Cr:3
	UChar* Rec[MAX_NUM_COMPONENT];
	Char* Ori_buffer[MAX_NUM_COMPONENT];
	UChar* Ori[MAX_NUM_COMPONENT];
	UChar* Filtered_buffer[MAX_NUM_COMPONENT];
	UChar** patchSet[MAX_NUM_COMPONENT];			//keep all possible patches in image, first dim represents patch index, second dim represents patch value			
	Double** curArray[MAX_NUM_COMPONENT];			//keep NumPch patches obtained by patchSearch

//Private function
	Void createBuffer();
	Void creatPatchSet();
	Void destroyBuffer();
	Void destroyPatchSet();
	Void clearPatchSet();
	Int getPos(Int x_pos, Int y_pos, Int uiWidth);
	vector<Int> sort_indexes(const vector<Double> &v);
	Void ImgToPatchSet(UChar* RecImg, ComponentID id);
	vector<Int> patchSearch(ComponentID id, Int curRowIdx, Int curColIdx, UChar** curPatchSet);
	Int pelToIndex(Int pelRow, Int pelCol, Int wNum);
	Void svdDecomp(Double** inputArray, Int Rows, Int Cols, Double Th);
	Void SANFProcess(ComponentID id);
	
public:
	SANF(){};
	~SANF(){};
	Void setParameter(ComponentID id, Int uiWidth, Int uiHeight, Double Th, Int TempWin, Int SearchWin, Int NumPatch, Int SlidDis, Int FrameNumber);
	Void init();
	Void destroy();
	Void clearBuffer();
	Void fillBuffer(ifstream& origFile, ifstream& recFile, Int cur_frame);
	Double calcPSNR(UChar* Ori, UChar* Result,ComponentID id);
	Void StructureFilter(ofstream& pFiltered, Int cur_frame);
};
#endif