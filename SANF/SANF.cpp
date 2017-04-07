#include "SANF.h"
#include <fstream>
#include <algorithm>
#include <Eigen/SVD>
#include <iomanip>
using namespace Eigen;

// convert Double into  UChar
UChar castToUChar(Double i){
	Double result;
	if(i < 0)
		result = 0;
	else if(i > 255)
		result = 255;
	else
		result = i;
	return (UChar)(result + 0.5);
}
//public function definition
Void SANF::init(){
	createBuffer();
	creatPatchSet();

	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		PSNR_before[chan] = 0;
		PSNR_after[chan] = 0;
	}
}
Void SANF::destroy(){
	destroyBuffer();
	destroyPatchSet();
}
Void SANF::StructureFilter(ofstream& pFiltered, Int cur_frame){
	// after processing, write output file and clear Filltered_buffer
	pFiltered.seekp(cur_frame * FrameSize, ios::beg);
	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){

		SANFProcess((ComponentID)chan);
		//SANFProcess((ComponentID)chan);
		PSNR_before[chan] += calcPSNR(Ori[chan], Rec[chan],(ComponentID)chan);
		PSNR_after[chan] += calcPSNR(Ori[chan], Filtered_buffer[chan], (ComponentID)chan);
		pFiltered.write((Char*)Filtered_buffer[chan], param[chan].Width * param[chan].Height);
		cout<<calcPSNR(Ori[chan], Rec[chan],(ComponentID)chan)<<"/"<<left<<setw(13)<<calcPSNR(Ori[chan], Filtered_buffer[chan], (ComponentID)chan);
	}
	cout<<endl;
	if(cur_frame == param[COMPONENT_Y].FrameNum - 1){
		cout<<endl;
		cout<<left<<setw(12)<<"Average";
		for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
			PSNR_before[chan] /= param[COMPONENT_Y].FrameNum;
			PSNR_after[chan] /= param[COMPONENT_Y].FrameNum;
			cout<<PSNR_before[chan]<<"/"<<left<<setw(13)<<PSNR_after[chan];
		}
		cout<<endl;
	}
	//printInfo();
	clearBuffer();
}
Void SANF::fillBuffer(ifstream& origFile, ifstream& recFile, Int cur_frame){

	origFile.seekg(cur_frame * FrameSize, ios::beg);
	recFile.seekg(cur_frame * FrameSize, ios::beg);

	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		origFile.read(Ori_buffer[chan], param[chan].Width * param[chan].Height);
		recFile.read(Rec_buffer[chan], param[chan].Width * param[chan].Height);
		Ori[chan] = (UChar*)Ori_buffer[chan];
		Rec[chan] = (UChar*)Rec_buffer[chan];
	}
}
Void SANF::clearBuffer(){
	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		memset(Filtered_buffer[chan], 0, param[chan].Width * param[chan].Height);
	}
}
Void SANF::setParameter(ComponentID id, Int uiWidth, Int uiHeight, Double Th, Int TempWin, Int SearchWin, Int NumPatch, Int SlidDis, Int FrameNumber){

	if(id == COMPONENT_Y){
		FrameSize = uiWidth * uiHeight * 3 / 2;
	}
	param[id].FrameNum = FrameNumber;
	param[id].Width = uiWidth;
	param[id].Height = uiHeight;
	param[id].Threshold = Th;
	param[id].TemplateWindowSize = TempWin;
	param[id].SearchWindowSize = SearchWin;
	param[id].NumberPatches = NumPatch;
	param[id].SlidingDis = SlidDis;
	param[id].HeiNum = uiHeight - TempWin + 1;
	param[id].WidNum = uiWidth - TempWin + 1;
	param[id].totalNum = param[id].HeiNum * param[id].WidNum;
}
Double SANF::calcPSNR(UChar* Ori, UChar* Result, ComponentID id){
	Int uiHeight = param[id].Height;
	Int uiWidth = param[id].Width;

	Double SSE = 0.0, MSE, PSNR;
	int offset;

	for(int rows = 0; rows < uiHeight; rows++)
		for(int cols = 0; cols < uiWidth; cols++){
			offset = rows * uiWidth + cols;
			SSE += ((Ori[offset] - Result[offset]) * (Ori[offset] - Result[offset]));
		}

		if(SSE == 0.0)
			return 0.0;
		else{
			MSE = SSE / (Double)(uiWidth * uiHeight);
			PSNR = 10.0 * log10((255 * 255) / MSE);
			return PSNR;
		}
}
Void SANF::SANFProcess(ComponentID id){
	UChar** curPatchSet = patchSet[id];
	Int uiWidth, uiHeight, NumPch, tempWin, uiSlidDis, uiWidNum;

	uiWidth = param[id].Width;
	uiHeight = param[id].Height;
	NumPch = param[id].NumberPatches;
	tempWin = param[id].TemplateWindowSize;
	uiSlidDis = param[id].SlidingDis;
	uiWidNum = param[id].WidNum;

	//construct patchSet & ImgToPatchSet
	ImgToPatchSet(Rec[id], id);
	vector<Int> indexArray(NumPch, -1);

	//define intermediate variables: tempImg & weight
	Double* tempImg = new Double[uiWidth * uiHeight];
	UInt* weight = new UInt[uiWidth * uiHeight];
	memset(tempImg, 0, sizeof(Double) * uiWidth * uiHeight);
	memset(weight, 0, sizeof(Int) * uiWidth * uiHeight);

	//for loops: 
	//	for every overlapped block, do patchSearch and obtain currArray,
	//	Output: indexArry in order to PatchToImg
	Int cycleWidNum = uiWidth / uiSlidDis - 1;
	Int cycleHeiNum = uiHeight / uiSlidDis - 1;

	for(Int row = 0; row < cycleHeiNum; row++){
		for(Int col = 0; col < cycleWidNum; col++){
			// 1. do patchSearch and return most similar patch's index array
			indexArray = patchSearch(id, row * uiSlidDis, col * uiSlidDis, curPatchSet);

			// 2. put process results into temp image, and update weights
			//define temporary variable: rowIndex-> patchIndex row index
			Int rowIndex, colIndex, ImRow, ImCol, Offest;
			for(Int idx = 0; idx < NumPch; idx++){
				rowIndex = indexArray[idx] / uiWidNum;
				colIndex = indexArray[idx] % uiWidNum;
				for(Int tmIdx = 0; tmIdx < tempWin * tempWin; tmIdx++){
					ImRow = rowIndex + tmIdx / tempWin;
					ImCol = colIndex  + tmIdx % tempWin;
					Offest = getPos(ImRow, ImCol, uiWidth);
					tempImg[Offest] += curArray[id][idx][tmIdx];
					weight[Offest] += 1;	
				}
			}
			//clear indexArray(),continue next iteration
			indexArray.clear();
		}
	}
	//every value in temporary image divide its weights, and get result
	for(Int i = 0; i < uiWidth * uiHeight; i++){
		Double tepValue = tempImg[i] / weight[i];
		Filtered_buffer[id][i] = castToUChar(tepValue);
	}
	if(tempImg){
		delete[] tempImg;
		tempImg = NULL;
	}
	if(weight){
		delete[] weight;
		weight = NULL;
	}	
}

/************************************************************************/
/* definition of private functions                                                               
/************************************************************************/
Void SANF::createBuffer(){
	Int uiWidth, uiHeight;
	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		uiWidth = param[chan].Width;
		uiHeight = param[chan].Height;
		Ori_buffer[chan] = new Char[uiWidth * uiHeight];
		Rec_buffer[chan] = new Char[uiWidth * uiHeight];
		Filtered_buffer[chan] = new UChar[uiWidth * uiHeight];
		memset(Ori_buffer[chan], 0, uiWidth * uiHeight);
		memset(Rec_buffer[chan], 0, uiWidth * uiHeight);
		memset(Filtered_buffer[chan], 0, uiWidth * uiHeight);

	}
}
Void SANF::creatPatchSet(){
	//create patchSet
	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		patchSet[chan] = new UChar*[param[chan].totalNum];
		for(Int index = 0; index < param[chan].totalNum; index++){
			patchSet[chan][index] = new UChar[param[chan].TemplateWindowSize * param[chan].TemplateWindowSize];
			memset(patchSet[chan][index], 0, param[chan].TemplateWindowSize * param[chan].TemplateWindowSize);
		}
	}

	//create curArry
	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		curArray[chan] = new Double*[param[chan].NumberPatches];
		for(Int index = 0; index < param[chan].NumberPatches; index++){
			curArray[chan][index] = new Double[param[chan].TemplateWindowSize * param[chan].TemplateWindowSize];
			memset(curArray[chan][index], 0, sizeof(Double) * param[chan].TemplateWindowSize * param[chan].TemplateWindowSize);
		}
	}
}
Void SANF::destroyBuffer(){
	for(UInt chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		if(Ori_buffer[chan]){
			delete[] Ori_buffer[chan];
			Ori_buffer[chan] = NULL;
			Ori[chan] = NULL;
		}
		if(Rec_buffer[chan]){
			delete[] Rec_buffer[chan];
			Rec_buffer[chan] = NULL;
			Rec[chan] = NULL;
		}
		if(Filtered_buffer){
			delete[] Filtered_buffer[chan];
			Filtered_buffer[chan] = NULL;
		}
	}
}
Void SANF::destroyPatchSet(){
	//destroy patchSet
	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		for(Int index = 0; index < param[chan].totalNum; index++){
			if(patchSet[chan][index]){
				delete[] patchSet[chan][index];
				patchSet[chan][index] = NULL;
			}
		}
		if(patchSet[chan]){
			delete[] patchSet[chan];
			patchSet[chan] = NULL;
		}
	}
	//destroy curArray
	for(Int chan = 0; chan < MAX_NUM_COMPONENT; chan++){
		for(Int index = 0; index < param[chan].NumberPatches; index++){
			if(curArray[chan][index]){
				delete[] curArray[chan][index];
				curArray[chan][index] = NULL;
			}
		}
		if(curArray[chan]){
			delete[] curArray[chan];
			curArray[chan] = NULL;
		}
	}
}
Int SANF::getPos(Int x_pos, Int y_pos, Int uiWidth){
	return x_pos * uiWidth + y_pos;
}
Void SANF::ImgToPatchSet(UChar* RecImg, ComponentID id){
	Int uiHeight, uiWidth, uiSlidDis, tempWin;
	Int uiHeiNum, uiWidNum;

	uiWidth = param[id].Width;
	uiHeight = param[id].Height;
	uiSlidDis = param[id].SlidingDis;
	tempWin = param[id].TemplateWindowSize;
	uiHeiNum = param[id].HeiNum;
	uiWidNum = param[id].WidNum;


	//imRow£ºrow index of current pixel
	//imCol£ºcol index of current pixel
	//patchIndex£ºpatch index of current patch
	//pixelPos£ºpixel index of current pixel in current patch
	Int imRow, imCol, patchIndex, pixelPos;
	for(Int indexRow = 0; indexRow < uiHeiNum; indexRow++)
		for(Int indexCol = 0; indexCol < uiWidNum; indexCol++){
			patchIndex = indexRow * uiWidNum + indexCol;
			//for loops: every template patch
			for(Int pixRow = 0; pixRow < tempWin; pixRow++)
				for(Int pixCol = 0; pixCol < tempWin; pixCol++){
					imRow = indexRow  + pixRow;
					imCol = indexCol  + pixCol;
					pixelPos = pixRow * tempWin+ pixCol;
					patchSet[id][patchIndex][pixelPos] = RecImg[getPos(imRow, imCol, uiWidth)];
				}
		}
}

/***************************************************************************************/
/* patchSearch£ºfind matching group, and do SVD decomposition and use threshold 
/*	            to process group matrix. 
/* return: index of Group patches in image, curArray store the processed results
/***************************************************************************************/
vector<Int> SANF::patchSearch(ComponentID id, Int curRowIdx, Int curColIdx, UChar** curPatchSet){

	//acquire parameter informations
	/*********************************************************************************/
	/* searchWin£ºa half of search window size, searchWindowSize = 2 * searchWin + 1 
	/* tempWin£ºblock length of side
	/*********************************************************************************/
	Int searchWin, tempWin, uiWidth, uiHeight, NumPch, uiWidNum; 

	searchWin = param[id].SearchWindowSize;
	tempWin = param[id].TemplateWindowSize;
	uiWidth = param[id].Width;
	uiHeight = param[id].Height;
	NumPch = param[id].NumberPatches;
	uiWidNum = param[id].WidNum;

	// according to the input curRowIndex and curColIndex, obtain the effective search range
	//(top left corner of patch's row and col informations in image), calculate the corresponding
	//index in patchSet, save it into patchIndex  (it can be defined as vector type)
	Int lboud = max(0, curColIdx - searchWin);
	Int rboud = min(uiWidth  - tempWin, curColIdx + searchWin);
	Int uboud = max(0, curRowIdx - searchWin);
	Int dboud = min(uiHeight - tempWin, curRowIdx + searchWin);		//all 4 bounds are available, don't forget >=  <=
	vector<Int> patchIndex;
	Int index;
	for(int row = uboud; row <= dboud; row++)
		for(int col = lboud; col<= rboud; col++){
			index = pelToIndex(row, col, uiWidNum);
			patchIndex.push_back(index);
		}
	Int curPchIdx = pelToIndex(curRowIdx, curColIdx, uiWidNum);

	//calculate MSE between current patch and available patchIndex in patchSet, 
	// and save it in MSE vector
	vector<Double> MSE;
	Int forPchIdx; // patchIndex found in loops
	Double SSE;
	for(Int i = 0; i < (Int)patchIndex.size();i++){
		SSE = 0;
		forPchIdx = patchIndex[i];

		//find corresponding block in patchSet, calculate MMS.
		// define intermediate variable SSE
		for(Int tmIdx = 0; tmIdx < tempWin * tempWin; tmIdx++){
			SSE += (curPatchSet[forPchIdx][tmIdx] - curPatchSet[curPchIdx][tmIdx]) * (curPatchSet[forPchIdx][tmIdx] - curPatchSet[curPchIdx][tmIdx]);
		}
		SSE = SSE / (tempWin * tempWin);
		MSE.push_back(SSE);
	}

	//sort MSE, and return sorted index of elements in MSE vector
	vector<Int> sortIndex;
	sortIndex = sort_indexes(MSE);
	vector<Int> result(NumPch);			//result saves curArray patches index in patchSet


	//choose first NumPch results from curPatchSet, save in curArray
	for(Int i = 0; i < NumPch; i++){
		result[i] = patchIndex[sortIndex[i]];
		for(Int j = 0; j < tempWin * tempWin; j++){
			curArray[id][i][j] = (Double)curPatchSet[result[i]][j];
		}
	}
	/************************************************************************/
	/* Important!!!!!
	/* This piece of code is very important, remove it can lead cavities in image
	/************************************************************************/
	result[0] = curPchIdx;
	for(Int j = 0; j < tempWin * tempWin; j++){
		curArray[id][0][j] = (Double)curPatchSet[curPchIdx][j];
	}

	//SVD decomposition to curArray, first convert curArray into Matrix Type, and convert back after process 
	svdDecomp(curArray[id], tempWin * tempWin, NumPch, param[id].Threshold);

	return result;
}
/************************************************************************/
/* Utility Functions                                                    */
/************************************************************************/
//sort elements according to index
vector<Int> SANF::sort_indexes(const vector<Double> &v){
	//initialize original index locations
	vector<Int> idx(v.size());
	for(Int i = 0; i < (Int)idx.size(); i++)	idx[i] = i;

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),[& v](Int i1, Int i2){ return v[i1] < v[i2]; });
	return idx;
}
Int SANF::pelToIndex(Int pelRow, Int pelCol, Int wNum){
	return pelRow * wNum + pelCol;
}
//incoming parameters: matrix rows and cols
//rows: TW * TW, cols: NumPch
Void SANF::svdDecomp(Double** inputArray, Int Rows, Int Cols, Double Th){
	//covert input arrays into matrix
	MatrixXd temp(Rows, Cols);
	for(Int col = 0; col < Cols; col++){
		for(Int row = 0; row < Rows; row++){
			temp(row, col) = inputArray[col][row];
		}
	}

	//SVD decomposition, obtain U, V matrix and singular vector
	JacobiSVD<MatrixXd> svd(temp, ComputeThinU | ComputeThinV);
	MatrixXd U = svd.matrixU();
	VectorXd singular = svd.singularValues();
	MatrixXd V = svd.matrixV();

	//process singular using threshold
	for(int i = 0; i < singular.size(); i++){
		if(singular(i) <= Th){
			singular(i) = 0;
		}
	}
	MatrixXd result = U * singular.asDiagonal() * V.transpose();


	//save resultant matrix into curArray
	for(Int col = 0; col < Cols; col++){
		for(Int row = 0; row < Rows; row++){
			inputArray[col][row] = result(row, col);
		}
	}
}