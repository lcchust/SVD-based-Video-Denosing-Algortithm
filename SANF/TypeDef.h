#ifndef __TYPEDEF__
#define __TYPEDEF__
#include <string>
#include <iomanip>
#define PATCHNUM	60
#define cfgParamNum 11

typedef       char                Char;
typedef       unsigned char       UChar;
typedef       short               Short;
typedef       unsigned short      UShort;
typedef       int                 Int;
typedef       unsigned int        UInt;
typedef       double              Double;
typedef       float               Float;
typedef		  void				  Void;
/*typedef		  string			  String;*/

/// chroma formats (according to semantics of chroma_format_idc)
enum ChromaFormat
{
	CHROMA_400        = 0,
	CHROMA_420        = 1,
	CHROMA_422        = 2,
	CHROMA_444        = 3,
	NUM_CHROMA_FORMAT = 4
};

enum ComponentID
{
	COMPONENT_Y       = 0,
	COMPONENT_Cb      = 1,
	COMPONENT_Cr      = 2,
	MAX_NUM_COMPONENT = 3
};

struct Parameter{
	Double Threshold;
	Int TemplateWindowSize;
	Int SearchWindowSize;
	Int NumberPatches;
	Int FrameNum;
	Int SlidingDis;
	Int Width;
	Int Height;
	Int HeiNum; //patchSet: patch numbers in height
	Int WidNum; //patchSet: patch numbers in width
	Int totalNum; //total patches in patchSet

};
#endif