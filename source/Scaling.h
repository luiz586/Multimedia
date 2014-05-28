#include "xBitstream.h"
#include "xEntropy.h"
#include "xHuffman.h"
#include "xSequence.h"
#include "xCfg.h"
#include "xTransform.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h> 


class Scaling{

public:
	int newWidth;
	int newHeight;
	int inputImageXSize;
	int inputImageYSize;
	float ratioX;
	float ratioY;
	
	static enum scalingTypes{
		simple = 0,
		bilinear = 1,
		bicubic = 2,
		lanczos = 3,
	};

public:
	Scaling();
	Scaling(int newWidth, int newHeight);
	void scaleImage(int scalingType);
	int simpleScaling(int16** inputPel);
	int Scaling::bilinearScaling(float samplingPointX, float samplingPointY, int16** inputPel);
	int Scaling::bicubicInterpolation();
	int16** scaleUsingType(int scaleType, int16** inputPel);
	FILE_HANDLE setOutputUri(int scalingType);

};
