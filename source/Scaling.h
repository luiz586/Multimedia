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
using namespace AVlib;


class Scaling{

public:
	int newWidth;
	int newHeight;
	int inputImageXSize;
	int inputImageYSize;
	float ratioX;
	float ratioY;
	xImg* mInputImage;
	xImg* mOutputImage;
	
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
	int bilinearScaling(float samplingPointX, float samplingPointY, int16** inputPel);
	int bicubicScaling(float samplingPointX, float samplingPointY, int16** inputPel);
	int lanczosScaling(float samplingPointX, float samplingPointY, int16** inputPel);
	int lanczos3_filter(int t);
	int clean(int t);
	double sinc(double x);
	int cubicInterpolate (int16 samples[4], float dx);
	int16** scaleUsingType(int scaleType, int16** inputPel);
	FILE_HANDLE setOutputUri(int scalingType);
	int secondCubic(float samplingPointX, float samplingPointY, int16** inputPel);
	float rFunction(float x);
	float pFunction(float x);

};
