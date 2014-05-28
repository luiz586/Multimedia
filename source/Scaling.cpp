#include "Scaling.h"
using namespace AVlib;
using namespace std;

Scaling::Scaling(int newWidth, int newHeight){
	this->newWidth = newWidth;
	this->newHeight = newHeight;
}


void Scaling::scaleImage(int scalingType){
	FILE_HANDLE inputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538.bmp"), TEXT("r"));
	FILE_HANDLE outputFile = setOutputUri(scalingType);

	xImg* inputImage = new xImg;
	xImg* outputImage = new xImg;
	inputImage->getBMP(inputFile);
    
	this->inputImageXSize = inputImage->m_Cmp[0]->m_SizeX;
	this->inputImageYSize = inputImage->m_Cmp[0]->m_SizeY;

	this->ratioX = (this->inputImageXSize*1.0)/(newWidth*1.0);
	this->ratioY = (this->inputImageYSize*1.0)/(newHeight*1.0);
	cout<<"ratioX: "<<ratioX<<", ratioY: "<<ratioY<<endl;

	// Tworzy taki sam obiekt xImg jak obrazu wejsciowego tylko ze zmieniony rozmiarami
	outputImage->createResizedCompatibile(inputImage, newWidth, newHeight);

	// tablica wskaŸników ze sk³adowymi obrazu , np. [R G B A]
	xCmp** m_cmp = outputImage->m_Cmp;  

	// Dla ka¿dej sk³adowej wykonaj skalowanie
	int16** outputPel;
	for(int i=0; i<inputImage->m_CmpTotal; i++){
		xCmp* inputComponent = inputImage->m_Cmp[i];
		int16** inputPel = inputComponent->m_Pel;

		outputPel = scaleUsingType(scalingType, inputPel);
		outputImage->m_Cmp[i]->m_Pel = outputPel;
	}

	outputImage->putBMP(outputFile);
	

	inputImage->destroy();
	outputImage->destroy();
	delete inputImage;
	delete outputImage;
}

int16** Scaling::scaleUsingType(int scaleType, int16** inputPel){	
	int16** outputPel = new int16*[newHeight]; 
	int samplingPointX, samplingPointY;

	switch(scaleType){
	case scalingTypes(simple):
		for(int j=0; j<newHeight; j++){
			outputPel[j] = new int16[newWidth];
			for(int k=0; k<newWidth; k++){
				samplingPointX = k*ratioX;
				samplingPointY = j*ratioY;	
				outputPel[j][k] = inputPel[samplingPointY][samplingPointX];
			}
		}
		break;
	case scalingTypes(bilinear):
		for(int j=0; j<newHeight; j++){
			outputPel[j] = new int16[newWidth];
			for(int k=0; k<newWidth; k++){
				samplingPointX = k*ratioX;
				samplingPointY = j*ratioY;					
				outputPel[j][k] = bilinearScaling(samplingPointX, samplingPointY, inputPel);
			}
		}
		break;
	case scalingTypes(bicubic):
		break;
	case scalingTypes(lanczos):
		break;
	default:
		break;

	}

	return outputPel;
}

int Scaling::bilinearScaling(float samplingPointX, float samplingPointY, int16** inputPel){
	int16 refPointsValues[4];
	int horizontalInterpolA, horizontalInterpolB, verticalInterpol;

	float flooredX;
	float flooredY;
	int ceiledX = ceil(samplingPointX);
	int ceiledY = ceil(samplingPointY);

	// Czêœci dziesiêtne wyliczonych punktów
	float a = modf(samplingPointX, &flooredX);
	float b = modf(samplingPointY, &flooredY);

	// Warunki brzegowe - interpolacja zwykla
	if(ceiledX >= inputImageXSize || ceiledY >= inputImageYSize){
		verticalInterpol = inputPel[(int)flooredY][(int)flooredX];
	}
	// Interpolacja biliniowa
	else{
		refPointsValues[0] = inputPel[(int)flooredY][(int)flooredX];					
		refPointsValues[1] = inputPel[(int)flooredY][ceiledX];				
		refPointsValues[2] = inputPel[ceiledY][(int)flooredX];					
		refPointsValues[3] = inputPel[ceiledY][ceiledX];				

		horizontalInterpolA = (1-a) * refPointsValues[0] + a * refPointsValues[1];
		horizontalInterpolB = (1-a) * refPointsValues[2] + a * refPointsValues[3];

		verticalInterpol = (1-b) * horizontalInterpolA + b * horizontalInterpolB;
	}
	
	return verticalInterpol;
}

FILE_HANDLE Scaling::setOutputUri(int scalingType){
	FILE_HANDLE outputFile;

	switch(scalingType){
	case scalingTypes(simple):
		outputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538_simple.bmp"), TEXT("w"));
		break;
	case scalingTypes(bilinear):
		outputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538_bilinear.bmp"), TEXT("w"));
		break;
	case scalingTypes(bicubic):
		outputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538_bisubic.bmp"), TEXT("w"));
		break;
	case scalingTypes(lanczos):
		outputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538_lanczos.bmp"), TEXT("w"));
		break;
	default:
		break;
	}

	return outputFile;
}