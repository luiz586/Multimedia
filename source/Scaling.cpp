#include "Scaling.h"
using namespace AVlib;
using namespace std;

#define M_PI 3.14159265358979323846
#define M_LANCZOS_FILTER_SIZE 3

Scaling::Scaling(int newWidth, int newHeight){
	this->newWidth = newWidth;
	this->newHeight = newHeight;
}


void Scaling::scaleImage(int scalingType){
	FILE_HANDLE inputFile = x_fopen(TEXT("img\\lena.bmp"), TEXT("r"));
	FILE_HANDLE outputFile = setOutputUri(scalingType);

	mInputImage = new xImg;
	mOutputImage = new xImg;
	mInputImage->getBMP(inputFile);
    
	this->inputImageXSize = mInputImage->m_Cmp[0]->m_SizeX;
	this->inputImageYSize = mInputImage->m_Cmp[0]->m_SizeY;

	this->ratioX = (this->inputImageXSize*1.0)/(newWidth*1.0);
	this->ratioY = (this->inputImageYSize*1.0)/(newHeight*1.0);
	cout<<"ratioX: "<<ratioX<<", ratioY: "<<ratioY<<endl;

	// Tworzy taki sam obiekt xImg jak obrazu wejsciowego tylko ze zmieniony rozmiarami
	mOutputImage->createResizedCompatibile(mInputImage, newWidth, newHeight);

	// tablica wskaŸników ze sk³adowymi obrazu , np. [R G B A]
	xCmp** m_cmp = mOutputImage->m_Cmp;  

	// Dla ka¿dej sk³adowej wykonaj skalowanie
	int16** outputPel;
	for(int i=0; i<mInputImage->m_CmpTotal; i++){
		xCmp* inputComponent = mInputImage->m_Cmp[i];
		int16** inputPel = inputComponent->m_Pel;

		outputPel = scaleUsingType(scalingType, inputPel);
		mOutputImage->m_Cmp[i]->m_Pel = outputPel;
	}

	mOutputImage->putBMP(outputFile);
	

	mInputImage->destroy();
	mOutputImage->destroy();
	delete mInputImage;
	delete mOutputImage;
}

int16** Scaling::scaleUsingType(int scaleType, int16** inputPel){	
	int16** outputPel = new int16*[newHeight]; 
	float samplingPointX, samplingPointY;

	switch(scaleType){
	case scalingTypes(simple):
		for(int j=0; j<newHeight; j++){
			outputPel[j] = new int16[newWidth];
			for(int k=0; k<newWidth; k++){
				samplingPointX = k*ratioX;
				samplingPointY = j*ratioY;	
				outputPel[j][k] = inputPel[(int)samplingPointY][(int)samplingPointX];
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
		for(int j=0; j<newHeight; j++){
			outputPel[j] = new int16[newWidth];
			for(int k=0; k<newWidth; k++){
				samplingPointX = k*ratioX;
				samplingPointY = j*ratioY;					
				outputPel[j][k] = bicubicScaling(samplingPointX, samplingPointY, inputPel);
				//outputPel[j][k] = secondCubic(samplingPointX, samplingPointY, inputPel);
			}
		}
		break;
	case scalingTypes(lanczos): 
		for(int j=0; j<newHeight; j++){
			outputPel[j] = new int16[newWidth];
			for(int k=0; k<newWidth; k++){
				samplingPointX = k*ratioX;
				samplingPointY = j*ratioY;					
				outputPel[j][k] = lanczosScaling(samplingPointX, samplingPointY, inputPel);
			}
		}
		break;
	default:
		break;

	}

	return outputPel;
}

int Scaling::lanczos3_filter(int t)
{
   if (t < 0.0f)
      t = -t;

   if (t < 3.0f)
      return clean(sinc(t) * sinc(t / 3.0f));
   else
      return (0.0f);
}


int Scaling::clean(int t)
{
   const int EPSILON = .0000125f;
   if (abs(t) < EPSILON)
      return 0.0f;
   return (int)t;
}

double Scaling::sinc(double x)
{
   x = (x * M_PI);

   if ((x < 0.01f) && (x > -0.01f))
      return 1.0f + x*x*(-1.0f/6.0f + x*x*1.0f/120.0f);

   return sin(x) / x;
}

int Scaling::lanczosScaling(float samplingPointX, float samplingPointY, int16** inputPel) {
	int outputValue = 0;

	for (int i = int(abs(samplingPointY) - M_LANCZOS_FILTER_SIZE + 1); 
		i < int(abs(samplingPointY) + M_LANCZOS_FILTER_SIZE); i++)  {
			for (int j = int(abs(samplingPointX) - M_LANCZOS_FILTER_SIZE + 1); 
				j < int(abs(samplingPointX) + M_LANCZOS_FILTER_SIZE); i++) {
					if(j == inputImageXSize || i == inputImageYSize || j < 0 || i < 0){
						return inputPel[(int)samplingPointY][(int)samplingPointX];
					}
					//outputValue += *(inputPel[(int)samplingPointY] + (int)samplngPointX) * lanczos3_filter(samplngPointX - i) * lanczos3_filter(samplingPointY - j);
					outputValue += inputPel[i][j] * lanczos3_filter(samplingPointY - i) * lanczos3_filter(samplingPointX - j);
			}

	}

	return outputValue;
}

int Scaling::bicubicScaling(float samplingPointX, float samplingPointY, int16** inputPel){	
	float flooredX;
	float flooredY;

	// Czêœci dziesiêtne wyliczonych punktów
	float dx = modf(samplingPointX, &flooredX);
	float dy = modf(samplingPointY, &flooredY);

	
	// Pierwszy punkt z szestanstu s¹siaduj¹cych
	int startPointX = flooredX - 1;
	int startPointY = flooredY - 1;

	// Warunki brzegowe
	if(startPointX < 0 || startPointY < 0 
		|| startPointX >= (this->inputImageXSize-3) 
		|| startPointY >= (this->inputImageYSize-3) ){			
			return *(inputPel[(int)samplingPointY] + (int)samplingPointX);
	}

	// Pobieram wartoœci 16 s¹siaduj¹cych próbek
	int16 columns[4][4];
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			columns[i][j] = inputPel[startPointY + j][startPointX + i];
		}
	}

	int16 rows[4];
	rows[0] = cubicInterpolate(columns[0], dx);
	rows[1] = cubicInterpolate(columns[1], dx);
	rows[2] = cubicInterpolate(columns[2], dx);
	rows[3] = cubicInterpolate(columns[3], dx);

	return cubicInterpolate(rows, dy);
}

int Scaling::cubicInterpolate(int16 cols[4], float dx) {
	//return p[1] + 0.5 * dx*(p[2] - p[0] + dx*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + dx*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	float a1 = 0.5;
	float a2 = 1.5;
	float a3 = 2.5;
	float x2 = dx*dx;
	float x3 = x2*dx;

	return ( -a1 * cols[0] + a2 * cols[1] - a2 * cols[2] + a1 * cols[3] ) * x3 
		+ (cols[0] - a3 * cols[1] + 2 * cols[2] - a1 * cols[3] ) * x2 
		+ (-a1 * cols[0] + a1 * cols[2]) * dx + cols[1];
}


/*int Scaling::secondCubic(float samplingPointX, float samplingPointY, int16** inputPel){
	float flooredX;
	float flooredY;

	// Czêœci dziesiêtne wyliczonych punktów
	float dx = modf(samplingPointX, &flooredX);
	float dy = modf(samplingPointY, &flooredY);

	int16 samples[4][4];
	// Pierwszy punkt z szestanstu s¹siaduj¹cych
	//int startPointX = flooredX - 1;
	//int startPointY = flooredY - 1;

	// Warunki brzegowe
	if(flooredX < 1 || flooredY < 1 
		|| flooredX >= (this->inputImageXSize-2) 
		|| flooredY >= (this->inputImageYSize-2) ){			
			return *(inputPel[(int)samplingPointY] + (int)samplingPointX);
	}
	// Pobieram wartoœci 16 s¹siaduj¹cych próbek
	float result = 0;
	for(int m=-1; m<3; m++){
		for(int n=-1; n<3; n++){
			result += inputPel[(int)flooredY + m][(int)flooredX + n] * rFunction(m-dx) * rFunction(dy-n);
		}
	}

	return result;

}

float Scaling::rFunction(float x){
	float x1 = pow(x+2,3);
	float x2 = pow(x+1,3);
	float x3 = pow(x,3);
	float x4 = pow(x-1,3);
	return 1/6 * ( pFunction(pow(x+2,3)) - 4*pFunction(pow(x+1,3)) + 6*pFunction(pow(x,3)) - 4*pFunction(pow(x-1,3)) );
}

float Scaling::pFunction(float x){
	return x > 0 ? x : 0;
} */


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
		outputFile = x_fopen(TEXT("img\\lena_simple.bmp"), TEXT("w"));
		break;
	case scalingTypes(bilinear):
		outputFile = x_fopen(TEXT("img\\lena_bilinear.bmp"), TEXT("w"));
		break;
	case scalingTypes(bicubic):
		outputFile = x_fopen(TEXT("img\\lena_bicubic.bmp"), TEXT("w"));
		break;
	case scalingTypes(lanczos):
		outputFile = x_fopen(TEXT("img\\lena_lanczos.bmp"), TEXT("w"));
		break;
	default:
		break;
	}

	return outputFile;
}