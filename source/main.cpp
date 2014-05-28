#include "xBitstream.h"
#include "xEntropy.h"
#include "xHuffman.h"
#include "xSequence.h"
#include "xCfg.h"
#include "xTransform.h"
#include "Scaling.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h> 


using namespace std;
using namespace AVlib;

void example_entropy();
void example_seq();
void example_bitmap();
void example_transform();
void example_huff();
void example_cfg();

void simpleResiezing(int width, int height);
void bilinearInterpolation(int width, int height);
int bilinear(float samplingPointX, float samplingPointY, int16** outputPel, int width, int height);

int wmain( int argc, wchar_t *argv[ ], wchar_t *envp[ ])
{
	int newWidth, newHeight, scalingType;
	char continueChar;

	do{
		cout<<"Podaj now¹ szerokoœæ obrazu"<<endl;
		cin>>newWidth;
		cout<<"Podaj now¹ wysokoœæ obrazu"<<endl;
		cin>>newHeight;
		cout<<"Wybierz rodzaj interpolacji:"<<endl<<"0 - szybka, 1 - dwuliniowa, 2 - dwuszeœcienna, 3 lanczosa"<<endl;
		cin>>scalingType;

		if(newWidth > 0 && newHeight > 0){
			Scaling* scaling = new Scaling(newWidth, newHeight);
			scaling->scaleImage(scalingType);
		}
		else{
			cout<<"Wysokoœæ i szerokoœæ musz¹ byæ dodatnie"<<endl;
		}

		cout<<"Aby kontynowaæ skalowanie wciœnij 'a'"<<endl;
		cin>>continueChar;

	} while (continueChar == 'a');






	//simpleResiezing(newWidth, newHeight);
	//bilinearInterpolation(newWidth, newHeight);


  system("pause");
  return EXIT_SUCCESS;
}

void Resize(int width, int height){
	FILE_HANDLE inputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538.bmp"), TEXT("r"));
	FILE_HANDLE outputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538_simple.bmp"), TEXT("w"));

	xImg* inputImage = new xImg;
	xImg* outputImage = new xImg;
	inputImage->getBMP(inputFile);
    
	int inputImageXSize = inputImage->m_Cmp[0]->m_SizeX;
	int inputImageYSize = inputImage->m_Cmp[0]->m_SizeY;

	int samplingPointX, samplingPointY;

	float ratioX = (inputImageXSize*1.0)/(width*1.0);
	float ratioY = (inputImageYSize*1.0)/(height*1.0);
	cout<<"ratioX: "<<ratioX<<", ratioY: "<<ratioY<<endl;

	outputImage->createResizedCompatibile(inputImage, width, height);

	xCmp** m_cmp = outputImage->m_Cmp;  // tablica wskaŸników z komponentami , np. R G B A

	xCmp* newPel = new xCmp[];
	for(int i=0; i<inputImage->m_CmpTotal; i++){
		xCmp* inputComponent = inputImage->m_Cmp[i];

		int16** inputPel = inputComponent->m_Pel;
		int16** outputPel = new	int16*[height]; 

		for(int j=0; j<height; j++){
			outputPel[j] = new int16[width];
			for(int k=0; k<width; k++){
				samplingPointX = k*ratioX;
				samplingPointY = j*ratioY;
				outputPel[j][k] = inputPel[samplingPointY][samplingPointX];
			}
			
		}
		outputImage->m_Cmp[i]->m_Pel = outputPel;
	}

	outputImage->putBMP(outputFile);
	

	inputImage->destroy();
	outputImage->destroy();
	delete inputImage;
	delete outputImage;
}

void simpleResiezing(int width, int height){
	FILE_HANDLE inputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538.bmp"), TEXT("r"));
	FILE_HANDLE outputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538_simple.bmp"), TEXT("w"));

	xImg* inputImage = new xImg;
	xImg* outputImage = new xImg;
	inputImage->getBMP(inputFile);
    
	int inputImageXSize = inputImage->m_Cmp[0]->m_SizeX;
	int inputImageYSize = inputImage->m_Cmp[0]->m_SizeY;

	int samplingPointX, samplingPointY;

	float ratioX = (inputImageXSize*1.0)/(width*1.0);
	float ratioY = (inputImageYSize*1.0)/(height*1.0);
	cout<<"ratioX: "<<ratioX<<", ratioY: "<<ratioY<<endl;

	outputImage->createResizedCompatibile(inputImage, width, height);

	xCmp** m_cmp = outputImage->m_Cmp;  // tablica wskaŸników z komponentami , np. R G B A

	for(int i=0; i<inputImage->m_CmpTotal; i++){
		xCmp* inputComponent = inputImage->m_Cmp[i];

		int16** inputPel = inputComponent->m_Pel;
		int16** outputPel = new	int16*[height]; 

		for(int j=0; j<height; j++){
			outputPel[j] = new int16[width];
			for(int k=0; k<width; k++){
				samplingPointX = k*ratioX;
				samplingPointY = j*ratioY;
				outputPel[j][k] = inputPel[samplingPointY][samplingPointX];
			}
			
		}
		outputImage->m_Cmp[i]->m_Pel = outputPel;
	}

	outputImage->putBMP(outputFile);
	

	inputImage->destroy();
	outputImage->destroy();
	delete inputImage;
	delete outputImage;
}


void bilinearInterpolation(int width, int height){
	FILE_HANDLE inputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538.bmp"), TEXT("r"));
	FILE_HANDLE outputFile = x_fopen(TEXT("C:\\Users\\Luiz\\Desktop\\samochody_1920x1080_0538_bilinear.bmp"), TEXT("w"));

	xImg* inputImage = new xImg;
	xImg* outputImage = new xImg;
	inputImage->getBMP(inputFile);

	float samplingPointX, samplingPointY;
    
	int inputImageXSize = inputImage->m_Cmp[0]->m_SizeX;
	int inputImageYSize = inputImage->m_Cmp[0]->m_SizeY;

	float ratioX = (inputImageXSize*1.0)/(width*1.0);
	float ratioY = (inputImageYSize*1.0)/(height*1.0);
	cout<<"ratioX: "<<ratioX<<", ratioY: "<<ratioY<<endl;

	outputImage->createResizedCompatibile(inputImage, width, height);

	xCmp** m_cmp = outputImage->m_Cmp;  // tablica wskaŸników z komponentami , np. R G B A

	xCmp* newPel = new xCmp[];
	for(int i=0; i<inputImage->m_CmpTotal; i++){
		xCmp* inputComponent = inputImage->m_Cmp[i];

		int16** inputPel = inputComponent->m_Pel;
		int16** outputPel = new	int16*[height]; 

		for(int j=0; j<height; j++){
			outputPel[j] = new int16[width];
			for(int k=0; k<width; k++){
					samplingPointX = k * ratioX;
					samplingPointY = j * ratioY;

					outputPel[j][k] = bilinear(samplingPointX, samplingPointY, inputPel, inputImageXSize, inputImageYSize);
			}	
		}
		// Zapisanie próbek komponentu
		outputImage->m_Cmp[i]->m_Pel = outputPel;
	}

	outputImage->putBMP(outputFile);
	

	inputImage->destroy();
	outputImage->destroy();
	delete inputImage;
	delete outputImage;
}

int bilinear(float samplingPointX, float samplingPointY, int16** inputPel, int width, int height){
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
	if(ceiledX >= width || ceiledY >= height){
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


void example_cfg()
{
  xCfgParser CfgParser;
  CfgParser.loadFromFile("example.cfg");
  string SampleField = CfgParser.getParamArg("Gamma_Selector", string("def"));
  CfgParser.storeToFile("example2.cfg");  
}

void example_huff()
{
  xSeq* InSeq = new xSeq;
  InSeq->createYUVRead(1920, 1088, 8, CrF_420, ClrSpc_BT709, TEXT("Poznan_Street_00_1920x1088_tex_cam03.yuv"), 8);
  InSeq->getFrame();

  xHuffman* Huff = new xHuffman;
  Huff->BuildCodeTable(InSeq->m_Img->m_Cmp[0]);
}

void example_transform()
{
  xSeq* InSeq = new xSeq;
  InSeq->createYUVRead(1920, 1088, 8, CrF_420, ClrSpc_BT709, TEXT("Poznan_Street_00_1920x1088_tex_cam03.yuv"), 8);
  xSeq* OutSeq = new xSeq;
  OutSeq->createYUVWrite(1920, 1088, 8, CrF_420, ClrSpc_BT709, TEXT("Poznan_Street_00_1920x1088_tex_cam03_b.yuv"), 8);

  uint64 Ticks = 0;

  for(int32 i=0; i<InSeq->m_NumOfFrames; i++)
  {    
    InSeq->getFrame();
    Ticks += OutSeq->m_Img->m_Cmp[0]->transform_4x4(InSeq->m_Img->m_Cmp[0]);
    OutSeq->m_Img->m_Cmp[1]->copy(InSeq->m_Img->m_Cmp[1]);
    OutSeq->m_Img->m_Cmp[2]->copy(InSeq->m_Img->m_Cmp[2]);
    OutSeq->putFrame();

    uint64 SAD = OutSeq->m_Img->m_Cmp[0]->calcSSD(InSeq->m_Img->m_Cmp[0]);

    fprintf(stdout, ".");
  }
  fprintf(stdout, "\n");

  InSeq->destroy();   delete InSeq;
  OutSeq->destroy();  delete OutSeq;

  fprintf(stdout, "\n10^9 TICKS: %f\n", (double)(Ticks/1000)/1000000.0);   
  fprintf(stdout, "END\n");
}

void example_bitmap()
{
  FILE_HANDLE InputFile = x_fopen(TEXT("samochody_1920x1080_0538.bmp"), TEXT("r"));
  FILE_HANDLE OutputFile = x_fopen(TEXT("samochody_1920x1080_0538_out.bmp"), TEXT("w"));

  xImg* InputImage = new xImg;
  xImg* OutputImage = new xImg;
  InputImage->getBMP(InputFile);
  OutputImage->createCompatible(InputImage);
  OutputImage->copy(InputImage);
  OutputImage->putBMP(OutputFile);

  InputImage->destroy();
  OutputImage->destroy();
  delete InputImage;
  delete OutputImage;

  x_fclose(InputFile);
  x_fclose(OutputFile);
}

void example_seq()
{
  {
    xSeq* InSeqYUV = new xSeq;
    InSeqYUV->createYUVRead(128, 128, 8, CrF_420, ClrSpc_BT709, TEXT("Crew_128x128.yuv"), 0);

    xSeq* OutSeqXSEQ = new xSeq;
    OutSeqXSEQ->createXSEQWrite(128, 128, 8, CrF_420, ClrSpc_BT709, TEXT("Crew_128x128_a.xseq"), SeqCmpr_LZ4, SeqPred_None, 0);

    uint64 Start = GetTickCount64();

    for(int32 i=0; i<InSeqYUV->m_NumOfFrames; i++)
    {
      fprintf(stdout, ".");
      InSeqYUV->getFrame();
      OutSeqXSEQ->m_Img->copy(InSeqYUV->m_Img);
      OutSeqXSEQ->putFrame();
    }
    fprintf(stdout, "\n");

    uint64 Time = GetTickCount64() - Start;
    double FPS  = 1000.0/((double)Time/(double)InSeqYUV->m_NumOfFrames);

    uint64 SizeYUV = x_fsize(InSeqYUV->m_File);
    uint64 SizeXSQE = x_fsize(OutSeqXSEQ->m_File);
    double Ratio = ((double)(SizeYUV-SizeXSQE)/(double)(SizeYUV))*100.0;
    
    fprintf(stdout, "Frames: %d, Time: %d ms, FPS: %.2f, Ratio: %.2f%%\n",InSeqYUV->m_NumOfFrames, Time, FPS, Ratio);

    InSeqYUV->destroy();   delete InSeqYUV;
    OutSeqXSEQ->destroy(); delete OutSeqXSEQ;
  }

  //-------------------------------------------------------------------------
  {
    xSeq* InSeqXSEQ = new xSeq;
    InSeqXSEQ->createXSEQRead(TEXT("Crew_128x128_a.xseq"), 0);

    xSeq* OutSeqYUV = new xSeq;
    OutSeqYUV->createYUVWrite(128, 128, 8, CrF_420, ClrSpc_BT709, TEXT("Crew_128x128_b.yuv"), 0);


    uint64 Start = GetTickCount64();

    for(int32 i=0; i<InSeqXSEQ->m_NumOfFrames; i++)
    {
      fprintf(stdout, ".");
      InSeqXSEQ->getFrame();
      OutSeqYUV->m_Img->copy(InSeqXSEQ->m_Img);
      OutSeqYUV->putFrame();
    }
    fprintf(stdout, "\n");

    uint64 Time = GetTickCount64() - Start;
    double FPS  = 1000.0/((double)Time/(double)InSeqXSEQ->m_NumOfFrames);
    uint64 SizeYUV = x_fsize(OutSeqYUV->m_File);
    uint64 SizeXSQE = x_fsize(InSeqXSEQ->m_File);
    double Ratio = ((double)(SizeYUV-SizeXSQE)/(double)(SizeYUV))*100.0;
    
    fprintf(stdout, "Frames: %d, Time: %d ms, FPS: %.2f, Ratio: %.2f%%\n",InSeqXSEQ->m_NumOfFrames, Time, FPS, Ratio);

    InSeqXSEQ->destroy(); delete InSeqXSEQ;
    OutSeqYUV->destroy(); delete OutSeqYUV;
  }
}

void example_entropy()
{
  //fstream output_file;
  //output_file.open(L"out.bin", ios_base::out | ios_base::binary);
  
  stringstream some_stream;

  xBitstreamWriter BitWriter;
  BitWriter.create(2048);
  BitWriter.init();
  xEntropy::WriteGolombRice(&BitWriter, 41, 1);
  BitWriter.writeAlignZero();
  BitWriter.flushToStream(&some_stream);

  some_stream.flush();
  some_stream.sync();
  
  xBitstreamReader BitReader;
  BitReader.create(2048);
  BitReader.init();
  BitReader.fetchFromStream(&some_stream);
  int32 abc = xEntropy::ReadGolombRice(&BitReader, 1);

  //output_file.flush();
  //output_file.close();
}