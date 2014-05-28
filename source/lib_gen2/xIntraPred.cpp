#include "xIntraPred.h"

namespace AVlib {

const int32 xIntraPred::m_AngTable[9]              = {0,    2,    5,   9,  13,  17,  21,  26,  32};
const int32 xIntraPred::m_InvAngTable[9]           = {0, 4096, 1638, 910, 630, 482, 390, 315, 256}; // (256 * 32) / Angle

void xIntraPred::PredictDC_NxN_STD(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, bool Filter)
{
  int32 DCVal = 0;
  for(int32 i=0; i<Size; i++) { DCVal += (Left[i] + Above[i]); }
  DCVal = (DCVal + Size) >> (xLog2(Size) + 1);

  if(Filter) //Filing samples with filtering
  {
    //top left sample
    Dst[0] = (Left[0] + 2*DCVal + Above[0] + 2) >> 2;
    //top line
    for (int32 x = 1; x < Size; x++) { Dst[x] = (3*DCVal + Above[x] + 2) >> 2; }
    Dst += DStride;
    //remaining lines
    for (int32 y = 1; y < Size; y++)
    {
      //left sample in line
      Dst[0] = (3*DCVal + Left[y] + 2) >> 2;
      //remaining samples
      for (int32 x = 1; x < Size; x++) { Dst[x] = DCVal; }
      Dst += DStride;
    }
  }
  else //Filing samples without filtering
  {    
    for (int32 y = 0; y < Size; y++)
    {
      for (int32 x = 0; x < Size; x++)
      {
        Dst[x] = DCVal;
      }
      Dst += DStride;
    }
  }
}
void xIntraPred::PredictDC_NxN_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, bool Filter)
{
  __m128i DC    = _mm_set1_epi16(0);
  __m128i Zero  = _mm_set1_epi16(0);
  __m128i Add   = _mm_set1_epi32(Size);
  int32   Shift = xLog2(Size)+1;

  //DC predictor evaluation  
  for(int32 i=0; i<Size; i+=8) { DC = _mm_add_epi16(DC, _mm_add_epi16(_mm_loadu_si128((__m128i*)(Left+i)), _mm_loadu_si128((__m128i*)(Above+i)))); }
  __m128i DC_line_0 = _mm_unpacklo_epi16(DC, Zero);
  __m128i	DC_line_1 = _mm_unpackhi_epi16(DC, Zero);
  DC = _mm_add_epi32(DC_line_0, DC_line_1);
  DC = _mm_hadd_epi32(DC, DC);
  DC = _mm_hadd_epi32(DC, DC);
  DC = _mm_add_epi32(DC, Add);
  DC = _mm_srai_epi32(DC, Shift);
  DC = _mm_packs_epi32(DC, DC);

  if(Filter) //Filing samples with filtering
  {
    __m128i C2 = _mm_set1_epi16(2);
    __m128i C3 = _mm_set1_epi16(3);
    __m128i C3DC = _mm_mullo_epi16(C3,DC);
    //top line - (3*DCVal + Above[x] + 2)/4
    for(int32 x=0; x<Size; x+=8) { _mm_storeu_si128((__m128i*)(Dst+x), _mm_srai_epi16(_mm_add_epi16(C3DC, _mm_add_epi16(_mm_loadu_si128((__m128i*)(Above+x)), C2)), 2)); }
    //top left sample
    Dst[0] = (Left[0] + 2*_mm_extract_epi16(DC,0) + Above[0] + 2) >> 2;
    Dst += DStride;
    
    //remaining lines
    //calculate left samples in lines: S[x]=(Left[x]+2+3*DC)/4
    __m128i S = _mm_srai_epi16(_mm_add_epi16(C3DC, _mm_add_epi16(_mm_loadu_si128((__m128i*)(Left)), C2)), 2);
    for (int32 y=1; y<8; y++)
    {
      S = _mm_srli_si128(S, 2);
      _mm_storeu_si128((__m128i*)(Dst), _mm_blend_epi16(DC, S, 0x1));
      for(int32 x=8; x<Size; x=x+8) { _mm_storeu_si128((__m128i*)(Dst+x),DC); }
      Dst += DStride;
    }

    for(int32 y=8; y<Size; y+=8)
    {
      //calculate left samples in lines: S[x]=(Left[x]+2+3*DC)/4
      S = _mm_srai_epi16(_mm_add_epi16(C3DC, _mm_add_epi16(_mm_loadu_si128((__m128i*)(Left+y)), C2)), 2);
      for (int32 yy=y; yy<y+8; yy++)
      {
        S = _mm_srli_si128(S, 2);
        _mm_storeu_si128((__m128i*)(Dst), _mm_blend_epi16(DC, S, 0x1));
        for(int32 x=8; x<Size; x=x+8) { _mm_storeu_si128((__m128i*)(Dst+x),DC); }
        Dst += DStride;
      }
    }
  }
  else //Filing samples without filtering
  {    
    for (int32 y = 0; y < Size; y++)
    {      
      for(int32 x=0; x<Size; x+=8) { _mm_storeu_si128((__m128i*)(Dst+x), DC); }  
      Dst += DStride;
    }
  }
}
void xIntraPred::PredictDC_4x4_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, bool Filter)
{
  //DC predictor evaluation
  __m128i C4 = _mm_set1_epi16(4);
  __m128i DC = _mm_add_epi16(_mm_loadl_epi64((__m128i*)Left), _mm_loadl_epi64((__m128i*)Above));
  DC = _mm_hadd_epi16(DC, DC);
  DC = _mm_hadd_epi16(DC, DC);
  DC = _mm_packs_epi32(DC, DC);
  DC = _mm_add_epi16(DC, C4);
  DC = _mm_srai_epi16(DC, 3);

  if(Filter) //Filing samples with filtering
  {
    __m128i C2   = _mm_set1_epi16(2);
    __m128i C3DC = _mm_mullo_epi16(DC, _mm_set1_epi16(3));
    //top line - (3*DCVal + Above[x] + 2)/4
    _mm_storel_epi64((__m128i*)(Dst), _mm_srai_epi16(_mm_add_epi16(C3DC, _mm_add_epi16(_mm_loadl_epi64((__m128i*)(Above)), C2)), 2));
    //top left sample
    Dst[0] = (Left[0] + 2*_mm_extract_epi16(DC, 0) + Above[0] + 2) >> 2;
    Dst += DStride;		

    //calculate left samples in lines: S[x]=(Left[x]+2+3*DC)/4    
    __m128i S = _mm_srai_epi16(_mm_add_epi16(C3DC, _mm_add_epi16(_mm_loadl_epi64((__m128i*)(Left)), C2)), 2);

    //remaining lines
    S = _mm_srli_si128(S, 2);
    _mm_storel_epi64((__m128i*)(Dst), _mm_blend_epi16(DC, S, 0x1));
    Dst += DStride;
    S = _mm_srli_si128(S, 2);
    _mm_storel_epi64((__m128i*)(Dst), _mm_blend_epi16(DC, S, 0x1));
    Dst += DStride;
    S = _mm_srli_si128(S, 2);
    _mm_storel_epi64((__m128i*)(Dst), _mm_blend_epi16(DC, S, 0x1));
  }
  else //Filing samples without filtering
  {    
    for (int32 y=0; y<4; y++)
    {
      _mm_storel_epi64((__m128i*)(Dst),DC);
      Dst += DStride;
    }
  } 
}

void xIntraPred::PredictPlanar_NxN_STD(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size)
{
  int32 BufferAbove[64];
  int32 BufferBelow[64];
  int32 BufferLeft [64];
  int32 BufferRight[64];

  int32  Offset2D = Size;
  int32  Shift1D  = xLog2(Size);
  int32  Shift2D  = Shift1D + 1;

  // Get left and above reference column and row
  // Prepare intermediate variables used in interpolation
  int32 BottomLeftSample = (int32)Left [Size];
  int32 TopRightSample   = (int32)Above[Size];
  for(int32 i=0; i<Size; i++)
  {
    int32 AboveSample = (int32)Above[i];
    int32 LeftSample  = (int32)Left [i];
    BufferAbove[i] = AboveSample<<Shift1D;
    BufferLeft [i] = LeftSample <<Shift1D;
    BufferBelow[i] = BottomLeftSample - AboveSample;
    BufferRight[i] = TopRightSample   - LeftSample;
  }

  // Generate prediction signal
  for(int32 y=0; y<Size; y++)
  {
    int32 HorPred = BufferLeft[y] + Offset2D;
    for(int32 x=0; x<Size; x++)
    {
      HorPred += BufferRight[y];
      BufferAbove[x] += BufferBelow[x];
      Dst[x] = (int16)( (HorPred + BufferAbove[x]) >> Shift2D );
    }
    Dst += DStride;
  }
}
void xIntraPred::PredictPlanar_NxN_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size)
{
  int32 BufferAbove[64];
  int32 BufferBelow[64];
  int32 BufferLeft [64];  
  int32 BufferRight[64];

  int32  Offset2D = Size;
  int32  Shift1D  = xLog2(Size);
  int32  Shift2D  = Shift1D + 1;

  __m128i Zero = _mm_setzero_si128();

  // Get left and above reference column and row
  // Prepare intermediate variables used in interpolation
  __m128i BottomLeftSampleV = _mm_set1_epi16((int32)Left [Size]);
  __m128i TopRightSampleV   = _mm_set1_epi16((int32)Above[Size]);

  // Prepare intermediate variables used in interpolation
  for(int32 i=0; i<Size; i+=8)
  {
    __m128i AboveV = _mm_loadu_si128((__m128i*)(Above+i));
    __m128i LeftV  = _mm_loadu_si128((__m128i*)(Left+i));
    __m128i BelowV = _mm_sub_epi16(BottomLeftSampleV, AboveV);
    __m128i Right  = _mm_sub_epi16(TopRightSampleV, LeftV);

    _mm_storeu_si128((__m128i*)(BufferAbove+i  ), _mm_slli_epi32(_mm_unpacklo_epi16(AboveV, Zero), Shift1D));
    _mm_storeu_si128((__m128i*)(BufferAbove+i+4), _mm_slli_epi32(_mm_unpackhi_epi16(AboveV, Zero), Shift1D));

    _mm_storeu_si128((__m128i*)(BufferLeft+i  ), _mm_slli_epi32(_mm_unpacklo_epi16(LeftV, Zero), Shift1D));
    _mm_storeu_si128((__m128i*)(BufferLeft+i+4), _mm_slli_epi32(_mm_unpackhi_epi16(LeftV, Zero), Shift1D));

    _mm_storeu_si128((__m128i*)(BufferBelow+i  ), _mm_cvtepi16_epi32(BelowV));
    _mm_storeu_si128((__m128i*)(BufferBelow+i+4), _mm_cvtepi16_epi32(_mm_srli_si128(BelowV, 8)));
    
    _mm_storeu_si128((__m128i*)(BufferRight+i  ),_mm_cvtepi16_epi32(Right));
    _mm_storeu_si128((__m128i*)(BufferRight+i+4),_mm_cvtepi16_epi32(_mm_srli_si128(Right, 8))); 
  }

  const __m128i Multiplier1V = _mm_setr_epi32(1, 2, 3, 4); 
  const __m128i Multiplier2V = _mm_setr_epi32(5, 6, 7, 8); 

  //Generate prediction signal
  for(int32 y=0; y<Size; y++)
  {
    __m128i HorPred1V = _mm_set1_epi32(BufferLeft[y] + Offset2D);
    __m128i HorPred2V = HorPred1V;
    __m128i RightV    = _mm_set1_epi32(BufferRight[y]);
    for(int32 x=0; x<Size; x+=8)
    { 
      __m128i Above1V = _mm_loadu_si128((__m128i*)&BufferAbove[x  ]);
      __m128i Above2V = _mm_loadu_si128((__m128i*)&BufferAbove[x+4]);
      __m128i Below1V = _mm_loadu_si128((__m128i*)&BufferBelow[x  ]);
      __m128i Below2V = _mm_loadu_si128((__m128i*)&BufferBelow[x+4]);
      Above1V = _mm_add_epi32(Above1V, Below1V);
      Above2V = _mm_add_epi32(Above2V, Below2V);
      _mm_storeu_si128((__m128i*)&BufferAbove[x  ], Above1V);
      _mm_storeu_si128((__m128i*)&BufferAbove[x+4], Above2V);
      HorPred1V = _mm_add_epi32(HorPred1V, _mm_mullo_epi32(Multiplier1V, RightV));
      HorPred2V = _mm_add_epi32(HorPred2V, _mm_mullo_epi32(Multiplier2V, RightV));
      __m128i Pred1V =  _mm_srai_epi32(_mm_add_epi32(HorPred1V, Above1V), Shift2D);
      __m128i Pred2V =  _mm_srai_epi32(_mm_add_epi32(HorPred2V, Above2V), Shift2D);
      __m128i Pred   = _mm_packs_epi32(Pred1V, Pred2V);
      _mm_storeu_si128((__m128i*)&Dst[x], Pred);
      HorPred1V = HorPred2V = _mm_shuffle_epi32(HorPred2V, 255);
    }
    Dst += DStride;
  }
}
void xIntraPred::PredictPlanar_4x4_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size)
{
  int32 BufferLeft [4];
  int32 BufferRight[4];

  int32  Offset2D = 4;
  __m128i Offset2DV = _mm_set1_epi32(Offset2D);
  int32  Shift1D  = 2;
  int32  Shift2D  = 3;

  // Get left and above reference column and row
  // Prepare intermediate variables used in interpolation
  __m128i BottomLeftSampleV = _mm_set1_epi32((int32)Left [4]);
  __m128i TopRightSampleV   = _mm_set1_epi32((int32)Above[4]);
  __m128i AboveV            = _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)(Above)));
  __m128i LeftV             = _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*)(Left)));
  __m128i Zero              = _mm_set1_epi16(0); 

  __m128i BufferBelowV = _mm_sub_epi32(BottomLeftSampleV, AboveV);
  _mm_storeu_si128((__m128i*)(BufferRight), _mm_sub_epi32(TopRightSampleV, LeftV));
  __m128i BufferAboveV = _mm_slli_epi32(AboveV, Shift1D);
  _mm_storeu_si128((__m128i*)(BufferLeft), _mm_add_epi32(_mm_slli_epi32(LeftV, Shift1D), Offset2DV));

  const __m128i MultiplierV = _mm_setr_epi32(1, 2, 3, 4); 

  // Generate prediction signal
  for (int32 y=0; y<4; y++)
  {
    __m128i HorPredV = _mm_set1_epi32(BufferLeft[y]);
    __m128i RightV   = _mm_set1_epi32(BufferRight[y]);
    BufferAboveV     = _mm_add_epi32(BufferAboveV, BufferBelowV);
    HorPredV         = _mm_add_epi32(HorPredV, _mm_mullo_epi32(MultiplierV, RightV));
    __m128i Pred     = _mm_packs_epi32(_mm_srai_epi32(_mm_add_epi32(HorPredV, BufferAboveV), Shift2D), Zero);
    _mm_storel_epi64((__m128i*)Dst, Pred);

    Dst += DStride;
  }
}

void xIntraPred::PredictHorizontal_NxN_STD(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter)
{
  if(Filter) //Filing samples with filtering
  {
    //top line
    int16 AboveLeftSample = Above[-1];
    for (int32 x=0; x<Size; x++) { Dst[x] = xClipU<int16>(Left[0] + ((Above[x] - AboveLeftSample)>>1), ClippingRange); } 
    Dst += DStride;
    //remaining lines
    for (int32 y=1; y<Size; y++)
    {
      int16 LeftSample = Left[y];
      for (int32 x = 0; x < Size; x++) { Dst[x] = LeftSample; }
      Dst += DStride;
    }   
  }
  else //Filing samples without filtering
  {    
    for(int32 y = 0; y < Size; y++)
    {
      int16 LeftSample = Left[y];
      for(int32 x=0; x<Size; x++) { Dst[x] = LeftSample; }
      Dst += DStride;
    }
  }
}
void xIntraPred::PredictHorizontal_NxN_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter)
{
  if(Filter)//Filing samples with filtering
  {
    //top line
    __m128i AboveLeftV     = _mm_set1_epi16(Above[-1]);
    __m128i LeftV          = _mm_set1_epi16(Left[0]);
    __m128i ClippingRangeV = _mm_set1_epi16(ClippingRange);
    __m128i ZeroV          = _mm_setzero_si128();
    for(int32 i=0; i<Size; i+=8)
    {
      __m128i AboveV = _mm_loadu_si128((__m128i*)(Above+i));
      __m128i Result  = _mm_max_epi16(ZeroV, _mm_min_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_subs_epi16(AboveV, AboveLeftV), 1), LeftV), ClippingRangeV));
      _mm_storeu_si128((__m128i*)(Dst+i), Result);
    }
    Dst += DStride;

    //remaining lines
    for (int32 y=1; y<Size; y++)
    {
      __m128i LeftV = _mm_set1_epi16(Left[y]);
      for(int32 x=0; x<Size; x+=8)
      {
        _mm_storeu_si128((__m128i*)(Dst+x), LeftV);
      }
      Dst += DStride;
    }   
  }
  else //Filing samples without filtering
  {    
    for (int32 y=0; y<Size; y++)
    {
      __m128i LeftV = _mm_set1_epi16(Left[y]);
      for(int32 x=0; x<Size; x+=8)
      {
        _mm_storeu_si128((__m128i*)(Dst+x), LeftV);
      }
      Dst += DStride;
    }
  }
}
void xIntraPred::PredictHorizontal_4x4_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter)
{
  if(Filter)//Filing samples with filtering
  {
    //top line
    __m128i AboveLeftV     = _mm_set1_epi16(Above[-1]);
    __m128i LeftV          = _mm_set1_epi16(Left[0]);
    __m128i ClippingRangeV = _mm_set1_epi16(ClippingRange);
    __m128i ZeroV          = _mm_setzero_si128();

    __m128i AboveV = _mm_loadl_epi64((__m128i*)Above);
    __m128i Result = _mm_max_epi16(ZeroV, _mm_min_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_subs_epi16(AboveV, AboveLeftV), 1), LeftV), ClippingRangeV));
    _mm_storel_epi64((__m128i*)Dst, Result);
    Dst += DStride;

    //remaining lines
    for (int32 y = 1; y < 4; y++)
    {
      __m128i LeftV = _mm_set1_epi16(Left[y]);
      _mm_storel_epi64((__m128i*)(Dst), LeftV);
      Dst += DStride;
    }   
  }
  else //Filing samples without filtering
  {    
    for (int32 y= 0; y<4; y++)
    {
      _mm_storel_epi64((__m128i*)(Dst),_mm_set1_epi16(Left[y]));
      Dst += DStride;
    }
  } 
}

void xIntraPred::PredictVertical_NxN_STD(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter)
{
  if(Filter) //Filing samples with filtering
  {
    int16 AboveLeftSample = Above[-1];

    for (int32 y=0; y<Size; y++)
    {          
      Dst[0] = xClipU<int16>(Above[0] + ((Left[y] - AboveLeftSample)>>1), ClippingRange);
      ::memcpy(Dst+1, Above+1, (Size-1)*sizeof(int16));
      Dst += DStride;
    }      
  }
  else //Filing samples without filtering
  {    
    for (int32 y=0; y<Size; y++)
    {
      ::memcpy(Dst, Above, Size*sizeof(int16));    
      Dst += DStride;
    }
  } 
}
void xIntraPred::PredictVertical_NxN_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter)
{
  if(Filter) //Filing samples with filtering (luma)
  {
    __m128i AboveLeftV = _mm_set1_epi16(Above[-1]);
    __m128i AboveSampleV = _mm_set1_epi16(Above[0]);
    __m128i ClippingRangeV = _mm_set1_epi16(ClippingRange);
    __m128i Zero = _mm_setzero_si128();     

    for(int32 y=0; y<Size; y+=8)
    {          
      __m128i LeftV = _mm_loadu_si128((__m128i*)(Left+y));
      __m128i Result = _mm_max_epi16(Zero, _mm_min_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_subs_epi16(LeftV, AboveLeftV), 1), AboveSampleV), ClippingRangeV));

      for(int32 yy=y; yy<y+8; yy++)
      {
        __m128i AboveV = _mm_loadu_si128((__m128i*)(Above));
        _mm_storeu_si128((__m128i*)(Dst), _mm_blend_epi16(AboveV, Result, 0x1));
        for(int32 x=8; x<Size; x+=8) { _mm_storeu_si128((__m128i*)(Dst+x), _mm_loadu_si128((__m128i*)(Above+x))); }
        Result = _mm_srli_si128(Result, 2);
        Dst += DStride;
      }
    }
  }
  else //Filing samples without filtering
  {    
    for(int32 y=0; y<Size; y++)
    {
      for (int32 x=0; x<Size; x+=8)
      {
        _mm_storeu_si128((__m128i*)(Dst+x),_mm_loadu_si128((__m128i*)(Above+x)));
      }
      Dst += DStride;
    }
  } 
}
void xIntraPred::PredictVertical_4x4_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter)
{
  if(Filter) //Filing samples with filtering (luma)
  {
    __m128i AboveLeftV = _mm_set1_epi16(Above[-1]);
    __m128i AboveSampleV = _mm_set1_epi16(Above[0]);
    __m128i ClippingRangeV = _mm_set1_epi16(ClippingRange);
    __m128i Zero = _mm_setzero_si128();     

    //calculate left samples
    __m128i LeftV  = _mm_loadl_epi64((__m128i*)Left );
    __m128i AboveV = _mm_loadl_epi64((__m128i*)Above); //1
    __m128i Result = _mm_max_epi16(Zero, _mm_min_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_subs_epi16(LeftV, AboveLeftV), 1), AboveSampleV), ClippingRangeV));
    
    //fill
    _mm_storel_epi64((__m128i*)Dst, _mm_blend_epi16(AboveV, Result, 0x1));
    Dst += DStride;
    _mm_storel_epi64((__m128i*)Dst, _mm_blend_epi16(AboveV, _mm_srli_si128(Result, 2), 0x1));
    Dst += DStride;
    _mm_storel_epi64((__m128i*)Dst, _mm_blend_epi16(AboveV, _mm_srli_si128(Result, 4), 0x1));
    Dst += DStride;
    _mm_storel_epi64((__m128i*)Dst, _mm_blend_epi16(AboveV, _mm_srli_si128(Result, 6), 0x1));    
  }
  else //Filing samples (chroma)
  {    
    for (int32 y=0; y<4; y++)
    {
      ::memcpy(Dst, Above, 4*sizeof(int16));
      Dst += DStride;
    }
  } 
}

void xIntraPred::PredictAngular_NxN_STD(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int32 PredMode)
{
  int16 RefArrTmp[130];

  // Map the mode index to main prediction direction and angle
  bool  ModeVer       = PredMode >= I_HOR_VER;
  int32 PredAngle     = ModeVer ? (int32)(PredMode - I_VER) : (int32)(I_HOR - PredMode);
  int32 AbsPredAngle  = xAbs(PredAngle);
  int32 SignPredAngle = PredAngle<0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32 InvPredAngle = m_InvAngTable[AbsPredAngle];
  AbsPredAngle       = m_AngTable[AbsPredAngle];
  PredAngle          = SignPredAngle * AbsPredAngle;

  // Initialise the reference array
  int16* RefArr;
  if (PredAngle < 0)
  {    
    int16* RefMain = (ModeVer ? Above : Left );
    int16* RefSide = (ModeVer ? Left  : Above);
    RefArr = RefArrTmp+Size+1;   
    ::memcpy(RefArr-1, RefMain-1, (Size+1)*sizeof(int16));

    int32 Range = ((int32)Size*PredAngle>>5)-1;
    int32 InvAngleSum  = 128; // rounding for (shift by 8)
    for (int32 i=-2; i>Range; i--) //max dist = 32 
    {
      InvAngleSum += InvPredAngle;
      int32 RefSideAddr = InvAngleSum>>8;
      RefArr[i] = RefSide[RefSideAddr-1];
    }    
  }
  else // PredAngle>=0
  {
    RefArr = (ModeVer ? Above : Left);
  }
    
  // Generate prediction signal
  int32 DeltaPos=0;
  if(ModeVer)
  {
    for (int32 y=0; y<Size; y++)
    {
      DeltaPos  += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        // Do linear filtering
        for (int32 x=0; x<Size; x++)
        {
          int32 RefMainIndex = x+DeltaInt;
          int32 A = RefArr[RefMainIndex  ];
          int32 B = RefArr[RefMainIndex+1];
          Dst[x] = (int16) (((DeltaFract*(B-A)+16)>>5)+A); //( (IvnDeltaFract*A + DeltaFract*B + 16) >> 5 )
        }
      }
      else
      {
        // Copy the integer samples
        ::memcpy(Dst, RefArr+DeltaInt, Size*sizeof(int16));
      }
      Dst += DStride;
    }
  }
  else //ModeHor
  {
    int32 DeltaPos=0;
    for (int32 x=0; x<Size; x++)
    {
      DeltaPos  += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        // Do linear filtering
        for (int32 y=0; y<Size; y++)
        {
          int32 RefMainIndex = y+DeltaInt;
          int32 A = RefArr[RefMainIndex  ];
          int32 B = RefArr[RefMainIndex+1];
          Dst[y*DStride+x]=(int16)(((DeltaFract*(B-A)+16)>>5)+A); //( (IvnDeltaFract*A + DeltaFract*B + 16) >> 5 )
        }
      }
      else
      {
        // Copy the integer samples
        for (int32 y=0; y<Size; y++)
        {
          Dst[y*DStride+x] = RefArr[y+DeltaInt];
        }
      }
    }
  } 
}
void xIntraPred::PredictAngular_NxN_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int32 PredMode)
{
  int16 RefArrTmp[140];  

  // Map the mode index to main prediction direction and angle
  bool  ModeVer       = PredMode >= I_HOR_VER;
  int32 PredAngle     = ModeVer ? (int32)(PredMode - I_VER) : (int32)(I_HOR - PredMode);
  int32 AbsPredAngle  = xAbs(PredAngle);
  int32 SignPredAngle = PredAngle<0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32 InvPredAngle = m_InvAngTable[AbsPredAngle];
  AbsPredAngle       = m_AngTable[AbsPredAngle];
  PredAngle          = SignPredAngle * AbsPredAngle;

  // Initialise the Main and Left reference array
  int16* RefArr;
  if(PredAngle<0)
  {      
    int16* RefMain = (ModeVer ? Above : Left );
    int16* RefSide = (ModeVer ? Left  : Above);
    RefArr = RefArrTmp+Size;
    ::memcpy(RefArr-1, RefMain-1, (Size+1)*sizeof(int16));

    // Extend the Main reference to the left.
    int32 Range = ((Size*PredAngle)>>5) - 1;
    int32 InvAngleSum = 128;       // rounding for (shift by 8)
    for (int32 i=-2; i>Range; i--)//Size*PredAngle>>5=8*PredAngle>>5=Predngle>>2
    {
      InvAngleSum += InvPredAngle;
      int32 RefSideAddr = InvAngleSum>>8;
      RefArr[i] = RefSide[RefSideAddr-1];
    }
  }
  else	// PredAngle>=0
  {
    RefArr = (ModeVer ? Above : Left); 
  }

  // Generate prediction signal
  int32 DeltaPos = 0;
  __m128i c16 = _mm_set1_epi16(16);
  if(ModeVer)
  {
    for (int32 k=0; k<Size; k++)
    {
      DeltaPos += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);

        for(int32 i=0; i<Size; i+=8)
        {
          __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
          __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i+1));
          __m128i PredV = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
          _mm_storeu_si128((__m128i*)(Dst+i), PredV);
        }
      }
      else
      {
        ::memcpy(Dst, RefArr+DeltaInt, Size*sizeof(int16));
      }
      Dst += DStride;
    }
  }
  else //(ModeHor)
  {
    int16 TmpCollumn[64];
    for (int32 k=0; k<Size; k++)
    {
      DeltaPos  += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);

        for(int32 i=0; i<Size; i+=8)
        {
          __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
          __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i+1));
          __m128i PredV = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
          _mm_storeu_si128((__m128i*)(TmpCollumn+i),PredV);
        }

        for (int l=0; l<Size; l++)
        {
          Dst[l*DStride+k]=TmpCollumn[l];
        }
      }
      else
      {
        for (int32 l=0; l<Size; l++)
        {
          Dst[l*DStride+k] = RefArr[l+DeltaInt];
        }
      }
    }
  }
}
void xIntraPred::PredictAngular_4x4_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int32 PredMode)
{
  int16 RefArrTmp[10];

  // Map the mode index to main prediction direction and angle
  bool  ModeVer       = PredMode >= I_HOR_VER;
  int32 PredAngle     = ModeVer ? (int32)(PredMode - I_VER) : (int32)(I_HOR - PredMode);
  int32 AbsPredAngle  = xAbs(PredAngle);
  int32 SignPredAngle = PredAngle<0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32 InvPredAngle = m_InvAngTable[AbsPredAngle];
  AbsPredAngle       = m_AngTable[AbsPredAngle];
  PredAngle          = SignPredAngle * AbsPredAngle;

  // Initialise the Main and Left reference array
  int16* RefArr;
  if(PredAngle<0)
  {
    int16* RefMain = (ModeVer ? Above : Left );
    int16* RefSide = (ModeVer ? Left  : Above);
    RefArr = RefArrTmp+4;    
    ::memcpy(RefArr-1, RefMain-1, (4+1)*sizeof(int16));

    // Extend the Main reference to the left.
    int32 Range = (PredAngle>>3)-1;
    int32 InvAngleSum = 128;     // rounding for (shift by 8)
    for (int i=-2; i>Range; i--) //Size*PredAngle>>5=4*PredAngle>>5=PredAngle>>3
    {
      InvAngleSum += InvPredAngle;
      int32 RefSideAddr = InvAngleSum>>8;
      RefArr[i] = RefSide[RefSideAddr-1];
    }
  }
  else	// PredAngle>=0
  {
    RefArr = (ModeVer ? Above : Left); 
  }

  // Generate prediction signal
  int32 DeltaPos=0;
  __m128i c16 = _mm_set1_epi16(16);
  if(ModeVer)
  {
    for(int32 k=0; k<2; k++)
    {
      DeltaPos += PredAngle;
      int32 DeltaInt1   = DeltaPos >> 5;
      int32 DeltaFract1 = DeltaPos & (32 - 1);

      DeltaPos += PredAngle;
      int32 DeltaInt2   = DeltaPos >> 5;
      int32 DeltaFract2 = DeltaPos & (32 -1);

      __m128i DeltaFracV=_mm_set_epi16(DeltaFract2, DeltaFract2, DeltaFract2, DeltaFract2, DeltaFract1, DeltaFract1, DeltaFract1, DeltaFract1);

      __m128i RefV0 = _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)(RefArr+DeltaInt1  )), _mm_loadl_epi64((__m128i*)(RefArr+DeltaInt2  ))); 
      __m128i RefV1 = _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)(RefArr+DeltaInt1+1)), _mm_loadl_epi64((__m128i*)(RefArr+DeltaInt2+1)));
      __m128i PredV = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFracV), c16), 5), RefV0);
      _mm_storel_epi64((__m128i*)Dst, PredV);
      Dst +=DStride;
      _mm_storel_epi64((__m128i*)Dst, _mm_srli_si128(PredV,8));
      Dst += DStride; 
    }
  }
  else //(ModeHor)
  {
    __m128i BuffPredV[2];
    for(int32 k=0; k<2; k++)
    {
      DeltaPos += PredAngle;
      int32 DeltaInt1   = DeltaPos >> 5;
      int32 DeltaFract1 = DeltaPos & (32 - 1);

      DeltaPos += PredAngle;
      int32 DeltaInt2   = DeltaPos >> 5;
      int32 DeltaFract2 = DeltaPos & (32 -1);

      __m128i DeltaFracV=_mm_set_epi16(DeltaFract2, DeltaFract2, DeltaFract2, DeltaFract2, DeltaFract1, DeltaFract1, DeltaFract1, DeltaFract1);

      __m128i RefV0 = _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)(RefArr+DeltaInt1  )), _mm_loadl_epi64((__m128i*)(RefArr+DeltaInt2  ))); 
      __m128i RefV1 = _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)(RefArr+DeltaInt1+1)), _mm_loadl_epi64((__m128i*)(RefArr+DeltaInt2+1)));
      BuffPredV[k] = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFracV), c16), 5), RefV0);
    }

    //transpose
    __m128i Transpose0 = _mm_unpacklo_epi16(BuffPredV[0], BuffPredV[1]);
    __m128i Transpose1 = _mm_unpackhi_epi16(BuffPredV[0], BuffPredV[1]);  
    BuffPredV[0] = _mm_unpacklo_epi16(Transpose0, Transpose1);
    BuffPredV[1] = _mm_unpackhi_epi16(Transpose0, Transpose1);

    //store in random area with DStride
    _mm_storel_epi64((__m128i*)Dst, BuffPredV[0]);
    BuffPredV[0] = _mm_srli_si128(BuffPredV[0], 8);
    Dst += DStride;
    _mm_storel_epi64((__m128i*)Dst, BuffPredV[0]);
    Dst += DStride;
    _mm_storel_epi64((__m128i*)Dst, BuffPredV[1]);
    BuffPredV[1] = _mm_srli_si128(BuffPredV[1], 8);
    Dst += DStride;
    _mm_storel_epi64((__m128i*)Dst, BuffPredV[1]); 
  }
}
void xIntraPred::PredictAngular_8x8_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int32 PredMode)
{
  int16 RefArrTmp[20];  

  // Map the mode index to main prediction direction and angle
  bool  ModeVer       = PredMode >= I_HOR_VER;
  int32 PredAngle     = ModeVer ? (int32)(PredMode - I_VER) : (int32)(I_HOR - PredMode);
  int32 AbsPredAngle  = xAbs(PredAngle);
  int32 SignPredAngle = PredAngle<0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32 InvPredAngle = m_InvAngTable[AbsPredAngle];
  AbsPredAngle       = m_AngTable[AbsPredAngle];
  PredAngle          = SignPredAngle * AbsPredAngle;

  // Initialise the Main and Left reference array.
  int16* RefArr;
  if(PredAngle<0)
  {      
    int16* RefMain = (ModeVer ? Above : Left );
    int16* RefSide = (ModeVer ? Left  : Above);
    RefArr = RefArrTmp+8;
    ::memcpy(RefArr-1, RefMain-1, (8+1)*sizeof(int16));

    // Extend the Main reference to the left.
    int32 Range = (PredAngle>>2)-1;
    int32 InvAngleSum = 128;       // rounding for (shift by 8)
    for (int32 i=-2; i>Range; i--)//Size*PredAngle>>5=8*PredAngle>>5=Predngle>>2
    {
      InvAngleSum += InvPredAngle;
      int32 RefSideAddr = InvAngleSum>>8;
      RefArr[i] = RefSide[RefSideAddr-1];
    }
  }
  else	// PredAngle>=0
  {
    RefArr = (ModeVer ? Above : Left); 
  }

  // Generate prediction signal
  int32 DeltaPos=0;
  __m128i c16=_mm_set1_epi16(16);
  if(ModeVer)
  {
    for (int32 k=0; k<8;k++)
    {
      DeltaPos += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);
        __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt  ));
        __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+1));
        __m128i PredV = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
        _mm_storeu_si128((__m128i*)(Dst), PredV);
      }
      else
      {
        _mm_storeu_si128((__m128i*)(Dst), _mm_loadu_si128((__m128i*)(RefArr+DeltaInt  )));
      }
      Dst += DStride;
    }
  }
  else //(ModeHor)
  {
    __m128i BuffPredV[8];
    for (int32 k=0; k<8; k++)
    {
      DeltaPos  += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);
        __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt  ));
        __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+1));
        BuffPredV[k] = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
      }
      else
      {
        BuffPredV[k] = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt  ));
      }
    }

    //transpose and write
    __m128i TransposeA0 = _mm_unpacklo_epi16(BuffPredV[0], BuffPredV[4]);
    __m128i TransposeA1 = _mm_unpackhi_epi16(BuffPredV[0], BuffPredV[4]);
    __m128i TransposeA2 = _mm_unpacklo_epi16(BuffPredV[1], BuffPredV[5]);
    __m128i TransposeA3 = _mm_unpackhi_epi16(BuffPredV[1], BuffPredV[5]);
    __m128i TransposeA4 = _mm_unpacklo_epi16(BuffPredV[2], BuffPredV[6]);
    __m128i TransposeA5 = _mm_unpackhi_epi16(BuffPredV[2], BuffPredV[6]);
    __m128i TransposeA6 = _mm_unpacklo_epi16(BuffPredV[3], BuffPredV[7]);
    __m128i TransposeA7 = _mm_unpackhi_epi16(BuffPredV[3], BuffPredV[7]);

    __m128i TransposeB0 = _mm_unpacklo_epi16(TransposeA0, TransposeA4);
    __m128i TransposeB1 = _mm_unpackhi_epi16(TransposeA0, TransposeA4);
    __m128i TransposeB2 = _mm_unpacklo_epi16(TransposeA1, TransposeA5);
    __m128i TransposeB3 = _mm_unpackhi_epi16(TransposeA1, TransposeA5);
    __m128i TransposeB4 = _mm_unpacklo_epi16(TransposeA2, TransposeA6);
    __m128i TransposeB5 = _mm_unpackhi_epi16(TransposeA2, TransposeA6);
    __m128i TransposeB6 = _mm_unpacklo_epi16(TransposeA3, TransposeA7);
    __m128i TransposeB7 = _mm_unpackhi_epi16(TransposeA3, TransposeA7);

    BuffPredV[0] = _mm_unpacklo_epi16(TransposeB0, TransposeB4);
    BuffPredV[1] = _mm_unpackhi_epi16(TransposeB0, TransposeB4);
    BuffPredV[2] = _mm_unpacklo_epi16(TransposeB1, TransposeB5);
    BuffPredV[3] = _mm_unpackhi_epi16(TransposeB1, TransposeB5);
    BuffPredV[4] = _mm_unpacklo_epi16(TransposeB2, TransposeB6);
    BuffPredV[5] = _mm_unpackhi_epi16(TransposeB2, TransposeB6);
    BuffPredV[6] = _mm_unpacklo_epi16(TransposeB3, TransposeB7);
    BuffPredV[7] = _mm_unpackhi_epi16(TransposeB3, TransposeB7);

    for (int32 k=0; k<8; k++)
    {
      _mm_storeu_si128((__m128i*)(Dst), BuffPredV[k]);
      Dst += DStride;
    }
  }
}
void xIntraPred::PredictAngular_16x16_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int32 PredMode)
{
  int16 RefArrTmp[40];  

  // Map the mode index to main prediction direction and angle
  bool  ModeVer       = PredMode >= I_HOR_VER;
  int32 PredAngle     = ModeVer ? (int32)(PredMode - I_VER) : (int32)(I_HOR - PredMode);
  int32 AbsPredAngle  = xAbs(PredAngle);
  int32 SignPredAngle = PredAngle<0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32 InvPredAngle = m_InvAngTable[AbsPredAngle];
  AbsPredAngle       = m_AngTable[AbsPredAngle];
  PredAngle          = SignPredAngle * AbsPredAngle;

  // Initialise the Main and Left reference array.
  int16* RefArr;
  if(PredAngle<0)
  {      
    int16* RefMain = (ModeVer ? Above : Left );
    int16* RefSide = (ModeVer ? Left  : Above);
    RefArr = RefArrTmp+16;
    ::memcpy(RefArr-1, RefMain-1, (16+1)*sizeof(int16));

    // Extend the Main reference to the left.
    int32 Range = (PredAngle>>1)-1;
    int32 InvAngleSum = 128;       // rounding for (shift by 8)
    for (int32 i=-2; i>Range; i--)//Size*PredAngle>>5=8*PredAngle>>5=Predngle>>2
    {
      InvAngleSum += InvPredAngle;
      int32 RefSideAddr = InvAngleSum>>8;
      RefArr[i] = RefSide[RefSideAddr-1];
    }
  }
  else	// PredAngle>=0
  {
    RefArr = (ModeVer ? Above : Left); 
  }

  // Generate prediction signal
  int32 DeltaPos=0;
  __m128i c16=_mm_set1_epi16(16);
  if(ModeVer)
  {
    for (int32 k=0; k<16 ;k++)
    {
      DeltaPos += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);

        for(int32 i=0; i<16; i+=8)
        {
          __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
          __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i+1));
          __m128i PredV = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
          _mm_storeu_si128((__m128i*)(Dst+i), PredV);
        }
      }
      else
      {
        ::memcpy(Dst, RefArr+DeltaInt, 32);
      }
      Dst += DStride;
    }
  }
  else //(ModeHor)
  {
    __m128i BuffPredV[32];
    __m128i BuffTrnsV[32];
    for (int32 k=0; k<32; k+=2)
    {
      DeltaPos  += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);

        for(int32 i=0; i<16; i+=8)
        {
          __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
          __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i+1));
          BuffPredV[k + (i>>3)] = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
        }
      }
      else
      {
        for(int32 i=0; i<16; i+=8)
        {
          BuffPredV[k + (i>>3)] = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
        }
      }
    }

    //transpose
    for(int32 i=0; i<16; i++)
    {
      BuffTrnsV[2*i  ] = _mm_unpacklo_epi16(BuffPredV[0+i], BuffPredV[16+i]);
      BuffTrnsV[2*i+1] = _mm_unpackhi_epi16(BuffPredV[0+i], BuffPredV[16+i]);
    }

    for(int32 i=0; i<16; i++)
    {
      BuffPredV[2*i  ] = _mm_unpacklo_epi16(BuffTrnsV[0+i], BuffTrnsV[16+i]);
      BuffPredV[2*i+1] = _mm_unpackhi_epi16(BuffTrnsV[0+i], BuffTrnsV[16+i]);
    }

    for(int32 i=0; i<16; i++)
    {
      BuffTrnsV[2*i  ] = _mm_unpacklo_epi16(BuffPredV[0+i], BuffPredV[16+i]);
      BuffTrnsV[2*i+1] = _mm_unpackhi_epi16(BuffPredV[0+i], BuffPredV[16+i]);
    }

    for(int32 i=0; i<16; i++)
    {
      BuffPredV[2*i  ] = _mm_unpacklo_epi16(BuffTrnsV[0+i], BuffTrnsV[16+i]);
      BuffPredV[2*i+1] = _mm_unpackhi_epi16(BuffTrnsV[0+i], BuffTrnsV[16+i]);
    }

    //store
    for(int32 j=0; j<32; j+=2)
    {
      _mm_storeu_si128((__m128i*)(Dst  ), BuffPredV[j  ]);
      _mm_storeu_si128((__m128i*)(Dst+8), BuffPredV[j+1]);
      Dst += DStride;
    }
  }
}
void xIntraPred::PredictAngular_32x32_SSE(int16* restrict Left, int16* restrict Above, int16* restrict Dst, int32 DStride, int32 Size, int32 PredMode)
{
  int16 RefArrTmp[80];  

  // Map the mode index to main prediction direction and angle
  bool  ModeVer       = PredMode >= I_HOR_VER;
  int32 PredAngle     = ModeVer ? (int32)(PredMode - I_VER) : (int32)(I_HOR - PredMode);
  int32 AbsPredAngle  = xAbs(PredAngle);
  int32 SignPredAngle = PredAngle<0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32 InvPredAngle = m_InvAngTable[AbsPredAngle];
  AbsPredAngle       = m_AngTable[AbsPredAngle];
  PredAngle          = SignPredAngle * AbsPredAngle;

  // Initialise the Main and Left reference array
  int16* RefArr;
  if(PredAngle<0)
  {      
    int16* RefMain = (ModeVer ? Above : Left );
    int16* RefSide = (ModeVer ? Left  : Above);
    RefArr = RefArrTmp+32;
    ::memcpy(RefArr-1, RefMain-1, (32+1)*sizeof(int16));

    // Extend the Main reference to the left.
    int32 Range = PredAngle-1;
    int32 InvAngleSum = 128;       // rounding for (shift by 8)
    for (int32 i=-2; i>Range; i--)//Size*PredAngle>>5=8*PredAngle>>5=Predngle>>2
    {
      InvAngleSum += InvPredAngle;
      int32 RefSideAddr = InvAngleSum>>8;
      RefArr[i] = RefSide[RefSideAddr-1];
    }
  }
  else	// PredAngle>=0
  {
    RefArr = (ModeVer ? Above : Left); 
  }

  // Generate prediction signal
  int32 DeltaPos = 0;
  __m128i c16 = _mm_set1_epi16(16);
  if(ModeVer)
  {
    for (int32 k=0; k<32; k++)
    {
      DeltaPos += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);

        for(int32 i=0; i<32; i+=8)
        {
          __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
          __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i+1));
          __m128i PredV = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
          _mm_storeu_si128((__m128i*)(Dst+i), PredV);
        }
      }
      else
      {
        ::memcpy(Dst, RefArr+DeltaInt, 64);
      }
      Dst += DStride;
    }
  }
  else //(ModeHor)
  //{
  //  int16 TmpCollumn[32];
  //  for (int32 k=0; k<32; k++)
  //  {
  //    DeltaPos  += PredAngle;
  //    int32 DeltaInt   = DeltaPos >> 5;
  //    int32 DeltaFract = DeltaPos & (32 - 1);

  //    if(DeltaFract)
  //    {
  //      __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);

  //      for(int32 i=0; i<32; i+=8)
  //      {
  //        __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
  //        __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i+1));
  //        __m128i PredV = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
  //        _mm_storeu_si128((__m128i*)(TmpCollumn+i),PredV);
  //      }

  //      for (int l=0; l<32; l++)
  //      {
  //        Dst[l*DStride+k]=TmpCollumn[l];
  //      }
  //    }
  //    else
  //    {
  //      for (int32 l=0; l<32; l++)
  //      {
  //        Dst[l*DStride+k] = RefArr[l+DeltaInt];
  //      }
  //    }
  //  }
  //}
  {
    __m128i BuffPredV[128];
    __m128i BuffTrnsV[128];
    for (int32 k=0; k<128; k+=4)
    {
      DeltaPos  += PredAngle;
      int32 DeltaInt   = DeltaPos >> 5;
      int32 DeltaFract = DeltaPos & (32 - 1);

      if(DeltaFract)
      {
        __m128i DeltaFractV = _mm_set1_epi16((int16)DeltaFract);

        for(int32 i=0; i<32; i+=8)
        {
          __m128i RefV0 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
          __m128i RefV1 = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i+1));
          BuffPredV[k + (i>>3)] = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(_mm_mullo_epi16(_mm_sub_epi16(RefV1, RefV0), DeltaFractV), c16), 5), RefV0);   // ((IvnDeltaFract*A + DeltaFract*B + 16) >> 5) -> ((deltaF*(R1-R)+16)>>5)+R 
        }
      }
      else
      {
        for(int32 i=0; i<32; i+=8)
        {
          BuffPredV[k + (i>>3)] = _mm_loadu_si128((__m128i*)(RefArr+DeltaInt+i  ));
        }
      }
    }

    //transpose
    for(int32 i=0; i<64; i++)
    {
      BuffTrnsV[2*i  ] = _mm_unpacklo_epi16(BuffPredV[0+i], BuffPredV[64+i]);
      BuffTrnsV[2*i+1] = _mm_unpackhi_epi16(BuffPredV[0+i], BuffPredV[64+i]);
    }

    for(int32 i=0; i<64; i++)
    {
      BuffPredV[2*i  ] = _mm_unpacklo_epi16(BuffTrnsV[0+i], BuffTrnsV[64+i]);
      BuffPredV[2*i+1] = _mm_unpackhi_epi16(BuffTrnsV[0+i], BuffTrnsV[64+i]);
    }

    for(int32 i=0; i<64; i++)
    {
      BuffTrnsV[2*i  ] = _mm_unpacklo_epi16(BuffPredV[0+i], BuffPredV[64+i]);
      BuffTrnsV[2*i+1] = _mm_unpackhi_epi16(BuffPredV[0+i], BuffPredV[64+i]);
    }

    for(int32 i=0; i<64; i++)
    {
      BuffPredV[2*i  ] = _mm_unpacklo_epi16(BuffTrnsV[0+i], BuffTrnsV[64+i]);
      BuffPredV[2*i+1] = _mm_unpackhi_epi16(BuffTrnsV[0+i], BuffTrnsV[64+i]);
    }

    for(int32 i=0; i<64; i++)
    {
      BuffTrnsV[2*i  ] = _mm_unpacklo_epi16(BuffPredV[0+i], BuffPredV[64+i]);
      BuffTrnsV[2*i+1] = _mm_unpackhi_epi16(BuffPredV[0+i], BuffPredV[64+i]);
    }

    //store
    for(int32 j=0; j<128; j+=4)
    {
      _mm_storeu_si128((__m128i*)(Dst   ), BuffPredV[j  ]);
      _mm_storeu_si128((__m128i*)(Dst+ 8), BuffPredV[j+1]);
      _mm_storeu_si128((__m128i*)(Dst+16), BuffPredV[j+2]);
      _mm_storeu_si128((__m128i*)(Dst+24), BuffPredV[j+3]);
      Dst += DStride;
    }
  }
}

} //end of namespace AVLib