#include "xPixel.h"

namespace AVlib {

void xPixel::xAddContinuous_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  for(int32 i=0; i<Area; i++)
  { 
    Dst[i] = Src0[i] + Src1[i];   
  }
}
void xPixel::xAddContinuous_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  int32 Area8 = (int32)((uint32)Area & (uint32)0xFFFFFFF8);
  for(int32 i=0; i<Area8; i+=8)
  {
    __m128i Src0Block  = _mm_load_si128((__m128i*)&Src0[i]);
    __m128i Src1Block  = _mm_load_si128((__m128i*)&Src1[i]);
    __m128i Add        = _mm_adds_epi16(Src0Block, Src1Block); //
    _mm_store_si128((__m128i*)&Dst[i], Add);
  }
  for(int32 i=Area8; i<Area; i++)
  {
    Dst[i] = Src0[i] + Src1[i];
  }
}
void xPixel::xAddContinuous_AVX(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);
  for(int32 i=0; i<Area16; i+=16)
  {
    __m256i Src0Block  = _mm256_load_si256((__m256i*)&Src0[i]);
    __m256i Src1Block  = _mm256_load_si256((__m256i*)&Src1[i]);
    __m256i Add        = _mm256_adds_epi16(Src0Block, Src1Block); //
    _mm256_store_si256((__m256i*)&Dst[i], Add);
  }
  for(int32 i=Area16; i<Area; i++)
  {
    Dst[i] = Src0[i] + Src1[i];
  }
}
void xPixel::xAdd_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {      
      Dst[x] = Src0[x] + Src1[x];
    }
    Src0 += Src0Stride;
    Src1 += Src1Stride;
    Dst  += DstStride;  
  }
}
void xPixel::xAdd_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height)
{
  if(((uint32)Width & (uint32)0x7)==0) //Width%8==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=8)
      {
        __m128i Src0Block  = _mm_load_si128((__m128i*)&Src0[x]);
        __m128i Src1Block  = _mm_load_si128((__m128i*)&Src1[x]);
        __m128i Add        = _mm_adds_epi16(Src0Block, Src1Block); //
        _mm_store_si128((__m128i*)&Dst[x], Add);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
  else if(((uint32)Width & (uint32)0x3)==0) //Width%4==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=4)
      {
        __m128i Src0Block  = _mm_loadl_epi64((__m128i*)&Src0[x]);
        __m128i Src1Block  = _mm_loadl_epi64((__m128i*)&Src1[x]);
        __m128i Add        = _mm_adds_epi16(Src0Block, Src1Block); //
        _mm_storel_epi64((__m128i*)&Dst[x], Add);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
  else
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x++)
      {      
        Dst[x] = Src0[x] + Src1[x];
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
}

void xPixel::xAddAndClipContinuous_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area, int16 ClippingRange)
{
  for(int32 i=0; i<Area; i++)
  { 
    Dst[i] = (int16)xClipU<int32>((int32)Src0[i] + (int32)Src1[i], (int32)ClippingRange);   
  }
}
void xPixel::xAddAndClipContinuous_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area, int16 ClippingRange)
{
  const __m128i Zero   = _mm_setzero_si128();
  const __m128i MaxVal = _mm_set1_epi16(ClippingRange);
  int32 Area8 = (int32)((uint32)Area & (uint32)0xFFFFFFF8);

  for(int32 i=0; i<Area8; i+=8)
  {
    __m128i Src0Block  = _mm_load_si128((__m128i*)&Src0[i]);
    __m128i Src1Block  = _mm_load_si128((__m128i*)&Src1[i]);
    __m128i AddClip    = _mm_adds_epi16(Src0Block, Src1Block); //
    AddClip = _mm_min_epi16(AddClip, MaxVal);
    AddClip = _mm_max_epi16(AddClip, Zero);
    _mm_store_si128((__m128i*)&Dst[i], AddClip);
  }
  for(int32 i=Area8; i<Area; i+=8)
  {
    Dst[i] = (int16)xClipU<int32>((int32)Src0[i] + (int32)Src1[i], (int32)ClippingRange);
  }
}
void xPixel::xAddAndClipContinuous_AVX(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area, int16 ClippingRange)
{
  const __m256i Zero   = _mm256_setzero_si256();
  const __m256i MaxVal = _mm256_set1_epi16(ClippingRange);
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);

  for(int32 i=0; i<Area16; i+=16)
  {
    __m256i Src0Block  = _mm256_load_si256((__m256i*)&Src0[i]);
    __m256i Src1Block  = _mm256_load_si256((__m256i*)&Src1[i]);
    __m256i AddClip    = _mm256_adds_epi16(Src0Block, Src1Block); //
    AddClip = _mm256_min_epi16(AddClip, MaxVal);
    AddClip = _mm256_max_epi16(AddClip, Zero);
    _mm256_store_si256((__m256i*)&Dst[i], AddClip);
  }
  for(int32 i=Area16; i<Area; i+=8)
  {
    Dst[i] = (int16)xClipU<int32>((int32)Src0[i] + (int32)Src1[i], (int32)ClippingRange);
  }
}
void xPixel::xAddAndClip_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height, int16 ClippingRange)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {      
      Dst[x] = xClipU<int16>(Src0[x] + Src1[x], ClippingRange);
    }
    Src0 += Src0Stride;
    Src1 += Src1Stride;
    Dst  += DstStride;  
  }
}
void xPixel::xAddAndClip_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height, int16 ClippingRange)
{
  const __m128i Zero   = _mm_setzero_si128();
  const __m128i MaxVal = _mm_set1_epi16(ClippingRange);

  if(((uint32)Width & (uint32)0x7)==0) //Width%8==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=8)
      {
        __m128i Src0Block  = _mm_load_si128((__m128i*)&Src0[x]);
        __m128i Src1Block  = _mm_load_si128((__m128i*)&Src1[x]);
        __m128i AddClip    = _mm_adds_epi16(Src0Block, Src1Block); //
        AddClip = _mm_min_epi16(AddClip, MaxVal);
        AddClip = _mm_max_epi16(AddClip, Zero);
        _mm_store_si128((__m128i*)&Dst[x], AddClip);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
  else if(((uint32)Width & (uint32)0x3)==0) //Width%4==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=4)
      {
        __m128i Src0Block  = _mm_loadl_epi64((__m128i*)&Src0[x]);
        __m128i Src1Block  = _mm_loadl_epi64((__m128i*)&Src1[x]);
        __m128i AddClip    = _mm_adds_epi16(Src0Block, Src1Block); //
        AddClip = _mm_min_epi16(AddClip, MaxVal);
        AddClip = _mm_max_epi16(AddClip, Zero);
        _mm_storel_epi64((__m128i*)&Dst[x], AddClip);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
  else
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x++)
      {      
        Dst[x] = xClipU<int16>(Src0[x] + Src1[x], ClippingRange);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
}

void xPixel::xSubContinuous_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  for(int32 i=0; i<Area; i++)
  {
    Dst[i] = Src0[i] - Src1[i];
  }
}
void xPixel::xSubContinuous_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  int32 Area8 = (int32)((uint32)Area & (uint32)0xFFFFFFF8);

  for(int32 i=0; i<Area8; i+=8)
  {
    __m128i Src0Block  = _mm_load_si128((__m128i*)&Src0[i]);
    __m128i Src1Block  = _mm_load_si128((__m128i*)&Src1[i]);
    __m128i Sub        = _mm_subs_epi16(Src0Block, Src1Block);
    _mm_store_si128((__m128i*)&Dst[i], Sub);
  }
  for(int32 i=Area8; i<Area; i++)
  {
    Dst[i] = Src0[i] - Src1[i];
  }
}
void xPixel::xSubContinuous_AVX(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);

  for(int32 i=0; i<Area16; i+=16)
  {
    __m256i Src0Block  = _mm256_load_si256((__m256i*)&Src0[i]);
    __m256i Src1Block  = _mm256_load_si256((__m256i*)&Src1[i]);
    __m256i Sub        = _mm256_subs_epi16(Src0Block, Src1Block);
    _mm256_store_si256((__m256i*)&Dst[i], Sub);
  }
  for(int32 i=Area16; i<Area; i++)
  {
    Dst[i] = Src0[i] - Src1[i];
  }
}
void xPixel::xSub_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {      
      Dst[x] = Src0[x] - Src1[x];
    }
    Src0 += Src0Stride;
    Src1 += Src1Stride;
    Dst  += DstStride;  
  }
}
void xPixel::xSub_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height)
{
  if(((uint32)Width & (uint32)0x7)==0) //Width%8==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=8)
      {
        __m128i Src0Block  = _mm_load_si128((__m128i*)&Src0[x]);
        __m128i Src1Block  = _mm_load_si128((__m128i*)&Src1[x]);
        __m128i Sub       = _mm_subs_epi16(Src0Block, Src1Block);
        _mm_store_si128((__m128i*)&Dst[x], Sub);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
  else if(((uint32)Width & (uint32)0x3)==0) //Width%4==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=4)
      {
        __m128i Src0Block  = _mm_loadl_epi64((__m128i*)&Src0[x]);
        __m128i Src1Block  = _mm_loadl_epi64((__m128i*)&Src1[x]);
        __m128i Sub        = _mm_subs_epi16(Src0Block, Src1Block);
        _mm_storel_epi64((__m128i*)&Dst[x], Sub);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
  else
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x++)
      {      
        Dst[x] = Src0[x] - Src1[x];
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;  
    }
  }
}

void xPixel::xBiPredRefEstContinuous_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  for(int32 i=0; i<Area; i++)
  {
    Dst[i] = (Src0[i]<<1) - Src1[i];
  }
}
void xPixel::xBiPredRefEstContinuous_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  int32 Area8 = (int32)((uint32)Area & (uint32)0xFFFFFFF8);

  for(int32 i=0; i<Area8; i+=8)
  {
    __m128i Src0Block  = _mm_load_si128((__m128i*)&Src0[i]);
    __m128i Src1Block  = _mm_load_si128((__m128i*)&Src1[i]);
    Src0Block          = _mm_slli_epi16(Src0Block, 1);
    __m128i Sub        = _mm_subs_epi16(Src0Block, Src1Block);
    _mm_store_si128((__m128i*)&Dst[i], Sub);
  }
  for(int32 i=Area8; i<Area; i++)
  {
    Dst[i] = (Src0[i]<<1) - Src1[i];
  }
}
void xPixel::xBiPredRefEstContinuous_AVX(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);

  for(int32 i=0; i<Area16; i+=16)
  {
    __m256i Src0Block  = _mm256_load_si256((__m256i*)&Src0[i]);
    __m256i Src1Block  = _mm256_load_si256((__m256i*)&Src1[i]);
    Src0Block          = _mm256_slli_epi16(Src0Block, 1);
    __m256i Sub        = _mm256_subs_epi16(Src0Block, Src1Block);
    _mm256_store_si256((__m256i*)&Dst[i], Sub);
  }
  for(int32 i=Area16; i<Area; i++)
  {
    Dst[i] = (Src0[i]<<1) - Src1[i];
  }
}
void xPixel::xBiPredRefEst_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {
      Dst[x] = (Src0[x]<<1) - Src1[x];
    }
    Src0 += Src0Stride;
    Src1 += Src1Stride;
    Dst  += DstStride;    
  }
}
void xPixel::xBiPredRefEst_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height)
{
  if(((uint32)Width & (uint32)0x7)==0) //Width%8==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=8)
      {
        __m128i Src0Block  = _mm_load_si128((__m128i*)&Src0[x]);
        __m128i Src1Block  = _mm_load_si128((__m128i*)&Src1[x]);
        Src0Block          = _mm_slli_epi16(Src0Block, 1);
        __m128i Sub        = _mm_subs_epi16(Src0Block, Src1Block);
        _mm_store_si128((__m128i*)&Dst[x], Sub);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
  else if(((uint32)Width & (uint32)0x3)==0) //Width%4==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=4)
      {
        __m128i Src0Block  = _mm_loadl_epi64((__m128i*)&Src0[x]);
        __m128i Src1Block  = _mm_loadl_epi64((__m128i*)&Src1[x]);
        Src0Block          = _mm_slli_epi16(Src0Block, 1);
        __m128i Sub        = _mm_subs_epi16(Src0Block, Src1Block);
        _mm_storel_epi64((__m128i*)&Dst[x], Sub);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
  else
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x++)
      {
        Dst[x] = (Src0[x]<<1) - Src1[x];
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
}

void xPixel::xAvgContinuous_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  for(int32 i=0; i<Area; i++)
  { 
    Dst[i] = (Src0[i] + Src1[i] + 1)>>1;
  }
}
void xPixel::xAvgContinuous_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  const __m128i ZeroV   = _mm_setzero_si128();
  const __m128i OffsetV = _mm_set1_epi16(1);
  int32 Area8 = (int32)((uint32)Area & (uint32)0xFFFFFFF8);

  for(int32 i=0; i<Area8; i+=8)
  {
    __m128i Src0BlockV  = _mm_load_si128((__m128i*)&Src0[i]);
    __m128i Src1BlockV  = _mm_load_si128((__m128i*)&Src1[i]);
    __m128i AddV        = _mm_adds_epi16(_mm_adds_epi16(Src0BlockV, Src1BlockV), OffsetV);
    __m128i ShiftV      = _mm_srai_epi16(AddV, 1);
    _mm_store_si128((__m128i*)&Dst[i], ShiftV);
  }
  for(int32 i=Area8; i<Area; i++)
  { 
    Dst[i] = (Src0[i] + Src1[i] + 1)>>1;
  }
}
void xPixel::xAvgContinuous_AVX(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area)
{
  const __m256i ZeroV   = _mm256_setzero_si256();
  const __m256i OffsetV = _mm256_set1_epi16(1);
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);

  for(int32 i=0; i<Area16; i+=16)
  {
    __m256i Src0BlockV  = _mm256_load_si256((__m256i*)&Src0[i]);
    __m256i Src1BlockV  = _mm256_load_si256((__m256i*)&Src1[i]);
    __m256i AddV        = _mm256_adds_epi16(_mm256_adds_epi16(Src0BlockV, Src1BlockV), OffsetV);
    __m256i ShiftV      = _mm256_srai_epi16(AddV, 1);
    _mm256_store_si256((__m256i*)&Dst[i], ShiftV);
  }
  for(int32 i=Area16; i<Area; i++)
  { 
    Dst[i] = (Src0[i] + Src1[i] + 1)>>1;
  }
}
void xPixel::xAvg_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {
      Dst[x] = (Src0[x] + Src1[x] + 1)>>1;
    }
    Src0 += Src0Stride;
    Src1 += Src1Stride;
    Dst  += DstStride;    
  }
}
void xPixel::xAvg_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height)
{
  const __m128i ZeroV   = _mm_setzero_si128();
  const __m128i OffsetV = _mm_set1_epi16(1);

  if(((uint32)Width & (uint32)0x7)==0) //Width%8==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=8)
      {
        __m128i Src0BlockV  = _mm_load_si128((__m128i*)&Src0[x]);
        __m128i Src1BlockV  = _mm_load_si128((__m128i*)&Src1[x]);
        __m128i AddV        = _mm_adds_epi16(_mm_adds_epi16(Src0BlockV, Src1BlockV), OffsetV);
        __m128i ShiftV      = _mm_srai_epi16(AddV, 1);
        _mm_store_si128((__m128i*)&Dst[x], ShiftV);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
  else if(((uint32)Width & (uint32)0x3)==0) //Width%4==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=4)
      {
        __m128i Src0BlockV  = _mm_loadl_epi64((__m128i*)&Src0[x]);
        __m128i Src1BlockV  = _mm_loadl_epi64((__m128i*)&Src1[x]);
        __m128i AddV        = _mm_adds_epi16(_mm_adds_epi16(Src0BlockV, Src1BlockV), OffsetV);
        __m128i ShiftV      = _mm_srai_epi16(AddV, 1);
        _mm_storel_epi64((__m128i*)&Dst[x], ShiftV);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
  else
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x++)
      {
        Dst[x] = (Src0[x] + Src1[x] + 1)>>1;
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
}

void xPixel::xAvgAndClipContinuous_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area, int16 ClippingRange)
{
  for(int32 i=0; i<Area; i++)
  { 
    Dst[i] = xClipU<int16>(((Src0[i] + Src1[i] + 1)>>1), ClippingRange);
  }
}
void xPixel::xAvgAndClipContinuous_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area, int16 ClippingRange)
{
  const __m128i ZeroV   = _mm_setzero_si128();
  const __m128i OffsetV = _mm_set1_epi16(1);
  const __m128i MaxValV = _mm_set1_epi16(ClippingRange);
  int32 Area8 = (int32)((uint32)Area & (uint32)0xFFFFFFF8);

  for(int32 i=0; i<Area8; i+=8)
  {
    __m128i Src0BlockV  = _mm_load_si128((__m128i*)&Src0[i]);
    __m128i Src1BlockV  = _mm_load_si128((__m128i*)&Src1[i]);
    __m128i AddV        = _mm_adds_epi16(_mm_adds_epi16(Src0BlockV, Src1BlockV), OffsetV);
    __m128i ShiftV      = _mm_srai_epi16(AddV, 1);
    __m128i ClipV       = _mm_max_epi16(_mm_min_epi16(ShiftV, MaxValV),ZeroV);
    _mm_store_si128((__m128i*)&Dst[i], ClipV);
  }
  for(int32 i=Area8; i<Area; i++)
  { 
    Dst[i] = xClipU<int16>(((Src0[i] + Src1[i] + 1)>>1), ClippingRange);
  }
}
void xPixel::xAvgAndClipContinuous_AVX(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 Area, int16 ClippingRange)
{
  const __m256i ZeroV   = _mm256_setzero_si256();
  const __m256i OffsetV = _mm256_set1_epi16(1);
  const __m256i MaxValV = _mm256_set1_epi16(ClippingRange);
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);

  for(int32 i=0; i<Area16; i+=16)
  {
    __m256i Src0BlockV  = _mm256_load_si256((__m256i*)&Src0[i]);
    __m256i Src1BlockV  = _mm256_load_si256((__m256i*)&Src1[i]);
    __m256i AddV        = _mm256_adds_epi16(_mm256_adds_epi16(Src0BlockV, Src1BlockV), OffsetV);
    __m256i ShiftV      = _mm256_srai_epi16(AddV, 1);
    __m256i ClipV       = _mm256_max_epi16(_mm256_min_epi16(ShiftV, MaxValV),ZeroV);
    _mm256_store_si256((__m256i*)&Dst[i], ClipV);
  }
  for(int32 i=Area16; i<Area; i++)
  { 
    Dst[i] = xClipU<int16>(((Src0[i] + Src1[i] + 1)>>1), ClippingRange);
  }
}
void xPixel::xAvgAndClip_STD(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height, int16 ClippingRange)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {
      Dst[x] = xClipU<int16>(((Src0[x] + Src1[x] + 1)>>1), ClippingRange);
    }
    Src0 += Src0Stride;
    Src1 += Src1Stride;
    Dst  += DstStride;    
  }
}
void xPixel::xAvgAndClip_SSE(int16* restrict Dst, int16* restrict Src0, int16* restrict Src1, int32 DstStride, int32 Src0Stride, int32 Src1Stride, int32 Width, int32 Height, int16 ClippingRange)
{
  const __m128i ZeroV   = _mm_setzero_si128();
  const __m128i OffsetV = _mm_set1_epi16(1);
  const __m128i MaxValV = _mm_set1_epi16(ClippingRange);

  if(((uint32)Width & (uint32)0x7)==0) //Width%8==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=8)
      {
        __m128i Src0BlockV  = _mm_load_si128((__m128i*)&Src0[x]);
        __m128i Src1BlockV  = _mm_load_si128((__m128i*)&Src1[x]);
        __m128i AddV        = _mm_adds_epi16(_mm_adds_epi16(Src0BlockV, Src1BlockV), OffsetV);
        __m128i ShiftV      = _mm_srai_epi16(AddV, 1);
        __m128i ClipV       = _mm_max_epi16(_mm_min_epi16(ShiftV, MaxValV), ZeroV);
        _mm_store_si128((__m128i*)&Dst[x], ClipV);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
  else if(((uint32)Width & (uint32)0x3)==0) //Width%4==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=4)
      {
        __m128i Src0BlockV  = _mm_loadl_epi64((__m128i*)&Src0[x]);
        __m128i Src1BlockV  = _mm_loadl_epi64((__m128i*)&Src1[x]);
        __m128i AddV        = _mm_adds_epi16(_mm_adds_epi16(Src0BlockV, Src1BlockV), OffsetV);
        __m128i ShiftV      = _mm_srai_epi16(AddV, 1);
        __m128i ClipV       = _mm_max_epi16(_mm_min_epi16(ShiftV, MaxValV),ZeroV);
        _mm_storel_epi64((__m128i*)&Dst[x], ClipV);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
  else
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x++)
      {
        Dst[x] = xClipU<int16>(((Src0[x] + Src1[x] + 1)>>1), ClippingRange);
      }
      Src0 += Src0Stride;
      Src1 += Src1Stride;
      Dst  += DstStride;    
    }
  }
}

void xPixel::xClipContinuous_STD(int16* restrict Dst, int16* restrict Src, int32 Area, int16 ClippingRangeLow, int16 ClippingRangeHigh)
{
  for(int32 i=0; i<Area; i++)
  { 
    Dst[i] = xClip<int16>(Src[i], ClippingRangeLow, ClippingRangeHigh);
  }
}
void xPixel::xClipContinuous_SSE(int16* restrict Dst, int16* restrict Src, int32 Area, int16 ClippingRangeLow, int16 ClippingRangeHigh)
{
  const __m128i RangeLowV  = _mm_set1_epi16(ClippingRangeLow );
  const __m128i RangeHighV = _mm_set1_epi16(ClippingRangeHigh);
  int32 Area8 = (int32)((uint32)Area & (uint32)0xFFFFFFF8);

  for(int32 i=0; i<Area8; i+=8)
  {
    __m128i SrcBlockV  = _mm_load_si128((__m128i*)&Src[i]);
    __m128i ClipV      = _mm_max_epi16(_mm_min_epi16(SrcBlockV, RangeHighV), RangeLowV);
    _mm_store_si128((__m128i*)&Dst[i], ClipV);
  }
  for(int32 i=Area8; i<Area; i++)
  { 
    Dst[i] = xClip<int16>(Src[i], ClippingRangeLow, ClippingRangeHigh);
  }
}
void xPixel::xClipContinuous_AVX(int16* restrict Dst, int16* restrict Src, int32 Area, int16 ClippingRangeLow, int16 ClippingRangeHigh)
{
  const __m256i RangeLowV  = _mm256_set1_epi16(ClippingRangeLow );
  const __m256i RangeHighV = _mm256_set1_epi16(ClippingRangeHigh);
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);

  for(int32 i=0; i<Area16; i+=16)
  {
    __m256i SrcBlockV  = _mm256_load_si256((__m256i*)&Src[i]);
    __m256i ClipV      = _mm256_max_epi16(_mm256_min_epi16(SrcBlockV, RangeHighV), RangeLowV);
    _mm256_store_si256((__m256i*)&Dst[i], ClipV);
  }
  for(int32 i=Area16; i<Area; i++)
  { 
    Dst[i] = xClip<int16>(Src[i], ClippingRangeLow, ClippingRangeHigh);
  }
}
void xPixel::xClip_STD(int16* restrict Dst, int16* restrict Src, int32 DstStride, int32 SrcStride, int32 Width, int32 Height, int16 ClippingRangeLow, int16 ClippingRangeHigh)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {
      Dst[x] = xClip<int16>(Src[x], ClippingRangeLow, ClippingRangeHigh);
    }
    Src += SrcStride;
    Dst += DstStride;    
  }
}
void xPixel::xClip_SSE(int16* restrict Dst, int16* restrict Src, int32 DstStride, int32 SrcStride, int32 Width, int32 Height, int16 ClippingRangeLow, int16 ClippingRangeHigh)
{
  const __m128i RangeLowV  = _mm_set1_epi16(ClippingRangeLow );
  const __m128i RangeHighV = _mm_set1_epi16(ClippingRangeHigh);

  if(((uint32)Width & (uint32)0x7)==0) //Width%8==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=8)
      {
        __m128i SrcBlockV  = _mm_load_si128((__m128i*)&Src[x]);
        __m128i ClipV      = _mm_max_epi16(_mm_min_epi16(SrcBlockV, RangeHighV), RangeLowV);
        _mm_store_si128((__m128i*)&Dst[x], ClipV);
      }
      Src += SrcStride;
      Dst += DstStride;    
    }
  }
  else if(((uint32)Width & (uint32)0x3)==0) //Width%4==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=4)
      {
        __m128i SrcBlockV  = _mm_loadl_epi64((__m128i*)&Src[x]);
        __m128i ClipV      = _mm_max_epi16(_mm_min_epi16(SrcBlockV, RangeHighV), RangeLowV);
        _mm_storel_epi64((__m128i*)&Dst[x], ClipV);
      }
      Src += SrcStride;
      Dst += DstStride;    
    }
  }
  else
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x++)
      {
        Dst[x] = xClip<int16>(Src[x], ClippingRangeLow, ClippingRangeHigh);
      }
      Src += SrcStride;
      Dst += DstStride;    
    }
  }
}

void xPixel::xPackI16ToU8Continuous_STD(uint8* restrict Dst, int16* restrict Src, int32 Area)
{
  for(int32 i=0; i<Area; i++)
  { 
    Dst[i] = (uint8)xClipU8(Src[i]);   
  }
}
void xPixel::xPackI16ToU8Continuous_SSE(uint8* restrict Dst, int16* restrict Src, int32 Area)
{
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);
  for(int32 i=0; i<Area16; i+=16)
  { 
    __m128i In1 = _mm_load_si128((__m128i*)&(Src[i  ]));
    __m128i In2 = _mm_load_si128((__m128i*)&(Src[i+8]));
    __m128i Out = _mm_packus_epi16(In1, In2);
    _mm_storeu_si128 ((__m128i*)(&(Dst[i])), Out);
  }
  for(int32 i=Area16; i<Area; i++)
  {
    Dst[i] = (uint8)xClipU8(Src[i]);
  }
}
void xPixel::xPackI16ToU8Continuous_AVX(uint8* restrict Dst, int16* restrict Src, int32 Area)
{ 
  int32 Area32 = (int32)((uint32)Area & (uint32)0xFFFFFFE0);
  for(int32 i=0; i<Area32; i+=32)
  { 
    __m256i In1 = _mm256_load_si256((__m256i*)&(Src[i   ]));
    __m256i In2 = _mm256_load_si256((__m256i*)&(Src[i+16]));
    __m256i Out = _mm256_packus_epi16(In1, In2);
    Out = _mm256_permute4x64_epi64(Out, 0xD8);
    _mm256_store_si256 ((__m256i*)(&(Dst[i])), Out);
  }
  for(int32 i=Area32; i<Area; i++)
  {
    Dst[i] = (uint8)xClipU8(Src[i]);
  }
}
void xPixel::xPackI16ToU8_STD(uint8* restrict Dst, int16* restrict Src, int32 DstStride, int32 SrcStride, int32 Width, int32 Height)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {      
      Dst[x] = (uint8)xClipU8(Src[x]);
    }
    Src += SrcStride;
    Dst += DstStride;  
  }
}
void xPixel::xPackI16ToU8_SSE(uint8* restrict Dst, int16* restrict Src, int32 DstStride, int32 SrcStride, int32 Width, int32 Height)
{
  if(((uint32)Width & (uint32)0xF)==0) //Width%16==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=16)
      {
        __m128i In1 = _mm_load_si128((__m128i*)&(Src[x  ]));
        __m128i In2 = _mm_load_si128((__m128i*)&(Src[x+8]));
        __m128i Out = _mm_packus_epi16(In1, In2);
        _mm_storeu_si128 ((__m128i*)(&(Dst[x])), Out);
      }
      Src += SrcStride;
      Dst += DstStride;    
    }
  }
  else
  {
    int32 Width16 = (int32)((uint32)Width & (uint32)0xFFFFFFF0);
    int32 Width8  = (int32)((uint32)Width & (uint32)0xFFFFFFF8);
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width16; x+=16)
      {
        __m128i In1 = _mm_load_si128((__m128i*)&(Src[x  ]));
        __m128i In2 = _mm_load_si128((__m128i*)&(Src[x+8]));
        __m128i Out = _mm_packus_epi16(In1, In2);
        _mm_storeu_si128 ((__m128i*)(&(Dst[x])), Out);
      }
      for(int32 x=Width16; x<Width8; x+=8)
      {
        __m128i In = _mm_load_si128((__m128i*)&(Src[x  ]));
        __m128i Out = _mm_packus_epi16(In, In);
        _mm_storel_epi64 ((__m128i*)(&(Dst[x])), Out);
      }
      for(int32 x=Width8 ; x<Width; x++)
      {
        Dst[x] = (uint8)xClipU8(Src[x]);
      }
      Src += SrcStride;
      Dst += DstStride;    
    }
  }
}

void xPixel::xUnpackU8ToI16Continuous_STD(int16* restrict Dst, uint8* restrict Src, int32 Area)
{
  for(int32 i=0; i<Area; i++)
  { 
    Dst[i] = (int16)(Src[i]);   
  }
}
void xPixel::xUnpackU8ToI16Continuous_SSE(int16* restrict Dst, uint8* restrict Src, int32 Area)
{
  int32 Area16 = (int32)((uint32)Area & (uint32)0xFFFFFFF0);
  __m128i Zero = _mm_setzero_si128();

  for(int32 i=0; i<Area16; i+=16)
  { 
    __m128i In = _mm_loadu_si128((__m128i*)&(Src[i]));
    __m128i Out1 = _mm_unpacklo_epi8(In, Zero);
    __m128i Out2 = _mm_unpackhi_epi8(In, Zero);
    _mm_store_si128 ((__m128i*)(&(Dst[i  ])), Out1);
    _mm_store_si128 ((__m128i*)(&(Dst[i+8])), Out2);
  }
  for(int32 i=Area16; i<Area; i++)
  {
    Dst[i] = (int16)(Src[i]);
  }
}
void xPixel::xUnpackU8ToI16Continuous_AVX(int16* restrict Dst, uint8* restrict Src, int32 Area)
{
  int32 Area32 = (int32)((uint32)Area & (uint32)0xFFFFFFE0);
  __m256i Zero = _mm256_setzero_si256();
  
  for(int32 i=0; i<Area32; i+=32)
  { 
    __m256i In = _mm256_load_si256((__m256i*)&(Src[i]));
    In = _mm256_permute4x64_epi64(In, 0xD8);
    __m256i Out1 = _mm256_unpacklo_epi8(In, Zero);
    __m256i Out2 = _mm256_unpackhi_epi8(In, Zero);
    _mm256_store_si256 ((__m256i*)(&(Dst[i   ])), Out1);
    _mm256_store_si256 ((__m256i*)(&(Dst[i+16])), Out2);
  }
  for(int32 i=Area32; i<Area; i++)
  {
    Dst[i] = (int16)(Src[i]);
  }
}
void xPixel::xUnpackU8ToI16_STD(int16* restrict Dst, uint8* restrict Src, int32 DstStride, int32 SrcStride, int32 Width, int32 Height)
{
  for(int32 y=0; y<Height; y++)
  {
    for(int32 x=0; x<Width; x++)
    {      
      Dst[x] = (int16)(Src[x]);
    }
    Src += SrcStride;
    Dst += DstStride;  
  }
}
void xPixel::xUnpackU8ToI16_SSE(int16* restrict Dst, uint8* restrict Src, int32 DstStride, int32 SrcStride, int32 Width, int32 Height)
{
  __m128i Zero = _mm_setzero_si128();

  if(((uint32)Width & (uint32)0xF)==0) //Width%16==0
  {
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width; x+=16)
      {
        __m128i In = _mm_loadu_si128((__m128i*)&(Src[x]));
        __m128i Out1 = _mm_unpacklo_epi8(In, Zero);
        __m128i Out2 = _mm_unpackhi_epi8(In, Zero);
        _mm_store_si128 ((__m128i*)(&(Dst[x  ])), Out1);
        _mm_store_si128 ((__m128i*)(&(Dst[x+8])), Out2);
      }
      Src += SrcStride;
      Dst += DstStride;    
    }
  }
  else
  {
    int32 Width16 = (int32)((uint32)Width & (uint32)0xFFFFFFF0);
    int32 Width8  = (int32)((uint32)Width & (uint32)0xFFFFFFF8);
    for(int32 y=0; y<Height; y++)
    {
      for(int32 x=0; x<Width16; x+=16)
      {
        __m128i In = _mm_loadu_si128((__m128i*)&(Src[x]));
        __m128i Out1 = _mm_unpacklo_epi8(In, Zero);
        __m128i Out2 = _mm_unpackhi_epi8(In, Zero);
        _mm_store_si128 ((__m128i*)(&(Dst[x  ])), Out1);
        _mm_store_si128 ((__m128i*)(&(Dst[x+8])), Out2);
      }
      for(int32 x=Width16; x<Width8; x+=8)
      {
        __m128i In = _mm_loadl_epi64((__m128i*)&(Src[x]));
        __m128i Out1 = _mm_unpacklo_epi8(In, Zero);
        _mm_store_si128 ((__m128i*)(&(Dst[x  ])), Out1);
      }
      for(int32 x=Width8 ; x<Width; x++)
      {
        Dst[x] = (int16)(Src[x]);
      }
      Src += SrcStride;
      Dst += DstStride;    
    }
  }
}


} //end of namespace AVLib


