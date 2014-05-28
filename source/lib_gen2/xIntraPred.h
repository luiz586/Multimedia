#pragma once
#ifndef _xIntrPredH_
#define _xIntrPredH_

#include "xLicense.h"
#include "xCommon.h"

namespace AVlib {

class xIntraPred
{
protected:
  static const int32 m_AngTable[9];              //Intra prediction angles
  static const int32 m_InvAngTable[9];           //Inverse intra prediction angles: (256 * 32) / Angle

public:
  static void PredictDC_NxN_STD        (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, bool Filter); //any N
  static void PredictDC_NxN_SSE        (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, bool Filter); //N = 8*k, integer k>=1
  static void PredictDC_4x4_SSE        (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, bool Filter);
  static void PredictPlanar_NxN_STD    (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size); //any N  
  static void PredictPlanar_NxN_SSE    (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size); //N = 8*k, integer k>=1
  static void PredictPlanar_4x4_SSE    (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size); 
  static void PredictHorizontal_NxN_STD(int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter); //any N
  static void PredictHorizontal_NxN_SSE(int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter); //N = 8*k, integer k>=1
  static void PredictHorizontal_4x4_SSE(int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter);
  static void PredictVertical_NxN_STD  (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter); //any N
  static void PredictVertical_NxN_SSE  (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter); //N = 8*k, integer k>=1
  static void PredictVertical_4x4_SSE  (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int16 ClippingRange, bool Filter);
  static void PredictAngular_NxN_STD   (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int32 PredMode); //any N  
  static void PredictAngular_NxN_SSE   (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int32 PredMode); //N = 8*k, integer k>=1
  static void PredictAngular_4x4_SSE   (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int32 PredMode); //uses SSE transposition
  static void PredictAngular_8x8_SSE   (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int32 PredMode); //uses SSE transposition
  static void PredictAngular_16x16_SSE (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int32 PredMode); //uses SSE transposition
  static void PredictAngular_32x32_SSE (int16* Left, int16* Above, int16* Dst, int32 DStride, int32 Size, int32 PredMode); //uses SSE transposition, slower than universal SSE

//=========================================================================================================================================



};

} //end of namespace AVLib

#endif //_xIntrPredH_




