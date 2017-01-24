#include "ipl.h"
#include "zbits.h"

// A class for typecasting ZBits images into Intel's IplImages
// This is trivially easy because the ZBits structure is built
// to be a near exact replica of IplImage.  Dividing this into an abstraction
// layer allows ZBits to be used without the Intel library for applications
// that don't need all the intel libs

struct Zipl {
	IplImage *iplImage;
		// This is the cached image. Cleaned up on destruct or clear()

	Zipl();
	Zipl( ZBits *zbits );
	Zipl( ZBits &zbits );
	~Zipl();
	operator IplImage *() { return iplImage; }
	virtual void set( ZBits *zbits );
	virtual void set( ZBits &zbits );
	virtual void clear();
	virtual void setROIxywh( int _x, int _y, int _w, int _h );
	virtual void setROIltrb( int _l, int _t, int _r, int _b );
};


////////////////////////////////////////////////////////////////////////////////////////////////////
// The following are all the IPL functions extracted from ipl.h just for quick reference
#if 0

void iplAllocateImage(IplImage* image, int doFill, int fillValue)
void iplAllocateImageFP(IplImage* image, int doFill, float fillValue)
IplImage* iplCreateImageJaehne ( int depth, int width, int height )
IplImage* iplCloneImage ( const IplImage* img ) 
void iplDeallocateHeader(IplImage* image)
void iplDeallocateImage(IplImage* image)
void iplDeallocate(IplImage* image, int flag)
IplROI *iplCreateROI(int coi,    int xOffset, int   yOffset, int width, int height )
void iplSetROI(IplROI*   roi,      int coi, int       xOffset,  int yOffset, int width,          int height)
void iplDeleteROI(IplROI* roi)
IplTileInfo* iplCreateTileInfo( IplCallBack  callBack, void* id, int width, int height )
void iplSetTileInfo ( IplTileInfo* tileInfo, IplCallBack  callBack, void* id, int width, int height )
void iplDeleteTileInfo (IplTileInfo* tileInfo)
IplImage* iplTranslateDIB(BITMAPINFOHEADER* dib, BOOL* cloneData)
void iplConvertFromDIB(BITMAPINFOHEADER* dib, IplImage* image)
void iplConvertToDIB(IplImage* image, BITMAPINFOHEADER* dib, int dither, int paletteConversion)
IPLStatus iplConvertFromDIBSep (BITMAPINFOHEADER* dibHeader, const char* dibData, IplImage* image) 
IPLStatus iplConvertToDIBSep (IplImage* image, BITMAPINFOHEADER* dibHeader, char* dibData, int dither, int paletteConversion) 

void iplCopy (IplImage* srcImage, IplImage* dstImage)
void iplExchange (IplImage* ImageA, IplImage* ImageB)
void iplSet (IplImage* image, int fillValue)
void iplSetFP (IplImage* image, float fillValue)
void iplPutPixel(IplImage* img, int x, int y, void* pixel)
void iplGetPixel(IplImage* img, int x, int y, void* pixel)
void iplConvert (IplImage *srcImage, IplImage *dstImage)

IPLStatus iplScale (const IplImage* srcImage, IplImage* dstImage )
IPLStatus iplScaleFP (const IplImage * srcImage,IplImage * dstImage, float  minVal, float  maxVal)
void iplAddS(IplImage* srcImage, IplImage* dstImage, int value)
void iplAddSFP(IplImage* srcImage, IplImage* dstImage, float value)
void iplSubtractS(IplImage* srcImage, IplImage* dstImage, int value, BOOL flip)
void iplSubtractSFP(IplImage* srcImage,IplImage* dstImage,float value, BOOL flip)
void iplMultiplyS(IplImage* srcImage, IplImage* dstImage, int value)
void iplMultiplySFP(IplImage* srcImage, IplImage* dstImage, float value)
void iplMultiplySScale(IplImage* srcImage, IplImage* dstImage, int value)
void iplAbs(IplImage* srcImage, IplImage* dstImage)
void iplSquare(IplImage* srcImage, IplImage* dstImage)
void iplAdd(IplImage* srcImageA, IplImage* srcImageB, IplImage* dstImage)
void iplSubtract(IplImage* srcImageA, IplImage* srcImageB, IplImage* dstImage)
void iplMultiply(IplImage* srcImageA, IplImage* srcImageB, IplImage* dstImage)
void iplMultiplyScale(IplImage* srcImageA, IplImage* srcImageB, IplImage* dstImage)
void iplNot(IplImage* srcImage, IplImage* dstImage)
void iplLShiftS(IplImage* srcImage, IplImage* dstImage, unsigned int nShift)
void iplRShiftS(IplImage* srcImage, IplImage* dstImage, unsigned int nShift)
void iplAndS(IplImage* srcImage, IplImage* dstImage, unsigned int value)
void iplOrS(IplImage* srcImage, IplImage* dstImage, unsigned int value)
void iplXorS(IplImage* srcImage, IplImage* dstImage, unsigned int value)
void iplAnd(IplImage* srcImageA, IplImage* srcImageB, IplImage* dstImage)
void iplOr(IplImage* srcImageA, IplImage* srcImageB, IplImage* dstImage)
void iplXor(IplImage* srcImageA, IplImage* srcImageB, IplImage* dstImage)
void iplAlphaComposite(IplImage* srcImageA, IplImage* srcImageB, IplImage* dstImage, int compositeType, IplImage* alphaImageA, IplImage* alphaImageB, IplImage* alphaImageDst,BOOL premulAlpha, BOOL divideMode)
void iplAlphaCompositeC(IplImage* srcImageA, IplImage* srcImageB,IplImage* dstImage, int compositeType, int aA, int aB,BOOL premulALpha, BOOL divideMode)
void iplPreMultiplyAlpha (IplImage* image, int alphaValue)

IplConvKernel* iplCreateConvKernel(int nCols, int nRows,int anchorX, int anchorY, int* values, int nShiftR)
IplConvKernelFP* iplCreateConvKernelFP(int nCols, int nRows,int anchorX, int anchorY, float* values)
IplConvKernel* iplCreateConvKernelChar(int nCols, int nRows,int anchorX, int anchorY, char* values, int nShiftR)
void iplGetConvKernel(IplConvKernel* kernel, int* nCols, int* nRows, int* anchorX, int* anchorY, int** values, int *nShiftR)
void iplGetConvKernelFP(IplConvKernelFP* kernel,int* nCols, int* nRows, int* anchorX, int* anchorY, float** values)
void iplGetConvKernelChar(IplConvKernel* kernel, int* nCols, int* nRows, int* anchorX, int* anchorY, char** values, int *nShiftR)
void iplDeleteConvKernel,(IplConvKernel* kernel)
void iplDeleteConvKernelFP(IplConvKernelFP* kernel)
void iplBlur(IplImage* srcImage, IplImage* dstImage, int nCols, int nRows, int anchorX, int anchorY)
void iplConvolve2D(IplImage* srcImage, IplImage* dstImage, IplConvKernel** kernel, int nKernels, int combineMethod)
void iplConvolve2DFP(IplImage* srcImage, IplImage* dstImage, IplConvKernelFP** kernel, int nKernels, int combineMethod)
IPLStatus iplFixedFilter(IplImage* srcImage, IplImage* dstImage, IplFilter filter)
void iplConvolveSep2D(IplImage* srcImage, IplImage* dstImage, IplConvKernel* xKernel, IplConvKernel* yKernel)
void iplConvolveSep2DFP(IplImage* srcImage, IplImage* dstImage, IplConvKernelFP* xKernel, IplConvKernelFP* yKernel)
void iplMedianFilter(IplImage* srcImage, IplImage* dstImage, int nCols, int nRows, int anchorX, int anchorY)
void iplColorMedianFilter(IplImage* srcImage, IplImage* dstImage,int nCols, int nRows, int anchorX, int anchorY)
void iplMaxFilter(IplImage* srcImage, IplImage* dstImage,int nCols, int nRows, int anchorX, int anchorY)
void iplMinFilter(IplImage* srcImage, IplImage* dstImage, int nCols,int nRows, int anchorX, int anchorY)
void iplRealFft2D(IplImage* srcImage, IplImage* dstImage, int flags)
void iplCcsFft2D(IplImage* srcImage, IplImage* dstImage, int flags)
void iplMpyRCPack2D(IplImage* srcA, IplImage* srcB, IplImage* dst)
void iplDCT2D(IplImage* src, IplImage* dst, int flags)
void iplErode(IplImage* srcImage, IplImage* dstImage,int nIterations)
void iplDilate(IplImage* srcImage, IplImage* dstImage,int nIterations)
void iplOpen(IplImage* srcImage, IplImage* dstImage,int nIterations)
void iplClose(IplImage* srcImage, IplImage* dstImage,int nIterations)

void iplReduceBits(IplImage* srcImage, IplImage* dstImage,int noise, int ditherType, int levels)
void iplBitonalToGray (IplImage *srcImage,IplImage *dstImage,int ZeroScale, int OneScale)
void iplColorToGray(IplImage* srcImage, IplImage* dstImage)
void iplGrayToColor (IplImage* srcImage, IplImage* dstImage,float FractR, float FractG, float FractB)
void iplRGB2HSV(IplImage* rgbImage, IplImage* hsvImage)
void iplHSV2RGB(IplImage* hsvImage, IplImage* rgbImage)
void iplRGB2HLS(IplImage* rgbImage, IplImage* hlsImage)
void iplHLS2RGB(IplImage* hlsImage, IplImage* rgbImage)
void iplRGB2XYZ(IplImage* srcImage, IplImage* dstImage)
void iplXYZ2RGB(IplImage* srcImage, IplImage* dstImage)
void iplRGB2LUV(IplImage* rgbImage, IplImage* LUVImage)
void iplLUV2RGB(IplImage* LUVImage, IplImage* rgbImage)
void iplRGB2YUV(IplImage* rgbImage, IplImage* YUVImage)
void iplYUV2RGB(IplImage* YUVImage, IplImage* rgbImage)
void iplRGB2YCrCb(IplImage* rgbImage, IplImage* YCrCbImage)
void iplYCrCb2RGB(IplImage* YCrCbImage, IplImage* rgbImage)
void iplYCC2RGB(IplImage* YCCImage, IplImage* rgbImage)
IplColorTwist* plCreateColorTwist(int data[16], int scalingValue)
void iplSetColorTwist (IplColorTwist *cTwist,int data[16],int scalingValue)
void iplDeleteColorTwist (IplColorTwist *cTwist)
void iplApplyColorTwist (IplImage* srcImage, IplImage* dstImage,IplColorTwist* cTwist, int offset)
IPLStatus iplColorTwistFP(const IplImage* srcImage, IplImage* dstImage, float * TwistFP)

void iplThreshold (IplImage* srcImage, IplImage* dstImage,int threshold)
void iplContrastStretch(IplImage* srcImage, IplImage* dstImage,IplLUT** lut)
void iplComputeHisto (IplImage* srcImage, IplLUT** lut)
void iplHistoEqualize(IplImage* srcImage, IplImage* dstImage,IplLUT** lut)

void iplWarpAffine (IplImage* srcImage, IplImage* dstImage, const double coeffs[2][3], int interpolate)
void iplRemap (IplImage* srcImage, IplImage* xMap, IplImage* yMap, IplImage* dstImage, int interpolate)
void iplShear (IplImage* srcImage, IplImage* dstImage,double xShear, double yShear, double xShift, double yShift, int interpolate)
void iplRotate (IplImage* srcImage, IplImage* dstImage, double angle,double xShift, double yShift, int interpolate)
void iplGetAffineQuad (IplImage* image, const double coeffs[2][3],double quad[4][2])
void iplGetAffineQuadROI (IplROI* roi, const double coeffs[2][3],double quad[4][2])
void iplGetAffineBound (IplImage* image, const double coeffs[2][3],double rect[2][2])
void iplGetAffineBoundROI (IplROI* roi, const double coeffs[2][3],double rect[2][2])
void iplGetAffineTransform (IplImage* image, double coeffs[2][3],const double quad[4][2])
void iplGetAffineTransformROI (IplROI* roi, double coeffs[2][3],const double quad[4][2])
void iplGetRotateShift (double xCenter, double yCenter, double angle,double *xShift, double *yShift)
void iplWarpBilinear  (IplImage* srcImage, IplImage* dstImage,const double coeffs[2][4], int warpFlag,int interpolate)
void iplWarpBilinearQ (IplImage* srcImage, IplImage* dstImage,const double quad[4][2], int warpFlag,int interpolate)
void iplWarpPerspective (IplImage* srcImage, IplImage* dstImage,const double coeffs[3][3], int warpFlag,int interpolate)
void iplWarpPerspectiveQ (IplImage* srcImage, IplImage* dstImage,const double quad[4][2], int warpFlag,int interpolate)
void iplGetBilinearQuad (IplImage* image, const double coeffs[2][4],double quadr[4][2])
void iplGetBilinearQuadROI (IplROI* roi, const double coeffs[2][4],double quadr[4][2])
void iplGetBilinearBound (IplImage* image, const double coeffs[2][4],double rect[2][2])
void iplGetBilinearBoundROI (IplROI* roi, const double coeffs[2][4],double rect[2][2])
void iplGetBilinearTransform (IplImage* image, double coeffs[2][4],const double quadr[4][2])
void iplGetBilinearTransformROI (IplROI* roi, double coeffs[2][4],const double quadr[4][2])
void iplGetPerspectiveQuad (IplImage* image, const double coeffs[3][3],double quadr[4][2])
void iplGetPerspectiveQuadROI (IplROI* roi, const double coeffs[3][3],double quadr[4][2])
void iplGetPerspectiveBound(IplImage* image, const double coeffs[3][3],double rect[2][2])
void iplGetPerspectiveBoundROI(IplROI* roi, const double coeffs[3][3],double rect[2][2])
void iplGetPerspectiveTransform (IplImage* image, double coeffs[3][3],const double quadr[4][2])
void iplGetPerspectiveTransformROI (IplROI* roi, double coeffs[3][3],const double quadr[4][2])
void iplResize (IplImage* srcImage, IplImage* dstImage,int xDst, int xSrc, int yDst, int ySrc, int interpolate)
void iplZoom(IplImage* srcImage, IplImage* dstImage, int xDst, int xSrc,int yDst, int ySrc, int interpolate)
void iplDecimate(IplImage* srcImage, IplImage* dstImage, int xDst,int xSrc, int yDst, int ySrc, int interpolate)
void iplDecimateBlur(IplImage* srcImage,IplImage* dstImage,int xDst,int xSrc, int yDst, int ySrc,int interpolate,int xMaskSize,int yMaskSize))
void iplMirror(IplImage* srcImage, IplImage* dstImage, int flipAxis)
double iplNorm(IplImage* srcImageA, IplImage* srcImageB, int normType)

void iplMoments(IplImage* img, IplMomentState stt)
double iplGetSpatialMoment(IplMomentState stt, int m, int n)
double iplGetNormalizedSpatialMoment(IplMomentState stt, int m, int n)
double iplGetCentralMoment(IplMomentState stt, int m, int n)
double iplGetNormalizedCentralMoment(IplMomentState stt, int m, int n)
double iplSpatialMoment(IplImage* img, int m, int n)
double iplNormalizedSpatialMoment(IplImage* img, int m, int n)
double iplCentralMoment(IplImage* img, int m, int n)
double iplNormalizedCentralMoment(IplImage* img, int m, int n)

IPLStatus iplGreater (IplImage* img1, IplImage* img2, IplImage* res)
IPLStatus iplLess    (IplImage* img1, IplImage* img2, IplImage* res)
IPLStatus iplEqual   (IplImage* img1, IplImage* img2, IplImage* res)
IPLStatus iplEqualFPEps (IplImage* img1, IplImage* img2, IplImage* res, float eps)
IPLStatus iplLessS(IplImage* img, int s, IplImage* res)
IPLStatus iplGreaterS(IplImage* img, int s, IplImage* res)
IPLStatus iplEqualS(IplImage* img, int s, IplImage* res)
IPLStatus iplLessSFP(IplImage* img, float s, IplImage* res)
IPLStatus iplGreaterSFP(IplImage* img, float s, IplImage* res)
IPLStatus iplEqualSFP(IplImage* img, float s, IplImage* res)
IPLStatus iplEqualSFPEps(IplImage* img, float s, IplImage* res, float eps)
IPLStatus iplMinMaxFP (const IplImage * srcImage,float * min, float * max)
IPLStatus iplNormCrossCorr(IplImage *tpl, IplImage *src, IplImage *dst)

IPLStatus iplWtInit (IplWtType type, int par1, int par2,IplWtKernel *wtKernel)
IPLStatus iplWtInitUserTaps (float *tap_filt[4],int len_filt[4], int ofs_filt[4], IplWtKernel *wtKernel)
IPLStatus iplWtInitUserFilter (const IplWtFilter *decLow,const IplWtFilter *decHigh, const IplWtFilter *recLow,const IplWtFilter *recHigh, IplWtKernel *wtKernel)
IPLStatus iplWtInitUserFilter4 (IplWtFilter filt[4],IplWtKernel *wtKernel)
IPLStatus iplWtDecompose (IplWtKernel *wtKernel, const IplImage *src,IplImage *approxDst, IplImage *xDetailDst, IplImage *yDetailDst,IplImage *dDetailDst))
IPLStatus iplWtReconstruct (IplWtKernel *wtKernel,const IplImage *approxSrc,  const IplImage *xDetailSrc,const IplImage *yDetailSrc, const IplImage *dDetailSrc, IplImage *dst)
IPLStatus iplWtFree (IplWtKernel *wtKernel)

void iplUserProcess(IplImage* srcImage, IplImage* dstImage, IplUserFunc cbFunc)
void iplUserProcessFP(IplImage* srcImage, IplImage* dstImage, IplUserFuncFP cbFunc)
void iplUserProcessPixel(IplImage* srcImage, IplImage* dstImage, IplUserFuncPixel cbFunc)
void iplNoiseUniformInit (IplNoiseParam* noiseParam, unsigned int seed, int low, int high)
void iplNoiseUniformInitFP (IplNoiseParam* noiseParam, unsigned int seed, float low, float high)
void iplNoiseGaussianInit (IplNoiseParam* noiseParam, unsigned int seed, int mean, int stDev)
void iplNoiseGaussianInitFP (IplNoiseParam* noiseParam, unsigned int seed, float mean, float stDev)
IPLStatus iplNoiseImage (IplImage* image, const IplNoiseParam* noiseParam)

#endif

