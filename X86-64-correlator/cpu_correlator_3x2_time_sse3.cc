#include <cstdlib>
#include <iostream>
#include <pmmintrin.h> // sse3

#include "X86-64_correlator.h"


using namespace std;

static unsigned char cellToStatX[MAX_CELLS], cellToStatY[MAX_CELLS];

static unsigned fillCellToStatTable(const unsigned nrStations)
{
    unsigned nrCells, stat0, stat2;

    for (stat2 = nrStations % 2 ? 1 : 2, nrCells = 0; stat2 < nrStations; stat2 += 2) {
	for (stat0 = 0; stat0 + 2 <= stat2; stat0 += 3, nrCells ++) {
	    cellToStatX[nrCells] = stat0;
	    cellToStatY[nrCells] = stat2;
	}
    }

    return nrCells;
}

static unsigned long long calcNrOps(const unsigned nrTimes, const unsigned nrStations, const unsigned nrChannels,
				    unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned nrCells, stat0, stat2;

    for (stat2 = nrStations % 2 ? 1 : 2, nrCells = 0; stat2 < nrStations; stat2 += 2) {
	for (stat0 = 0; stat0 + 2 <= stat2; stat0 += 3, nrCells ++);
    }

    unsigned long long ops = nrChannels * nrCells * nrTimes * 32L * 6L;

    *bytesLoaded = nrChannels * nrCells * nrTimes * 4L * 4L * 5L;
    *bytesStored = nrChannels * nrCells * 2L * 4L * 5L;

    return ops;
}

unsigned long long cpuCorrelator_3x2_time_sse3(const float* __restrict__ samples, float* __restrict__ visibilities, 
						      const unsigned nrTimes, const unsigned nrTimesWidth, 
						      const unsigned nrStations, const unsigned nrChannels,
						      unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned nrCells = fillCellToStatTable(nrStations);

    for (unsigned channel = 0; channel < nrChannels; channel++) {
	for (unsigned cell = 0; cell < nrCells; cell++) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat3 = cellToStatY[cell];
	    unsigned index0 = SAMPLE_INDEX(stat0+0, channel, 0, 0, 0);
	    unsigned index1 = SAMPLE_INDEX(stat0+1, channel, 0, 0, 0);
	    unsigned index2 = SAMPLE_INDEX(stat0+2, channel, 0, 0, 0);
	    unsigned index3 = SAMPLE_INDEX(stat3+0, channel, 0, 0, 0);
	    unsigned index4 = SAMPLE_INDEX(stat3+1, channel, 0, 0, 0);

	    __m128 xxr03 = _mm_setzero_ps ();
	    __m128 xxi03 = _mm_setzero_ps ();
	    __m128 xyr03 = _mm_setzero_ps ();
	    __m128 xyi03 = _mm_setzero_ps ();
	    __m128 yxr03 = _mm_setzero_ps ();
	    __m128 yxi03 = _mm_setzero_ps ();
	    __m128 yyr03 = _mm_setzero_ps ();
	    __m128 yyi03 = _mm_setzero_ps ();

	    __m128 xxr13 = _mm_setzero_ps ();
	    __m128 xxi13 = _mm_setzero_ps ();
	    __m128 xyr13 = _mm_setzero_ps ();
	    __m128 xyi13 = _mm_setzero_ps ();
	    __m128 yxr13 = _mm_setzero_ps ();
	    __m128 yxi13 = _mm_setzero_ps ();
	    __m128 yyr13 = _mm_setzero_ps ();
	    __m128 yyi13 = _mm_setzero_ps ();

	    __m128 xxr23 = _mm_setzero_ps ();
	    __m128 xxi23 = _mm_setzero_ps ();
	    __m128 xyr23 = _mm_setzero_ps ();
	    __m128 xyi23 = _mm_setzero_ps ();
	    __m128 yxr23 = _mm_setzero_ps ();
	    __m128 yxi23 = _mm_setzero_ps ();
	    __m128 yyr23 = _mm_setzero_ps ();
	    __m128 yyi23 = _mm_setzero_ps ();


	    __m128 xxr04 = _mm_setzero_ps ();
	    __m128 xxi04 = _mm_setzero_ps ();
	    __m128 xyr04 = _mm_setzero_ps ();
	    __m128 xyi04 = _mm_setzero_ps ();
	    __m128 yxr04 = _mm_setzero_ps ();
	    __m128 yxi04 = _mm_setzero_ps ();
	    __m128 yyr04 = _mm_setzero_ps ();
	    __m128 yyi04 = _mm_setzero_ps ();

	    __m128 xxr14 = _mm_setzero_ps ();
	    __m128 xxi14 = _mm_setzero_ps ();
	    __m128 xyr14 = _mm_setzero_ps ();
	    __m128 xyi14 = _mm_setzero_ps ();
	    __m128 yxr14 = _mm_setzero_ps ();
	    __m128 yxi14 = _mm_setzero_ps ();
	    __m128 yyr14 = _mm_setzero_ps ();
	    __m128 yyi14 = _mm_setzero_ps ();

	    __m128 xxr24 = _mm_setzero_ps ();
	    __m128 xxi24 = _mm_setzero_ps ();
	    __m128 xyr24 = _mm_setzero_ps ();
	    __m128 xyi24 = _mm_setzero_ps ();
	    __m128 yxr24 = _mm_setzero_ps ();
	    __m128 yxi24 = _mm_setzero_ps ();
	    __m128 yyr24 = _mm_setzero_ps ();
	    __m128 yyi24 = _mm_setzero_ps ();




	    // assume nrTimes is divisable by 4
	    for (unsigned time = 0; time < nrTimes; time += 4) {
	      // This assumes a different memory layout with four samples in time behind each other:
	      // (xr1 xr2 xr3 xr4) (xi1 xi2 xi3 xi4) (yr1 yr2 yr3 yr4) (yi1 yi2 yi3 yi4), etc.
		__m128 sample0xr = _mm_load_ps(samples + time*4 + index0+0);
		__m128 sample0xi = _mm_load_ps(samples + time*4 + index0+4);
		__m128 sample0yr = _mm_load_ps(samples + time*4 + index0+8);
		__m128 sample0yi = _mm_load_ps(samples + time*4 + index0+12);

		__m128 sample1xr = _mm_load_ps(samples + time*4 + index1+0);
		__m128 sample1xi = _mm_load_ps(samples + time*4 + index1+4);
		__m128 sample1yr = _mm_load_ps(samples + time*4 + index1+8);
		__m128 sample1yi = _mm_load_ps(samples + time*4 + index1+12);

		__m128 sample2xr = _mm_load_ps(samples + time*4 + index2+0);
		__m128 sample2xi = _mm_load_ps(samples + time*4 + index2+4);
		__m128 sample2yr = _mm_load_ps(samples + time*4 + index2+8);
		__m128 sample2yi = _mm_load_ps(samples + time*4 + index2+12);

		__m128 sample3xr = _mm_load_ps(samples + time*4 + index3+0);
		__m128 sample3xi = _mm_load_ps(samples + time*4 + index3+4);
		__m128 sample3yr = _mm_load_ps(samples + time*4 + index3+8);
		__m128 sample3yi = _mm_load_ps(samples + time*4 + index3+12);

		__m128 sample4xr = _mm_load_ps(samples + time*4 + index4+0);
		__m128 sample4xi = _mm_load_ps(samples + time*4 + index4+4);
		__m128 sample4yr = _mm_load_ps(samples + time*4 + index4+8);
		__m128 sample4yi = _mm_load_ps(samples + time*4 + index4+12);

                xxr03 += sample0xr * sample3xr;
                xxi03 += sample0xi * sample3xr;
                xyr03 += sample0xr * sample3yr;
                xyi03 += sample0xi * sample3yr;
                yxr03 += sample0yr * sample3xr;
                yxi03 += sample0yi * sample3xr;
                yyr03 += sample0yr * sample3yr;
                yyi03 += sample0yi * sample3yr;
                xxr03 += sample0xi * sample3xi;
                xxi03 -= sample0xr * sample3xi;
                xyr03 += sample0xi * sample3yi;
                xyi03 -= sample0xr * sample3yi;
                yxr03 += sample0yi * sample3xi;
                yxi03 -= sample0yr * sample3xi;
                yyr03 += sample0yi * sample3yi;
                yyi03 -= sample0yr * sample3yi;

                xxr13 += sample1xr * sample3xr;
                xxi13 += sample1xi * sample3xr;
                xyr13 += sample1xr * sample3yr;
                xyi13 += sample1xi * sample3yr;
                yxr13 += sample1yr * sample3xr;
                yxi13 += sample1yi * sample3xr;
                yyr13 += sample1yr * sample3yr;
                yyi13 += sample1yi * sample3yr;
                xxr13 += sample1xi * sample3xi;
                xxi13 -= sample1xr * sample3xi;
                xyr13 += sample1xi * sample3yi;
                xyi13 -= sample1xr * sample3yi;
                yxr13 += sample1yi * sample3xi;
                yxi13 -= sample1yr * sample3xi;
                yyr13 += sample1yi * sample3yi;
                yyi13 -= sample1yr * sample3yi;

                xxr23 += sample2xr * sample3xr;
                xxi23 += sample2xi * sample3xr;
                xyr23 += sample2xr * sample3yr;
                xyi23 += sample2xi * sample3yr;
                yxr23 += sample2yr * sample3xr;
                yxi23 += sample2yi * sample3xr;
                yyr23 += sample2yr * sample3yr;
                yyi23 += sample2yi * sample3yr;
                xxr23 += sample2xi * sample3xi;
                xxi23 -= sample2xr * sample3xi;
                xyr23 += sample2xi * sample3yi;
                xyi23 -= sample2xr * sample3yi;
                yxr23 += sample2yi * sample3xi;
                yxi23 -= sample2yr * sample3xi;
                yyr23 += sample2yi * sample3yi;
                yyi23 -= sample2yr * sample3yi;


                xxr04 += sample0xr * sample4xr;
                xxi04 += sample0xi * sample4xr;
                xyr04 += sample0xr * sample4yr;
                xyi04 += sample0xi * sample4yr;
                yxr04 += sample0yr * sample4xr;
                yxi04 += sample0yi * sample4xr;
                yyr04 += sample0yr * sample4yr;
                yyi04 += sample0yi * sample4yr;
                xxr04 += sample0xi * sample4xi;
                xxi04 -= sample0xr * sample4xi;
                xyr04 += sample0xi * sample4yi;
                xyi04 -= sample0xr * sample4yi;
                yxr04 += sample0yi * sample4xi;
                yxi04 -= sample0yr * sample4xi;
                yyr04 += sample0yi * sample4yi;
                yyi04 -= sample0yr * sample4yi;

                xxr14 += sample1xr * sample4xr;
                xxi14 += sample1xi * sample4xr;
                xyr14 += sample1xr * sample4yr;
                xyi14 += sample1xi * sample4yr;
                yxr14 += sample1yr * sample4xr;
                yxi14 += sample1yi * sample4xr;
                yyr14 += sample1yr * sample4yr;
                yyi14 += sample1yi * sample4yr;
                xxr14 += sample1xi * sample4xi;
                xxi14 -= sample1xr * sample4xi;
                xyr14 += sample1xi * sample4yi;
                xyi14 -= sample1xr * sample4yi;
                yxr14 += sample1yi * sample4xi;
                yxi14 -= sample1yr * sample4xi;
                yyr14 += sample1yi * sample4yi;
                yyi14 -= sample1yr * sample4yi;

                xxr24 += sample2xr * sample4xr;
                xxi24 += sample2xi * sample4xr;
                xyr24 += sample2xr * sample4yr;
                xyi24 += sample2xi * sample4yr;
                yxr24 += sample2yr * sample4xr;
                yxi24 += sample2yi * sample4xr;
                yyr24 += sample2yr * sample4yr;
                yyi24 += sample2yi * sample4yr;
                xxr24 += sample2xi * sample4xi;
                xxi24 -= sample2xr * sample4xi;
                xyr24 += sample2xi * sample4yi;
                xyi24 -= sample2xr * sample4yi;
                yxr24 += sample2yi * sample4xi;
                yxi24 -= sample2yr * sample4xi;
                yyr24 += sample2yi * sample4yi;
                yyi24 -= sample2yr * sample4yi;
	    }

	    // store 8 floats per position, for a 3x2 cell that means 48 floats
	    unsigned vis_index = cell * 48;

	    // now, we have to sum the 4 values inside the regs.
	    __m128 tmp0, tmp1, xx_xy, tmp2, tmp3, yx_yy;

	    tmp0  = _mm_hadd_ps(xxr03, xxi03);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
	    tmp1  = _mm_hadd_ps(xyr03, xyi03);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
	    xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
	    _mm_store_ps(visibilities + vis_index + 0, xx_xy);
	    tmp2  = _mm_hadd_ps(yxr03, yxi03);
	    tmp3  = _mm_hadd_ps(yyr03, yyi03);
	    yx_yy = _mm_hadd_ps(tmp2, tmp3);
	    _mm_store_ps(visibilities + vis_index + 4, yx_yy);

	    tmp0  = _mm_hadd_ps(xxr13, xxi13);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
	    tmp1  = _mm_hadd_ps(xyr13, xyi13);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
	    xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
	    _mm_store_ps(visibilities + vis_index + 8, xx_xy);
	    tmp2  = _mm_hadd_ps(yxr13, yxi13);
	    tmp3  = _mm_hadd_ps(yyr13, yyi13);
	    yx_yy = _mm_hadd_ps(tmp2, tmp3);
	    _mm_store_ps(visibilities + vis_index + 12, yx_yy);

	    tmp0  = _mm_hadd_ps(xxr23, xxi23);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
	    tmp1  = _mm_hadd_ps(xyr23, xyi23);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
	    xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
	    _mm_store_ps(visibilities + vis_index + 16, xx_xy);
	    tmp2  = _mm_hadd_ps(yxr23, yxi23);
	    tmp3  = _mm_hadd_ps(yyr23, yyi23);
	    yx_yy = _mm_hadd_ps(tmp2, tmp3);
	    _mm_store_ps(visibilities + vis_index + 20, yx_yy);


	    tmp0  = _mm_hadd_ps(xxr04, xxi04);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
	    tmp1  = _mm_hadd_ps(xyr04, xyi04);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
	    xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
	    _mm_store_ps(visibilities + vis_index + 24, xx_xy);
	    tmp2  = _mm_hadd_ps(yxr04, yxi04);
	    tmp3  = _mm_hadd_ps(yyr04, yyi04);
	    yx_yy = _mm_hadd_ps(tmp2, tmp3);
	    _mm_store_ps(visibilities + vis_index + 28, yx_yy);

	    tmp0  = _mm_hadd_ps(xxr14, xxi14);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
	    tmp1  = _mm_hadd_ps(xyr14, xyi14);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
	    xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
	    _mm_store_ps(visibilities + vis_index + 32, xx_xy);
	    tmp2  = _mm_hadd_ps(yxr14, yxi14);
	    tmp3  = _mm_hadd_ps(yyr14, yyi14);
	    yx_yy = _mm_hadd_ps(tmp2, tmp3);
	    _mm_store_ps(visibilities + vis_index + 36, yx_yy);

	    tmp0  = _mm_hadd_ps(xxr24, xxi24);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
	    tmp1  = _mm_hadd_ps(xyr24, xyi24);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
	    xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
	    _mm_store_ps(visibilities + vis_index + 40, xx_xy);
	    tmp2  = _mm_hadd_ps(yxr24, yxi24);
	    tmp3  = _mm_hadd_ps(yyr24, yyi24);
	    yx_yy = _mm_hadd_ps(tmp2, tmp3);
	    _mm_store_ps(visibilities + vis_index + 44, yx_yy);
	}
    }

    return calcNrOps(nrTimes, nrStations, nrChannels, bytesLoaded, bytesStored);
}
