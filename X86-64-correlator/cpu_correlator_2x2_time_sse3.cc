#include <cstdlib>
#include <iostream>
#include <pmmintrin.h> // sse3

#include "X86-64_correlator.h"

using namespace std;

static unsigned char cellToStatX[MAX_CELLS], cellToStatY[MAX_CELLS];

static unsigned fillCellToStatTable(unsigned nrStations)
{
    unsigned nrCells, stat0, stat2;

    for (stat2 = nrStations % 2 ? 1 : 2, nrCells = 0; stat2 < nrStations; stat2 += 2) {
	for (stat0 = 0; stat0 + 2 <= stat2; stat0 += 2, nrCells ++) {
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
	for (stat0 = 0; stat0 + 2 <= stat2; stat0 += 2, nrCells ++);
    }

    unsigned long long ops = nrChannels * nrCells * nrTimes * 32L * 4L;

    *bytesLoaded = nrChannels * nrCells * nrTimes * 4L * 4L * 4L;
    *bytesStored = nrChannels * nrCells * 2L * 4L * 4L;

    return ops;
}

unsigned long long cpuCorrelator_2x2_time_sse3(const float* __restrict__ samples, float* __restrict__ visibilities, 
						      const unsigned nrTimes, const unsigned nrTimesWidth, 
						      const unsigned nrStations, const unsigned nrChannels,
						      unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    const unsigned nrBaselines = NR_BASELINES(nrStations);
    bool missedBaselines[nrBaselines*2]; // we sometimes go beyond the #baslines, but then we don't store the visibilities, so that is ok.
    for(unsigned baseline=0; baseline < nrBaselines; baseline++) {
	missedBaselines[baseline] = true;
    }

    const unsigned nrCells = fillCellToStatTable(nrStations);

    for (unsigned channel = 0; channel < nrChannels; channel++) {
	for (unsigned cell = 0; cell < nrCells; cell++) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat2 = cellToStatY[cell];

	    missedBaselines[BASELINE(stat0, stat2)] = false;
	    missedBaselines[BASELINE(stat0+1, stat2)] = false;
	    missedBaselines[BASELINE(stat0, stat2+1)] = false;
	    missedBaselines[BASELINE(stat0+1, stat2+1)] = false;

	    unsigned index0 = SAMPLE_INDEX(stat0+0, channel, 0, 0, 0);
	    unsigned index1 = SAMPLE_INDEX(stat0+1, channel, 0, 0, 0);
	    unsigned index2 = SAMPLE_INDEX(stat2+0, channel, 0, 0, 0);
	    unsigned index3 = SAMPLE_INDEX(stat2+1, channel, 0, 0, 0);

	    __m128 xxr02 = _mm_setzero_ps ();
	    __m128 xxi02 = _mm_setzero_ps ();
	    __m128 xyr02 = _mm_setzero_ps ();
	    __m128 xyi02 = _mm_setzero_ps ();
	    __m128 yxr02 = _mm_setzero_ps ();
	    __m128 yxi02 = _mm_setzero_ps ();
	    __m128 yyr02 = _mm_setzero_ps ();
	    __m128 yyi02 = _mm_setzero_ps ();

	    __m128 xxr12 = _mm_setzero_ps ();
	    __m128 xxi12 = _mm_setzero_ps ();
	    __m128 xyr12 = _mm_setzero_ps ();
	    __m128 xyi12 = _mm_setzero_ps ();
	    __m128 yxr12 = _mm_setzero_ps ();
	    __m128 yxi12 = _mm_setzero_ps ();
	    __m128 yyr12 = _mm_setzero_ps ();
	    __m128 yyi12 = _mm_setzero_ps ();

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

                xxr02 += sample0xr * sample2xr;
                xxi02 += sample0xi * sample2xr;
                xyr02 += sample0xr * sample2yr;
                xyi02 += sample0xi * sample2yr;
                yxr02 += sample0yr * sample2xr;
                yxi02 += sample0yi * sample2xr;
                yyr02 += sample0yr * sample2yr;
                yyi02 += sample0yi * sample2yr;
                xxr02 += sample0xi * sample2xi;
                xxi02 -= sample0xr * sample2xi;
                xyr02 += sample0xi * sample2yi;
                xyi02 -= sample0xr * sample2yi;
                yxr02 += sample0yi * sample2xi;
                yxi02 -= sample0yr * sample2xi;
                yyr02 += sample0yi * sample2yi;
                yyi02 -= sample0yr * sample2yi;

                xxr12 += sample1xr * sample2xr;
                xxi12 += sample1xi * sample2xr;
                xyr12 += sample1xr * sample2yr;
                xyi12 += sample1xi * sample2yr;
                yxr12 += sample1yr * sample2xr;
                yxi12 += sample1yi * sample2xr;
                yyr12 += sample1yr * sample2yr;
                yyi12 += sample1yi * sample2yr;
                xxr12 += sample1xi * sample2xi;
                xxi12 -= sample1xr * sample2xi;
                xyr12 += sample1xi * sample2yi;
                xyi12 -= sample1xr * sample2yi;
                yxr12 += sample1yi * sample2xi;
                yxi12 -= sample1yr * sample2xi;
                yyr12 += sample1yi * sample2yi;
                yyi12 -= sample1yr * sample2yi;

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
	    }

	    if (cell < nrCells) {
		// store 8 floats per position, for a 2x2 cell that means 32 floats
		unsigned vis_index = cell * 32;

		// now, we have to sum the 4 values inside the regs.
		__m128 tmp0, tmp1, xx_xy, tmp2, tmp3, yx_yy;

		tmp0  = _mm_hadd_ps(xxr02, xxi02);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
		tmp1  = _mm_hadd_ps(xyr02, xyi02);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
		xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
		_mm_store_ps(visibilities + vis_index + 0, xx_xy);
		
		tmp2  = _mm_hadd_ps(yxr02, yxi02);
		tmp3  = _mm_hadd_ps(yyr02, yyi02);
		yx_yy = _mm_hadd_ps(tmp2, tmp3);
		_mm_store_ps(visibilities + vis_index + 4, yx_yy);
		
		tmp0  = _mm_hadd_ps(xxr12, xxi12);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
		tmp1  = _mm_hadd_ps(xyr12, xyi12);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
		xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
		_mm_store_ps(visibilities + vis_index + 8, xx_xy);
		
		tmp2  = _mm_hadd_ps(yxr12, yxi12);
		tmp3  = _mm_hadd_ps(yyr12, yyi12);
		yx_yy = _mm_hadd_ps(tmp2, tmp3);
		_mm_store_ps(visibilities + vis_index + 12, yx_yy);
		
		
		tmp0  = _mm_hadd_ps(xxr03, xxi03);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
		tmp1  = _mm_hadd_ps(xyr03, xyi03);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
		xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
		_mm_store_ps(visibilities + vis_index + 16, xx_xy);
		
		tmp2  = _mm_hadd_ps(yxr03, yxi03);
		tmp3  = _mm_hadd_ps(yyr03, yyi03);
		yx_yy = _mm_hadd_ps(tmp2, tmp3);
		_mm_store_ps(visibilities + vis_index + 20, yx_yy);
		
		
		tmp0  = _mm_hadd_ps(xxr13, xxi13);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
		tmp1  = _mm_hadd_ps(xyr13, xyi13);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
		xx_xy = _mm_hadd_ps(tmp0, tmp1);     // xxr, xxi, xyr,xyi
		_mm_store_ps(visibilities + vis_index + 24, xx_xy);
		
		tmp2  = _mm_hadd_ps(yxr13, yxi13);
		tmp3  = _mm_hadd_ps(yyr13, yyi13);
		yx_yy = _mm_hadd_ps(tmp2, tmp3);
		_mm_store_ps(visibilities + vis_index + 28, yx_yy);
	    }
	}
    }
    
    unsigned long long missedBytesLoaded, missedBytesStored;
    unsigned long long missedOps = computeMissedBaselines(samples, visibilities, missedBaselines,
							  nrTimes, nrTimesWidth, nrStations, nrChannels,
							  &missedBytesLoaded, &missedBytesStored);

    unsigned long long ops = missedOps + calcNrOps(nrTimes, nrStations, nrChannels, bytesLoaded, bytesStored);
    *bytesLoaded += missedBytesLoaded;
    *bytesLoaded += missedBytesStored;
    return ops;
}
