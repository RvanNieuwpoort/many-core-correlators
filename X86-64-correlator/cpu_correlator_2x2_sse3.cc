#include <cstdlib>
#include <iostream>
#include <pmmintrin.h> // sse3

#include "X86-64_correlator.h"


using namespace std;

static unsigned char cellToStatX[MAX_NR_CELLS], cellToStatY[MAX_NR_CELLS];

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

static unsigned long long calcNrOps(unsigned nrCells, unsigned nrTimes, unsigned nrStations, unsigned nrChannels,
				    unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned long long ops = (unsigned long long)nrChannels * nrCells * nrTimes * 32L * 4L;
    *bytesLoaded = (unsigned long long)nrChannels * nrCells * nrTimes * 4L * 4L * 4L;
    *bytesStored = (unsigned long long)nrChannels * nrCells * 2L * 4L * 4L;
    return ops;
}

static inline void store(float* visibilities, unsigned stat0, unsigned stat1, unsigned channel, unsigned nrChannels, __m128 real, __m128 imag)
{
    unsigned baseline = BASELINE(stat0, stat1);
    size_t vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);

    __m128 vis1tmp = _mm_shuffle_ps(real, imag, _MM_SHUFFLE(3, 2, 3, 2));
    __m128 vis1    = _mm_shuffle_ps(vis1tmp, vis1tmp, _MM_SHUFFLE(2, 0, 3, 1));
    _mm_store_ps(&(visibilities[vis_index+0]), vis1);

    __m128 vis2tmp = _mm_shuffle_ps(real, imag, _MM_SHUFFLE(1, 0, 1, 0));
    __m128 vis2    = _mm_shuffle_ps(vis2tmp, vis2tmp, _MM_SHUFFLE(2, 0, 3, 1));
    _mm_store_ps(&(visibilities[vis_index+4]), vis2);
}

unsigned long long cpuCorrelator_2x2_sse3(const float* __restrict__ samples, float* __restrict__ visibilities,
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

	    size_t index0 = SAMPLE_INDEX(stat0, channel, 0, 0, 0);
	    size_t index2 = SAMPLE_INDEX(stat2, channel, 0, 0, 0);

	    __m128 xxr_xyr_yxr_yyr_02 = _mm_setzero_ps();
	    __m128 xxi_xyi_yxi_yyi_02 = _mm_setzero_ps();

	    __m128 xxr_xyr_yxr_yyr_12 = _mm_setzero_ps();
	    __m128 xxi_xyi_yxi_yyi_12 = _mm_setzero_ps();

	    __m128 xxr_xyr_yxr_yyr_03 = _mm_setzero_ps();
	    __m128 xxi_xyi_yxi_yyi_03 = _mm_setzero_ps();

	    __m128 xxr_xyr_yxr_yyr_13 = _mm_setzero_ps();
	    __m128 xxi_xyi_yxi_yyi_13 = _mm_setzero_ps();

	    for (unsigned time = 0; time < nrTimes; time ++) {
		__m128 sample0 = _mm_load_ps(samples + time*4 + index0); // real_pol1, imag_pol1, real_pol2, imag_pol2
		__m128 sample1 = _mm_load_ps(samples + time*4 + index0 + nrTimesWidth); // real_pol1, imag_pol1, real_pol2, imag_pol2
		__m128 sample2 = _mm_load_ps(samples + time*4 + index2); // real_pol1, imag_pol1, real_pol2, imag_pol2
		__m128 sample3 = _mm_load_ps(samples + time*4 + index2 + nrTimesWidth); // real_pol1, imag_pol1, real_pol2, imag_pol2

                // _MM_SHUFFLE(z,y,x,w) selects x&w 32 bit double words from m1 and z&y from m2
		__m128 sample0_xr_xr_yr_yr = _mm_shuffle_ps(sample0, sample0, _MM_SHUFFLE(0, 0, 2, 2));
		__m128 sample1_xr_xr_yr_yr = _mm_shuffle_ps(sample1, sample1, _MM_SHUFFLE(0, 0, 2, 2));
		__m128 sample2_xr_yr_xr_yr = _mm_shuffle_ps(sample2, sample2, _MM_SHUFFLE(0, 2, 0, 2));
		__m128 sample3_xr_yr_xr_yr = _mm_shuffle_ps(sample3, sample3, _MM_SHUFFLE(0, 2, 0, 2));
		__m128 sample0_xi_xi_yi_yi = _mm_shuffle_ps(sample0, sample0, _MM_SHUFFLE(1, 1, 3, 3));
		__m128 sample1_xi_xi_yi_yi = _mm_shuffle_ps(sample1, sample1, _MM_SHUFFLE(1, 1, 3, 3));
		__m128 sample2_xi_yi_xi_yi = _mm_shuffle_ps(sample2, sample2, _MM_SHUFFLE(1, 3, 1, 3));
		__m128 sample3_xi_yi_xi_yi = _mm_shuffle_ps(sample3, sample3, _MM_SHUFFLE(1, 3, 1, 3));

		xxr_xyr_yxr_yyr_02 = _mm_add_ps(xxr_xyr_yxr_yyr_02, _mm_add_ps(_mm_mul_ps(sample0_xr_xr_yr_yr, sample2_xr_yr_xr_yr), 
									       _mm_mul_ps(sample0_xi_xi_yi_yi, sample2_xi_yi_xi_yi)));
		
		xxi_xyi_yxi_yyi_02 = _mm_add_ps(xxi_xyi_yxi_yyi_02, _mm_sub_ps(_mm_mul_ps(sample0_xi_xi_yi_yi, sample2_xr_yr_xr_yr),
									       _mm_mul_ps(sample0_xr_xr_yr_yr, sample2_xi_yi_xi_yi)));


		xxr_xyr_yxr_yyr_12 = _mm_add_ps(xxr_xyr_yxr_yyr_12, _mm_add_ps(_mm_mul_ps(sample1_xr_xr_yr_yr, sample2_xr_yr_xr_yr), 
									       _mm_mul_ps(sample1_xi_xi_yi_yi, sample2_xi_yi_xi_yi)));
		
		xxi_xyi_yxi_yyi_12 = _mm_add_ps(xxi_xyi_yxi_yyi_12, _mm_sub_ps(_mm_mul_ps(sample1_xi_xi_yi_yi, sample2_xr_yr_xr_yr),
									       _mm_mul_ps(sample1_xr_xr_yr_yr, sample2_xi_yi_xi_yi)));


		xxr_xyr_yxr_yyr_03 = _mm_add_ps(xxr_xyr_yxr_yyr_03, _mm_add_ps(_mm_mul_ps(sample0_xr_xr_yr_yr, sample3_xr_yr_xr_yr), 
									       _mm_mul_ps(sample0_xi_xi_yi_yi, sample3_xi_yi_xi_yi)));
		
		xxi_xyi_yxi_yyi_03 = _mm_add_ps(xxi_xyi_yxi_yyi_03, _mm_sub_ps(_mm_mul_ps(sample0_xi_xi_yi_yi, sample3_xr_yr_xr_yr),
									       _mm_mul_ps(sample0_xr_xr_yr_yr, sample3_xi_yi_xi_yi)));
		
		
		xxr_xyr_yxr_yyr_13 = _mm_add_ps(xxr_xyr_yxr_yyr_13, _mm_add_ps(_mm_mul_ps(sample1_xr_xr_yr_yr, sample3_xr_yr_xr_yr), 
									       _mm_mul_ps(sample1_xi_xi_yi_yi, sample3_xi_yi_xi_yi)));
		
		xxi_xyi_yxi_yyi_13 = _mm_add_ps(xxi_xyi_yxi_yyi_13, _mm_sub_ps(_mm_mul_ps(sample1_xi_xi_yi_yi, sample3_xr_yr_xr_yr),
									       _mm_mul_ps(sample1_xr_xr_yr_yr, sample3_xi_yi_xi_yi)));
	    }

	    if (cell < nrCells) {
		missedBaselines[BASELINE(stat0, stat2)] = false;
		missedBaselines[BASELINE(stat0+1, stat2)] = false;
		missedBaselines[BASELINE(stat0, stat2+1)] = false;
		missedBaselines[BASELINE(stat0+1, stat2+1)] = false;
   
		store(visibilities, stat0,   stat2,   channel, nrChannels, xxr_xyr_yxr_yyr_02, xxi_xyi_yxi_yyi_02);
		store(visibilities, stat0+1, stat2,   channel, nrChannels, xxr_xyr_yxr_yyr_12, xxi_xyi_yxi_yyi_12);
		store(visibilities, stat0,   stat2+1, channel, nrChannels, xxr_xyr_yxr_yyr_03, xxi_xyi_yxi_yyi_03);
		store(visibilities, stat0+1, stat2+1, channel, nrChannels, xxr_xyr_yxr_yyr_13, xxi_xyi_yxi_yyi_13);
	    }
	}
    }

    unsigned long long missedBytesLoaded, missedBytesStored;
    unsigned long long missedOps = computeMissedBaselines(samples, visibilities, missedBaselines,
							  nrTimes, nrTimesWidth, nrStations, nrChannels,
							  &missedBytesLoaded, &missedBytesStored);

    unsigned long long ops = missedOps + calcNrOps(nrCells, nrTimes, nrStations, nrChannels, bytesLoaded, bytesStored);
    *bytesLoaded += missedBytesLoaded;
    *bytesLoaded += missedBytesStored;
    return ops;
}
