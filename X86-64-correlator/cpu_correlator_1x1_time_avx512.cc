#include <cstdlib>
#include <iostream>
#include <immintrin.h> // avx512

#include "X86-64_correlator.h"

#if defined(__AVX512F__)


using namespace std;

// used for 1x1 correlator only
static unsigned char baselineToStat1[MAX_NR_CELLS], baselineToStat2[MAX_NR_CELLS];

static unsigned fillCellToStatTable(const unsigned nrStations)
{
    unsigned baseline;

    for (unsigned stat2 = baseline = 0; stat2 < nrStations; stat2 ++) {
	for (unsigned stat1 = 0; stat1 <= stat2; stat1 ++, baseline ++) {
	    baselineToStat1[baseline] = stat1;
	    baselineToStat2[baseline] = stat2;
	}
    }

    return baseline;
}

unsigned long long cpuCorrelator_1x1_time_avx512(const float* __restrict__ samples, float* __restrict__ visibilities,
						 const unsigned nrTimes, const unsigned nrTimesWidth,
						 const unsigned nrStations, const unsigned nrChannels,
						 unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned nrBaselines = fillCellToStatTable(nrStations);

//    __mmask16 mask = _mm512_int2mask(0b01010101);

    
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned baseline = 0; baseline < nrBaselines; baseline++) {
	    unsigned stat1 = baselineToStat1[baseline];
	    unsigned stat2 = baselineToStat2[baseline];

	    size_t index1 = SAMPLE_INDEX(stat1, channel, 0, 0, 0);
	    size_t index2 = SAMPLE_INDEX(stat2, channel, 0, 0, 0);

	    __m512 xxr = _mm512_setzero_ps ();
	    __m512 xxi = _mm512_setzero_ps ();
	    __m512 xyr = _mm512_setzero_ps ();
	    __m512 xyi = _mm512_setzero_ps ();
	    __m512 yxr = _mm512_setzero_ps ();
	    __m512 yxi = _mm512_setzero_ps ();
	    __m512 yyr = _mm512_setzero_ps ();
	    __m512 yyi = _mm512_setzero_ps ();

	    // assume nrTimes is divisable by vector width
	    for (unsigned time = 0; time < nrTimes; time += VECTOR_WIDTH_IN_FLOATS) {
	      // This assumes a different memory layout with four samples in time behind each other:
	      // (xr1 xr2 xr3 xr4) (xi1 xi2 xi3 xi4) (yr1 yr2 yr3 yr4) (yi1 yi2 yi3 yi4), etc.
		__m512 sample1xr = _mm512_load_ps(samples + time*4 + index1+0*VECTOR_WIDTH_IN_FLOATS);
		__m512 sample1xi = _mm512_load_ps(samples + time*4 + index1+1*VECTOR_WIDTH_IN_FLOATS);
		__m512 sample1yr = _mm512_load_ps(samples + time*4 + index1+2*VECTOR_WIDTH_IN_FLOATS);
		__m512 sample1yi = _mm512_load_ps(samples + time*4 + index1+3*VECTOR_WIDTH_IN_FLOATS);

		__m512 sample2xr = _mm512_load_ps(samples + time*4 + index2+0*VECTOR_WIDTH_IN_FLOATS);
		__m512 sample2xi = _mm512_load_ps(samples + time*4 + index2+1*VECTOR_WIDTH_IN_FLOATS);
		__m512 sample2yr = _mm512_load_ps(samples + time*4 + index2+2*VECTOR_WIDTH_IN_FLOATS);
		__m512 sample2yi = _mm512_load_ps(samples + time*4 + index2+3*VECTOR_WIDTH_IN_FLOATS);

                xxr += sample1xr * sample2xr;
                xxi += sample1xi * sample2xr;
                xyr += sample1xr * sample2yr;
                xyi += sample1xi * sample2yr;
                yxr += sample1yr * sample2xr;
                yxi += sample1yi * sample2xr;
                yyr += sample1yr * sample2yr;
                yyi += sample1yi * sample2yr;
                xxr += sample1xi * sample2xi;
                xxi -= sample1xr * sample2xi;
                xyr += sample1xi * sample2yi;
                xyi -= sample1xr * sample2yi;
                yxr += sample1yi * sample2xi;
                yxi -= sample1yr * sample2xi;
                yyr += sample1yi * sample2yi;
                yyi -= sample1yr * sample2yi;
	    }

	    size_t vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);

	    float res_xxr = _mm512_reduce_add_ps(xxr);
	    float res_xxi = _mm512_reduce_add_ps(xxi);
	    float res_xyr = _mm512_reduce_add_ps(xyr);
	    float res_xyi = _mm512_reduce_add_ps(xyi);
	    float res_yxr = _mm512_reduce_add_ps(yxr);
	    float res_yxi = _mm512_reduce_add_ps(yxi);
	    float res_yyr = _mm512_reduce_add_ps(yyr);
	    float res_yyi = _mm512_reduce_add_ps(yyi);
	    
	    visibilities[vis_index+0] = res_xxr;
	    visibilities[vis_index+1] = res_xxi;
	    visibilities[vis_index+2] = res_xyr;
	    visibilities[vis_index+3] = res_xyi;
	    visibilities[vis_index+4] = res_yxr;
	    visibilities[vis_index+5] = res_yxi;
	    visibilities[vis_index+6] = res_yyr;
	    visibilities[vis_index+7] = res_yyi;
	}
    }

    *bytesLoaded = nrChannels * nrBaselines * nrTimes * 8 * sizeof(float); // samples
    *bytesStored = nrChannels * nrBaselines * 8 * sizeof(float); // vis

    unsigned long long ops = nrChannels * nrBaselines * nrTimes * 16L * 2L;
    return ops;
}

#endif // defined(__AVX512F__)
