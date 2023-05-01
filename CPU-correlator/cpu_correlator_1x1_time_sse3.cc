/*
Copyright (C) 2009 Rob van Nieuwpoort & John Romein
Astron
P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, nieuwpoort@astron.nl

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <cstdlib>
#include <iostream>
#include <pmmintrin.h> // sse3

#include "cpu_correlator.h"


using namespace std;

// used for 1x1 correlator only
static unsigned char baselineToStat1[64 * 65 / 2], baselineToStat2[64 * 65 / 2];

static unsigned fillCellToStatTable()
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

// assume samples are in the correct order: xr0, xr1, xr2, xr3; xi0, xi1, xi2, xi3; yr0, yr1, yr2, yr3; ui0, yi1, yi2, yi3

unsigned long long cpuCorrelator_1x1_time_sse3(float* samples, float* visibilities, 
					       unsigned nrTimes, unsigned nrTimesWidth,
					       unsigned nrStations, unsigned nrChannels,
					       unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned nrBaselines = fillCellToStatTable();

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned baseline = 0; baseline < nrBaselines; baseline++) {
	    unsigned stat1 = baselineToStat1[baseline];
	    unsigned stat2 = baselineToStat2[baseline];

	    unsigned index1 = SAMPLE_INDEX(stat1, channel, 0, 0, 0);
	    unsigned index2 = SAMPLE_INDEX(stat2, channel, 0, 0, 0);

	    __m128 xxr = _mm_setzero_ps ();
	    __m128 xxi = _mm_setzero_ps ();
	    __m128 xyr = _mm_setzero_ps ();
	    __m128 xyi = _mm_setzero_ps ();
	    __m128 yxr = _mm_setzero_ps ();
	    __m128 yxi = _mm_setzero_ps ();
	    __m128 yyr = _mm_setzero_ps ();
	    __m128 yyi = _mm_setzero_ps ();

	    // assume nrTimes is divisable by 4
	    for (unsigned time = 0; time < nrTimes; time += 4) {
	      // This assumes a different memory layout with four samples in time behind each other:
	      // (xr1 xr2 xr3 xr4) (xi1 xi2 xi3 xi4) (yr1 yr2 yr3 yr4) (yi1 yi2 yi3 yi4), etc.
		__m128 sample1xr = _mm_load_ps(samples + time*4 + index1+0);
		__m128 sample1xi = _mm_load_ps(samples + time*4 + index1+4);
		__m128 sample1yr = _mm_load_ps(samples + time*4 + index1+8);
		__m128 sample1yi = _mm_load_ps(samples + time*4 + index1+12);

		__m128 sample2xr = _mm_load_ps(samples + time*4 + index2+0);
		__m128 sample2xi = _mm_load_ps(samples + time*4 + index2+4);
		__m128 sample2yr = _mm_load_ps(samples + time*4 + index2+8);
		__m128 sample2yi = _mm_load_ps(samples + time*4 + index2+12);

#if 1
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
#else
		xxr = _mm_add_ps(xxr, _mm_add_ps(_mm_mul_ps(sample1xr, sample2xr), _mm_mul_ps(sample1xi, sample2xi)));
		xxi = _mm_add_ps(xxi, _mm_sub_ps(_mm_mul_ps(sample1xi, sample2xr), _mm_mul_ps(sample1xr, sample2xi)));

		xyr = _mm_add_ps(xyr, _mm_add_ps(_mm_mul_ps(sample1xr, sample2yr), _mm_mul_ps(sample1xi, sample2yi)));
		xyi = _mm_add_ps(xyi, _mm_sub_ps(_mm_mul_ps(sample1xi, sample2yr), _mm_mul_ps(sample1xr, sample2yi))); 

		yxr = _mm_add_ps(yxr, _mm_add_ps(_mm_mul_ps(sample1yr, sample2xr), _mm_mul_ps(sample1yi, sample2xi)));
		yxi = _mm_add_ps(yxi, _mm_sub_ps(_mm_mul_ps(sample1yi, sample2xr), _mm_mul_ps(sample1yr, sample2xi))); 

		yyr = _mm_add_ps(yyr, _mm_add_ps(_mm_mul_ps(sample1yr, sample2yr), _mm_mul_ps(sample1yi, sample2yi)));
		yyi = _mm_add_ps(yyi, _mm_sub_ps(_mm_mul_ps(sample1yi, sample2yr), _mm_mul_ps(sample1yr, sample2yi))); 
#endif
	    }

	    if (baseline < nrBaselines) {
		unsigned vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);

		// now, we have to sum the 4 values inside the regs.
		__m128 tmp0  = _mm_hadd_ps(xxr, xxi);   // xxr0+xxr1, xxr2+xxr3, xxi0+xxi1, xxi2+xxi3
		__m128 tmp1  = _mm_hadd_ps(xyr, xyi);   // xyr0+xyr1, xyr2+xyr3, xyi0+xyi1, xyi2+xyi3
		__m128 xx_xy = _mm_hadd_ps(tmp0, tmp1); // xxr, xxi, xyr,xyi

		__m128 tmp2  = _mm_hadd_ps(yxr, yxi);
		__m128 tmp3  = _mm_hadd_ps(yyr, yyi);
		__m128 yx_yy = _mm_hadd_ps(tmp2, tmp3);

		_mm_store_ps(visibilities + vis_index + 0, xx_xy);
		_mm_store_ps(visibilities + vis_index + 4, yx_yy);
	    } else {
		    cout << "X" << endl;
	    }
	}
    }

    *bytesLoaded = nrChannels * nrBaselines * nrTimes * 8 * sizeof(float); // samples
    *bytesStored = nrChannels * nrBaselines * 8 * sizeof(float); // vis

    unsigned long long ops = nrChannels * nrBaselines * nrTimes * 16L * 2L;
    return ops;
}
