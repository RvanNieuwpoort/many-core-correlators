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

unsigned long long cpuCorrelator_1x1_sse3(float* samples, float* visibilities, 
					  unsigned nrTimes, unsigned nrTimesWidth, unsigned nrStations, unsigned nrChannels,
				     unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned nrBaselines = fillCellToStatTable();

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned baseline = 0; baseline < nrBaselines; baseline++) {
	    unsigned stat1 = baselineToStat1[baseline];
	    unsigned stat2 = baselineToStat2[baseline];

	    unsigned index1 = SAMPLE_INDEX(stat1, channel, 0, 0, 0);
	    unsigned index2 = SAMPLE_INDEX(stat2, channel, 0, 0, 0);

	    __m128 xxr_xyr_yxr_yyr = _mm_setzero_ps();
	    __m128 xxi_xyi_yxi_yyi = _mm_setzero_ps();

	    for (unsigned time = 0; time < nrTimes; time ++) {
		__m128 sample1 = _mm_load_ps(&(samples[index1])); // real_pol1, imag_pol1, real_pol2, imag_pol2
		__m128 sample2 = _mm_load_ps(&(samples[index2]));

                // _MM_SHUFFLE(z,y,x,w) selects x&w 32 bit double words from m1 and z&y from m2
		__m128 sample1_xr_xr_yr_yr = _mm_shuffle_ps(sample1, sample1, _MM_SHUFFLE(0, 0, 2, 2));
		__m128 sample2_xr_yr_xr_yr = _mm_shuffle_ps(sample2, sample2, _MM_SHUFFLE(0, 2, 0, 2));
		__m128 sample1_xi_xi_yi_yi = _mm_shuffle_ps(sample1, sample1, _MM_SHUFFLE(1, 1, 3, 3));
		__m128 sample2_xi_yi_xi_yi = _mm_shuffle_ps(sample2, sample2, _MM_SHUFFLE(1, 3, 1, 3));

		xxr_xyr_yxr_yyr = _mm_add_ps(xxr_xyr_yxr_yyr, _mm_add_ps(_mm_mul_ps(sample1_xr_xr_yr_yr, sample2_xr_yr_xr_yr), 
									 _mm_mul_ps(sample1_xi_xi_yi_yi, sample2_xi_yi_xi_yi)));
		
		xxi_xyi_yxi_yyi = _mm_add_ps(xxi_xyi_yxi_yyi, _mm_sub_ps(_mm_mul_ps(sample1_xi_xi_yi_yi, sample2_xr_yr_xr_yr),
									 _mm_mul_ps(sample1_xr_xr_yr_yr, sample2_xi_yi_xi_yi)));
		index1 += 4;
		index2 += 4;
	    }
	    if (baseline < nrBaselines) {
		unsigned vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);

		__m128 vis1tmp = _mm_shuffle_ps(xxr_xyr_yxr_yyr, xxi_xyi_yxi_yyi, _MM_SHUFFLE(3, 2, 3, 2));
		__m128 vis1    = _mm_shuffle_ps(vis1tmp, vis1tmp, _MM_SHUFFLE(2, 0, 3, 1));
		_mm_store_ps(&(visibilities[vis_index+0]), vis1);

		__m128 vis2tmp = _mm_shuffle_ps(xxr_xyr_yxr_yyr, xxi_xyi_yxi_yyi, _MM_SHUFFLE(1, 0, 1, 0));
		__m128 vis2    = _mm_shuffle_ps(vis2tmp, vis2tmp, _MM_SHUFFLE(2, 0, 3, 1));
		_mm_store_ps(&(visibilities[vis_index+4]), vis2);
/*
		fprintf(stderr, "vis: %f %f %f %f %f %f %f %f\n", 
			visibilities[vis_index+0], visibilities[vis_index+1], visibilities[vis_index+2], visibilities[vis_index+3],
			visibilities[vis_index+4], visibilities[vis_index+5], visibilities[vis_index+6], visibilities[vis_index+7]);
*/
	    }
	}
    }

    unsigned long long ops = nrChannels * nrBaselines * nrTimes * 16L * 2L;
    *bytesLoaded = nrChannels * nrBaselines * nrTimes * 8 * sizeof(float); // samples
    *bytesStored = nrChannels * nrBaselines * 8 * sizeof(float); // vis
    
    return ops;
}
