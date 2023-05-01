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

unsigned long long cpuCorrelator_1x1(float* samples, float* visibilities, 
				     unsigned nrTimes, unsigned nrTimesWidth, unsigned nrStations, unsigned nrChannels,
				     unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned nrBaselines = fillCellToStatTable();

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned baseline = 0; baseline < nrBaselines; baseline++) {
	    unsigned stat1 = baselineToStat1[baseline];
	    unsigned stat2 = baselineToStat2[baseline];

	    float xxr = 0, xxi = 0, xyr = 0, xyi = 0, yxr = 0, yxi = 0, yyr = 0, yyi = 0;
	    unsigned index1 = SAMPLE_INDEX(stat1, channel, 0, 0, 0);
	    unsigned index2 = SAMPLE_INDEX(stat2, channel, 0, 0, 0);
	    for (unsigned time = 0; time < nrTimes; time ++) {

		float sample1xr = samples[index1+0];
		float sample1xi = samples[index1+1];
		float sample1yr = samples[index1+2];
		float sample1yi = samples[index1+3];
		float sample2xr = samples[index2+0];
		float sample2xi = samples[index2+1];
		float sample2yr = samples[index2+2];
		float sample2yi = samples[index2+3];

		xxr += sample1xr * sample2xr + sample1xi * sample2xi;
		xxi += sample1xi * sample2xr - sample1xr * sample2xi;

		xyr += sample1xr * sample2yr + sample1xi * sample2yi;
		xyi += sample1xi * sample2yr - sample1xr * sample2yi;

		yxr += sample1yr * sample2xr + sample1yi * sample2xi;
		yxi += sample1yi * sample2xr - sample1yr * sample2xi;

		yyr += sample1yr * sample2yr + sample1yi * sample2yi;
		yyi += sample1yi * sample2yr - sample1yr * sample2yi;

		index1 += 4;
		index2 += 4;
	    }
	    if (baseline < nrBaselines) {
		unsigned vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
		visibilities[vis_index+0] = xxr;
		visibilities[vis_index+1] = xxi;
		visibilities[vis_index+2] = xyr;
		visibilities[vis_index+3] = xyi;
		visibilities[vis_index+4] = yxr;
		visibilities[vis_index+5] = yxi;
		visibilities[vis_index+6] = yyr;
		visibilities[vis_index+7] = yyi;
	    }
	}
    }

    unsigned long long ops = nrChannels * nrBaselines * nrTimes * 16L * 2L;

    *bytesLoaded = nrChannels * nrBaselines * nrTimes * 8 * sizeof(float); // samples
    *bytesStored = nrChannels * nrBaselines * 8 * sizeof(float); // vis

    return ops;
}
