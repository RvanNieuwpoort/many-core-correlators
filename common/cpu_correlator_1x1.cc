#include <cstdlib>
#include <iostream>

#include "common.h"

using namespace std;

unsigned long long cpuCorrelator_1x1(const float* __restrict__ samples, float* __restrict__ visibilities, 
				     const unsigned nrTimes, const unsigned nrTimesWidth,
				     const unsigned nrStations, const unsigned nrChannels,
				     unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    const unsigned nrBaselines = nrStations * (nrStations + 1) / 2;

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for(unsigned int stationY = 0; stationY < nrStations; stationY++) {
	    for(unsigned int stationX = 0; stationX <= stationY; stationX++) {

		float xxr = 0, xxi = 0, xyr = 0, xyi = 0, yxr = 0, yxi = 0, yyr = 0, yyi = 0;
		unsigned index1 = SAMPLE_INDEX(stationX, channel, 0, 0, 0);
		unsigned index2 = SAMPLE_INDEX(stationY, channel, 0, 0, 0);

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
		unsigned baseline = BASELINE(stationX, stationY);
		if (baseline < nrBaselines) {
		    const unsigned vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
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
    }

    *bytesLoaded = nrChannels * nrBaselines * nrTimes * 8L * sizeof(float); // samples
    *bytesStored = nrChannels * nrBaselines * 8L * sizeof(float); // vis
    
    return nrChannels * nrBaselines * nrTimes * 16L * 2L;
}
