#include "common.h"
#include "cpu_correlator.h"
#include <iostream>

using namespace std;

unsigned long long cpuCorrelator_2x1(const float* __restrict__ samples, float* __restrict__ visibilities,
				       const unsigned nrTimes, const unsigned nrTimesWidth,
				       const unsigned nrStations, const unsigned nrChannels,
				       unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    for (unsigned channel = 0; channel < nrChannels; channel++) {
	for (unsigned stationY = 1; stationY < nrStations; stationY++) {
	    for (unsigned stationX = 0; stationX + 2 <= stationY; stationX += 2) {
		float s0xxr = 0, s0xxi = 0;
		float s0xyr = 0, s0xyi = 0;
		float s0yxr = 0, s0yxi = 0;
		float s0yyr = 0, s0yyi = 0;

		float s1xxr = 0, s1xxi = 0;
		float s1xyr = 0, s1xyi = 0;
		float s1yxr = 0, s1yxi = 0;
		float s1yyr = 0, s1yyi = 0;

		unsigned index0 = SAMPLE_INDEX(stationX,   channel, 0, 0, 0);
		unsigned index1 = SAMPLE_INDEX(stationX+1, channel, 0, 0, 0);

		for (unsigned time = 0; time < nrTimes; time ++) {
		    float sample0xr = samples[index0+0];
		    float sample0xi = samples[index0+1];
		    float sample0yr = samples[index0+2];
		    float sample0yi = samples[index0+3];

		    float sample1xr = samples[index1+0];
		    float sample1xi = samples[index1+1];
		    float sample1yr = samples[index1+2];
		    float sample1yi = samples[index1+3];
		
		    s0xxr  += sample0xr * sample1xr;
		    s0xxi  += sample0xi * sample1xr;
		    s0xyr  += sample0xr * sample1yr;
		    s0xyi  += sample0xi * sample1yr;
		    s0yxr  += sample0yr * sample1xr;
		    s0yxi  += sample0yi * sample1xr;
		    s0yyr  += sample0yr * sample1yr;
		    s0yyi  += sample0yi * sample1yr;
		    s0xxr  += sample0xi * sample1xi;
		    s0xxi  -= sample0xr * sample1xi;
		    s0xyr  += sample0xi * sample1yi;
		    s0xyi  -= sample0xr * sample1yi;
		    s0yxr  += sample0yi * sample1xi;
		    s0yxi  -= sample0yr * sample1xi;
		    s0yyr  += sample0yi * sample1yi;
		    s0yyi  -= sample0yr * sample1yi;

		    s1xxr  += sample1xr * sample1xr;
		    s1xxi  += sample1xi * sample1xr;
		    s1xyr  += sample1xr * sample1yr;
		    s1xyi  += sample1xi * sample1yr;
		    s1yxr  += sample1yr * sample1xr;
		    s1yxi  += sample1yi * sample1xr;
		    s1yyr  += sample1yr * sample1yr;
		    s1yyi  += sample1yi * sample1yr;
		    s1xxr  += sample1xi * sample1xi;
		    s1xxi  -= sample1xr * sample1xi;
		    s1xyr  += sample1xi * sample1yi;
		    s1xyi  -= sample1xr * sample1yi;
		    s1yxr  += sample1yi * sample1xi;
		    s1yxi  -= sample1yr * sample1xi;
		    s1yyr  += sample1yi * sample1yi;
		    s1yyi  -= sample1yr * sample1yi;

		    index0 += 4;
		    index1 += 4;
		}

		unsigned baseline = BASELINE(stationX, stationY);
		unsigned vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
		visibilities[vis_index+0] = s0xxr; 
		visibilities[vis_index+1] = s0xxi;
		visibilities[vis_index+2] = s0xyr;
		visibilities[vis_index+3] = s0xyi;
		visibilities[vis_index+4] = s0yxr;
		visibilities[vis_index+5] = s0yxi;
		visibilities[vis_index+6] = s0yyr;
		visibilities[vis_index+7] = s0yyi;
		    
		baseline = BASELINE(stationX+1, stationY);
		vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
		visibilities[vis_index+0] = s1xxr; 
		visibilities[vis_index+1] = s1xxi;
		visibilities[vis_index+2] = s1xyr;
		visibilities[vis_index+3] = s1xyi;
		visibilities[vis_index+4] = s1yxr;
		visibilities[vis_index+5] = s1yxi;
		visibilities[vis_index+6] = s1yyr;
		visibilities[vis_index+7] = s1yyi;
	    }
	}
    }
    
    unsigned nrCells = calcNrCells(2, 1, nrStations);    
    *bytesLoaded = nrChannels * nrCells * nrTimes * 8L * sizeof(float); // samples
    *bytesStored = nrChannels * nrCells * 8L * 2L * sizeof(float); // visibilities
    return nrChannels * nrCells * nrTimes * 16L * 2L * 2L;
}
