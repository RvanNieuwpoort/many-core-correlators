#include "riscv_correlator.h"
#include <iostream>

static unsigned char cellToStatX[MAX_CELLS], cellToStatY[MAX_CELLS];

static unsigned fillCellToStatTable(unsigned nrStations)
{
    unsigned nrCells, stat0, stat2;

    for (stat2 = nrStations % 2 ? 1 : 2, nrCells = 0; stat2 < nrStations; stat2 += 2) {
	for (stat0 = 0; stat0 + 2 <= stat2; stat0 += 2, nrCells++) {
	    cellToStatX[nrCells] = stat0;
	    cellToStatY[nrCells] = stat2;
	}
    }

    return nrCells;
}

unsigned long long riscvCorrelator_2x2(const float* __restrict__ samples, float* __restrict__ visibilities,
				       const unsigned nrTimes, const unsigned nrTimesWidth,
				       const unsigned nrStations, const unsigned nrChannels,
				       unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    const unsigned nrCells = fillCellToStatTable(nrStations);

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned cell = 0; cell < nrCells; cell ++) {
	    const unsigned stat0 = cellToStatX[cell];
	    const unsigned stat2 = cellToStatY[cell];

	    float v0x2xr = 0, v0x2xi = 0;
	    float v0x2yr = 0, v0x2yi = 0;
	    float v0y2xr = 0, v0y2xi = 0;
	    float v0y2yr = 0, v0y2yi = 0;

	    float v1x2xr = 0, v1x2xi = 0;
	    float v1x2yr = 0, v1x2yi = 0;
	    float v1y2xr = 0, v1y2xi = 0;
	    float v1y2yr = 0, v1y2yi = 0;

	    float v0x3xr = 0, v0x3xi = 0;
	    float v0x3yr = 0, v0x3yi = 0;
	    float v0y3xr = 0, v0y3xi = 0;
	    float v0y3yr = 0, v0y3yi = 0;

	    float v1x3xr = 0, v1x3xi = 0;
	    float v1x3yr = 0, v1x3yi = 0;
	    float v1y3xr = 0, v1y3xi = 0;
	    float v1y3yr = 0, v1y3yi = 0;

	    unsigned index0 = SAMPLE_INDEX(stat0,   channel, 0, 0, 0);
	    unsigned index1 = SAMPLE_INDEX(stat0+1, channel, 0, 0, 0);
	    unsigned index2 = SAMPLE_INDEX(stat2,   channel, 0, 0, 0);
	    unsigned index3 = SAMPLE_INDEX(stat2+1, channel, 0, 0, 0);

	    for (unsigned time = 0; time < nrTimes; time ++) {
		float sample0xr = samples[index0+0];
		float sample0xi = samples[index0+1];
		float sample0yr = samples[index0+2];
		float sample0yi = samples[index0+3];

		float sample1xr = samples[index1+0];
		float sample1xi = samples[index1+1];
		float sample1yr = samples[index1+2];
		float sample1yi = samples[index1+3];

		float sample2xr = samples[index2+0];
		float sample2xi = samples[index2+1];
		float sample2yr = samples[index2+2];
		float sample2yi = samples[index2+3];

		float sample3xr = samples[index3+0];
		float sample3xi = samples[index3+1];
		float sample3yr = samples[index3+2];
		float sample3yi = samples[index3+3];
		
		v0x2xr  += sample0xr * sample2xr;
		v0x2xi  += sample0xi * sample2xr;
		v0x2yr  += sample0xr * sample2yr;
		v0x2yi  += sample0xi * sample2yr;
		v0y2xr  += sample0yr * sample2xr;
		v0y2xi  += sample0yi * sample2xr;
		v0y2yr  += sample0yr * sample2yr;
		v0y2yi  += sample0yi * sample2yr;
		v0x2xr  += sample0xi * sample2xi;
		v0x2xi  -= sample0xr * sample2xi;
		v0x2yr  += sample0xi * sample2yi;
		v0x2yi  -= sample0xr * sample2yi;
		v0y2xr  += sample0yi * sample2xi;
		v0y2xi  -= sample0yr * sample2xi;
		v0y2yr  += sample0yi * sample2yi;
		v0y2yi  -= sample0yr * sample2yi;

		v1x2xr  += sample1xr * sample2xr;
		v1x2xi  += sample1xi * sample2xr;
		v1x2yr  += sample1xr * sample2yr;
		v1x2yi  += sample1xi * sample2yr;
		v1y2xr  += sample1yr * sample2xr;
		v1y2xi  += sample1yi * sample2xr;
		v1y2yr  += sample1yr * sample2yr;
		v1y2yi  += sample1yi * sample2yr;
		v1x2xr  += sample1xi * sample2xi;
		v1x2xi  -= sample1xr * sample2xi;
		v1x2yr  += sample1xi * sample2yi;
		v1x2yi  -= sample1xr * sample2yi;
		v1y2xr  += sample1yi * sample2xi;
		v1y2xi  -= sample1yr * sample2xi;
		v1y2yr  += sample1yi * sample2yi;
		v1y2yi  -= sample1yr * sample2yi;

		v0x3xr  += sample0xr * sample3xr;
		v0x3xi  += sample0xi * sample3xr;
		v0x3yr  += sample0xr * sample3yr;
		v0x3yi  += sample0xi * sample3yr;
		v0y3xr  += sample0yr * sample3xr;
		v0y3xi  += sample0yi * sample3xr;
		v0y3yr  += sample0yr * sample3yr;
		v0y3yi  += sample0yi * sample3yr;
		v0x3xr  += sample0xi * sample3xi;
		v0x3xi  -= sample0xr * sample3xi;
		v0x3yr  += sample0xi * sample3yi;
		v0x3yi  -= sample0xr * sample3yi;
		v0y3xr  += sample0yi * sample3xi;
		v0y3xi  -= sample0yr * sample3xi;
		v0y3yr  += sample0yi * sample3yi;
		v0y3yi  -= sample0yr * sample3yi;

		v1x3xr  += sample1xr * sample3xr;
		v1x3xi  += sample1xi * sample3xr;
		v1x3yr  += sample1xr * sample3yr;
		v1x3yi  += sample1xi * sample3yr;
		v1y3xr  += sample1yr * sample3xr;
		v1y3xi  += sample1yi * sample3xr;
		v1y3yr  += sample1yr * sample3yr;
		v1y3yi  += sample1yi * sample3yr;
		v1x3xr  += sample1xi * sample3xi;
		v1x3xi  -= sample1xr * sample3xi;
		v1x3yr  += sample1xi * sample3yi;
		v1x3yi  -= sample1xr * sample3yi;
		v1y3xr  += sample1yi * sample3xi;
		v1y3xi  -= sample1yr * sample3xi;
		v1y3yr  += sample1yi * sample3yi;
		v1y3yi  -= sample1yr * sample3yi;

		index0 += 4;
		index1 += 4;
		index2 += 4;
		index3 += 4;
	    }

	    if (cell < nrCells) {
		unsigned baseline = BASELINE(stat0, stat2);
		unsigned vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
		visibilities[vis_index+0] = v0x2xr; 
		visibilities[vis_index+1] = v0x2xi;
		visibilities[vis_index+2] = v0x2yr;
		visibilities[vis_index+3] = v0x2yi;
		visibilities[vis_index+4] = v0y2xr;
		visibilities[vis_index+5] = v0y2xi;
		visibilities[vis_index+6] = v0y2yr;
		visibilities[vis_index+7] = v0y2yi;

		baseline = BASELINE(stat0+1, stat2);
		vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
		visibilities[vis_index+0] = v1x2xr; 
		visibilities[vis_index+1] = v1x2xi;
		visibilities[vis_index+2] = v1x2yr;
		visibilities[vis_index+3] = v1x2yi;
		visibilities[vis_index+4] = v1y2xr;
		visibilities[vis_index+5] = v1y2xi;
		visibilities[vis_index+6] = v1y2yr;
		visibilities[vis_index+7] = v1y2yi;

		baseline = BASELINE(stat0, stat2+1);
		vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
		visibilities[vis_index+0] = v0x3xr; 
		visibilities[vis_index+1] = v0x3xi;
		visibilities[vis_index+2] = v0x3yr;
		visibilities[vis_index+3] = v0x3yi;
		visibilities[vis_index+4] = v0y3xr;
		visibilities[vis_index+5] = v0y3xi;
		visibilities[vis_index+6] = v0y3yr;
		visibilities[vis_index+7] = v0y3yi;

		baseline = BASELINE(stat0+1, stat2+1);
		vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
		visibilities[vis_index+0] = v1x3xr; 
		visibilities[vis_index+1] = v1x3xi;
		visibilities[vis_index+2] = v1x3yr;
		visibilities[vis_index+3] = v1x3yi;
		visibilities[vis_index+4] = v1y3xr;
		visibilities[vis_index+5] = v1y3xi;
		visibilities[vis_index+6] = v1y3yr;
		visibilities[vis_index+7] = v1y3yi;
	    }
	}
    }

    *bytesLoaded = nrChannels * nrCells * nrTimes * 16L * sizeof(float); // samples
    *bytesStored = nrChannels * nrCells * 8L * 4L * sizeof(float); // vis

    return nrChannels * nrCells * 4L * nrTimes * 16L * 2L;
}
