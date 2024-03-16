#include "cpu_correlator.h"

static unsigned char cellToStatX[MAX_CELLS], cellToStatY[MAX_CELLS];

static unsigned fillCellToStatTable()
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


static unsigned long long calcNrOps()
{
    unsigned nrCells, stat0, stat2;

    for (stat2 = nrStations % 2 ? 1 : 2, nrCells = 0; stat2 < nrStations; stat2 += 2) {
	for (stat0 = 0; stat0 + 2 <= stat2; stat0 += 2, nrCells ++);
    }

    unsigned long long ops = nrChannels * nrCells * nrTimes * 16L * 4L * 2L;
    return ops;
}

void correlate_2x2(complex<float>* samples, complex<float>* visibilities, 
		   unsigned nrTimes, unsigned nrTimesWidth, unsigned nrStations, unsigned nrChannels)
{
    unsigned nrCells = fillCellToStatTable();
    unsigned loopCount = LOOP_COUNT(nrCells, 1);

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned cell = 0; cell < loopCount; cell ++) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat2 = cellToStatY[cell];
	    unsigned index0 = SAMPLE_INDEX(stat0, channel, 0, 0) / 2;
	    unsigned index2 = SAMPLE_INDEX(stat2, channel, 0, 0) / 2;

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

	    for (unsigned time = 0; time < nrTimes; time ++) {
		float4 sample0 = tex1Dfetch(tex0, index0 + time);
		float4 sample1 = tex1Dfetch(tex0, index0 + time + nrTimesWidth);
		float4 sample2 = tex1Dfetch(tex0, index2 + time);
		float4 sample3 = tex1Dfetch(tex0, index2 + time + nrTimesWidth);

		v0x2xr  += sample0.x * sample2.x;
		v0x2xi  += sample0.y * sample2.x;
		v0x2yr  += sample0.x * sample2.z;
		v0x2yi  += sample0.y * sample2.z;
		v0y2xr  += sample0.z * sample2.x;
		v0y2xi  += sample0.w * sample2.x;
		v0y2yr  += sample0.z * sample2.z;
		v0y2yi  += sample0.w * sample2.z;
		v0x2xr  += sample0.y * sample2.y;
		v0x2xi  -= sample0.x * sample2.y;
		v0x2yr  += sample0.y * sample2.w;
		v0x2yi  -= sample0.x * sample2.w;
		v0y2xr  += sample0.w * sample2.y;
		v0y2xi  -= sample0.z * sample2.y;
		v0y2yr  += sample0.w * sample2.w;
		v0y2yi  -= sample0.z * sample2.w;

		// the load of sample 1 can be done here, but for
		// some reason it was slower...

		v1x2xr  += sample1.x * sample2.x;
		v1x2xi  += sample1.y * sample2.x;
		v1x2yr  += sample1.x * sample2.z;
		v1x2yi  += sample1.y * sample2.z;
		v1y2xr  += sample1.z * sample2.x;
		v1y2xi  += sample1.w * sample2.x;
		v1y2yr  += sample1.z * sample2.z;
		v1y2yi  += sample1.w * sample2.z;
		v1x2xr  += sample1.y * sample2.y;
		v1x2xi  -= sample1.x * sample2.y;
		v1x2yr  += sample1.y * sample2.w;
		v1x2yi  -= sample1.x * sample2.w;
		v1y2xr  += sample1.w * sample2.y;
		v1y2xi  -= sample1.z * sample2.y;
		v1y2yr  += sample1.w * sample2.w;
		v1y2yi  -= sample1.z * sample2.w;

		// load of sample3 can be done here, reusing sample2's registers

		v0x3xr  += sample0.x * sample3.x;
		v0x3xi  += sample0.y * sample3.x;
		v0x3yr  += sample0.x * sample3.z;
		v0x3yi  += sample0.y * sample3.z;
		v0y3xr  += sample0.z * sample3.x;
		v0y3xi  += sample0.w * sample3.x;
		v0y3yr  += sample0.z * sample3.z;
		v0y3yi  += sample0.w * sample3.z;
		v0x3xr  += sample0.y * sample3.y;
		v0x3xi  -= sample0.x * sample3.y;
		v0x3yr  += sample0.y * sample3.w;
		v0x3yi  -= sample0.x * sample3.w;
		v0y3xr  += sample0.w * sample3.y;
		v0y3xi  -= sample0.z * sample3.y;
		v0y3yr  += sample0.w * sample3.w;
		v0y3yi  -= sample0.z * sample3.w;

		v1x3xr  += sample1.x * sample3.x;
		v1x3xi  += sample1.y * sample3.x;
		v1x3yr  += sample1.x * sample3.z;
		v1x3yi  += sample1.y * sample3.z;
		v1y3xr  += sample1.z * sample3.x;
		v1y3xi  += sample1.w * sample3.x;
		v1y3yr  += sample1.z * sample3.z;
		v1y3yi  += sample1.w * sample3.z;
		v1x3xr  += sample1.y * sample3.y;
		v1x3xi  -= sample1.x * sample3.y;
		v1x3yr  += sample1.y * sample3.w;
		v1x3yi  -= sample1.x * sample3.w;
		v1y3xr  += sample1.w * sample3.y;
		v1y3xi  -= sample1.z * sample3.y;
		v1y3yr  += sample1.w * sample3.w;
		v1y3yi  -= sample1.z * sample3.w;
	    }

	    if (cell < nrCells) {
		unsigned baseline = BASELINE(stat0, stat2);
		float4 *dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		dst[0].x = v0x2xr;
		dst[0].y = v0x2xi;
		dst[0].z = v0x2yr;
		dst[0].w = v0x2yi;
		dst[1].x = v0y2xr;
		dst[1].y = v0y2xi;
		dst[1].z = v0y2yr;
		dst[1].w = v0y2yi;

//		baseline ++;
//		dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));
		dst += nrChannels;

		dst[0].x = v1x2xr;
		dst[0].y = v1x2xi;
		dst[0].z = v1x2yr;
		dst[0].w = v1x2yi;
		dst[1].x = v1y2xr;
		dst[1].y = v1y2xi;
		dst[1].z = v1y2yr;
		dst[1].w = v1y2yi;

		baseline += stat2 + 1;
		dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		dst[0].x = v0x3xr;
		dst[0].y = v0x3xi;
		dst[0].z = v0x3yr;
		dst[0].w = v0x3yi;
		dst[1].x = v0y3xr;
		dst[1].y = v0y3xi;
		dst[1].z = v0y3yr;
		dst[1].w = v0y3yi;

//		baseline ++;
//		dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));
		dst += nrChannels;

		dst[0].x = v1x3xr;
		dst[0].y = v1x3xi;
		dst[0].z = v1x3yr;
		dst[0].w = v1x3yi;
		dst[1].x = v1y3xr;
		dst[1].y = v1y3xi;
		dst[1].z = v1y3yr;
		dst[1].w = v1y3yi;
	    }
	}
    }
}
