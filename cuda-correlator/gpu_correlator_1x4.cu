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
__global__ void correlate_1x4(float *devVisibilities, unsigned nrTimes, unsigned nrTimesWidth, 
			      unsigned nrStations, unsigned nrChannels, unsigned nrCells,
			      unsigned nrBlocks, unsigned nrThreads, unsigned loopCount)
{
    unsigned myBlock   = blockIdx.x + gridDim.x * blockIdx.y;
    unsigned myThread  = threadIdx.x + blockDim.x * threadIdx.y;

    for (unsigned channel = myBlock; channel < nrChannels; channel += nrBlocks) {
	PREFETCH();

#if USE_ALTERNATE_MEMORY_LAYOUT
	unsigned index = SAMPLE_INDEX(0, channel, 0, 0) / 2;
#endif
	for (unsigned cell = myThread; cell < loopCount; cell += nrThreads) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat1 = cellToStatY[cell];

#if !USE_ALTERNATE_MEMORY_LAYOUT
	    unsigned index0 = SAMPLE_INDEX(stat0, channel, 0, 0) / 2;
	    unsigned index1 = SAMPLE_INDEX(stat1, channel, 0, 0) / 2;
#endif
	    float v0x1xr = 0, v0x1xi = 0;
	    float v0x1yr = 0, v0x1yi = 0;
	    float v0y1xr = 0, v0y1xi = 0;
	    float v0y1yr = 0, v0y1yi = 0;
	    float v0x2xr = 0, v0x2xi = 0;
	    float v0x2yr = 0, v0x2yi = 0;
	    float v0y2xr = 0, v0y2xi = 0;
	    float v0y2yr = 0, v0y2yi = 0;
	    float v0x3xr = 0, v0x3xi = 0;
	    float v0x3yr = 0, v0x3yi = 0;
	    float v0y3xr = 0, v0y3xi = 0;
	    float v0y3yr = 0, v0y3yi = 0;
	    float v0x4xr = 0, v0x4xi = 0;
	    float v0x4yr = 0, v0x4yi = 0;
	    float v0y4xr = 0, v0y4xi = 0;
	    float v0y4yr = 0, v0y4yi = 0;

	    for (unsigned time = 0; time < nrTimes; time ++) {
		float4 sample0 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample1 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample2 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample3 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample4 = {0.0f, 0.0f, 0.0f, 0.0f};

#if !DO_LOADS
		// A dummy condition that always evaluates to false.
		// Needed to trick the compiler
		if(nrTimes > 100000) {
#endif

#if USE_ALTERNATE_MEMORY_LAYOUT
		    sample0 = tex1Dfetch(tex0, index + time * nrStations + stat0+0);

		    sample1 = tex1Dfetch(tex0, index + time * nrStations + stat1+0);
		    sample2 = tex1Dfetch(tex0, index + time * nrStations + stat1+1);
		    sample3 = tex1Dfetch(tex0, index + time * nrStations + stat1+2);
		    sample4 = tex1Dfetch(tex0, index + time * nrStations + stat1+3);
#else
		    // Just load the samples, even though cell >= nrCells.
		    // We are loading junk in that case, but still this is faster than testing.
		    sample0 = tex1Dfetch(tex0, index0 + time + 0 * nrTimesWidth);

		    sample1 = tex1Dfetch(tex0, index1 + time + 0 * nrTimesWidth);
		    sample2 = tex1Dfetch(tex0, index1 + time + 1 * nrTimesWidth);
		    sample3 = tex1Dfetch(tex0, index1 + time + 2 * nrTimesWidth);
		    sample4 = tex1Dfetch(tex0, index1 + time + 3 * nrTimesWidth);
#endif // USE_ALTERNATE_MEMORY_LAYOUT

#if !DO_LOADS
		}
#endif

#if SYNC_THREADS
		__syncthreads();
#endif

		v0x1xr  += sample0.x * sample1.x;
		v0x1xi  += sample0.y * sample1.x;
		v0x1yr  += sample0.x * sample1.z;
		v0x1yi  += sample0.y * sample1.z;
		v0y1xr  += sample0.z * sample1.x;
		v0y1xi  += sample0.w * sample1.x;
		v0y1yr  += sample0.z * sample1.z;
		v0y1yi  += sample0.w * sample1.z;

		v0x1xr  += sample0.y * sample1.y;
		v0x1xi  -= sample0.x * sample1.y;
		v0x1yr  += sample0.y * sample1.w;
		v0x1yi  -= sample0.x * sample1.w;
		v0y1xr  += sample0.w * sample1.y;
		v0y1xi  -= sample0.z * sample1.y;
		v0y1yr  += sample0.w * sample1.w;
		v0y1yi  -= sample0.z * sample1.w;


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

		v0x4xr  += sample0.x * sample4.x;
		v0x4xi  += sample0.y * sample4.x;
		v0x4yr  += sample0.x * sample4.z;
		v0x4yi  += sample0.y * sample4.z;
		v0y4xr  += sample0.z * sample4.x;
		v0y4xi  += sample0.w * sample4.x;
		v0y4yr  += sample0.z * sample4.z;
		v0y4yi  += sample0.w * sample4.z;

		v0x4xr  += sample0.y * sample4.y;
		v0x4xi  -= sample0.x * sample4.y;
		v0x4yr  += sample0.y * sample4.w;
		v0x4yi  -= sample0.x * sample4.w;
		v0y4xr  += sample0.w * sample4.y;
		v0y4xi  -= sample0.z * sample4.y;
		v0y4yr  += sample0.w * sample4.w;
		v0y4yi  -= sample0.z * sample4.w;
	    }

#if !DO_STORES
	    if(nrTimes > 100000) { // Dummy condition that always evaluates to false; trick the compiler
#endif
		if (cell < nrCells) {
		    unsigned baseline = BASELINE(stat0, stat1);
		    float4 *dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x1xr;
		    dst[0].y = v0x1xi;
		    dst[0].z = v0x1yr;
		    dst[0].w = v0x1yi;
		    dst[1].x = v0y1xr;
		    dst[1].y = v0y1xi;
		    dst[1].z = v0y1yr;
		    dst[1].w = v0y1yi;

		    baseline = BASELINE(stat0, stat1+1);
		    dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x2xr;
		    dst[0].y = v0x2xi;
		    dst[0].z = v0x2yr;
		    dst[0].w = v0x2yi;
		    dst[1].x = v0y2xr;
		    dst[1].y = v0y2xi;
		    dst[1].z = v0y2yr;
		    dst[1].w = v0y2yi;

		    baseline = BASELINE(stat0, stat1+2);
		    dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x3xr;
		    dst[0].y = v0x3xi;
		    dst[0].z = v0x3yr;
		    dst[0].w = v0x3yi;
		    dst[1].x = v0y3xr;
		    dst[1].y = v0y3xi;
		    dst[1].z = v0y3yr;
		    dst[1].w = v0y3yi;

		    baseline = BASELINE(stat0, stat1+3);
		    dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x4xr;
		    dst[0].y = v0x4xi;
		    dst[0].z = v0x4yr;
		    dst[0].w = v0x4yi;
		    dst[1].x = v0y4xr;
		    dst[1].y = v0y4xi;
		    dst[1].z = v0y4yr;
		    dst[1].w = v0y4yi;
		}
#if !DO_STORES
	    }
#endif
	}
    }
}
