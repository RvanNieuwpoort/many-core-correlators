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
__global__ void correlate_4x3(float *devVisibilities, unsigned nrTimes, unsigned nrTimesWidth, 
			      unsigned nrStations, unsigned nrChannels, unsigned nrCells,
			      unsigned nrBlocks, unsigned nrThreads, unsigned loopCount)
{
    unsigned myBlock   = blockIdx.x + gridDim.x * blockIdx.y;
    unsigned myThread  = threadIdx.x + blockDim.x * threadIdx.y;

    for (unsigned channel = myBlock; channel < nrChannels; channel += nrBlocks) {
	PREFETCH();
	for (unsigned cell = myThread; cell < loopCount; cell += nrThreads) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat4 = cellToStatY[cell];

#if !USE_ALTERNATE_MEMORY_LAYOUT
	    unsigned index0 = SAMPLE_INDEX(stat0, channel, 0, 0) / 2;
	    unsigned index4 = SAMPLE_INDEX(stat4, channel, 0, 0) / 2;
#else
	    unsigned index = SAMPLE_INDEX(0, channel, 0, 0) / 2;
#endif
	    float v0x4xr = 0, v0x4xi = 0;
	    float v0x4yr = 0, v0x4yi = 0;
	    float v0y4xr = 0, v0y4xi = 0;
	    float v0y4yr = 0, v0y4yi = 0;
	    float v1x4xr = 0, v1x4xi = 0;
	    float v1x4yr = 0, v1x4yi = 0;
	    float v1y4xr = 0, v1y4xi = 0;
	    float v1y4yr = 0, v1y4yi = 0;
	    float v2x4xr = 0, v2x4xi = 0;
	    float v2x4yr = 0, v2x4yi = 0;
	    float v2y4xr = 0, v2y4xi = 0;
	    float v2y4yr = 0, v2y4yi = 0;
	    float v3x4xr = 0, v3x4xi = 0;
	    float v3x4yr = 0, v3x4yi = 0;
	    float v3y4xr = 0, v3y4xi = 0;
	    float v3y4yr = 0, v3y4yi = 0;

	    float v0x5xr = 0, v0x5xi = 0;
	    float v0x5yr = 0, v0x5yi = 0;
	    float v0y5xr = 0, v0y5xi = 0;
	    float v0y5yr = 0, v0y5yi = 0;
	    float v1x5xr = 0, v1x5xi = 0;
	    float v1x5yr = 0, v1x5yi = 0;
	    float v1y5xr = 0, v1y5xi = 0;
	    float v1y5yr = 0, v1y5yi = 0;
	    float v2x5xr = 0, v2x5xi = 0;
	    float v2x5yr = 0, v2x5yi = 0;
	    float v2y5xr = 0, v2y5xi = 0;
	    float v2y5yr = 0, v2y5yi = 0;
	    float v3x5xr = 0, v3x5xi = 0;
	    float v3x5yr = 0, v3x5yi = 0;
	    float v3y5xr = 0, v3y5xi = 0;
	    float v3y5yr = 0, v3y5yi = 0;

	    float v0x6xr = 0, v0x6xi = 0;
	    float v0x6yr = 0, v0x6yi = 0;
	    float v0y6xr = 0, v0y6xi = 0;
	    float v0y6yr = 0, v0y6yi = 0;
	    float v1x6xr = 0, v1x6xi = 0;
	    float v1x6yr = 0, v1x6yi = 0;
	    float v1y6xr = 0, v1y6xi = 0;
	    float v1y6yr = 0, v1y6yi = 0;
	    float v2x6xr = 0, v2x6xi = 0;
	    float v2x6yr = 0, v2x6yi = 0;
	    float v2y6xr = 0, v2y6xi = 0;
	    float v2y6yr = 0, v2y6yi = 0;
	    float v3x6xr = 0, v3x6xi = 0;
	    float v3x6yr = 0, v3x6yi = 0;
	    float v3y6xr = 0, v3y6xi = 0;
	    float v3y6yr = 0, v3y6yi = 0;

	    for (unsigned time = 0; time < nrTimes; time ++) {
		float4 sample0 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample1 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample2 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample3 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample4 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample5 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample6 = {0.0f, 0.0f, 0.0f, 0.0f};

#if !DO_LOADS
		if(nrTimes > 100000) { // Dummy condition that always evaluates to false; trick the compiler
#endif
#if USE_ALTERNATE_MEMORY_LAYOUT
		    sample0 = tex1Dfetch(tex0, index + stat0+0);
		    sample1 = tex1Dfetch(tex0, index + stat0+1);
		    sample2 = tex1Dfetch(tex0, index + stat0+2);
		    sample3 = tex1Dfetch(tex0, index + stat0+3);

		    sample4 = tex1Dfetch(tex0, index + stat4+0);
		    sample5 = tex1Dfetch(tex0, index + stat4+1);
		    sample6 = tex1Dfetch(tex0, index + stat4+2);
		    index += nrStations; // increase time
#else
		    sample0 = tex1Dfetch(tex0, index0 + time + 0 * nrTimesWidth);
		    sample1 = tex1Dfetch(tex0, index0 + time + 1 * nrTimesWidth);
		    sample2 = tex1Dfetch(tex0, index0 + time + 2 * nrTimesWidth);
		    sample3 = tex1Dfetch(tex0, index0 + time + 3 * nrTimesWidth);

		    sample4 = tex1Dfetch(tex0, index4 + time + 0 * nrTimesWidth);
		    sample5 = tex1Dfetch(tex0, index4 + time + 1 * nrTimesWidth);
		    sample6 = tex1Dfetch(tex0, index4 + time + 2 * nrTimesWidth);
#endif // USE_ALTERNATE_MEMORY_LAYOUT

#if !DO_LOADS
		}
#endif

#if SYNC_THREADS
		__syncthreads();
#endif

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

		v1x4xr  += sample1.x * sample4.x;
		v1x4xi  += sample1.y * sample4.x;
		v1x4yr  += sample1.x * sample4.z;
		v1x4yi  += sample1.y * sample4.z;
		v1y4xr  += sample1.z * sample4.x;
		v1y4xi  += sample1.w * sample4.x;
		v1y4yr  += sample1.z * sample4.z;
		v1y4yi  += sample1.w * sample4.z;
		v1x4xr  += sample1.y * sample4.y;
		v1x4xi  -= sample1.x * sample4.y;
		v1x4yr  += sample1.y * sample4.w;
		v1x4yi  -= sample1.x * sample4.w;
		v1y4xr  += sample1.w * sample4.y;
		v1y4xi  -= sample1.z * sample4.y;
		v1y4yr  += sample1.w * sample4.w;
		v1y4yi  -= sample1.z * sample4.w;

		v2x4xr  += sample2.x * sample4.x;
		v2x4xi  += sample2.y * sample4.x;
		v2x4yr  += sample2.x * sample4.z;
		v2x4yi  += sample2.y * sample4.z;
		v2y4xr  += sample2.z * sample4.x;
		v2y4xi  += sample2.w * sample4.x;
		v2y4yr  += sample2.z * sample4.z;
		v2y4yi  += sample2.w * sample4.z;
		v2x4xr  += sample2.y * sample4.y;
		v2x4xi  -= sample2.x * sample4.y;
		v2x4yr  += sample2.y * sample4.w;
		v2x4yi  -= sample2.x * sample4.w;
		v2y4xr  += sample2.w * sample4.y;
		v2y4xi  -= sample2.z * sample4.y;
		v2y4yr  += sample2.w * sample4.w;
		v2y4yi  -= sample2.z * sample4.w;

		v3x4xr  += sample3.x * sample4.x;
		v3x4xi  += sample3.y * sample4.x;
		v3x4yr  += sample3.x * sample4.z;
		v3x4yi  += sample3.y * sample4.z;
		v3y4xr  += sample3.z * sample4.x;
		v3y4xi  += sample3.w * sample4.x;
		v3y4yr  += sample3.z * sample4.z;
		v3y4yi  += sample3.w * sample4.z;
		v3x4xr  += sample3.y * sample4.y;
		v3x4xi  -= sample3.x * sample4.y;
		v3x4yr  += sample3.y * sample4.w;
		v3x4yi  -= sample3.x * sample4.w;
		v3y4xr  += sample3.w * sample4.y;
		v3y4xi  -= sample3.z * sample4.y;
		v3y4yr  += sample3.w * sample4.w;
		v3y4yi  -= sample3.z * sample4.w;



		v0x5xr  += sample0.x * sample5.x;
		v0x5xi  += sample0.y * sample5.x;
		v0x5yr  += sample0.x * sample5.z;
		v0x5yi  += sample0.y * sample5.z;
		v0y5xr  += sample0.z * sample5.x;
		v0y5xi  += sample0.w * sample5.x;
		v0y5yr  += sample0.z * sample5.z;
		v0y5yi  += sample0.w * sample5.z;
		v0x5xr  += sample0.y * sample5.y;
		v0x5xi  -= sample0.x * sample5.y;
		v0x5yr  += sample0.y * sample5.w;
		v0x5yi  -= sample0.x * sample5.w;
		v0y5xr  += sample0.w * sample5.y;
		v0y5xi  -= sample0.z * sample5.y;
		v0y5yr  += sample0.w * sample5.w;
		v0y5yi  -= sample0.z * sample5.w;

		v1x5xr  += sample1.x * sample5.x;
		v1x5xi  += sample1.y * sample5.x;
		v1x5yr  += sample1.x * sample5.z;
		v1x5yi  += sample1.y * sample5.z;
		v1y5xr  += sample1.z * sample5.x;
		v1y5xi  += sample1.w * sample5.x;
		v1y5yr  += sample1.z * sample5.z;
		v1y5yi  += sample1.w * sample5.z;
		v1x5xr  += sample1.y * sample5.y;
		v1x5xi  -= sample1.x * sample5.y;
		v1x5yr  += sample1.y * sample5.w;
		v1x5yi  -= sample1.x * sample5.w;
		v1y5xr  += sample1.w * sample5.y;
		v1y5xi  -= sample1.z * sample5.y;
		v1y5yr  += sample1.w * sample5.w;
		v1y5yi  -= sample1.z * sample5.w;

		v2x5xr  += sample2.x * sample5.x;
		v2x5xi  += sample2.y * sample5.x;
		v2x5yr  += sample2.x * sample5.z;
		v2x5yi  += sample2.y * sample5.z;
		v2y5xr  += sample2.z * sample5.x;
		v2y5xi  += sample2.w * sample5.x;
		v2y5yr  += sample2.z * sample5.z;
		v2y5yi  += sample2.w * sample5.z;
		v2x5xr  += sample2.y * sample5.y;
		v2x5xi  -= sample2.x * sample5.y;
		v2x5yr  += sample2.y * sample5.w;
		v2x5yi  -= sample2.x * sample5.w;
		v2y5xr  += sample2.w * sample5.y;
		v2y5xi  -= sample2.z * sample5.y;
		v2y5yr  += sample2.w * sample5.w;
		v2y5yi  -= sample2.z * sample5.w;

		v3x5xr  += sample3.x * sample5.x;
		v3x5xi  += sample3.y * sample5.x;
		v3x5yr  += sample3.x * sample5.z;
		v3x5yi  += sample3.y * sample5.z;
		v3y5xr  += sample3.z * sample5.x;
		v3y5xi  += sample3.w * sample5.x;
		v3y5yr  += sample3.z * sample5.z;
		v3y5yi  += sample3.w * sample5.z;
		v3x5xr  += sample3.y * sample5.y;
		v3x5xi  -= sample3.x * sample5.y;
		v3x5yr  += sample3.y * sample5.w;
		v3x5yi  -= sample3.x * sample5.w;
		v3y5xr  += sample3.w * sample5.y;
		v3y5xi  -= sample3.z * sample5.y;
		v3y5yr  += sample3.w * sample5.w;
		v3y5yi  -= sample3.z * sample5.w;



		v0x6xr  += sample0.x * sample6.x;
		v0x6xi  += sample0.y * sample6.x;
		v0x6yr  += sample0.x * sample6.z;
		v0x6yi  += sample0.y * sample6.z;
		v0y6xr  += sample0.z * sample6.x;
		v0y6xi  += sample0.w * sample6.x;
		v0y6yr  += sample0.z * sample6.z;
		v0y6yi  += sample0.w * sample6.z;
		v0x6xr  += sample0.y * sample6.y;
		v0x6xi  -= sample0.x * sample6.y;
		v0x6yr  += sample0.y * sample6.w;
		v0x6yi  -= sample0.x * sample6.w;
		v0y6xr  += sample0.w * sample6.y;
		v0y6xi  -= sample0.z * sample6.y;
		v0y6yr  += sample0.w * sample6.w;
		v0y6yi  -= sample0.z * sample6.w;

		v1x6xr  += sample1.x * sample6.x;
		v1x6xi  += sample1.y * sample6.x;
		v1x6yr  += sample1.x * sample6.z;
		v1x6yi  += sample1.y * sample6.z;
		v1y6xr  += sample1.z * sample6.x;
		v1y6xi  += sample1.w * sample6.x;
		v1y6yr  += sample1.z * sample6.z;
		v1y6yi  += sample1.w * sample6.z;
		v1x6xr  += sample1.y * sample6.y;
		v1x6xi  -= sample1.x * sample6.y;
		v1x6yr  += sample1.y * sample6.w;
		v1x6yi  -= sample1.x * sample6.w;
		v1y6xr  += sample1.w * sample6.y;
		v1y6xi  -= sample1.z * sample6.y;
		v1y6yr  += sample1.w * sample6.w;
		v1y6yi  -= sample1.z * sample6.w;

		v2x6xr  += sample2.x * sample6.x;
		v2x6xi  += sample2.y * sample6.x;
		v2x6yr  += sample2.x * sample6.z;
		v2x6yi  += sample2.y * sample6.z;
		v2y6xr  += sample2.z * sample6.x;
		v2y6xi  += sample2.w * sample6.x;
		v2y6yr  += sample2.z * sample6.z;
		v2y6yi  += sample2.w * sample6.z;
		v2x6xr  += sample2.y * sample6.y;
		v2x6xi  -= sample2.x * sample6.y;
		v2x6yr  += sample2.y * sample6.w;
		v2x6yi  -= sample2.x * sample6.w;
		v2y6xr  += sample2.w * sample6.y;
		v2y6xi  -= sample2.z * sample6.y;
		v2y6yr  += sample2.w * sample6.w;
		v2y6yi  -= sample2.z * sample6.w;

		v3x6xr  += sample3.x * sample6.x;
		v3x6xi  += sample3.y * sample6.x;
		v3x6yr  += sample3.x * sample6.z;
		v3x6yi  += sample3.y * sample6.z;
		v3y6xr  += sample3.z * sample6.x;
		v3y6xi  += sample3.w * sample6.x;
		v3y6yr  += sample3.z * sample6.z;
		v3y6yi  += sample3.w * sample6.z;
		v3x6xr  += sample3.y * sample6.y;
		v3x6xi  -= sample3.x * sample6.y;
		v3x6yr  += sample3.y * sample6.w;
		v3x6yi  -= sample3.x * sample6.w;
		v3y6xr  += sample3.w * sample6.y;
		v3y6xi  -= sample3.z * sample6.y;
		v3y6yr  += sample3.w * sample6.w;
		v3y6yi  -= sample3.z * sample6.w;
	    }

#if !DO_STORES
	    if(nrTimes > 100000) { // Dummy condition that always evaluates to false; trick the compiler
#endif
		if (cell < nrCells) {
		    unsigned baseline = BASELINE(stat0, stat4);
		    float4 *dst = (float4*) (devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x4xr;
		    dst[0].y = v0x4xi;
		    dst[0].z = v0x4yr;
		    dst[0].w = v0x4yi;
		    dst[1].x = v0y4xr;
		    dst[1].y = v0y4xi;
		    dst[1].z = v0y4yr;
		    dst[1].w = v0y4yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v1x4xr;
		    dst[0].y = v1x4xi;
		    dst[0].z = v1x4yr;
		    dst[0].w = v1x4yi;
		    dst[1].x = v1y4xr;
		    dst[1].y = v1y4xi;
		    dst[1].z = v1y4yr;
		    dst[1].w = v1y4yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v2x4xr;
		    dst[0].y = v2x4xi;
		    dst[0].z = v2x4yr;
		    dst[0].w = v2x4yi;
		    dst[1].x = v2y4xr;
		    dst[1].y = v2y4xi;
		    dst[1].z = v2y4yr;
		    dst[1].w = v2y4yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v3x4xr;
		    dst[0].y = v3x4xi;
		    dst[0].z = v3x4yr;
		    dst[0].w = v3x4yi;
		    dst[1].x = v3y4xr;
		    dst[1].y = v3y4xi;
		    dst[1].z = v3y4yr;
		    dst[1].w = v3y4yi;


		    baseline = BASELINE(stat0, stat4+1);
		    dst = (float4*) (devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x5xr;
		    dst[0].y = v0x5xi;
		    dst[0].z = v0x5yr;
		    dst[0].w = v0x5yi;
		    dst[1].x = v0y5xr;
		    dst[1].y = v0y5xi;
		    dst[1].z = v0y5yr;
		    dst[1].w = v0y5yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v1x5xr;
		    dst[0].y = v1x5xi;
		    dst[0].z = v1x5yr;
		    dst[0].w = v1x5yi;
		    dst[1].x = v1y5xr;
		    dst[1].y = v1y5xi;
		    dst[1].z = v1y5yr;
		    dst[1].w = v1y5yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v2x5xr;
		    dst[0].y = v2x5xi;
		    dst[0].z = v2x5yr;
		    dst[0].w = v2x5yi;
		    dst[1].x = v2y5xr;
		    dst[1].y = v2y5xi;
		    dst[1].z = v2y5yr;
		    dst[1].w = v2y5yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v3x5xr;
		    dst[0].y = v3x5xi;
		    dst[0].z = v3x5yr;
		    dst[0].w = v3x5yi;
		    dst[1].x = v3y5xr;
		    dst[1].y = v3y5xi;
		    dst[1].z = v3y5yr;
		    dst[1].w = v3y5yi;


		    baseline = BASELINE(stat0, stat4+2);
		    dst = (float4*) (devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x6xr;
		    dst[0].y = v0x6xi;
		    dst[0].z = v0x6yr;
		    dst[0].w = v0x6yi;
		    dst[1].x = v0y6xr;
		    dst[1].y = v0y6xi;
		    dst[1].z = v0y6yr;
		    dst[1].w = v0y6yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v1x6xr;
		    dst[0].y = v1x6xi;
		    dst[0].z = v1x6yr;
		    dst[0].w = v1x6yi;
		    dst[1].x = v1y6xr;
		    dst[1].y = v1y6xi;
		    dst[1].z = v1y6yr;
		    dst[1].w = v1y6yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v2x6xr;
		    dst[0].y = v2x6xi;
		    dst[0].z = v2x6yr;
		    dst[0].w = v2x6yi;
		    dst[1].x = v2y6xr;
		    dst[1].y = v2y6xi;
		    dst[1].z = v2y6yr;
		    dst[1].w = v2y6yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v3x6xr;
		    dst[0].y = v3x6xi;
		    dst[0].z = v3x6yr;
		    dst[0].w = v3x6yi;
		    dst[1].x = v3y6xr;
		    dst[1].y = v3y6xi;
		    dst[1].z = v3y6yr;
		    dst[1].w = v3y6yi;
		}
#if !DO_STORES
	    }
#endif
	}
    }
}
