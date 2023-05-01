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
__global__ void correlate_1x1_sharedmem(float* devSamples, float *devVisibilities, unsigned nrTimes, 
					unsigned nrTimesWidth, 
					unsigned nrStations, unsigned nrChannels,
					unsigned nrCells, unsigned nrBlocks, unsigned nrThreads, 
					unsigned loopCount)
{
    unsigned myBlock   = blockIdx.x + gridDim.x * blockIdx.y;
    unsigned myThread  = threadIdx.x + blockDim.x * threadIdx.y;

    for (unsigned channel = myBlock; channel < nrChannels; channel += nrBlocks) {
	for (unsigned cell = myThread; cell < loopCount; cell += nrThreads) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat1 = cellToStatY[cell];

	    float xxr = 0, xxi = 0, xyr = 0, xyi = 0, yxr = 0, yxi = 0, yyr = 0, yyi = 0;

	    for (unsigned time = 0; time < nrTimes; time ++) {
		float4 sample0;
		float4 sample1;

		// Pull in all samples in one go.
		// This code is the same, regardless of the memory layout
		for(unsigned i=myThread; i<nrStations; i+=nrThreads) {
		    unsigned index = SAMPLE_INDEX(i, channel, time, 0) / 2;
		    samples[i] = tex1Dfetch(tex0, index);
		}

		// In this case, we have to do a syncthreads for correctness.
		__syncthreads();

		sample0 = samples[stat0];
		sample1 = samples[stat1];

		xxr += sample0.x * sample1.x;
		xxi += sample0.y * sample1.x;
		xyr += sample0.x * sample1.z;
		xyi += sample0.y * sample1.z;
		yxr += sample0.z * sample1.x;
		yxi += sample0.w * sample1.x;
		yyr += sample0.z * sample1.z;
		yyi += sample0.w * sample1.z;
		xxr += sample0.y * sample1.y;
		xxi -= sample0.x * sample1.y;
		xyr += sample0.y * sample1.w;
		xyi -= sample0.x * sample1.w;
		yxr += sample0.w * sample1.y;
		yxi -= sample0.z * sample1.y;
		yyr += sample0.w * sample1.w;
		yyi -= sample0.z * sample1.w;
	    }

#if !DO_STORES
	    if(nrTimes > 100000) { // Dummy condition that always evaluates to false; trick the compiler
#endif
		if (cell < nrCells) {
		    unsigned baseline = BASELINE(stat0, stat1);
		    float *dst = devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0);
		    dst[0] = xxr;
		    dst[1] = xxi;
		    dst[2] = xyr;
		    dst[3] = xyi;
		    dst[4] = yxr;
		    dst[5] = yxi;
		    dst[6] = yyr;
		    dst[7] = yyi;
		}
#if !DO_STORES
	    }
#endif
	}
    }
}
