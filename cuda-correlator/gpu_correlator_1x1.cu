__global__ void correlate_1x1(float *devVisibilities, unsigned nrTimes, unsigned nrTimesWidth, 
			      unsigned nrStations, unsigned nrChannels,
			      unsigned nrCells, unsigned nrBlocks, unsigned nrThreads, unsigned loopCount)
{
    unsigned myBlock   = blockIdx.x + gridDim.x * blockIdx.y;
    unsigned myThread  = threadIdx.x + blockDim.x * threadIdx.y;

    for (unsigned channel = myBlock; channel < nrChannels; channel += nrBlocks) {	
	PREFETCH();

#if USE_ALTERNATE_MEMORY_LAYOUT
	unsigned index = SAMPLE_INDEX(0, channel, 0, 0) / 2;
#endif // USE_ALTERNATE_MEMORY_LAYOUT

	for (unsigned cell = myThread; cell < loopCount; cell += nrThreads) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat1 = cellToStatY[cell];

#if !USE_ALTERNATE_MEMORY_LAYOUT
	    unsigned index0 = SAMPLE_INDEX(stat0, channel, 0, 0) / 2;
	    unsigned index1 = SAMPLE_INDEX(stat1, channel, 0, 0) / 2;
#endif // USE_ALTERNATE_MEMORY_LAYOUT

	    float xxr = 0, xxi = 0, xyr = 0, xyi = 0, yxr = 0, yxi = 0, yyr = 0, yyi = 0;

	    for (unsigned time = 0; time < nrTimes; time ++) {
		float4 sample0;
		float4 sample1;

#if !DO_LOADS
		// A dummy condition that always evaluates to false.
		// Needed to trick the compiler
		if(nrTimes > 100000) {
#endif

#if USE_ALTERNATE_MEMORY_LAYOUT
		    sample0 = tex1Dfetch(tex0, index + stat0+0);
		    sample1 = tex1Dfetch(tex0, index + stat1+0);
		    index += nrStations; // increase time
#else
		    sample0 = tex1Dfetch(tex0, index0 + time);
		    sample1 = tex1Dfetch(tex0, index1 + time);
#endif // USE_ALTERNATE_MEMORY_LAYOUT

#if !DO_LOADS
		}
#endif

#if SYNC_THREADS
		__syncthreads();
#endif

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
