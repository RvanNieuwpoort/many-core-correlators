__global__ void correlate_3x2(float* devSamplesIn, float* devSamples, float *devVisibilities, 
			      unsigned nrTimes, unsigned nrTimesWidth, 
			      unsigned nrStations, unsigned nrChannels, unsigned nrCells,
			      unsigned nrBlocks, unsigned nrThreads, unsigned loopCount, unsigned nrMissedCells, unsigned missingLoopCount)
{
    unsigned myBlock   = blockIdx.x + gridDim.x * blockIdx.y;
    unsigned myThread  = threadIdx.x + blockDim.x * threadIdx.y;

#if USE_ALTERNATE_MEMORY_LAYOUT && DO_MEMORY_CONVERSION && CONVERT_ON_GPU
    for(unsigned channel=myThread; channel<nrChannels; channel+=nrBlocks) {
	unsigned destIndex = ALTERNATE_SAMPLE_INDEX(0, channel, 0, 0);
	for(unsigned time=0; time<nrTimes; time+=nrThreads) {
	    unsigned srcIndex = NORMAL_SAMPLE_INDEX(0, channel, time, 0);
	    for(unsigned station=0; station<nrStations; station++) {
		float4 sample = tex1Dfetch(unconvertedTexture, srcIndex);
		*((float4*)(&devSamples[destIndex])) = sample;
		destIndex += 2;
		srcIndex += nrTimesWidth * 2;
	    }
	}
    }
    __syncthreads();
#endif

    for (unsigned channel = myBlock; channel < nrChannels; channel += nrBlocks) {
	PREFETCH();
	for (unsigned cell = myThread; cell < loopCount; cell += nrThreads) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat3 = cellToStatY[cell];

#if !USE_ALTERNATE_MEMORY_LAYOUT
	    unsigned index0 = SAMPLE_INDEX(stat0, channel, 0, 0) / 2;
	    unsigned index3 = SAMPLE_INDEX(stat3, channel, 0, 0) / 2;
#else
	    unsigned index = SAMPLE_INDEX(0, channel, 0, 0) / 2;
#endif
	    // 48 registers
	    float v0x3xr = 0.0f, v0x3xi = 0.0f;
	    float v0x3yr = 0.0f, v0x3yi = 0.0f;
	    float v0y3xr = 0.0f, v0y3xi = 0.0f;
	    float v0y3yr = 0.0f, v0y3yi = 0.0f;
	    float v1x3xr = 0.0f, v1x3xi = 0.0f;
	    float v1x3yr = 0.0f, v1x3yi = 0.0f;
	    float v1y3xr = 0.0f, v1y3xi = 0.0f;
	    float v1y3yr = 0.0f, v1y3yi = 0.0f;
	    float v2x3xr = 0.0f, v2x3xi = 0.0f;
	    float v2x3yr = 0.0f, v2x3yi = 0.0f;
	    float v2y3xr = 0.0f, v2y3xi = 0.0f;
	    float v2y3yr = 0.0f, v2y3yi = 0.0f;
	    float v0x4xr = 0.0f, v0x4xi = 0.0f;
	    float v0x4yr = 0.0f, v0x4yi = 0.0f;
	    float v0y4xr = 0.0f, v0y4xi = 0.0f;
	    float v0y4yr = 0.0f, v0y4yi = 0.0f;
	    float v1x4xr = 0.0f, v1x4xi = 0.0f;
	    float v1x4yr = 0.0f, v1x4yi = 0.0f;
	    float v1y4xr = 0.0f, v1y4xi = 0.0f;
	    float v1y4yr = 0.0f, v1y4yi = 0.0f;
	    float v2x4xr = 0.0f, v2x4xi = 0.0f;
	    float v2x4yr = 0.0f, v2x4yi = 0.0f;
	    float v2y4xr = 0.0f, v2y4xi = 0.0f;
	    float v2y4yr = 0.0f, v2y4yi = 0.0f;

	    for (unsigned time = 0; time < nrTimes; time ++) {
		// 20 registers, 68 total
		float4 sample0 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample1 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample2 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample3 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample4 = {0.0f, 0.0f, 0.0f, 0.0f};

#if !DO_LOADS
		if(nrTimes > 100000) { // Dummy condition that always evaluates to false; trick the compiler
#endif
#if USE_ALTERNATE_MEMORY_LAYOUT
		    sample0 = tex1Dfetch(tex0, index + time * nrStations + stat0+0);
		    sample1 = tex1Dfetch(tex0, index + time * nrStations + stat0+1);
		    sample2 = tex1Dfetch(tex0, index + time * nrStations + stat0+2);

		    sample3 = tex1Dfetch(tex0, index + time * nrStations + stat3+0);
		    sample4 = tex1Dfetch(tex0, index + time * nrStations + stat3+1);
#else
		    // 20 registers, 68 total
		    sample0 = tex1Dfetch(tex0, index0 + time + 0 * nrTimesWidth);
		    sample1 = tex1Dfetch(tex0, index0 + time + 1 * nrTimesWidth);
		    sample2 = tex1Dfetch(tex0, index0 + time + 2 * nrTimesWidth);

		    sample3 = tex1Dfetch(tex0, index3 + time + 0 * nrTimesWidth);
		    sample4 = tex1Dfetch(tex0, index3 + time + 1 * nrTimesWidth);
#endif // USE_ALTERNATE_MEMORY_LAYOUT
#if !DO_LOADS
		}
#endif

#if SYNC_THREADS
//		if ((time & 63) == 0)
		__syncthreads();
#endif

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

		v2x3xr  += sample2.x * sample3.x;
		v2x3xi  += sample2.y * sample3.x;
		v2x3yr  += sample2.x * sample3.z;
		v2x3yi  += sample2.y * sample3.z;
		v2y3xr  += sample2.z * sample3.x;
		v2y3xi  += sample2.w * sample3.x;
		v2y3yr  += sample2.z * sample3.z;
		v2y3yi  += sample2.w * sample3.z;
		v2x3xr  += sample2.y * sample3.y;
		v2x3xi  -= sample2.x * sample3.y;
		v2x3yr  += sample2.y * sample3.w;
		v2x3yi  -= sample2.x * sample3.w;
		v2y3xr  += sample2.w * sample3.y;
		v2y3xi  -= sample2.z * sample3.y;
		v2y3yr  += sample2.w * sample3.w;
		v2y3yi  -= sample2.z * sample3.w;

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
	    }

#if !DO_STORES
	    if(nrTimes > 100000) { // Dummy condition that always evaluates to false; trick the compiler
#endif
		if (cell < nrCells) {

		    unsigned baseline = BASELINE(stat0, stat3);
		    float4 *dst = (float4*) (devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x3xr;
		    dst[0].y = v0x3xi;
		    dst[0].z = v0x3yr;
		    dst[0].w = v0x3yi;
		    dst[1].x = v0y3xr;
		    dst[1].y = v0y3xi;
		    dst[1].z = v0y3yr;
		    dst[1].w = v0y3yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v1x3xr;
		    dst[0].y = v1x3xi;
		    dst[0].z = v1x3yr;
		    dst[0].w = v1x3yi;
		    dst[1].x = v1y3xr;
		    dst[1].y = v1y3xi;
		    dst[1].z = v1y3yr;
		    dst[1].w = v1y3yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v2x3xr;
		    dst[0].y = v2x3xi;
		    dst[0].z = v2x3yr;
		    dst[0].w = v2x3yi;
		    dst[1].x = v2y3xr;
		    dst[1].y = v2y3xi;
		    dst[1].z = v2y3yr;
		    dst[1].w = v2y3yi;

		    baseline = BASELINE(stat0, stat3+1);
		    dst = (float4*) (devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

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
		}
#if !DO_STORES
	    }
#endif
	}
    }

#if CALCULATE_MISSING
    // calculate everything we missed by tiling
    for (unsigned channel = myBlock; channel < nrChannels; channel += nrBlocks) {
	for (unsigned cell = myThread; cell < nrMissedCells; cell += nrThreads) {
	    unsigned stat0 = missedStatX[cell];
	    unsigned stat1 = missedStatY[cell];
	    unsigned index0 = SAMPLE_INDEX(stat0, channel, 0, 0) / 2;
	    unsigned index1 = SAMPLE_INDEX(stat1, channel, 0, 0) / 2;
	    
	    float v0x1xr = 0.0f, v0x1xi = 0.0f;
	    float v0x1yr = 0.0f, v0x1yi = 0.0f;
	    float v0y1xr = 0.0f, v0y1xi = 0.0f;
	    float v0y1yr = 0.0f, v0y1yi = 0.0f;

	    for (unsigned time = 0; time < nrTimes; time ++) {
#if USE_ALTERNATE_MEMORY_LAYOUT
		float4 sample0 = tex1Dfetch(tex0, index0 + time * nrStations);
		float4 sample1 = tex1Dfetch(tex0, index1 + time * nrStations);
#else
		float4 sample0 = tex1Dfetch(tex0, index0 + time);
		float4 sample1 = tex1Dfetch(tex0, index1 + time);
#endif // USE_ALTERNATE_MEMORY_LAYOUT

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
	    }

//	    if(cell < nrMissedCells) {
		unsigned baseline = BASELINE(stat0, stat1);
		float4 *dst = (float4*) (devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		dst[0].x = v0x1xr;
		dst[0].y = v0x1xi;
		dst[0].z = v0x1yr;
		dst[0].w = v0x1yi;
		dst[1].x = v0y1xr;
		dst[1].y = v0y1xi;
		dst[1].z = v0y1yr;
		dst[1].w = v0y1yi;
//	    }
	}
    }
#endif // CALCULATE_MISSING
}
