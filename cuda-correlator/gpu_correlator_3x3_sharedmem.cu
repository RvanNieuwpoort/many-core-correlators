__global__ void correlate_3x3_sharedmem(float* devSamples, float *devVisibilities, 
					unsigned nrTimes, unsigned nrTimesWidth, 
					unsigned nrStations, unsigned nrChannels, unsigned nrCells,
					unsigned nrBlocks, unsigned nrThreads, unsigned loopCount)
{
    unsigned myBlock   = blockIdx.x + gridDim.x * blockIdx.y;
    unsigned myThread  = threadIdx.x + blockDim.x * threadIdx.y;

    for (unsigned channel = myBlock; channel < nrChannels; channel += nrBlocks) {
	for (unsigned cell = myThread; cell < loopCount; cell += nrThreads) {
	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat3 = cellToStatY[cell];

	    float v0x3xr = 0, v0x3xi = 0;
	    float v0x3yr = 0, v0x3yi = 0;
	    float v0y3xr = 0, v0y3xi = 0;
	    float v0y3yr = 0, v0y3yi = 0;
	    float v1x3xr = 0, v1x3xi = 0;
	    float v1x3yr = 0, v1x3yi = 0;
	    float v1y3xr = 0, v1y3xi = 0;
	    float v1y3yr = 0, v1y3yi = 0;
	    float v2x3xr = 0, v2x3xi = 0;
	    float v2x3yr = 0, v2x3yi = 0;
	    float v2y3xr = 0, v2y3xi = 0;
	    float v2y3yr = 0, v2y3yi = 0;
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

	    for (unsigned time = 0; time < nrTimes; time ++) {
		float4 sample0 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample1 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample2 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample3 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample4 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample5 = {0.0f, 0.0f, 0.0f, 0.0f};

#if 0
		// pull in all samples in one go
		for(unsigned i = myThread; i < nrStations; i += nrThreads) {
		    unsigned index = SAMPLE_INDEX(i, channel, time, 0) / 2;
//		    samples[i] = ((float4*) devSamples)[index];

		    samples[i] = tex1Dfetch(tex0, index);
		}
#endif


		// In this case, we have to do a syncthreads for correctness.
		__syncthreads();

#if !DO_LOADS
		if(nrTimes > 100000) { // Dummy condition that always evaluates to false; trick the compiler
#endif
		    sample0 = samples[stat0+0];
		    sample1 = samples[stat0+1];
		    sample2 = samples[stat0+2];
		    
		    sample3 = samples[stat3+0];
		    sample4 = samples[stat3+1];
		    sample5 = samples[stat3+2];
#if !DO_LOADS
		}
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


		    baseline = BASELINE(stat0, (stat3+1));
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

		    baseline = BASELINE(stat0, (stat3+2));
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
		}
#if !DO_STORES
	    }
#endif

	}
    }
}
