typedef struct __align__(16) {
    unsigned channel[MAX_THREADS];
    unsigned cell[MAX_THREADS];
    unsigned stat0[MAX_THREADS];
    unsigned stat3[MAX_THREADS];
    unsigned index0[MAX_THREADS];
    unsigned index3[MAX_THREADS];
} spilled_t;

__shared__ spilled_t spilled;

typedef struct __align__(16) {
    unsigned index_table[MAX_STATIONS][MAX_CHANNELS];
    unsigned vis_table[MAX_STATIONS][MAX_STATIONS][MAX_CHANNELS];
} device_tables_t;


__device__ device_tables_t tables;

unsigned index_tableHost[MAX_STATIONS][MAX_CHANNELS];
unsigned vis_tableHost[MAX_STATIONS][MAX_STATIONS][MAX_CHANNELS];

void initTables()
{
    for(unsigned stat = 0; stat<nrStations; stat++) {
	for(unsigned channel = 0; channel < nrChannels; channel++) {
	    index_tableHost[stat][channel] = SAMPLE_INDEX(stat, channel, 0, 0) / 2;
	}
    }
    cudaMemcpyToSymbol(tables.index_table, index_tableHost, sizeof(tables.index_table));

    for(unsigned stat0 = 0; stat0 < nrStations; stat0++) {
	for(unsigned stat3 = 0; stat3 < nrStations; stat3++) {
	    for(unsigned channel = 0; channel < nrChannels; channel++) {
		unsigned baseline = BASELINE(stat0, stat3);
		vis_tableHost[stat0][stat3][channel] = 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0);
	    }
	}
    }
    cudaMemcpyToSymbol(tables.vis_table, vis_tableHost, sizeof(tables.vis_table));
}

__global__ void correlate_3x2_spill(float *devVisibilities, unsigned nrTimes, unsigned nrTimesWidth, unsigned nrStations, unsigned nrChannels, unsigned nrCells,
			      unsigned nrBlocks, unsigned nrThreads, unsigned loopCount) 
{
    // 3x2 implementation
    unsigned myBlock   = blockIdx.x + gridDim.x * blockIdx.y;
    unsigned myThread  = threadIdx.x + blockDim.x * threadIdx.y;

    for (spilled.channel[myThread] = myBlock; spilled.channel[myThread] < nrChannels; spilled.channel[myThread] += nrBlocks) {
	for (spilled.cell[myThread] = myThread; spilled.cell[myThread] < loopCount; spilled.cell[myThread] += nrThreads) {
	    spilled.stat0[myThread] = cellToStatX[spilled.cell[myThread]];
	    spilled.stat3[myThread] = cellToStatY[spilled.cell[myThread]];

	    unsigned index0 = tables.index_table[spilled.stat0[myThread]][spilled.channel[myThread]];
	    spilled.index3[myThread] = tables.index_table[spilled.stat3[myThread]][spilled.channel[myThread]];

	    // 48 registers
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

	    for (unsigned time = 0; time < nrTimes; time ++) {
		// 12 registers, 60 total
		float4 sample0 = tex1Dfetch(tex0, index0 + time);
		float4 sample3 = tex1Dfetch(tex0, spilled.index3[myThread] + time);
		float4 sample4 = tex1Dfetch(tex0, spilled.index3[myThread] + time + nrTimesWidth);

//		if ((time & 63) == 0)
#if SYNC_THREADS
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

		// really sample1
		sample0 = tex1Dfetch(tex0, index0 + time + nrTimesWidth);

		v1x3xr  += sample0.x * sample3.x;
		v1x3xi  += sample0.y * sample3.x;
		v1x3yr  += sample0.x * sample3.z;
		v1x3yi  += sample0.y * sample3.z;
		v1y3xr  += sample0.z * sample3.x;
		v1y3xi  += sample0.w * sample3.x;
		v1y3yr  += sample0.z * sample3.z;
		v1y3yi  += sample0.w * sample3.z;
		v1x3xr  += sample0.y * sample3.y;
		v1x3xi  -= sample0.x * sample3.y;
		v1x3yr  += sample0.y * sample3.w;
		v1x3yi  -= sample0.x * sample3.w;
		v1y3xr  += sample0.w * sample3.y;
		v1y3xi  -= sample0.z * sample3.y;
		v1y3yr  += sample0.w * sample3.w;
		v1y3yi  -= sample0.z * sample3.w;

		v1x4xr  += sample0.x * sample4.x;
		v1x4xi  += sample0.y * sample4.x;
		v1x4yr  += sample0.x * sample4.z;
		v1x4yi  += sample0.y * sample4.z;
		v1y4xr  += sample0.z * sample4.x;
		v1y4xi  += sample0.w * sample4.x;
		v1y4yr  += sample0.z * sample4.z;
		v1y4yi  += sample0.w * sample4.z;
		v1x4xr  += sample0.y * sample4.y;
		v1x4xi  -= sample0.x * sample4.y;
		v1x4yr  += sample0.y * sample4.w;
		v1x4yi  -= sample0.x * sample4.w;
		v1y4xr  += sample0.w * sample4.y;
		v1y4xi  -= sample0.z * sample4.y;
		v1y4yr  += sample0.w * sample4.w;
		v1y4yi  -= sample0.z * sample4.w;

		// really sample 2
		sample0 = tex1Dfetch(tex0, index0 + time + 2 * nrTimesWidth);
		v2x3xr  += sample0.x * sample3.x;
		v2x3xi  += sample0.y * sample3.x;
		v2x3yr  += sample0.x * sample3.z;
		v2x3yi  += sample0.y * sample3.z;
		v2y3xr  += sample0.z * sample3.x;
		v2y3xi  += sample0.w * sample3.x;
		v2y3yr  += sample0.z * sample3.z;
		v2y3yi  += sample0.w * sample3.z;
		v2x3xr  += sample0.y * sample3.y;
		v2x3xi  -= sample0.x * sample3.y;
		v2x3yr  += sample0.y * sample3.w;
		v2x3yi  -= sample0.x * sample3.w;
		v2y3xr  += sample0.w * sample3.y;
		v2y3xi  -= sample0.z * sample3.y;
		v2y3yr  += sample0.w * sample3.w;
		v2y3yi  -= sample0.z * sample3.w;

		v2x4xr  += sample0.x * sample4.x;
		v2x4xi  += sample0.y * sample4.x;
		v2x4yr  += sample0.x * sample4.z;
		v2x4yi  += sample0.y * sample4.z;
		v2y4xr  += sample0.z * sample4.x;
		v2y4xi  += sample0.w * sample4.x;
		v2y4yr  += sample0.z * sample4.z;
		v2y4yi  += sample0.w * sample4.z;
		v2x4xr  += sample0.y * sample4.y;
		v2x4xi  -= sample0.x * sample4.y;
		v2x4yr  += sample0.y * sample4.w;
		v2x4yi  -= sample0.x * sample4.w;
		v2y4xr  += sample0.w * sample4.y;
		v2y4xi  -= sample0.z * sample4.y;
		v2y4yr  += sample0.w * sample4.w;
		v2y4yi  -= sample0.z * sample4.w;
	    }

	    if (spilled.cell[myThread] < nrCells) {
		unsigned baseline = BASELINE(spilled.stat0[myThread], spilled.stat3[myThread]);
		float4 *dst = (float4*) (devVisibilities + 2 * VISIBILITIES_INDEX(baseline, spilled.channel[myThread], 0, 0));

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

		baseline = BASELINE(spilled.stat0[myThread], spilled.stat3[myThread]+1);
		dst = (float4*) (devVisibilities + 2 * VISIBILITIES_INDEX(baseline, spilled.channel[myThread], 0, 0));

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
	}
    }
}
