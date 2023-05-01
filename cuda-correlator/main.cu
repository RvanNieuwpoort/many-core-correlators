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
#include "gpu_correlator.h"
#include "timer.h"
#include "common.h"
#include "gpu_complex.h" // #include <complex> does not work in -deviceemu mode, so we include our own simple version
#include <cassert>
#include <cuda_runtime.h>
#include <cstdlib>
#include <iostream>

bool verify = false;
bool printResult = false;
unsigned nrStations	 = 64;
unsigned nrTimes	 = 768, nrTimesWidth = 770;
unsigned cudaDevice      = 0;
#if USE_STREAMING
unsigned nrChannels	 = 256; // 2048; // 512; 
#else
unsigned nrChannels	 = 1200; // 512;
#endif
unsigned nrPolarizations = 2;

int iters = 10;

unsigned correlator = CORRELATOR_1x1;
unsigned cellWidth = 1;
unsigned cellHeight = 1;

unsigned blockX = 128, blockY = 1;
unsigned gridX = 0, gridY = 1; // gridX is initialized later to the number of multiprocessors.

unsigned nrMultiprocessors = 30;

// if major = 1 and minor = 0, we cannot overlap communication and kernels, so nrStreams should be 1
#if USE_STREAMING
    int nrStreams = 1;
#else 
    int nrStreams = 1;
#endif

texture<float4> tex0;
texture<float4> unconvertedTexture;

__constant__ unsigned char cellToStatX[MAX_CELLS], cellToStatY[MAX_CELLS];
unsigned char cellToStatXHost[MAX_CELLS], cellToStatYHost[MAX_CELLS];

// table is much to large, but who cares...
__constant__ unsigned char missedStatX[4*MAX_STATIONS], missedStatY[4*MAX_STATIONS];
unsigned char missedStatXHost[4*MAX_STATIONS], missedStatYHost[4*MAX_STATIONS];

// for the shared mem versions
__shared__ float4 samples[MAX_STATIONS];

// is the combination of x and y correlated by the w by h cell-version?
// table is much to large, but who cares...
unsigned char correlated[MAX_STATIONS][MAX_STATIONS];

#include "gpu_correlator_1x1.cu"
#include "gpu_correlator_1x1_sharedmem.cu"
#include "gpu_correlator_1x2.cu"
#include "gpu_correlator_1x3.cu"
#include "gpu_correlator_1x4.cu"
#include "gpu_correlator_1x6.cu"
#include "gpu_correlator_2x2.cu"
#include "gpu_correlator_3x2.cu"
#include "gpu_correlator_3x2_spill.cu"
#include "gpu_correlator_3x3.cu"
#include "gpu_correlator_3x3_sharedmem.cu"
#include "gpu_correlator_4x3.cu"

unsigned fillCellToStatTable(const unsigned w, const unsigned h, const unsigned nrStations) 
{
    memset(correlated, 0, MAX_STATIONS * MAX_STATIONS);

    unsigned nrCells = 0;
    for(int statY = nrStations - h; statY >= 0; statY -= h) {
	for(int statX = 0; statX + w - 1 <= statY; statX += w) {
	    cellToStatXHost[nrCells] = statX;
	    cellToStatYHost[nrCells] = statY;
	    fprintf(stdout, "cellToStat[%u] = %u, %u\n", nrCells, statX, statY);

	    // fill the correlatedTable
	    for(int y=0; y<h; y++) {
		for(int x=0; x<w; x++) {
		    correlated[statX + x][statY + y] = 1;
		}
	    }

	    nrCells++;
	}
    }

    cudaMemcpyToSymbol(cellToStatX, cellToStatXHost, sizeof(cellToStatX));
    cudaMemcpyToSymbol(cellToStatY, cellToStatYHost, sizeof(cellToStatY));

    return nrCells;
}

#if CALCULATE_MISSING
unsigned fillMissedCellToStatTable(const unsigned w, const unsigned h, const unsigned nrStations) 
{
    unsigned nrMissed = 0;
    for(unsigned y=0; y<nrStations; y++) {
	for(unsigned x=0; x<=y; x++) {
	    if(!correlated[x][y]) {
		missedStatXHost[nrMissed] = x;
		missedStatYHost[nrMissed] = y;
		fprintf(stdout, "missedStat[%u] = %u, %u\n", nrMissed, x, y);
		nrMissed++;
	    }
	}
    }

    cudaMemcpyToSymbol(missedStatX, missedStatXHost, sizeof(missedStatX));
    cudaMemcpyToSymbol(missedStatY, missedStatYHost, sizeof(missedStatY));

    for(unsigned y=0; y<nrStations; y++) {
	std::cout.width(2);
	std::cout << y << " ";
	for(unsigned x=0; x<nrStations; x++) {
	    if(x>y) {
		std::cout << "  ";
	    } else {
		if(correlated[x][y]) {
		    std::cout << ".. ";
		} else {
		    std::cout << "** ";
		}
	    }
	}
	std::cout << std::endl;
    }

    std::cout << "   ";
    for(unsigned x=0; x<nrStations; x++) {
	std::cout.width(2);
	std::cout << x << " ";
    }
    std::cout << std::endl;

    return nrMissed;
}
#endif // CALCULATE_MISSING

void checkCudaCall(cudaError_t result)
{
    if (result != cudaSuccess) {
	std::cerr << "cuda error: " << cudaGetErrorString(result) << std::endl;
	exit(1);
    }
}

void gpu_correlate(float* devSamplesIn, float* devSamples, float *devVisibilities, 
		   unsigned nrTimes, unsigned nrTimesWidth, 
		   unsigned nrStations, unsigned nrChannels, unsigned nrCells, unsigned nrMissedCells,
		   cudaStream_t stream)
{
    unsigned nrBlocks = gridX * gridY;
    unsigned nrThreads = blockX * blockY;
    unsigned loopCount = LOOP_COUNT(nrCells, nrThreads);
    unsigned missingLoopCount = LOOP_COUNT(nrMissedCells, nrThreads);
    dim3 grid(gridX, gridY), block(blockX, blockY);

    // reset the error code
    cudaError_t error = cudaGetLastError();
    cudaGetLastError();
    if(error != cudaSuccess) {
	std::cout << "LAUNCH ERROR: " << cudaGetErrorString(error) << std::endl;
	exit(1);
    }

    if (correlator == CORRELATOR_1x1) {
	correlate_1x1<<<grid, block, 0, stream>>>(devVisibilities, nrTimes, nrTimesWidth, nrStations, 
				       nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_1x1_SHAREDMEM) {
	std::cout << "using 1x1 correlator with shared memory" << std::endl;
	correlate_1x1_sharedmem<<<grid, block, 0, stream>>>(devSamples, devVisibilities, nrTimes, nrTimesWidth, nrStations, 
						 nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_1x2) {
	correlate_1x2<<<grid, block, 0, stream>>>(devVisibilities, nrTimes, nrTimesWidth, nrStations, 
				       nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_1x3) {
	correlate_1x3<<<grid, block, 0, stream>>>(devVisibilities, nrTimes, nrTimesWidth, nrStations, 
				       nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_1x4) {
	correlate_1x4<<<grid, block, 0, stream>>>(devVisibilities, nrTimes, nrTimesWidth, nrStations, 
				       nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_1x6) {
	correlate_1x6<<<grid, block, 0, stream>>>(devVisibilities, nrTimes, nrTimesWidth, nrStations, 
				       nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_2x2) {
	correlate_2x2<<<grid, block, 0, stream>>>((float4*)devSamples, devVisibilities, nrTimes, nrTimesWidth, nrStations, 
				       nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_3x2) {
	correlate_3x2<<<grid, block, 0, stream>>>(devSamplesIn, devSamplesIn, devVisibilities, nrTimes, nrTimesWidth, nrStations,
						  nrChannels, nrCells, nrBlocks, nrThreads, loopCount, nrMissedCells, missingLoopCount);
    } else if (correlator == CORRELATOR_3x2_SPILL) {
	std::cout << "using 3x2 correlator with manual spilling" << std::endl;
	correlate_3x2_spill<<<grid, block, 0, stream>>>(devVisibilities, nrTimes, nrTimesWidth, nrStations,
					     nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_3x3) {
	correlate_3x3<<<grid, block, 0, stream>>>(devVisibilities, nrTimes, nrTimesWidth, nrStations, 
						  nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_3x3_SHAREDMEM) {
	std::cout << "using 3x3 correlator with shared memory" << std::endl;
	correlate_3x3_sharedmem<<<grid, block, 0, stream>>>(devSamples, devVisibilities, nrTimes, nrTimesWidth, nrStations,
				       nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else if (correlator == CORRELATOR_4x3) {
	correlate_4x3<<<grid, block, 0, stream>>>(devVisibilities, nrTimes, nrTimesWidth, nrStations, 
						  nrChannels, nrCells, nrBlocks, nrThreads, loopCount);
    } else {
	std::cerr << "error illegal correlator" << std::endl;
	exit(1);   
    }

    error = cudaGetLastError();
    if(error != cudaSuccess) {
	std::cout << "LAUNCH ERROR: " << cudaGetErrorString(error) << std::endl;
	exit(1);
    }
}

int main(int argc, char** argv)
{
    char* size = NULL;

    for(unsigned i=1; i<argc; i++) {
        if(strcmp(argv[i], "-correlator") == 0) {
            i++;
            size = argv[i];
        } else if(strcmp(argv[i], "-blockX") == 0) {
            i++;
            blockX = atoi(argv[i]);
        } else if(strcmp(argv[i], "-blockY") == 0) {
            i++;
            blockY = atoi(argv[i]);
        } else if(strcmp(argv[i], "-gridX") == 0) {
            i++;
            gridX = atoi(argv[i]);
        } else if(strcmp(argv[i], "-gridY") == 0) {
            i++;
            gridY = atoi(argv[i]);
        } else if(strcmp(argv[i], "-streams") == 0) {
            i++;
            nrStreams = atoi(argv[i]);
        } else if(strcmp(argv[i], "-iters") == 0) {
            i++;
            iters = atoi(argv[i]);
        } else if(strcmp(argv[i], "-device") == 0) {
            i++;
            cudaDevice = atoi(argv[i]);
        } else if(strcmp(argv[i], "-stations") == 0) {
            i++;
            nrStations = atoi(argv[i]);
        } else if(strcmp(argv[i], "-channels") == 0) {
            i++;
            nrChannels = atoi(argv[i]);
        } else if(strcmp(argv[i], "-verify") == 0) {
            verify = true;
        } else if(strcmp(argv[i], "-print") == 0) {
            printResult = true;
	} else {
	    std::cout << "illegal parameter: " << argv[i] << std::endl;
	    exit(1);
	}
    }

    if(size == NULL) {
        size = "1x1";
    }

    if(!strcmp(size, "1x1")) {
	correlator = CORRELATOR_1x1;
	cellWidth = 1;
	cellHeight = 1;
    } else if(!strcmp(size, "1x1-sharedmem")) {
	correlator = CORRELATOR_1x1_SHAREDMEM;
	cellWidth = 1;
	cellHeight = 1;
    } else if(!strcmp(size, "1x2")) {
	correlator = CORRELATOR_1x2;
	cellWidth = 1;
	cellHeight = 2;
    } else if(!strcmp(size, "1x3")) {
	correlator = CORRELATOR_1x3;
	cellWidth = 1;
	cellHeight = 3;
    } else if(!strcmp(size, "1x4")) {
	correlator = CORRELATOR_1x4;
	cellWidth = 1;
	cellHeight = 4;
    } else if(!strcmp(size, "1x6")) {
	correlator = CORRELATOR_1x6;
	cellWidth = 1;
	cellHeight = 6;
    } else if(!strcmp(size, "2x2")) {
	correlator = CORRELATOR_2x2;
	cellWidth = 2;
	cellHeight = 2;
    } else if(!strcmp(size, "3x2")) {
	correlator = CORRELATOR_3x2;
	cellWidth = 3;
	cellHeight = 2;
    } else if(!strcmp(size, "3x2-spill")) {
	correlator = CORRELATOR_3x2_SPILL;
	cellWidth = 3;
	cellHeight = 2;
    } else if(!strcmp(size, "3x3")) {
	correlator = CORRELATOR_3x3;
	cellWidth = 3;
	cellHeight = 3;
    } else if(!strcmp(size, "3x3-sharedmem")) {
	correlator = CORRELATOR_3x3_SHAREDMEM;
	cellWidth = 3;
	cellHeight = 3;
    } else if(!strcmp(size, "4x3")) {
	correlator = CORRELATOR_4x3;
	cellWidth = 4;
	cellHeight = 3;
    } else {
	std::cout << "illegal cell size" << std::endl;
	exit(1);
    }

    if(iters % nrStreams) {
	std::cout << "illegal argument, nr iterations must be dividable by the number of streams." << std:: endl;
	exit(1);
    }

    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    cudaDeviceProp p;
    for (int device=0; device<deviceCount; device++) {
	cudaGetDeviceProperties(&p, device);

	std::cout << "device " << device << " name: " << p.name << std::endl 
		  << "  global mem: " << p.totalGlobalMem
		  << " sharedMemPerBlock: " << p.sharedMemPerBlock << " regsPerBlock: " << p.regsPerBlock
		  << " warpSize: " << p.warpSize << " memPitch: " << p.memPitch << std::endl 
		  << "  maxThreadsPerBlock: " << p.maxThreadsPerBlock 
		  << " maxBlockSize: (" << p.maxThreadsDim[0] << ", "<< p.maxThreadsDim[1] << ", "<< p.maxThreadsDim[2]
		  << ") maxGridSize: (" << p.maxGridSize[0] << ", "<< p.maxGridSize[1] << ", "<< p.maxGridSize[2]
		  << ") const mem: " << p.totalConstMem << std::endl
		  << "  version: " << p.major << "." << p.minor 
		  << " clock rate: " << p.clockRate << " tex alignment: " << p.textureAlignment 
		  << " multiprocessors: " << p.multiProcessorCount << std::endl;

	if(device == cudaDevice) {
	    nrMultiprocessors = p.multiProcessorCount;
	    if(gridX == 0) {
		gridX = nrMultiprocessors;
	    }
	}
    }

    std::cout << "using device " << cudaDevice << std::endl;

    cudaSetDevice(cudaDevice);

    if(nrChannels % gridX) {
	std::cout << "WARNING, the number of channels is not divisable by gridX (the number of stream processors)." 
		  << std::endl << "  This leads to suboptimal performance" << std::endl;
    }

    unsigned nrBaselines = nrStations * (nrStations + 1) / 2;
    size_t samplesSize	  = nrStations * nrChannels * nrTimesWidth * nrPolarizations * sizeof(complex<float>);
    size_t visibilitiesSize = nrBaselines * nrChannels * nrPolarizations * nrPolarizations * sizeof(complex<float>);

    std::cout << "sample buf: " << (samplesSize / (1204.0 * 1024.0)) << " MB "
	      << "vis buf: " << (visibilitiesSize / (1024.0*1024.0)) << " MB " << std::endl;

    complex<float> *devSamples[nrStreams], *hostSamples[nrStreams];
    complex<float> *devVisibilities[nrStreams], *hostVisibilities[nrStreams];
    complex<float> *devSamplesIn[nrStreams];

#if USE_ALTERNATE_MEMORY_LAYOUT
    complex<float> *hostSamplesIn[nrStreams];
#endif

    // the 3.0 times is because we do three operations per cycle (multiply-add) + 1 multiply or add (on gtx280)
    double maxFlops = p.clockRate * 3.0 * (double)(nrMultiprocessors * NR_PROCESSORS_PER_MULTIPROCESSOR) / 1000000.0;
    unsigned nrCells = fillCellToStatTable(cellWidth, cellHeight, nrStations);

#if CALCULATE_MISSING
    unsigned nrMissedCells = fillMissedCellToStatTable(cellWidth, cellHeight, nrStations);
    std::cout << "The optimized correlator misses " << nrMissedCells << " cells" << std::endl;
#else
    unsigned nrMissedCells = 0;
#endif

    if(correlator == CORRELATOR_3x2_SPILL) {
	initTables();
    }

    unsigned long long loadedDeviceBytes = (unsigned long long)nrChannels * nrCells * (cellWidth + cellHeight) * nrTimes * sizeof(float4);
    double gb = (double) loadedDeviceBytes / (1024.0 * 1024.0 * 1024.0);
    
    // 4 mads x 2 ops (mul + add) x 4 regs (float4)
    // 4 * 2 * 4 * cellW * cellH
    unsigned long long ops =  (unsigned long long) nrChannels * nrCells 
	* nrTimes * 4L * 2L * 4L * cellWidth * cellHeight;

    std::cout << "using " << cellWidth << "X" << cellHeight << " correlator" << std::endl;
    std::cout << "correlating with " << nrStations << " stations, " << nrBaselines << " baselines, " 
	      << nrCells << " cells, " << nrChannels << " channels" << std::endl;
    std::cout << "ops: " << ops << ", loaded device memory bytes = " << gb << " GB, " << std::endl;

#if USE_STREAMING
    std::cout << "streaming version, with " << nrStreams << " streams" << std::endl;
#else
    std::cout << "non-streaming version" << std::endl;
#endif

#if DO_LOADS
    std::cout << "doing normal memory loads" << std::endl;
#else
    std::cout << "doing NO memory loads" << std::endl;
#endif

#if DO_STORES
    std::cout << "doing normal device memory stores" << std::endl;
#else
    std::cout << "doing NO device memory stores" << std::endl;
#endif

#if SYNC_THREADS
    std::cout << "doing syncthreads after memory loads" << std::endl;
#else
    std::cout << "NOT doing syncthreads after memory loads" << std::endl;
#endif

#if USE_ALTERNATE_MEMORY_LAYOUT
    std::cout << "using alternate memory layout, grouping all station data in time" << std::endl;

#if DO_MEMORY_CONVERSION
    std::cout << "performing memory conversion" << std::endl;    

#if CONVERT_ON_GPU
    std::cout << "performing memory conversion on the GPU" << std::endl;    
#else
#if USE_MULTITHREADED_MEMORY_CONVERSION
    std::cout << "using multi-treaded memory conversion, with " << MULTITHREADED_MEMORY_CONVERSION_THREADS << " threads" << std::endl;    
#else // !USE_MULTITHREADED_MEMORY_CONVERSION
    std::cout << "using single-treaded memory conversion" << std::endl;    
#endif // USE_MULTITHREADED_MEMORY_CONVERSION
#endif // CONVERT_ON_GPU

#else // !DO_MEMORY_CONVERSION
    std::cout << "NOT performing memory conversion" << std::endl;    
#endif

#else // !USE_ALTERNATE_MEMORY_LAYOUT
    std::cout << "using normal memory layout" << std::endl;
#endif // USE_ALTERNATE_MEMORY_LAYOUT


#if USE_WRITE_COMBINING
    std::cout << "storing host samples in pinned and write-combined memory" << std::endl;
#else
    std::cout << "storing host samples in pinned, but NOT write-combined memory" << std::endl;
#endif

    std::cout << "allocating memory..." << std::endl;

    for(int i=0; i<nrStreams; i++) {	
	checkCudaCall(cudaMalloc((void **) &devSamples[i], samplesSize));
	assert(devSamples[i]);
	checkCudaCall(cudaMemset(devSamples[i], 0, samplesSize));

#if USE_ALTERNATE_MEMORY_LAYOUT
	// also pin the in buffers, to avoid swapping
	checkCudaCall(cudaHostAlloc((void **) &hostSamplesIn[i], samplesSize, 0));
	assert(hostSamplesIn[i]);
	devSamplesIn[i] = 0;
#if CONVERT_ON_GPU
	checkCudaCall(cudaMalloc((void **) &devSamplesIn[i], samplesSize));
	assert(devSamplesIn[i]);
	checkCudaCall(cudaMemset(devSamplesIn[i], 0, samplesSize));
#endif  // CONVERT_ON_GPU
#endif // USE_ALTERNATE_MEMORY_LAYOUT

	// pinned memory, write combining
#if USE_WRITE_COMBINING
	checkCudaCall(cudaHostAlloc((void **) &hostSamples[i], samplesSize, cudaHostAllocWriteCombined));
#else
	checkCudaCall(cudaHostAlloc((void **) &hostSamples[i], samplesSize, 0));
        // for unpinned transfer, we could use malloc
//	hostSamples[i] = (complex<float>*) malloc(samplesSize);
#endif
	assert(hostSamples[i]);
	memset(hostSamples[i], 0, samplesSize);

	checkCudaCall(cudaMalloc((void **) &devVisibilities[i], visibilitiesSize));
	assert(devVisibilities[i]);
	checkCudaCall(cudaMemset(devVisibilities[i], 0, visibilitiesSize));
	checkCudaCall(cudaHostAlloc((void **) &hostVisibilities[i], visibilitiesSize, 0));
	assert(hostVisibilities[i]);
	memset(hostVisibilities[i], 0, visibilitiesSize);
    }

    std::cout << "allocating memory DONE" << std::endl;

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();

#if USE_ALTERNATE_MEMORY_LAYOUT && DO_MEMORY_CONVERSION && CONVERT_ON_GPU
    // create a texture for the unconverted samples. 
    size_t offset;
    checkCudaCall(cudaBindTexture(&offset, unconvertedTexture, devSamplesIn[0], channelDesc, samplesSize));
    assert(offset == 0);

    checkCudaCall(cudaBindTexture(&offset, tex0, devSamples[0], channelDesc, samplesSize));
    assert(offset == 0);
#else
    // create a texture. 
    size_t offset;
    checkCudaCall(cudaBindTexture(&offset, tex0, devSamples[0], channelDesc, samplesSize));
    assert(offset == 0);
#endif

    std::cout << "initializing samples..." << std::endl;

    // Initialize the samples to some pseudo-random value.
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned stat = 0; stat < nrStations; stat ++) {
	    for (unsigned time = 0; time < nrTimes; time ++) {
		for(int i=0; i<nrStreams; i++) {
#if USE_ALTERNATE_MEMORY_LAYOUT
#if DO_MEMORY_CONVERSION
		    hostSamplesIn[i][NORMAL_SAMPLE_INDEX(stat, channel, time, 0)] = complex<float>(time % 8, stat);
		    hostSamplesIn[i][NORMAL_SAMPLE_INDEX(stat, channel, time, 1)] = complex<float>(time % 8 + 1, channel%50);
#else
		    hostSamplesIn[i][ALTERNATE_SAMPLE_INDEX(stat, channel, time, 0)] = complex<float>(time % 8, stat);
		    hostSamplesIn[i][ALTERNATE_SAMPLE_INDEX(stat, channel, time, 1)] = complex<float>(time % 8 + 1, channel%50);
#endif // DO_MEMORY_CONVERSION

#else
		    hostSamples[i][NORMAL_SAMPLE_INDEX(stat, channel, time, 0)] = complex<float>(time % 8, stat);
		    hostSamples[i][NORMAL_SAMPLE_INDEX(stat, channel, time, 1)] = complex<float>(time % 8 + 1, channel%50);
#endif
		}
	    }
	}
    }
    std::cout << "initializing samples DONE" << std::endl;

    unsigned nrThreads = blockX * blockY;
    unsigned loopCount = LOOP_COUNT(nrCells, nrThreads);
    unsigned missingLoopCount = LOOP_COUNT(nrMissedCells, nrThreads);

    std::cout << "using a " << gridX << "x" << gridY << " grid, with " 
	 << blockX << "x" << blockY << " blocks, nr of threads per block: " 
	      << nrThreads << " loopcount = " << loopCount << ", missing loopcount = " << missingLoopCount << std::endl;

#if ! USE_STREAMING

    // create cuda event handles
    cudaEvent_t start, stop;
    checkCudaCall(cudaEventCreate(&start));
    checkCudaCall(cudaEventCreate(&stop));

    timer totalTimer("total");
    totalTimer.start();
    for (int i = 0; i < iters; i ++) {
	timer loadTimer("load transfer");

#if DO_IO

#if USE_ALTERNATE_MEMORY_LAYOUT && DO_MEMORY_CONVERSION && !CONVERT_ON_GPU
	convertBufferLayout(hostSamplesIn[0], hostSamples[0]);
#endif // USE_ALTERNATE_MEMORY_LAYOUT && DO_MEMORY_CONVERSION && !CONVERT_ON_GPU

	loadTimer.start();
	checkCudaCall(cudaMemcpy(devSamples[0], hostSamples[0], samplesSize, cudaMemcpyHostToDevice));
	loadTimer.stop();

	std::cout << "PCI-e host-to-device load time of " << samplesSize << " bytes was " << loadTimer.getTimeInSeconds() << " achieved " 
	     << (samplesSize / (1024.0 * 1024.0 * 1024.0)) / loadTimer.getTimeInSeconds() << " GB/s" << std::endl;
#endif // DO_IO

	checkCudaCall(cudaEventRecord(start, 0));
	gpu_correlate(reinterpret_cast<float *>(devSamplesIn[0]), reinterpret_cast<float *>(devSamples[0]), reinterpret_cast<float *>(devVisibilities[0]), 
		      nrTimes, nrTimesWidth, nrStations, nrChannels, nrCells, nrMissedCells, 0);
	checkCudaCall(cudaEventRecord(stop, 0));
	checkCudaCall(cudaEventSynchronize(stop));

#if DO_IO
	checkCudaCall(cudaMemcpy(hostVisibilities[0], devVisibilities[0], visibilitiesSize, cudaMemcpyDeviceToHost));
#endif

#if 0
	timer missingTimer("missing");
	missingTimer.start();
	correlateMissingOnHost(nrMissedCells, hostSamples[0], hostVisibilities[0]);
	missingTimer.stop();
	std::cout << "correlation of missed baselines took: " << missingTimer.getTimeInSeconds() << std::endl;
#endif

	float elapsed;
	checkCudaCall(cudaEventElapsedTime(&elapsed, start, stop)); // time in ms
	elapsed /= 1000.0;
	double flops = (ops / elapsed) / 1000000000.0;
	double efficiency = (flops / maxFlops) * 100.0;
	double memEfficiency = ((gb/elapsed) / MAX_DEVICE_MEM_BANDWIDTH) * 100.0;

	std::cout << "correlate took " << elapsed << "s, max Gflops = " << maxFlops 
		  << ", achieved " << flops << " Gflops, " << efficiency << " % efficiency, loaded "
		  << gb << " GB, " << (gb / elapsed) << " GB/s, " << memEfficiency << "%" << std::endl;
	
	double gbAssumingReuse = (samplesSize + visibilitiesSize) / (1024.0*1024.0*1024.0);

	std::cout << "memory throughput, assuming perfect reuse: " << gbAssumingReuse << " GB, "
		  << (gbAssumingReuse / elapsed) << " GB/s, " << std::endl;

	double totalElapsed = loadTimer.getTimeInSeconds() + elapsed;
	double totalFlops =  (ops / totalElapsed) / 1000000000.0;
	double totalEfficiency = (totalFlops / maxFlops) * 100.0;
	double totalMemEfficiency = ((gb/totalElapsed) / MAX_DEVICE_MEM_BANDWIDTH) * 100.0;

	std::cout.precision(10);
	std::cout << "including host-to-device copy: " << totalFlops << " gflops, " << totalEfficiency 
		  << " %, bandwidth = " << (gb/totalElapsed) << " GB/s, " 
		  << totalMemEfficiency << " %" << std::endl;
    }
    totalTimer.stop();
    std::cout << totalTimer;

#else // USE_STREAMING

    cudaStream_t streamID[nrStreams];

    // create the streams
    for(int i=0; i<nrStreams; i++) {
	checkCudaCall(cudaStreamCreate(&streamID[i]));
    }

    // create cuda event handles
    cudaEvent_t start, stop;
    checkCudaCall(cudaEventCreate(&start));
    checkCudaCall(cudaEventCreate(&stop));

    timer totalTimer("total");
    totalTimer.start();

#if USE_ALTERNATE_MEMORY_LAYOUT && DO_MEMORY_CONVERSION && !CONVERT_ON_GPU
    // convert the first buffer
    convertBufferLayout(hostSamplesIn[0], hostSamples[0]);
#endif // USE_ALTERNATE_MEMORY_LAYOUT && DO_MEMORY_CONVERSION && !CONVERT_ON_GPU

    for (int i = 0; i <iters/nrStreams; i++) {

	cudaEventRecord(start, 0);
	for(int stream = 0; stream < nrStreams; stream++) {
#if DO_IO
	    checkCudaCall(cudaMemcpyAsync(devSamples[stream], hostSamples[stream], 
					  samplesSize, cudaMemcpyHostToDevice, streamID[stream]));
#endif
	    gpu_correlate(reinterpret_cast<float *>(devSamplesIn[stream]), reinterpret_cast<float *>(devSamples[stream]), reinterpret_cast<float *>(devVisibilities[stream]), 
			  nrTimes, nrTimesWidth, nrStations, nrChannels, nrCells, nrMissedCells, streamID[stream]);

#if USE_ALTERNATE_MEMORY_LAYOUT && DO_MEMORY_CONVERSION && !CONVERT_ON_GPU
	    // convert the next buffer while correlating
	    if( stream != nrStreams-1) {
		convertBufferLayout(hostSamplesIn[stream+1], hostSamples[stream+1]);
	    }
#endif // USE_ALTERNATE_MEMORY_LAYOUT && DO_MEMORY_CONVERSION && !CONVERT_ON_GPU

	    // for some reason, this kills all async behaviour. Just copy everything back in one loop later, using sync copy.
//	    checkCudaCall(cudaMemcpyAsync(hostVisibilities[stream], devVisibilities[stream], visibilitiesSize, cudaMemcpyDeviceToHost, stream));
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

#if DO_IO
	// copy visibilities back to host mem
	for(int stream = 0; stream < nrStreams; stream++) {
	    checkCudaCall(cudaMemcpy(hostVisibilities[stream], devVisibilities[stream], visibilitiesSize, cudaMemcpyDeviceToHost));
	}
#endif
    }

    totalTimer.stop();
    std::cout << totalTimer;

#endif // USE_STREAMING

    ops *= iters;
    gb *= iters;

    double elapsed = totalTimer.getTimeInSeconds();
    double gflops = (ops / elapsed) / 1000000000.0;
    double efficiency = (gflops / maxFlops) * 100.0;
    double memEfficiency = ((gb/elapsed) / MAX_DEVICE_MEM_BANDWIDTH) * 100.0;

    std::cout.precision(5);
    std::cout << "total time = " << elapsed << " s, achieved " << gflops << " gflops" 
	      << ", " << efficiency << "%, bandwidth = " 
	      << (gb / elapsed) << " GB/s, " << memEfficiency << "%" << std::endl;

    if(verify) {
	checkResult(visibilitiesSize, hostSamples[0], hostVisibilities[0], nrBaselines, nrCells);
    }


    checkCudaCall(cudaUnbindTexture(tex0));

    for(int i=0; i<nrStreams; i++) {
	checkCudaCall(cudaFree(devSamples[i]));    
	checkCudaCall(cudaFreeHost(hostSamples[i]));
	checkCudaCall(cudaFree(devVisibilities[i]));
	checkCudaCall(cudaFreeHost(hostVisibilities[i]));
    }

    return 0;
}
