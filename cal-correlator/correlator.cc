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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <cassert>

#include "correlator.h"
#include "timer.h"

#include "correlator-1x1.h"
#include "correlator-2x2.h"
#include "correlator-1x1-mem.h"
#include "correlator-2x2-mem.h"
#include "correlator-3x2-mem.h"
#include "correlator-3x3-mem.h"
#include "correlator-4x3-mem.h"

int main(int argc, char** argv)
{
    unsigned nrStations      = 64;
    unsigned nrTimes	     = 768;
    unsigned nrChannels      = 100;
    unsigned nrBaselines     = nrStations * (nrStations + 1) / 2;
    unsigned iters           = 100;

    unsigned nrPolarizations = 2;
    unsigned nrSamples      = nrChannels * nrStations * nrTimes;
    unsigned nrVisibilities = nrBaselines * nrChannels * nrPolarizations;

    unsigned cellSize = CELL_1X1;
    unsigned cellWidth = 1;
    unsigned cellHeight = 1;
    bool useOutputRegisters = false;
    
    disassemble = false;
    print = false;
    bool sync = false;
    const char* size = 0;
    const CALchar* kernel = NULL;

    for(int i=0; i<argc; i++) {
	if(strcmp(argv[i], "-d") == 0) {
	    disassemble = true;
	}
	if(strcmp(argv[i], "-p") == 0) {
	    print = true;
	}
	if(strcmp(argv[i], "-sync") == 0) {
	    sync = true;
	}
	if(strcmp(argv[i], "-regs") == 0) {
	    useOutputRegisters = true;
	}
	if(strcmp(argv[i], "-size") == 0) {
	    i++;
	    size = argv[i];
	}
    }

    if(size == NULL) {
	size = "1x1";
    }

    if(!strcmp(size, "1x1")) {
	cellSize = CELL_1X1;
	cellWidth = 1;
	cellHeight = 1;
	if(useOutputRegisters) {
	    kernel = Kernel_1x1;
	} else {
	    kernel = Kernel_1x1_mem;
	}
    } else if(!strcmp(size, "2x2")) {
	cellSize = CELL_2X2;
	cellWidth = 2;
	cellHeight = 2;
	if(useOutputRegisters) {
	    kernel = Kernel_2x2;
	} else {
	    kernel = Kernel_2x2_mem;
	}
    } else if(!strcmp(size, "3x2")) {
	cellSize = CELL_3X2;
	cellWidth = 3;
	cellHeight = 2;
	if(useOutputRegisters) {
	    std::cout << "illegal cell size, only 8 output registers available, cannot do more than 2x2" << std::endl;
	} else {
	    kernel = Kernel_3x2_mem;
	}
    } else if(!strcmp(size, "3x3")) {
	cellSize = CELL_3X3;
	cellWidth = 3;
	cellHeight = 3;
	if(useOutputRegisters) {
	    std::cout << "illegal cell size, only 8 output registers available, cannot do more than 2x2" << std::endl;
	} else {
	    kernel = Kernel_3x3_mem;
	}
    } else if(!strcmp(size, "4x3")) {
	cellSize = CELL_4X3;
	cellWidth = 4;
	cellHeight = 3;
	if(useOutputRegisters) {
	    std::cout << "illegal cell size, only 8 output registers available, cannot do more than 2x2" << std::endl;
	} else {
	    kernel = Kernel_4x3_mem;
	}
    } else {
	std::cout << "illegal cell size" << std::endl;
    }
   
    unsigned nrCells = calcNrCells(cellWidth, cellHeight, nrStations);
    unsigned nrOutputsPerCell = calcNrOutputsPerCell(cellWidth, cellHeight, nrPolarizations);

    cal_check(calInit());

    CALuint numDevices = 0;
    calDeviceGetCount(&numDevices);

    std::cout << "number of devices = " << numDevices << std::endl;

    if(numDevices < 1) {
	    std::cout << "no CAL device found" << std::endl;
	    exit(1);
    }

    // Opening device
    CALdevice device = 0;
    cal_check(calDeviceOpen(&device, 0));

    // Querying device attribs
    CALdeviceattribs attribs;
    attribs.struct_size = sizeof(CALdeviceattribs);
    calDeviceGetAttribs(&attribs, 0);
    CALuint version[3];
    calGetVersion(&version[0], &version[1], &version[2]);

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "CAL Runtime version: " << version[0] << "." << version[1] << "." << version[2] << std::endl;
    std::cout << "target: " << attribs.target << ", physical RAM: " << attribs.localRAM 
	      << " MB, uncached remote RAM: " << attribs.uncachedRemoteRAM 
	      << " MB, cached remote RAM: " << attribs.cachedRemoteRAM 
	      << " MB, engine clock: " << attribs.engineClock 
	      << " MHz, memory clock: " << attribs.memoryClock << " MHZ" << std::endl;
    std::cout << "wavefront size: " << attribs.wavefrontSize 
	      << ", number of SIMDs: " << attribs.numberOfSIMD << std::endl;

    CALdeviceinfo info;
    cal_check(calDeviceGetInfo(&info, 0));
    std::cout << "max 1D elts: " << info.maxResource1DWidth
	      << ", max 2D elts: " << info.maxResource2DWidth
	      << " X " << info.maxResource2DHeight << std::endl;


    // Creating context w.r.t. to opened device
    CALcontext ctx = 0;
    calCtxCreate(&ctx, device);

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "running " << cellWidth << "X" << cellHeight << " correlator with " << nrStations << " stations, " 
	      << nrBaselines << " baselines, " << nrCells << " cells, " << nrChannels << " channels, nrTimes = " << nrTimes << std::endl;
    std::cout << "samples: " << nrSamples << " (" << nrTimes << " X " << (nrChannels * nrStations) << " x 4) elements, " 
	      << (4.0 * sizeof(float) * nrSamples) / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "visibilities: " << nrVisibilities << " (" << nrBaselines << " x " << nrChannels << " x " 
	      << nrOutputsPerCell << ") elements, "
	      << (nrVisibilities * sizeof(float) * 2.0) / (1024.0 * 1024.0) << " MB" << std::endl;

    if(useOutputRegisters) {
	std::cout << "using output registers" << std::endl;
    } else {
	std::cout << "using direct memory writes" << std::endl;
    }

    initCounters();

    // Input and output resources
    CALresource hostSamplesRes[2];
    CALresource deviceSamplesRes[2];
    CALresource deviceCellToStatXRes = 0;
    CALresource deviceCellToStatYRes = 0;

    // use 2D sample buffer:
    // - elts are float4, representing both polarizations X real/imag
    // - layout is nrTimes x channel x station
    // in 2D, we see nrTimes as one dimension, (channel, station) as the other.
    //
    // vis buffer: nrCells * nrChannels * nrPolarizations * nrPolarizations * real/imag;
    // use 2 2d buffers, one for reals, one for imags.
    // - real and imag are nrCells * nrChannels * float4, where each float is a combination of polarizations.

    // allocate device memory, format is width, height
    cal_check(calResAllocRemote2D(&hostSamplesRes[0], &device, 1, nrTimes, nrChannels * nrStations, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocRemote2D(&hostSamplesRes[1], &device, 1, nrTimes, nrChannels * nrStations, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(&deviceSamplesRes[0], device, nrTimes, nrChannels * nrStations, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(&deviceSamplesRes[1], device, nrTimes, nrChannels * nrStations, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal1D(&deviceCellToStatXRes, device, nrCells, CAL_FORMAT_FLOAT_1, 0));
    cal_check(calResAllocLocal1D(&deviceCellToStatYRes, device, nrCells, CAL_FORMAT_FLOAT_1, 0));

    CALmem hostSamplesMem[2];
    CALmem deviceSamplesMem[2];
    cal_check(calCtxGetMem(&hostSamplesMem[0], ctx, hostSamplesRes[0]));
    cal_check(calCtxGetMem(&hostSamplesMem[1], ctx, hostSamplesRes[1]));
    cal_check(calCtxGetMem(&deviceSamplesMem[0], ctx, deviceSamplesRes[0]));
    cal_check(calCtxGetMem(&deviceSamplesMem[1], ctx, deviceSamplesRes[1]));

    CALresource deviceVisReal02 = 0;
    CALresource deviceVisImag02 = 0;
    CALresource deviceVisReal12 = 0;
    CALresource deviceVisImag12 = 0;
    CALresource deviceVisReal03 = 0;
    CALresource deviceVisImag03 = 0;
    CALresource deviceVisReal13 = 0;
    CALresource deviceVisImag13 = 0;
    CALresource hostVisibilitiesRes = 0;

    if(useOutputRegisters) {
	if(cellSize == CELL_1X1) {
	    allocateVis1X1Reg(&deviceVisReal02, &deviceVisImag02,
			      ctx, device, nrCells, nrChannels);
	} else if(cellSize == CELL_2X2) {
	    allocateVis2X2Reg(&deviceVisReal02, &deviceVisImag02, &deviceVisReal12, &deviceVisImag12,
			      &deviceVisReal03, &deviceVisImag03, &deviceVisReal13, &deviceVisImag13,
			      ctx, device, nrCells, nrChannels);
	} else {
	    std::cout << "illegal cell size" << std::endl;
	}
    } else {    
	std::cout << "elts: " << (nrCells * nrChannels) << " x " << nrOutputsPerCell 
		  << " = " << (nrCells * nrChannels * nrOutputsPerCell) << std::endl;

	// Allocate host memory, we cannot allocate more than 8192 elements.
	cal_check(calResAllocRemote2D(&hostVisibilitiesRes, &device, 1, nrCells, nrOutputsPerCell * nrChannels, CAL_FORMAT_FLOAT_4, 0));

	// wipe output buffers
	wipeOutputMem(&ctx, hostVisibilitiesRes, nrChannels);
    }

    // Constants
    CALresource constRes = 0;
    cal_check(calResAllocRemote1D(&constRes, &device, 1, 1, CAL_FORMAT_FLOAT_4, 0));
    CALfloat* constPtr = NULL;
    CALuint constPitch = 0;
    CALmem constMem = 0;
    cal_check(calCtxGetMem(&constMem, ctx, constRes));
    cal_check(calResMap((CALvoid**)&constPtr, &constPitch, constRes, 0));
    constPtr[0] = (float) nrStations;
    constPtr[1] = (float) nrTimes;
    constPtr[2] = (float) nrChannels;
    constPtr[3] = (float) nrCells;
    cal_check(calResUnmap(constRes));
    calCtxReleaseMem(ctx, constMem);

#if 0
    // init samples
    CALfloat* fdata = NULL;
    CALuint pitch = 0;
    CALmem deviceSamplesMem = 0;
    cal_check(calCtxGetMem(&deviceSamplesMem, ctx, deviceSamplesRes0));
    cal_check(calResMap((CALvoid**)&fdata, &pitch, deviceSamplesRes0, 0));
    std::cout << "sample pitch: " << pitch << " nrTimes = " << nrTimes << " width = " << (nrTimes * sizeof(float4)) << std::endl;
    initSamples(nrStations, nrChannels, nrTimes, pitch, (float4*)fdata);
    timer unmapTimer("unmap");
    unmapTimer.start();
    cal_check(calResUnmap(deviceSamplesRes0));
    unmapTimer.stop();
    double elapsedUnmap = unmapTimer.getTimeInSeconds();
    double realSampleGb = (pitch * sizeof(float4)) * nrChannels * nrStations / (1024.0 * 1024.0 * 1024.0);
    std::cout << "unmap took: " << elapsedUnmap << " s, throughput = " << (realSampleGb / elapsedUnmap)<< " GB/s" << std::endl;
    calCtxReleaseMem(ctx, deviceSamplesMem);
#endif

    // init samples
    CALfloat* fdata = NULL;
    CALuint pitch = 0;

    cal_check(calResMap((CALvoid**)&fdata, &pitch, hostSamplesRes[0], 0)); // does not copy, both are remote
    initSamples(nrStations, nrChannels, nrTimes, pitch, (float4*)fdata);
    cal_check(calResUnmap(hostSamplesRes[0])); // does not copy, both are remote

    cal_check(calResMap((CALvoid**)&fdata, &pitch, hostSamplesRes[1], 0)); // does not copy, both are remote
    initSamples(nrStations, nrChannels, nrTimes, pitch, (float4*)fdata);
    cal_check(calResUnmap(hostSamplesRes[1])); // does not copy, both are remote

    double realSampleGb = (pitch * sizeof(float4)) * nrChannels * nrStations / (1024.0 * 1024.0 * 1024.0);
    std::cout << "sample pitch: " << pitch << " nrTimes = " << nrTimes << " width = " << (nrTimes * sizeof(float4)) << std::endl;

    // init cell to station tables
    CALfloat* fdataX = NULL, *fdataY = NULL;
    CALuint pitchX = 0, pitchY = 0;
    CALmem hostCellToStatX = 0, hostCellToStatY = 0;
    cal_check(calCtxGetMem(&hostCellToStatX, ctx, deviceCellToStatXRes));
    cal_check(calCtxGetMem(&hostCellToStatY, ctx, deviceCellToStatYRes));
    cal_check(calResMap((CALvoid**)&fdataY, &pitchY, deviceCellToStatYRes, 0));
    cal_check(calResMap((CALvoid**)&fdataX, &pitchX, deviceCellToStatXRes, 0));
    fillCellToStatTable(cellWidth, cellHeight, nrStations, fdataX, fdataY);
    cal_check(calResUnmap(deviceCellToStatXRes));
    cal_check(calResUnmap(deviceCellToStatYRes));

    //-------------------------------------------------------------------------
    // Loading module and setting domain
    //-------------------------------------------------------------------------

    // Compiling Device Program
    CALobject obj = NULL;
    CALimage image = NULL;


    if (calclCompile(&obj, CAL_LANGUAGE_IL, kernel, attribs.target) != CAL_RESULT_OK) {
        fprintf(stdout, "Program compilation failed. Exiting.\n");
        return 1;
    }

    if (calclLink(&image, &obj, 1) != CAL_RESULT_OK) {
        fprintf(stdout, "Program linking failed. Exiting.\n");
        return 1;
    }

    if(disassemble) {
	calclDisassembleImage(image, __logger);
    }

    // Creating module using compiled image
    CALmodule module = 0;
    if (calModuleLoad(&module, ctx, image) != CAL_RESULT_OK) {
        fprintf(stdout, "Program loading failed. Exiting.\n");
        return 1;
    }

    // Defining symbols in module
    CALfunc func = 0;

    // Defining entry point for the module
    if (calModuleGetEntry(&func, ctx, module, "main") != CAL_RESULT_OK) {
        fprintf(stdout, "Could not locate 'main'. Exiting.\n");
        return 1;
    }

    bindVariable(&ctx, &module, deviceCellToStatXRes, 1, VAR_INPUT);
    bindVariable(&ctx, &module, deviceCellToStatYRes, 2, VAR_INPUT);
    bindVariable(&ctx, &module, constRes,          0, VAR_CONSTANT);

    if(useOutputRegisters) {
	if(cellSize == CELL_1X1) {
	    bindVariable(&ctx, &module, deviceVisReal02,   0, VAR_OUTPUT);
	    bindVariable(&ctx, &module, deviceVisImag02,   1, VAR_OUTPUT);
	} else if (cellSize == CELL_2X2) {
	    bindVariable(&ctx, &module, deviceVisReal02,   0, VAR_OUTPUT);
	    bindVariable(&ctx, &module, deviceVisImag02,   1, VAR_OUTPUT);
	    bindVariable(&ctx, &module, deviceVisReal12,   2, VAR_OUTPUT);
	    bindVariable(&ctx, &module, deviceVisImag12,   3, VAR_OUTPUT);
	    bindVariable(&ctx, &module, deviceVisReal03,   4, VAR_OUTPUT);
	    bindVariable(&ctx, &module, deviceVisImag03,   5, VAR_OUTPUT);
	    bindVariable(&ctx, &module, deviceVisReal13,   6, VAR_OUTPUT);
	    bindVariable(&ctx, &module, deviceVisImag13,   7, VAR_OUTPUT);
	} else {
	    std::cout << "illegal cell size, only 8 output registers available, cannot do more than 2x2" << std::endl;
	}
    } else {
	bindVariable(&ctx, &module, hostVisibilitiesRes,   0, VAR_GLOBAL);
    }

#if 0 // does not work :-(
    CALfuncInfo funcInfo;
    cal_check(calModuleGetFuncInfo(&funcInfo, ctx, module, func));
    std::cout << "func info: maxScratchRegs: " << funcInfo.maxScratchRegsNeeded
	      << ", shared GPRs: " << funcInfo.numSharedGPRUser
	      << ", shared GPRs total: " << funcInfo.numSharedGPRTotal
	      << ", numThreadPerGroup: " << funcInfo.numThreadPerGroup
	      << ", totalNumThreadGroups: " << funcInfo.totalNumThreadGroup
	      << ", waveFrontsPerSIMD: " << funcInfo.wavefrontPerSIMD
	      << ", waveFrontsPerSIMD2: " << funcInfo.numWavefrontPerSIMD
	      << ", isMaxNumWavePerSIMD: " << funcInfo.isMaxNumWavePerSIMD << std::endl;
#endif

    //-------------------------------------------------------------------------
    // Executing program and waiting for program to terminate
    //-------------------------------------------------------------------------
    // Setting domain
    // TODO: use a value that is tuned to the hardware.
    CALdomain domain = {0, 0, nrCells, nrChannels};
    std::cout << "domain is " << nrCells << " X " << nrChannels << ", " << nrCells * nrChannels << " threads in total" << std::endl;

    CALcounter idleCounter;
    CALcounter cacheCounter;
    if(calCtxCreateCounterExt(&idleCounter, ctx, CAL_COUNTER_IDLE) != CAL_RESULT_OK) {
	fprintf(stdout, "Counter creation failed: %s. Exiting.\n", calGetErrorString());
	return 1;
    }

    if(calCtxCreateCounterExt(&cacheCounter, ctx, CAL_COUNTER_INPUT_CACHE_HIT_RATE) != CAL_RESULT_OK) {
	fprintf(stdout, "Counter creation failed: %s. Exiting.\n", calGetErrorString());
	return 1;
    }

    CALevent e = 0; // Event to check completion of the program
    CALevent copyEvent[2]; // Events to check completion of the copies
    unsigned currBuf = 0;
    unsigned nextBuf = 1;


    // measure host -> GPU throughput
    unsigned memIters = 50;
    timer copyTimer("copy");
    copyTimer.start();
    for(unsigned i = 0; i< memIters; i++) {
	copyMem(hostSamplesMem[currBuf], deviceSamplesMem[currBuf], ctx);
    }
    copyTimer.stop();
    double elapsedCopy = copyTimer.getTimeInSeconds();
    std::cout << "copy of " << (realSampleGb*memIters) << " GB took "
	      << elapsedCopy << ", throughput = " 
	      << ((realSampleGb * memIters) / elapsedCopy) 
	      << " GB/s" << std::endl;


    calCtxBeginCounterExt(ctx, idleCounter);
    calCtxBeginCounterExt(ctx, cacheCounter);

    timer runTimer("run");
    runTimer.start();

    bool doCopy = true;
    timer waitTimer("wait");
    timer kernelTimer("kernel");

    if(doCopy) {
	// copy the first buffer by hand
	copyEvent[currBuf] = copyMemAsync(hostSamplesMem[currBuf], deviceSamplesMem[currBuf], ctx);
	calCtxIsEventDone(ctx, copyEvent[currBuf]);

	if(sync) {
	    waitTimer.start();
	    while (calCtxIsEventDone(ctx, copyEvent[currBuf]) == CAL_RESULT_PENDING);
	    waitTimer.stop();
	}
    }

    for(unsigned i=0; i<iters; i++) {
	if(doCopy) {
	    if( i != iters-1) {
		// start copy of next buffer
		copyEvent[nextBuf] = copyMemAsync(hostSamplesMem[nextBuf], deviceSamplesMem[nextBuf], ctx);
		calCtxIsEventDone(ctx, copyEvent[nextBuf]);
		    if(sync) {
		    waitTimer.start();
		    while (calCtxIsEventDone(ctx, copyEvent[nextBuf]) == CAL_RESULT_PENDING);
		    waitTimer.stop();
		}
	    }
	}

	bindVariable(&ctx, &module, deviceSamplesRes[currBuf], 0, VAR_INPUT); // TODO, use optimized version

	if(doCopy && !sync) {
	    waitTimer.start();
	    // wait for the copy of the current buffer to complete
	    while (calCtxIsEventDone(ctx, copyEvent[currBuf]) == CAL_RESULT_PENDING);
	    waitTimer.stop();
	}
	
	// run the kernel
	kernelTimer.start();
	cal_check(calCtxRunProgram(&e, ctx, func, &domain));
	
	// Checking whether the execution of the program is complete or not
	while (calCtxIsEventDone(ctx, e) == CAL_RESULT_PENDING);
	kernelTimer.stop();

	if(useOutputRegisters) {
	    if(cellSize == CELL_1X1) {
		downloadVis1x1(deviceVisReal02, deviceVisImag02,
			       ctx, nrCells, nrChannels, nrPolarizations);
	    } else if(cellSize == CELL_2X2) {
		downloadVis2x2(deviceVisReal02, deviceVisImag02, deviceVisReal12, deviceVisImag12,
			       deviceVisReal03, deviceVisImag03, deviceVisReal13, deviceVisImag13, 
			       ctx, nrCells, nrChannels, nrPolarizations);
	    } else {
		std::cout << "illegal cell size, only 8 output registers available, cannot do more than 2x2" << std::endl;
	    }
	}
	currBuf = nextBuf;
	nextBuf = !currBuf;
    }

    runTimer.stop();
	
    calCtxEndCounterExt(ctx, idleCounter);
    calCtxEndCounterExt(ctx, cacheCounter);

    CALfloat idleResult;
    if(calCtxGetCounterExt(&idleResult, ctx, idleCounter) != CAL_RESULT_OK) {
	std::cout << "eek" << std::endl;
	return 1;
    }
    CALfloat cacheResult;
    if(calCtxGetCounterExt(&cacheResult, ctx, cacheCounter) != CAL_RESULT_OK) {
	std::cout << "eek" << std::endl;
	return 1;
    }
    
    unsigned long long loadedBytes = (unsigned long long)iters * nrChannels * nrCells * (cellWidth + cellHeight) * nrTimes * sizeof(float4);
    double gb = (double) loadedBytes / (1024.0 * 1024.0 * 1024.0);
    
    // 4 mads x 2 ops (mul + add) x 4 regs (float4)
    // 4 * 2 * 4 * cellW * cellH
    
    unsigned long long ops =  (unsigned long long) iters * nrChannels * nrCells * nrTimes * 4L * 2L * 4L * cellWidth * cellHeight;

    std::cout << "ops: " << ops << std::endl;

    double elapsed = runTimer.getTimeInSeconds();
    double flops = (ops / elapsed) / 1000000000.0;
    double kernelFlops = (ops / kernelTimer.getTimeInSeconds()) / 1000000000.0;

    double opsPerByte = (double) ops / (double)loadedBytes;
    realSampleGb *= iters;
    if(sync) {
	std::cout << "sync version" << std::endl;
    } else {
	std::cout << "async version" << std::endl;
    }
    std::cout << "Used " << cellWidth << "X" << cellHeight << " correlator, which has " << opsPerByte << " ops per byte" << std::endl;
    std::cout << "correlate kernel took " << kernelTimer.getTimeInSeconds() << " s, achieved " << kernelFlops << " Gflops" <<std::endl;
    std::cout << "entire correlate including mem I/O took " << elapsed << " s, achieved " << flops << " Gflops" <<std::endl;
    std::cout << "loaded " << realSampleGb << " GB sample data from host to GPU, " << (realSampleGb / elapsed) << "GB/s" << std::endl;
    std::cout << "loaded " << gb << " GB from GPU memory to regs, " << (gb / kernelTimer.getTimeInSeconds()) << " GB/s" << ", reload factor = " << (gb / realSampleGb) << std::endl;
    std::cout << "kernel time: " << kernelTimer.getTimeInSeconds() << ", " << (kernelTimer.getTimeInSeconds() / elapsed * 100.0) << " %" << std::endl;
    std::cout << "total wait time: " << waitTimer.getTimeInSeconds() << ", " << (waitTimer.getTimeInSeconds() / elapsed * 100.0) << " %" << std::endl;
    std::cout << "processor was " << (idleResult*100.0) << " % idle" << std::endl;
    std::cout << "cache hit rate was " << (cacheResult*100.0) << " %" << std::endl;

    if(!useOutputRegisters && print) {
	CALmem hostVisibilitiesMem = 0;
	float4* fdataVis = NULL;
	cal_check(calCtxGetMem(&hostVisibilitiesMem, ctx, hostVisibilitiesRes));
	cal_check(calResMap((CALvoid**)&fdataVis, &pitch, hostVisibilitiesRes, 0));
	
	for(unsigned i=0; i< nrCells * nrOutputsPerCell * nrChannels; i++) {
	    float4 val = fdataVis[i];
	    fprintf(stderr, "index %5d: %8.1fr %8.1fi %8.1fr %8.1fi\n", 
		    i, val.x, val.y, val.z, val.w);
	}
	cal_check(calResUnmap(hostVisibilitiesRes));
	cal_check(calCtxReleaseMem(ctx, hostVisibilitiesMem));
    }


    //-------------------------------------------------------------------------
    // Cleaning up
    //-------------------------------------------------------------------------

    calCtxDestroyCounterExt(ctx, idleCounter);
    calCtxDestroyCounterExt(ctx, cacheCounter);

    // Unloading the module
    calModuleUnload(ctx, module);

    // Freeing compiled program binary
    calclFreeImage(image);
    calclFreeObject(obj);

    // Releasing resource from context
    calCtxReleaseMem(ctx, hostCellToStatX);
    calCtxReleaseMem(ctx, hostCellToStatY);
    calCtxReleaseMem(ctx, hostSamplesMem[0]);
    calCtxReleaseMem(ctx, hostSamplesMem[1]);
    calCtxReleaseMem(ctx, deviceSamplesMem[0]);
    calCtxReleaseMem(ctx, deviceSamplesMem[1]);

    // Deallocating resources
    calResFree(hostSamplesRes[0]);
    calResFree(hostSamplesRes[1]);
    calResFree(deviceSamplesRes[0]);
    calResFree(deviceSamplesRes[1]);
    calResFree(deviceCellToStatXRes);
    calResFree(deviceCellToStatYRes);


    if(useOutputRegisters) {
	if(cellSize == CELL_1X1) {
	    calResFree(deviceVisReal02);
	    calResFree(deviceVisImag02);
	} else if(cellSize == CELL_2X2) {
	    calResFree(deviceVisReal02);
	    calResFree(deviceVisImag02);
	    calResFree(deviceVisReal12);
	    calResFree(deviceVisImag12);
	    calResFree(deviceVisReal03);
	    calResFree(deviceVisImag03);
	    calResFree(deviceVisReal13);
	    calResFree(deviceVisImag13);
	} else {
	    std::cout << "illegal cell size, only 8 output registers available, cannot do more than 2x2" << std::endl;
	}
    }

    // Destroying context
    calCtxDestroy(ctx);

    // Closing device
    calDeviceClose(device);

    // Shutting down CAL
    calShutdown();
}
