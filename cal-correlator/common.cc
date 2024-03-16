#include <stdio.h>
#include <iostream>
#include <cassert>
#include <stdlib.h>
#include <string.h>

#include "correlator.h"

bool disassemble;
bool print;

PFNCALCTXCREATECOUNTER  calCtxCreateCounterExt;
PFNCALCTXDESTROYCOUNTER calCtxDestroyCounterExt;
PFNCALCTXBEGINCOUNTER   calCtxBeginCounterExt;
PFNCALCTXENDCOUNTER     calCtxEndCounterExt;
PFNCALCTXGETCOUNTER     calCtxGetCounterExt;

void cal_check(CALresult r)
{
    if (r != CAL_RESULT_OK) {
        fprintf(stdout, "Error detected: %s.\n", calGetErrorString());
        abort();
    }
}

void __logger(const CALchar *msg)
{
    fprintf(stderr, "%s\n", msg);
}

unsigned fillCellToStatTable(unsigned w, unsigned h, unsigned nrStations, float* cellToStatX, float* cellToStatY)
{
    unsigned nrCells;
    for (unsigned statY = nrStations % 2 ? 1 : 2, nrCells = 0; statY < nrStations; statY += h) {
	for (unsigned statX = 0; statX + w <= statY; statX += w, nrCells ++) {
	    cellToStatX[nrCells] = statX;
	    cellToStatY[nrCells] = statY;
	}
    }
    return nrCells;
}

unsigned calcNrCells(unsigned w, unsigned h, unsigned nrStations)
{
    unsigned nrCells = 0;
    for (unsigned statY = nrStations % 2 ? 1 : 2; statY < nrStations; statY += h) {
	for (unsigned statX = 0; statX + w <= statY; statX += w) {
	    nrCells++;
	}
    }
    return nrCells;
}

unsigned power (unsigned base, unsigned n) {
    unsigned p = 1;
    for (unsigned i = 1; i <= n; ++i)
	p *= base;
    return p;
}

unsigned calcNrOutputsPerCell(unsigned w, unsigned h, unsigned nrPolarizations)
{
    return w * h * nrPolarizations;
}

void printResult(int error, unsigned pitch, unsigned nrCells, unsigned nrChannels, unsigned nrPolarizations, 
		 float4* visReal02, float4* visImag02,
		 float4* visReal12, float4* visImag12,
		 float4* visReal03, float4* visImag03,
		 float4* visReal13, float4* visImag13)
{
// vis: nrCells x nrChannels x float4
#if 1
    fprintf(stderr, "vis mode\n");
    for(unsigned channel = 0; channel < nrChannels; channel++) {
      float4* tmpReal02 = &visReal02[channel * pitch];
      float4* tmpImag02 = &visImag02[channel * pitch];
      float4* tmpReal12 = &visReal12[channel * pitch];
      float4* tmpImag12 = &visImag12[channel * pitch];
      float4* tmpReal03 = &visReal03[channel * pitch];
      float4* tmpImag03 = &visImag03[channel * pitch];
      float4* tmpReal13 = &visReal13[channel * pitch];
      float4* tmpImag13 = &visImag13[channel * pitch];

  for(unsigned cell = 0; cell < nrCells; cell++) {
	fprintf(stderr, "chan %d, cell %d 02(%8.1fr %8.1fi %8.1fr %8.1fi %8.1fr %8.1fi %8.1fr %8.1fi)\n", 
		channel, cell,
		tmpReal02[cell].x, tmpImag02[cell].x, 
		tmpReal02[cell].y, tmpImag02[cell].y,
		tmpReal02[cell].z, tmpImag02[cell].z, 
		tmpReal02[cell].w, tmpImag02[cell].w);
	fprintf(stderr, "chan %d, cell %d 12(%8.1fr %8.1fi %8.1fr %8.1fi %8.1fr %8.1fi %8.1fr %8.1fi)\n", 
		channel, cell,
		tmpReal12[cell].x, tmpImag12[cell].x, 
		tmpReal12[cell].y, tmpImag12[cell].y,
		tmpReal12[cell].z, tmpImag12[cell].z, 
		tmpReal12[cell].w, tmpImag12[cell].w);
	fprintf(stderr, "chan %d, cell %d 03(%8.1fr %8.1fi %8.1fr %8.1fi %8.1fr %8.1fi %8.1fr %8.1fi)\n", 
		channel, cell,
		tmpReal03[cell].x, tmpImag03[cell].x, 
		tmpReal03[cell].y, tmpImag03[cell].y,
		tmpReal03[cell].z, tmpImag03[cell].z, 
		tmpReal03[cell].w, tmpImag03[cell].w);
	fprintf(stderr, "chan %d, cell %d 13(%8.1fr %8.1fi %8.1fr %8.1fi %8.1fr %8.1fi %8.1fr %8.1fi)\n", 
		channel, cell,
		tmpReal13[cell].x, tmpImag13[cell].x, 
		tmpReal13[cell].y, tmpImag13[cell].y,
		tmpReal13[cell].z, tmpImag13[cell].z, 
		tmpReal13[cell].w, tmpImag13[cell].w);
    }
  }
#endif
#if 1
    fprintf(stderr, "reg mode\n");
    for(unsigned channel = 0; channel < nrChannels; channel++) {
      float4* tmpReal02 = &visReal02[channel * pitch];
      float4* tmpImag02 = &visImag02[channel * pitch];
      float4* tmpReal12 = &visReal12[channel * pitch];
      float4* tmpImag12 = &visImag12[channel * pitch];
      float4* tmpReal03 = &visReal03[channel * pitch];
      float4* tmpImag03 = &visImag03[channel * pitch];
      float4* tmpReal13 = &visReal13[channel * pitch];
      float4* tmpImag13 = &visImag13[channel * pitch];

      for(unsigned cell = 0; cell < nrCells; cell++) {
	fprintf(stderr, "chan %d, cell %d o0(%8.1f %8.1f %8.1f %8.1f)\n", 
		channel, cell,
		tmpReal02[cell].x, tmpReal02[cell].y, 
		tmpReal02[cell].z, tmpReal02[cell].w);
	fprintf(stderr, "chan %d, cell %d o1(%8.1f %8.1f %8.1f %8.1f)\n", 
		channel, cell,
		tmpImag02[cell].x, tmpImag02[cell].y, 
		tmpImag02[cell].z, tmpImag02[cell].w);
	fprintf(stderr, "chan %d, cell %d o2(%8.1f %8.1f %8.1f %8.1f)\n", 
		channel, cell,
		tmpReal12[cell].x, tmpReal12[cell].y, 
		tmpReal12[cell].z, tmpReal12[cell].w);
	fprintf(stderr, "chan %d, cell %d o3(%8.1f %8.1f %8.1f %8.1f)\n", 
		channel, cell,
		tmpImag12[cell].x, tmpImag12[cell].y, 
		tmpImag12[cell].z, tmpImag12[cell].w);
	fprintf(stderr, "chan %d, cell %d o4(%8.1f %8.1f %8.1f %8.1f)\n", 
		channel, cell,
		tmpReal03[cell].x, tmpReal03[cell].y, 
		tmpReal03[cell].z, tmpReal03[cell].w);
	fprintf(stderr, "chan %d, cell %d o5(%8.1f %8.1f %8.1f %8.1f)\n", 
		channel, cell,
		tmpImag03[cell].x, tmpImag03[cell].y, 
		tmpImag03[cell].z, tmpImag03[cell].w);
	fprintf(stderr, "chan %d, cell %d o6(%8.1f %8.1f %8.1f %8.1f)\n", 
		channel, cell,
		tmpReal13[cell].x, tmpReal13[cell].y, 
		tmpReal13[cell].z, tmpReal13[cell].w);
	fprintf(stderr, "chan %d, cell %d o7(%8.1f %8.1f %8.1f %8.1f)\n", 
		channel, cell,
		tmpImag13[cell].x, tmpImag13[cell].y, 
		tmpImag13[cell].z, tmpImag13[cell].w);
    }
  }
#endif
}

void initSamples(unsigned nrStations, unsigned nrChannels, unsigned nrTimes, unsigned pitch, float4* samples)
{
    // use 2D sample buffer:
    // - elts are float4, representing both polarizations X real/imag
    // - layout is station x channel x nrTimes
    // in 2D, we see (channel, station) as one dimension, nrTimes as the other.

    unsigned yMax = nrChannels * nrStations;

    for(unsigned y = 0; y<yMax; y++) {
	float4* tmpSamples = &samples[y * pitch];
	for(unsigned time=0; time<nrTimes; time++) {
	    tmpSamples[time].x = time % 8;
	    tmpSamples[time].y = y;
	    tmpSamples[time].z = time % 8 + 1;
	    tmpSamples[time].w = y+1;
 	}
    }
}

void bindVariable(CALcontext* ctx, CALmodule* module, CALresource res, int varNumber, int varType)
{
    char buffer[100];

    switch(varType) {
    case VAR_INPUT:
	sprintf(buffer, "i%d", varNumber);
	break;
    case VAR_OUTPUT:
	sprintf(buffer, "o%d", varNumber);
	break;
    case VAR_CONSTANT:
	sprintf(buffer, "cb%d", varNumber);
	break;
    case VAR_GLOBAL:
	sprintf(buffer, "g[]");
	break;
    default:
        fprintf(stdout, "Unknown variable type %d\n", varType);
	exit(1);	
    }

    CALname name;
    if (calModuleGetName(&name, *ctx, *module, buffer) != CAL_RESULT_OK) {
        fprintf(stdout, "Could not locate variable %s: %s Exiting.\n", buffer, calGetErrorString());
	exit(1);
    }
    CALmem mem = 0;
    if(calCtxGetMem(&mem, *ctx, res) != CAL_RESULT_OK) {
	fprintf(stderr, "Could not get mem for %s: %s\n", buffer, calGetErrorString());
	exit(1);
    }

    if(calCtxSetMem(*ctx, name, mem) != CAL_RESULT_OK) {
	fprintf(stderr, "Could not set mem for %s: %s\n", buffer, calGetErrorString());
	exit(1);
    }
}

void initCounters()
{
    if (calExtSupported((CALextid)CAL_EXT_COUNTERS) != CAL_RESULT_OK) {
        exit(1);
    }

    if (calExtGetProc((CALextproc*)&calCtxCreateCounterExt, (CALextid)CAL_EXT_COUNTERS, "calCtxCreateCounter")) {
        exit(1);
    }

    if (calExtGetProc((CALextproc*)&calCtxDestroyCounterExt, (CALextid)CAL_EXT_COUNTERS, "calCtxDestroyCounter")) {
        exit(1);
    }

    if (calExtGetProc((CALextproc*)&calCtxBeginCounterExt, (CALextid)CAL_EXT_COUNTERS, "calCtxBeginCounter")) {
        exit(1);
    }

    if (calExtGetProc((CALextproc*)&calCtxEndCounterExt, (CALextid)CAL_EXT_COUNTERS, "calCtxEndCounter")) {
        exit(1);
    }

    if (calExtGetProc((CALextproc*)&calCtxGetCounterExt, (CALextid)CAL_EXT_COUNTERS, "calCtxGetCounter")) {
        exit(1);
    }
}

void wipeOutputMem(CALcontext* ctx, CALresource res, unsigned heigth)
{
    CALmem mem = 0;
    CALfloat* data = NULL;
    CALuint pitch = 0;

    calCtxGetMem(&mem, *ctx, res);
    calResMap((CALvoid**)&data, &pitch, res, 0);
    memset(data, 0, pitch * heigth * sizeof(float));
    calResUnmap(res);
}

void allocateVis2X2Reg(CALresource* deviceVisReal02, CALresource* deviceVisImag02,
		       CALresource* deviceVisReal12, CALresource* deviceVisImag12,
		       CALresource* deviceVisReal03, CALresource* deviceVisImag03,
		       CALresource* deviceVisReal13, CALresource* deviceVisImag13,
		       CALcontext ctx, CALdevice device, unsigned nrCells, unsigned nrChannels)
{
    cal_check(calResAllocLocal2D(deviceVisReal02, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(deviceVisImag02, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(deviceVisReal12, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(deviceVisImag12, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(deviceVisReal03, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(deviceVisImag03, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(deviceVisReal13, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(deviceVisImag13, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));

    // wipe output buffers
    wipeOutputMem(&ctx, *deviceVisReal02, nrChannels);
    wipeOutputMem(&ctx, *deviceVisImag02, nrChannels);
    wipeOutputMem(&ctx, *deviceVisReal12, nrChannels);
    wipeOutputMem(&ctx, *deviceVisImag12, nrChannels);
    wipeOutputMem(&ctx, *deviceVisReal03, nrChannels);
    wipeOutputMem(&ctx, *deviceVisImag03, nrChannels);
    wipeOutputMem(&ctx, *deviceVisReal13, nrChannels);
    wipeOutputMem(&ctx, *deviceVisImag13, nrChannels);
}

void allocateVis1X1Reg(CALresource* deviceVisReal02, CALresource* deviceVisImag02,
		       CALcontext ctx, CALdevice device, unsigned nrCells, unsigned nrChannels)
{
    cal_check(calResAllocLocal2D(deviceVisReal02, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));
    cal_check(calResAllocLocal2D(deviceVisImag02, device, nrCells, nrChannels, CAL_FORMAT_FLOAT_4, 0));

    // wipe output buffers
    wipeOutputMem(&ctx, *deviceVisReal02, nrChannels);
    wipeOutputMem(&ctx, *deviceVisImag02, nrChannels);
}

void downloadVis1x1(CALresource deviceVisReal02, CALresource deviceVisImag02,
		    CALcontext ctx, unsigned nrCells, unsigned nrChannels, unsigned nrPolarizations)
{
    float4* fdataReal02 = NULL, *fdataImag02 = NULL;
    CALuint pitchReal = 0, pitchImag = 0;
    CALmem hostVisReal02 = 0, hostVisImag02 = 0;
    cal_check(calCtxGetMem(&hostVisReal02, ctx, deviceVisReal02));
    cal_check(calCtxGetMem(&hostVisImag02, ctx, deviceVisImag02));

    cal_check(calResMap((CALvoid**)&fdataReal02, &pitchReal, deviceVisReal02, 0));
    cal_check(calResMap((CALvoid**)&fdataImag02, &pitchImag, deviceVisImag02, 0));

    cal_check(calResUnmap(deviceVisReal02));
    cal_check(calResUnmap(deviceVisImag02));

    cal_check(calCtxReleaseMem(ctx, hostVisReal02));
    cal_check(calCtxReleaseMem(ctx, hostVisImag02));
}

void downloadVis2x2(CALresource deviceVisReal02, CALresource deviceVisImag02,
		    CALresource deviceVisReal12, CALresource deviceVisImag12,
		    CALresource deviceVisReal03, CALresource deviceVisImag03,
		    CALresource deviceVisReal13, CALresource deviceVisImag13,
		    CALcontext ctx, unsigned nrCells, unsigned nrChannels, unsigned nrPolarizations)
{
    float4* fdataReal02 = NULL, *fdataImag02 = NULL;
    float4* fdataReal12 = NULL, *fdataImag12 = NULL;
    float4* fdataReal03 = NULL, *fdataImag03 = NULL;
    float4* fdataReal13 = NULL, *fdataImag13 = NULL;
    CALuint pitchReal = 0, pitchImag = 0;
    CALmem hostVisReal02 = 0, hostVisImag02 = 0;
    CALmem hostVisReal12 = 0, hostVisImag12 = 0;
    CALmem hostVisReal03 = 0, hostVisImag03 = 0;
    CALmem hostVisReal13 = 0, hostVisImag13 = 0;
    cal_check(calCtxGetMem(&hostVisReal02, ctx, deviceVisReal02));
    cal_check(calCtxGetMem(&hostVisImag02, ctx, deviceVisImag02));
    cal_check(calCtxGetMem(&hostVisReal12, ctx, deviceVisReal12));
    cal_check(calCtxGetMem(&hostVisImag12, ctx, deviceVisImag12));
    cal_check(calCtxGetMem(&hostVisReal03, ctx, deviceVisReal03));
    cal_check(calCtxGetMem(&hostVisImag03, ctx, deviceVisImag03));
    cal_check(calCtxGetMem(&hostVisReal13, ctx, deviceVisReal13));
    cal_check(calCtxGetMem(&hostVisImag13, ctx, deviceVisImag13));

    cal_check(calResMap((CALvoid**)&fdataReal02, &pitchReal, deviceVisReal02, 0));
    cal_check(calResMap((CALvoid**)&fdataImag02, &pitchImag, deviceVisImag02, 0));
    cal_check(calResMap((CALvoid**)&fdataReal12, &pitchReal, deviceVisReal12, 0));
    cal_check(calResMap((CALvoid**)&fdataImag12, &pitchImag, deviceVisImag12, 0));
    cal_check(calResMap((CALvoid**)&fdataReal03, &pitchReal, deviceVisReal03, 0));
    cal_check(calResMap((CALvoid**)&fdataImag03, &pitchImag, deviceVisImag03, 0));
    cal_check(calResMap((CALvoid**)&fdataReal13, &pitchReal, deviceVisReal13, 0));
    cal_check(calResMap((CALvoid**)&fdataImag13, &pitchImag, deviceVisImag13, 0));

    if(print) {
	printResult(0, pitchReal, nrCells, nrChannels, nrPolarizations, 
		    fdataReal02, fdataImag02, fdataReal12, fdataImag12,
		    fdataReal03, fdataImag03, fdataReal13, fdataImag13);
    }

    cal_check(calResUnmap(deviceVisReal02));
    cal_check(calResUnmap(deviceVisImag02));
    cal_check(calResUnmap(deviceVisReal12));
    cal_check(calResUnmap(deviceVisImag12));
    cal_check(calResUnmap(deviceVisReal03));
    cal_check(calResUnmap(deviceVisImag03));
    cal_check(calResUnmap(deviceVisReal13));
    cal_check(calResUnmap(deviceVisImag13));

    cal_check(calCtxReleaseMem(ctx, hostVisReal02));
    cal_check(calCtxReleaseMem(ctx, hostVisImag02));
    cal_check(calCtxReleaseMem(ctx, hostVisReal12));
    cal_check(calCtxReleaseMem(ctx, hostVisImag12));
    cal_check(calCtxReleaseMem(ctx, hostVisReal03));
    cal_check(calCtxReleaseMem(ctx, hostVisImag03));
    cal_check(calCtxReleaseMem(ctx, hostVisReal13));
    cal_check(calCtxReleaseMem(ctx, hostVisImag13));
}

void copy(CALresource src, CALresource dst, CALcontext ctx)
{
    CALmem srcMem = 0;
    CALmem dstMem = 0;

    cal_check(calCtxGetMem(&srcMem, ctx, src));
    cal_check(calCtxGetMem(&dstMem, ctx, dst));

    CALevent e = 0;
    cal_check(calMemCopy(&e, ctx, srcMem, dstMem, 0));

    while(calCtxIsEventDone(ctx, e) == CAL_RESULT_PENDING);

    cal_check(calCtxReleaseMem(ctx, srcMem));
    cal_check(calCtxReleaseMem(ctx, dstMem));
}
