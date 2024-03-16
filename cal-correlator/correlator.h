#ifndef _CORRELATOR_H_
#define _CORRELATOR_H_

#include "cal.h"
#include "calcl.h"
#include "cal_ext.h"
#include "cal_ext_counter.h"

#define CELL_1X1 1
#define CELL_2X2 2
#define CELL_3X2 3
#define CELL_3X3 4
#define CELL_4X3 5

#define DO_LOADS 1
#define DO_COMPUTATION 1

typedef enum VarTypeEmum {
    VAR_INPUT,
    VAR_OUTPUT,
    VAR_CONSTANT,
    VAR_GLOBAL
} VarType;

typedef struct float4Rec {
    CALfloat x; 
    CALfloat y; 
    CALfloat z; 
    CALfloat w; 
} float4;

extern bool disassemble;
extern bool print;

extern PFNCALCTXCREATECOUNTER  calCtxCreateCounterExt;
extern PFNCALCTXDESTROYCOUNTER calCtxDestroyCounterExt;
extern PFNCALCTXBEGINCOUNTER   calCtxBeginCounterExt;
extern PFNCALCTXENDCOUNTER     calCtxEndCounterExt;
extern PFNCALCTXGETCOUNTER     calCtxGetCounterExt;

extern void cal_check(CALresult r);
extern void __logger(const CALchar *msg);
extern unsigned fillCellToStatTable(unsigned w, unsigned h, unsigned nrStations, float* cellToStatX, float* cellToStatY);
extern unsigned power (unsigned base, unsigned n);
extern unsigned calcNrCells(unsigned w, unsigned h, unsigned nrStations);
extern unsigned calcNrOutputsPerCell(unsigned w, unsigned h, unsigned nrPolarizations);
extern void printResult(int error, unsigned pitch, unsigned nrCells, unsigned nrChannels, unsigned nrPolarizations, 
		 float4* visReal02, float4* visImag02,
		 float4* visReal12, float4* visImag12,
		 float4* visReal03, float4* visImag03,
		 float4* visReal13, float4* visImag13);
extern void initSamples(unsigned nrStations, unsigned nrChannels, unsigned nrTimes, unsigned pitch, float4* samples);
extern void bindVariable(CALcontext* ctx, CALmodule* module, CALresource res, int varNumber, int varType);
extern void initCounters();
extern void wipeOutputMem(CALcontext* ctx, CALresource res, unsigned heigth);

extern void allocateVis1X1Reg(CALresource* deviceVisReal02, CALresource* deviceVisImag02,
			      CALcontext ctx, CALdevice device, unsigned nrCells, unsigned nrChannels);

extern void downloadVis1x1(CALresource deviceVisReal02, CALresource deviceVisImag02,
			   CALcontext ctx, unsigned nrCells, unsigned nrChannels, unsigned nrPolarizations);

extern void allocateVis2X2Reg(CALresource* deviceVisReal02, CALresource* deviceVisImag02,
			      CALresource* deviceVisReal12, CALresource* deviceVisImag12,
			      CALresource* deviceVisReal03, CALresource* deviceVisImag03,
			      CALresource* deviceVisReal13, CALresource* deviceVisImag13,
			      CALcontext ctx, CALdevice device, unsigned nrCells, unsigned nrChannels);

extern void downloadVis2x2(CALresource deviceVisReal02, CALresource deviceVisImag02,
			   CALresource deviceVisReal12, CALresource deviceVisImag12,
			   CALresource deviceVisReal03, CALresource deviceVisImag03,
			   CALresource deviceVisReal13, CALresource deviceVisImag13,
			   CALcontext ctx, unsigned nrCells, unsigned nrChannels, unsigned nrPolarizations);

extern void copy(CALresource src, CALresource dst, CALcontext ctx);

inline void copyMem(CALmem srcMem, CALmem dstMem, CALcontext ctx)
{
    CALevent e = 0;
    cal_check(calMemCopy(&e, ctx, srcMem, dstMem, 0));
    while(calCtxIsEventDone(ctx, e) == CAL_RESULT_PENDING);
}

inline CALevent copyMemAsync(CALmem srcMem, CALmem dstMem, CALcontext ctx)
{
    CALevent e = 0;
    cal_check(calMemCopy(&e, ctx, srcMem, dstMem, 0));
    return e;
}


#endif // _CORRELATOR_H_
