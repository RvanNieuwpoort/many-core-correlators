#include "../common/common.h"
#include "riscv_correlator.h"

#include <cstring> 
#include <complex>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <riscv_vector.h>

#define PRINT_RESULT 0
#define COMPUTE_GFLOPS 0

#define VECTOR_WIDTH_IN_FLOATS 8
#define ALIGNMENT 256 // in bytes
#define REAL_PEAK_GFLOPS_PER_CORE (1.6 /*GHz*/ * 8 /*vector width*/ * 2 /*fma*/)


using namespace std;

static const bool validateResults = true;
static const bool verbose = true;
static const unsigned nrStations = 64;
static const unsigned nrBaselines = NR_BASELINES(nrStations);
static const unsigned nrTimes = 768, nrTimesWidth = 768; // 770
static const unsigned nrChannels = 256;
static const unsigned iter = 1;
static const unsigned nrThreads = 8;

vfloat32m8_t init_vfloat32m8(const float val)
{
    float result[VECTOR_WIDTH_IN_FLOATS];
    for(int i=0; i<VECTOR_WIDTH_IN_FLOATS; i++) {
	result[i] = val;
    }

    return __riscv_vle32_v_f32m8(result, VECTOR_WIDTH_IN_FLOATS);
}

void print_vfloat32m8(const vfloat32m8_t f, const char* name)
{
    float result[VECTOR_WIDTH_IN_FLOATS];
    for(int i=0; i<VECTOR_WIDTH_IN_FLOATS; i++) {
	result[i] = -1.0f;
    }

    __riscv_vse32_v_f32m8(result, f, VECTOR_WIDTH_IN_FLOATS);
#if 0
    cout << "float vec " << name << " = { ";
    for(int i=0; i<VECTOR_WIDTH_IN_FLOATS; i++) {
	cout << result[i] << " ";
    }
    cout << " }" << endl;
#endif
}

void* calcMaxFlops(void* data)
{
    // Do a simple vectorized computation here to compute the maximum peak performance.
    // Vectors are 256 bit, so 8 32 bit floats.

    vfloat32m8_t va = init_vfloat32m8(0.2E-9f);
    vfloat32m8_t vb = init_vfloat32m8(0.3E-9f);
    vfloat32m8_t vc = init_vfloat32m8(0.4E-9f);
    vfloat32m8_t vd = init_vfloat32m8(0.5E-9f);
    vfloat32m8_t ve = init_vfloat32m8(0.6E-9f);
    vfloat32m8_t vf = init_vfloat32m8(0.7E-9f);
    vfloat32m8_t vg = init_vfloat32m8(0.8E-9f);
    vfloat32m8_t vh = init_vfloat32m8(0.9E-9f);
    vfloat32m8_t vi = init_vfloat32m8(0.10E-9f);
    vfloat32m8_t vj = init_vfloat32m8(0.11E-9f);
    vfloat32m8_t vk = init_vfloat32m8(0.12E-9f);
    vfloat32m8_t vl = init_vfloat32m8(0.13E-9f);
    vfloat32m8_t vm = init_vfloat32m8(0.14E-9f);
    vfloat32m8_t vn = init_vfloat32m8(0.15E-9f);
    vfloat32m8_t vo = init_vfloat32m8(0.16E-9f);
    vfloat32m8_t vp = init_vfloat32m8(0.17E-9f);

    // We do 16 ops per line in the loop. 1 FMA * 8 floats.
    for (unsigned long long x = 0; x < 1E9L / (VECTOR_WIDTH_IN_FLOATS*2*16); x++) {
	va = __riscv_vfmacc_vv_f32m8(va, va, va, VECTOR_WIDTH_IN_FLOATS);
	vb = __riscv_vfmacc_vv_f32m8(vb, vb, vb, VECTOR_WIDTH_IN_FLOATS);
	vc = __riscv_vfmacc_vv_f32m8(vc, vc, vc, VECTOR_WIDTH_IN_FLOATS);
	vd = __riscv_vfmacc_vv_f32m8(vd, vd, vd, VECTOR_WIDTH_IN_FLOATS);
	ve = __riscv_vfmacc_vv_f32m8(ve, ve, ve, VECTOR_WIDTH_IN_FLOATS);
	vf = __riscv_vfmacc_vv_f32m8(vf, vf, vf, VECTOR_WIDTH_IN_FLOATS);
	vg = __riscv_vfmacc_vv_f32m8(vg, vg, vg, VECTOR_WIDTH_IN_FLOATS);
	vh = __riscv_vfmacc_vv_f32m8(vh, vh, vh, VECTOR_WIDTH_IN_FLOATS);
	vi = __riscv_vfmacc_vv_f32m8(vi, vi, vi, VECTOR_WIDTH_IN_FLOATS);
	vj = __riscv_vfmacc_vv_f32m8(vj, vj, vj, VECTOR_WIDTH_IN_FLOATS);
	vk = __riscv_vfmacc_vv_f32m8(vk, vk, vk, VECTOR_WIDTH_IN_FLOATS);
	vl = __riscv_vfmacc_vv_f32m8(vl, vl, vl, VECTOR_WIDTH_IN_FLOATS);
	vm = __riscv_vfmacc_vv_f32m8(vm, vm, vm, VECTOR_WIDTH_IN_FLOATS);
	vn = __riscv_vfmacc_vv_f32m8(vn, vn, vn, VECTOR_WIDTH_IN_FLOATS);
	vo = __riscv_vfmacc_vv_f32m8(vo, vo, vo, VECTOR_WIDTH_IN_FLOATS);
	vp = __riscv_vfmacc_vv_f32m8(vp, vp, vp, VECTOR_WIDTH_IN_FLOATS);
    }


    FlopTestParams* p = (FlopTestParams*) data;
    p->gFlops = 1; 
    
    // Print results just to make sure the compiler doesn't optimize everything away.
    print_vfloat32m8(va, "va");
    print_vfloat32m8(vb, "vb");
    print_vfloat32m8(vc, "vc");
    print_vfloat32m8(vd, "vd");
    print_vfloat32m8(ve, "ve");
    print_vfloat32m8(vf, "vf");
    print_vfloat32m8(vg, "vg");
    print_vfloat32m8(vh, "vh");
    print_vfloat32m8(vi, "vi");
    print_vfloat32m8(vj, "vj");
    print_vfloat32m8(vk, "vk");
    print_vfloat32m8(vl, "vl");
    print_vfloat32m8(vm, "vm");
    print_vfloat32m8(vn, "vn");
    print_vfloat32m8(vo, "vo");
    print_vfloat32m8(vp, "vp");

    return 0;
}

void printCorrelatorType(int correlatorType)
{
    switch (correlatorType) {
    case CORRELATOR_REFERENCE:
	cout << "reference";
	break;
    case CORRELATOR_1X1:
	cout << "1x1";
	break;
    case CORRELATOR_2X1:
	cout << "2x1";
	break;
    case CORRELATOR_2X2:
	cout << "2x2";
	break;
    default:
	cout << "illegal correlator" << endl;
	exit(66);
    }
}

void* runCorrelator(void* data)
{
    CorrelatorParams* p = (CorrelatorParams*) data;

    for(unsigned i=0; i<iter; i++) {
	switch (p->correlatorType) {
	case CORRELATOR_REFERENCE:
	    p->ops = referenceCorrelator(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_1X1:
	    p->ops = cpuCorrelator_1x1(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_2X1:
	    p->ops = cpuCorrelator_2x1(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_2X2:
	    p->ops = cpuCorrelator_2x2(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	default:
	    cout << "illegal correlator" << endl;
	    exit(66);
	}
    }

    p->ops *= iter;
    p->bytesLoaded *= iter;
    p->bytesStored *= iter;

#if PRINT_RESULT
    printResult(visibilities, nrChannels, nrBaselines);
#endif

    return 0;
}

int main()
{
    // set the vector length
    size_t vl = __riscv_vsetvl_e32m8(VECTOR_WIDTH_IN_FLOATS);
    cout << "vector length = " << vl << " 32 bit floats, ";
    cout << "timer resolution is " << chrono::high_resolution_clock::period::den << " ticks per second" << endl;
    cout << "Configuration: " << nrStations << " stations, " << nrBaselines <<
	" baselines, " << nrChannels << " channels, " << nrTimes << " time samples" << endl;

    const unsigned arraySize = nrStations*nrChannels*nrTimesWidth*NR_POLARIZATIONS*2;
    const unsigned visArraySize = nrBaselines*nrChannels*NR_POLARIZATIONS* NR_POLARIZATIONS*2;

    cout << "sample array size = " << ((arraySize * nrThreads * sizeof(float))/(1024*1024)) << " mbytes, "
	 << "vis array size = " << ((visArraySize * nrThreads * sizeof(float))/(1024*1024)) << " mbytes" << endl;

#if COMPUTE_GFLOPS
    double maxGflops = computeMaxGflops(nrThreads, calcMaxFlops, NULL);
#else 
    double maxGflops = REAL_PEAK_GFLOPS_PER_CORE * nrThreads;
#endif // COMPUTE_GFLOPS
    cout << "peak flops with " << nrThreads << " threads is: " << maxGflops << " gflops" << std::endl;
    
    float* samples = new (align_val_t{ALIGNMENT}) float[nrThreads*arraySize];
    float* visibilities= new (align_val_t{ALIGNMENT}) float[nrThreads*visArraySize];
    memset(visibilities, 0, nrThreads*visArraySize*sizeof(float));

    initSamples(samples, nrThreads, nrTimes, nrTimesWidth, nrStations, nrChannels, arraySize);

    spawnCorrelatorThreads(CORRELATOR_1X1, runCorrelator, samples, arraySize,
			   visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			   nrThreads, iter, maxGflops, verbose, validateResults);

    spawnCorrelatorThreads(CORRELATOR_2X1, runCorrelator, samples, arraySize,
			   visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			   nrThreads, iter, maxGflops, verbose, validateResults);

    spawnCorrelatorThreads(CORRELATOR_2X2, runCorrelator, samples, arraySize,
			   visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			   nrThreads, iter, maxGflops, verbose, validateResults);

    endCommon();
    
    delete[] samples;
    delete[] visibilities;

    return 0;
}
