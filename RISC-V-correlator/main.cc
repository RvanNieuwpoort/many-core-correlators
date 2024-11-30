#include "riscv_correlator.h"

#include <cstring> // for memset
#include <complex>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <riscv_vector.h>

#define VECTOR_WIDTH_IN_FLOATS 8
#define PRINT_RESULT 0
#define ALIGNMENT 256 // in bytes


using namespace std;

const unsigned nrStations = 32;
const unsigned nrBaselines = nrStations * (nrStations + 1) / 2;
const unsigned nrTimes = 768, nrTimesWidth = 768; // 770
const unsigned nrChannels = 256;
const unsigned nrPolarizations = 2;
const unsigned iter = 1;
const unsigned nrThreads = 8;

class params {
public:
    int correlatorType;
    float* samples; 
    float* visibilities;
    unsigned long long ops, bytesLoaded, bytesStored;
};


vfloat32m8_t init_vfloat32m8(float val)
{
    float result[VECTOR_WIDTH_IN_FLOATS];
    for(int i=0; i<VECTOR_WIDTH_IN_FLOATS; i++) {
	result[i] = val;
    }

    return __riscv_vle32_v_f32m8(result, VECTOR_WIDTH_IN_FLOATS);
}

void print_vfloat32m8(vfloat32m8_t f, const char* name)
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
        
    // We do 16 ops in the loop. 1 FMA * 8 floats.
    for (unsigned long long x = 0; x < (1E9L/(16*16)); x++) {
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

    // print results to make sure the compiler doesn't optimize everything away
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
    case CORRELATOR_2X2:
	cout << "2x2";
	break;
    default:
	cout << "illegal correlator" << endl;
	exit(66);
    }
}

void printResult(float* visibilities)
{
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned baseline = 0; baseline < nrBaselines; baseline ++) {
	    for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
		for (unsigned pol1 = 0; pol1 < 2; pol1 ++) {
		    std::cout.precision(15);
		    std::cout << visibilities[VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 0)] << ", " 
			      << visibilities[VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 1)]
			      << std::endl;
		    
		}
	    }
	}
    }
}

void* runCorrelator(void* data)
{
    params* p = (params*) data;

    for(unsigned i=0; i<iter; i++) {
	switch (p->correlatorType) {
	case CORRELATOR_REFERENCE:
	    p->ops = referenceCorrelator(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_1X1:
	    p->ops = riscvCorrelator_1x1(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_2X2:
	    p->ops = riscvCorrelator_2x2(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
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
    printResult(visibilities);
#endif

    return 0;
}

double computeMaxFlops()
{
    pthread_t threads[nrThreads];
  
    auto start_time = chrono::high_resolution_clock::now();

    for(unsigned i=0; i<nrThreads; i++) {
	if (pthread_create(&threads[i], 0, calcMaxFlops, 0) != 0) {
	    std::cout << "could not create thread" << std::endl;
	    exit(1);
	}
    }

    for(unsigned i=0; i<nrThreads; i++) {
	if (pthread_join(threads[i], 0) != 0) {
	    std::cout << "could not join thread" << std::endl;
	    exit(1);
	}
    }
  
    auto end_time = chrono::high_resolution_clock::now();
    long long nanos = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count();
    double elapsed = (double) nanos / 1.0E9;
    unsigned long long gigaOps = 16L * nrThreads; 
    double maxGflops = (double) gigaOps / elapsed; 
    cout << "elapsed time in nanoseconds = " << nanos << ", in seconds = " << elapsed << endl;
    cout << "peak flops with " << nrThreads << " threads is: " << maxGflops << " gflops" << std::endl;
    return maxGflops;
}

void spawnCorrelatorThreads(int correlatorType, float* samples, const unsigned arraySize,
			    float* visibilities, const unsigned visArraySize, double maxFlops)
{
    cout << "Running correlator \"";
    printCorrelatorType(correlatorType);
    cout << "\" with " << nrThreads << " threads, for " << iter << " iterations" << endl;

    pthread_t threads[nrThreads];
    params p[nrThreads];
    for(unsigned t=0; t<nrThreads; t++) {
	p[t].ops = 0;
	p[t].bytesLoaded = 0;
	p[t].bytesStored = 0;
	p[t].samples = &samples[t*arraySize];
	p[t].visibilities = &visibilities[t*visArraySize];
	p[t].correlatorType = correlatorType;
    }

    auto start_time = chrono::high_resolution_clock::now();

    for(unsigned t=0; t<nrThreads; t++) {
	if (pthread_create(&threads[t], 0, runCorrelator, &p[t]) != 0) {
	    std::cout << "could not create thread" << std::endl;
	    exit(1);
	}
    }

    for(unsigned t=0; t<nrThreads; t++) {
	if (pthread_join(threads[t], 0) != 0) {
	    std::cout << "could not join thread" << std::endl;
	    exit(1);
	}
    }

    auto end_time = chrono::high_resolution_clock::now();

    unsigned long long ops=0, bytesLoaded=0, bytesStored=0;
    for(unsigned t=0; t<nrThreads; t++) {
	ops += p[t].ops;
	bytesLoaded += p[t].bytesLoaded;
	bytesStored += p[t].bytesStored;
    }

    long long nanos = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count();
    double elapsed = (double) nanos / 1.0E9;
    double flops = (ops / elapsed) / 1000000000.0;
    double efficiency = (flops / maxFlops) * 100.0;
    double gbsLoad = (double) (bytesLoaded / (1024.0 * 1024.0 * 1024.0)) / elapsed;
    double mbsStore = (double) (bytesStored / (1024.0 * 1024.0)) / elapsed;

    cout << "correlate took " << elapsed << " s, achieved " << flops << " Gflops, " << efficiency << " % efficiency" << endl;
    cout << "throughput: " << gbsLoad << " GB/s load, " << mbsStore << " MB/s store" << endl;
}

void initSamples(float* samples, const unsigned arraySize)
{
    for(unsigned t = 0; t<nrThreads; t++) {
	for (unsigned channel = 0; channel < nrChannels; channel ++) {
	    for (unsigned stat = 0; stat < nrStations; stat ++) {
		for (unsigned time = 0; time < nrTimes; time ++) {
		    samples[t*arraySize + SAMPLE_INDEX(stat, channel, time, 0, 0)] = 1.0f; // time % 8;
		    samples[t*arraySize + SAMPLE_INDEX(stat, channel, time, 0, 1)] = 2.0f; // stat;
		    
		    samples[t*arraySize + SAMPLE_INDEX(stat, channel, time, 1, 0)] = 1.0f; // channel;
		    samples[t*arraySize + SAMPLE_INDEX(stat, channel, time, 1, 1)] = 2.0f; // 0;
		}
	    }
	}
    }
}

void checkResult(float* samples, const unsigned arraySize,
		 float* visibilities, const unsigned visArraySize, double maxFlops)
{
    float* checkVisibilities = new (align_val_t{ALIGNMENT}) float[nrThreads*visArraySize];
    memset(checkVisibilities, 0, nrThreads*visArraySize*sizeof(float));
    
    spawnCorrelatorThreads(CORRELATOR_REFERENCE, samples, arraySize, checkVisibilities, visArraySize, maxFlops);

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned stat1 = 0; stat1 < nrStations; stat1 ++) {
	    for (unsigned stat0 = 0; stat0 <= stat1; stat0 ++) {
		for (unsigned time = 0; time < nrTimes; time ++) {
		    for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
			for (unsigned pol1 = 0; pol1 < 2; pol1 ++) { 
			    unsigned baseline = BASELINE(stat0, stat1);
			    unsigned vis_index_real = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 0);
			    unsigned vis_index_imag = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 1);
			    
			    if(visibilities[vis_index_real] != checkVisibilities[vis_index_real]) {
				std::cout << "ERROR: channel = " << channel << ", baseline = " << baseline <<
				    ", pol = " << pol0 << '/' << pol1 << ": " << visibilities[vis_index_real] <<
				    " != " << checkVisibilities[vis_index_real] << std::endl;
				return;
			    }
			    if(visibilities[vis_index_imag] != checkVisibilities[vis_index_imag]) {
				std::cout << "ERROR: channel = " << channel << ", baseline = " << baseline <<
				    ", pol = " << pol0 << '/' << pol1 << ": " << visibilities[vis_index_imag] <<
				    " != " << checkVisibilities[vis_index_imag] << std::endl;
			    }
			}
		    }
		}
	    }
	}
    }

    cout << "result validated OK" << endl;
    delete [] checkVisibilities;
}

int main()
{
    // set the vector length
    size_t vl = __riscv_vsetvl_e32m8(VECTOR_WIDTH_IN_FLOATS);
    cout << "vector length = " << vl << " 32 bit floats" << endl;
    cout << "Timer resolution is " << chrono::high_resolution_clock::period::den << " ticks per second" << endl;
    cout << "Configuration: " << nrStations << " stations, " << nrBaselines <<
	" baselines, " << nrChannels << " channels, " << nrTimes << " time samples, "
	 << nrPolarizations << " polarizations"  << endl;

    const unsigned arraySize = nrStations*nrChannels*nrTimesWidth*nrPolarizations*2;
    const unsigned visArraySize = nrBaselines*nrChannels*nrPolarizations*nrPolarizations*2;

    cout << "sample array size = " << ((arraySize * nrThreads * sizeof(float))/(1024*1024)) << " mbytes, "
	 << "vis array size = " << ((visArraySize * nrThreads * sizeof(float))/(1024*1024)) << " mbytes, " << endl;

    double maxFlops = computeMaxFlops();
   
    float* samples = new (align_val_t{ALIGNMENT}) float[nrThreads*arraySize];
    float* visibilities= new (align_val_t{ALIGNMENT}) float[nrThreads*visArraySize];
    memset(visibilities, 0, nrThreads*visArraySize*sizeof(float));

    initSamples(samples, arraySize);

    spawnCorrelatorThreads(CORRELATOR_1X1, samples, arraySize, visibilities, visArraySize, maxFlops);
    checkResult(samples, arraySize, visibilities, visArraySize, maxFlops);

    spawnCorrelatorThreads(CORRELATOR_2X2, samples, arraySize, visibilities, visArraySize, maxFlops);
    checkResult(samples, arraySize, visibilities, visArraySize, maxFlops);

    delete[] samples;
    delete[] visibilities;

    return 0;
}
