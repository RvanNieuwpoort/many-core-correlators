#include "common.h"

#include <iostream>
#include <chrono>
#include <cstring>

using namespace std;

static float* checkVisibilities = NULL;


unsigned calcNrCells(const unsigned w, const unsigned h, const unsigned nrStations) 
{
    unsigned nrCells = 0;
    for (unsigned stationY = h; stationY < nrStations - h + 1; stationY += h) {
	for (unsigned stationX = 0; stationX + w <= stationY; stationX += w) {
	    nrCells++;
	}
    }
    return nrCells;
}

void printResult(const float* visibilities, const unsigned nrChannels, const unsigned nrBaselines)
{
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned baseline = 0; baseline < nrBaselines; baseline ++) {
	    for (unsigned pol0 = 0; pol0 < NR_POLARIZATIONS; pol0 ++) {
		for (unsigned pol1 = 0; pol1 < NR_POLARIZATIONS; pol1 ++) {
		    std::cout.precision(15);
		    std::cout << visibilities[VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 0)] << ", " 
			      << visibilities[VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 1)]
			      << std::endl;
		    
		}
	    }
	}
    }
}

double computeMaxGflops(const unsigned nrThreads, void *(*runMaxFlopsTest) (void *), void* /* unused */)
{
    pthread_t threads[nrThreads];

    FlopTestParams p[nrThreads];
    for(unsigned t=0; t<nrThreads; t++) {
	p[t].gFlops = 0;
    }
    
    auto start_time = chrono::high_resolution_clock::now();

    for(unsigned i=0; i<nrThreads; i++) {
	if (pthread_create(&threads[i], 0, runMaxFlopsTest, &p[i]) != 0) {
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

    unsigned long long gFlops = 0L;
    for(unsigned t=0; t<nrThreads; t++) {
	gFlops += p[t].gFlops;
    }
    double maxGflops = (double) gFlops / elapsed; 
//    cout << "elapsed time in nanoseconds = " << nanos << ", in seconds = " << elapsed << endl;

    return maxGflops;
}

void initSamples(float* __restrict__ samples,
		 const unsigned nrThreads, const unsigned nrTimes, const unsigned nrTimesWidth,
			const unsigned nrStations, const unsigned nrChannels, const unsigned arraySize)
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

void checkResult(const float* __restrict__ samples, void *(*runCorrelator) (void *),
		 const float* __restrict__ visibilities, 
		 const unsigned nrThreads,
		 const unsigned nrTimes, const unsigned nrStations,
		 const unsigned nrChannels, const unsigned arraySize, const unsigned visArraySize)
{
	if(checkVisibilities == NULL) {
		checkVisibilities = new float[nrThreads*visArraySize];
		memset(checkVisibilities, 0, nrThreads*visArraySize*sizeof(float));
	
		spawnCorrelatorThreads(CORRELATOR_REFERENCE, runCorrelator, samples,
				       arraySize, checkVisibilities, visArraySize,
				       nrTimes, nrStations, nrChannels,
				       nrThreads, 1, 0.0, false, false);
	}
    
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned stat1 = 0; stat1 < nrStations; stat1 ++) {
	    for (unsigned stat0 = 0; stat0 <= stat1; stat0 ++) {
		for (unsigned time = 0; time < nrTimes; time ++) {
		    for (unsigned pol0 = 0; pol0 < NR_POLARIZATIONS; pol0 ++) {
			for (unsigned pol1 = 0; pol1 < NR_POLARIZATIONS; pol1 ++) { 
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
}

void endCommon()
{
	if(checkVisibilities != NULL) {
		delete [] checkVisibilities;
	}
}

void spawnCorrelatorThreads(int correlatorType, void *(*runCorrelator) (void *),
			    const float* __restrict__ samples, const unsigned arraySize,
			    float* __restrict__ visibilities, const unsigned visArraySize,
			    const unsigned nrTimes, const unsigned nrStations, const unsigned nrChannels,
			    const unsigned nrThreads, const unsigned iter, double maxFlops,
			    const bool verbose, const bool validateResults)
{
    if(verbose) {
	cout << "Running correlator \"";
	printCorrelatorType(correlatorType);
	cout << "\" with " << nrThreads << " threads, for " << iter << " iterations" << endl;
    }

    pthread_t threads[nrThreads];
    CorrelatorParams p[nrThreads];
    for(unsigned t=0; t<nrThreads; t++) {
	p[t].ops = 0;
	p[t].bytesLoaded = 0;
	p[t].bytesStored = 0;
	p[t].samples = &samples[t*arraySize];
	p[t].visibilities = &visibilities[t*visArraySize];
	p[t].correlatorType = correlatorType;
	p[t].verbose = verbose;
	p[t].iter = iter;
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

    if(verbose) {
	cout << "correlate ";
	printCorrelatorType(correlatorType);
	cout << " took " << elapsed << " s, achieved " << flops << " Gflops, " << efficiency << " % efficiency" << endl;
	cout << "throughput ";
	printCorrelatorType(correlatorType);
	cout << " is "<< gbsLoad << " GB/s load, " << mbsStore << " MB/s store" << endl;
    }

    if(validateResults) {
	    checkResult(samples, runCorrelator, visibilities, nrThreads, nrTimes, nrStations,
			nrChannels, arraySize, visArraySize);
    }
}
