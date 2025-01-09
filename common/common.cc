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
		    std::cout << "BL " << baseline << " ";
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

    return (double) gFlops / elapsed; 
}

void initSamples(float* __restrict__ samples,
		 const unsigned nrThreads, const unsigned nrTimes, const unsigned nrTimesWidth,
		 const unsigned nrStations, const unsigned nrChannels, const size_t arraySize)
{
    for(unsigned t = 0; t < nrThreads; t++) {
	for (unsigned stat = 0; stat < nrStations; stat++) {
	    for (unsigned channel = 0; channel < nrChannels; channel++) {
		for (unsigned time = 0; time < nrTimes; time++) {
		    
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
		 const unsigned nrChannels, const size_t arraySize, const size_t visArraySize) {

    if(checkVisibilities == NULL) {
	checkVisibilities = new float[nrThreads*visArraySize];
	memset(checkVisibilities, 0, nrThreads*visArraySize*sizeof(float));
	
	spawnCorrelatorThreads(CORRELATOR_REFERENCE, runCorrelator, samples,
			       arraySize, checkVisibilities, visArraySize,
			       nrTimes, nrStations, nrChannels,
			       nrThreads, 0.0, false, false);
    }

    int result = memcmp(visibilities, checkVisibilities, nrThreads*visArraySize);

    /*
    */

    if(result == 0) {
	cout << "result validated OK" << endl;
    } else {
	cout << "ERROR validating result" << endl;
	detailedComparison(visibilities, nrThreads, nrTimes, nrStations, nrChannels, visArraySize);
    }
}

void detailedComparison(const float* __restrict__ visibilities, 
			const unsigned nrThreads, const unsigned nrTimes,
			const unsigned nrStations, const unsigned nrChannels, const size_t visArraySize) {

    for(unsigned t = 0; t < nrThreads; t++) {
	for (unsigned channel = 0; channel < nrChannels; channel ++) {
	    for (unsigned baseline = 0; baseline < NR_BASELINES(nrStations); baseline++) {
		for (unsigned pol0 = 0; pol0 < NR_POLARIZATIONS; pol0 ++) {
		    for (unsigned pol1 = 0; pol1 < NR_POLARIZATIONS; pol1 ++) { 
			size_t vis_index_real = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 0);
			size_t vis_index_imag = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 1);
			
			if(visibilities[t*visArraySize + vis_index_real] != checkVisibilities[vis_index_real]) {
			    std::cout << "ERROR: channel = " << channel << ", baseline = " << baseline <<
				", pol = " << pol0 << '/' << pol1 << ": " << visibilities[vis_index_real] <<
				" != " << checkVisibilities[vis_index_real] << std::endl;
			    return;
			} 
			if(visibilities[t*visArraySize + vis_index_imag] != checkVisibilities[vis_index_imag]) {
			    std::cout << "ERROR: channel = " << channel << ", baseline = " << baseline <<
				", pol = " << pol0 << '/' << pol1 << ": " << visibilities[vis_index_imag] <<
				" != " << checkVisibilities[vis_index_imag] << std::endl;
			    return;
			} 
		    }
		}
	    }
	}
    }
}

void endCommon()
{
	if(checkVisibilities != NULL) {
		delete [] checkVisibilities;
	}
}

void spawnCorrelatorThreads(int correlatorType, void *(*runCorrelator) (void *),
			    const float* __restrict__ samples, const size_t arraySize,
			    float* __restrict__ visibilities, const size_t visArraySize,
			    const unsigned nrTimes, const unsigned nrStations, const unsigned nrChannels,
			    const unsigned nrThreads, double maxFlops,
			    const bool verbose, const bool validateResults)
{
    if(verbose) {
	cout << "Running correlator \"";
	printCorrelatorType(correlatorType);
	cout << "\" with " << nrThreads << " threads" << endl;
    }

    memset(visibilities, 0, nrThreads*visArraySize*sizeof(float));

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


unsigned long long computeMissedBaselines(const float* __restrict__ samples, float* __restrict__ visibilities, 
					  const bool* __restrict__ missedBaselines,
					  const unsigned nrTimes, const unsigned nrTimesWidth,
					  const unsigned nrStations, const unsigned nrChannels,
					  unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned computedBaselines = 0;
    for(unsigned stationY = 0; stationY < nrStations; stationY++) {
	for(unsigned stationX = 0; stationX <= stationY; stationX++) {
	    const unsigned baseline = BASELINE(stationX, stationY);
	    if(missedBaselines[baseline]) computedBaselines++;
	}
    }
    
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for(unsigned stationY = 0; stationY < nrStations; stationY++) {
	    for(unsigned stationX = 0; stationX <= stationY; stationX++) {

		const unsigned baseline = BASELINE(stationX, stationY);
		if(!missedBaselines[baseline]) continue;

//		cout << "computing missed bl " << baseline << endl;
		
		float xxr = 0, xxi = 0, xyr = 0, xyi = 0, yxr = 0, yxi = 0, yyr = 0, yyi = 0;
		size_t index1 = SAMPLE_INDEX(stationX, channel, 0, 0, 0);
		size_t index2 = SAMPLE_INDEX(stationY, channel, 0, 0, 0);

		for (unsigned time = 0; time < nrTimes; time ++) {
		    float sample1xr = samples[index1+0];
		    float sample1xi = samples[index1+1];
		    float sample1yr = samples[index1+2];
		    float sample1yi = samples[index1+3];
		    float sample2xr = samples[index2+0];
		    float sample2xi = samples[index2+1];
		    float sample2yr = samples[index2+2];
		    float sample2yi = samples[index2+3];
		    
		    xxr += sample1xr * sample2xr + sample1xi * sample2xi;
		    xxi += sample1xi * sample2xr - sample1xr * sample2xi;
		    
		    xyr += sample1xr * sample2yr + sample1xi * sample2yi;
		    xyi += sample1xi * sample2yr - sample1xr * sample2yi;
		    
		    yxr += sample1yr * sample2xr + sample1yi * sample2xi;
		    yxi += sample1yi * sample2xr - sample1yr * sample2xi;
		    
		    yyr += sample1yr * sample2yr + sample1yi * sample2yi;
		    yyi += sample1yi * sample2yr - sample1yr * sample2yi;
		    
		    index1 += 4;
		    index2 += 4;
		}
		const size_t vis_index = VISIBILITIES_INDEX(baseline, channel, 0, 0, 0);
		visibilities[vis_index+0] = xxr;
		visibilities[vis_index+1] = xxi;
		visibilities[vis_index+2] = xyr;
		visibilities[vis_index+3] = xyi;
		visibilities[vis_index+4] = yxr;
		visibilities[vis_index+5] = yxi;
		visibilities[vis_index+6] = yyr;
		visibilities[vis_index+7] = yyi;
	    }
	}
    }
    
    *bytesLoaded = nrChannels * computedBaselines * nrTimes * 8L * sizeof(float); // samples
    *bytesStored = nrChannels * computedBaselines * 8L * sizeof(float); // vis
    
    return nrChannels * computedBaselines * nrTimes * 16L * 2L;
}
