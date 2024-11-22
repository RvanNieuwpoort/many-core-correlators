#include "cpu_correlator.h"
//#include "timer.h"

#include <complex>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <chrono>

#define CORRELATOR_VERSION CORRELATOR_1X1

using namespace std;

const unsigned nrStations	 = 32;
const unsigned nrTimes	 = 768, nrTimesWidth = 768; // 770
const unsigned nrChannels	 = 128;
const unsigned nrPolarizations = 2;
const unsigned iter = 2;
const unsigned nrThreads = 8;

class params {
public:
    int version;
    float* samples; 
    float* visibilities;
    unsigned long long ops, bytesLoaded, bytesStored;
};

void* calcMaxFlops(void* data)
{
  // Do some (preferably) vectorized computation here to compute the maximum peak performance.
  return 0;
}

void printCorrelatorType()
{
      switch (CORRELATOR_VERSION) {
	case CORRELATOR_REFERENCE:
	  cout << "reference";
	  break;
	case CORRELATOR_1X1:
	  cout << "1x1";
	    break;
	default:
	    cout << "illegal correlator" << endl;
	    exit(66);
	}
}

void* runCorrelator(void* data)
{
    params* p = (params*) data;

    for(unsigned i=0; i<iter; i++) {
      switch (p->version) {
	case CORRELATOR_REFERENCE:
	    p->ops = referenceCorrelator(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_1X1:
	    p->ops = cpuCorrelator_1x1(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	default:
	    cout << "illegal correlator" << endl;
	    exit(66);
	}
    }

    p->ops *=iter;
    p->bytesLoaded *= iter;
    p->bytesStored *=iter;

#if 0
    unsigned nrBaselines = nrStations * (nrStations + 1) / 2;
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
#endif

    return 0;
}

int main()
{
  pthread_t threads[nrThreads];
  
  //  timer calcTimer("calc");
  //  calcTimer.start();
  /*
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
  //  calcTimer.stop();

  double time = 1.0; //calcTimer.getTimeInSeconds();

  unsigned long long tmpOps = 16L * 4L * nrThreads;

  //  double maxFlops = (double)tmpOps / time; // gigaflops
  
  cout << "total maxFlops with " << nrThreads << " threads is: " << maxFlops << std::endl;
  */
  double maxFlops = 100;
  const unsigned nrBaselines = nrStations * (nrStations + 1) / 2;

  const unsigned arraySize = nrStations*nrChannels*nrTimesWidth*nrPolarizations*2;
  const unsigned visArraySize = nrBaselines*nrChannels*nrPolarizations*nrPolarizations*2;

  float* samples = new float[nrThreads*nrStations*nrChannels*nrTimesWidth*nrPolarizations*2];
  float* visibilities= new float[nrThreads*nrBaselines*nrChannels*nrPolarizations*nrPolarizations*2];

#if 0
    for (unsigned time = 0; time < nrTimes; time ++) {
	samples[0][7][time][0] = 1;
	samples[0][7][time][1] = complex<float>(3, 4);
	samples[5][7][time][0] = 1;
    }
#else
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
#endif

    cout << "Running correlator \"";
    printCorrelatorType();
    cout << "\" with " << nrThreads << " threads, for " << iter << " iterations" << endl;
    cout << "Configuration: " << nrStations << " stations, " << nrChannels << " channels, " << nrTimes << " time samples, " << nrPolarizations << " polarizations"  << endl;

    cout << "Timer resolution is " << chrono::high_resolution_clock::period::den << " ticks per second" << endl;

    params p[nrThreads];

    auto start_time = chrono::high_resolution_clock::now();
    
    for(unsigned t=0; t<nrThreads; t++) {
	p[t].ops = 0;
	p[t].bytesLoaded = 0;
	p[t].bytesStored = 0;
	p[t].samples = &samples[t*arraySize];
	p[t].visibilities = &visibilities[t*visArraySize];
	p[t].version = CORRELATOR_VERSION;

	if (pthread_create(&threads[t], 0, runCorrelator, &p[t]) != 0) {
	    std::cout << "could not create thread" << std::endl;
	    exit(1);
	}
    }

    unsigned long long ops=0, bytesLoaded=0, bytesStored=0;
    for(unsigned t=0; t<nrThreads; t++) {
	if (pthread_join(threads[t], 0) != 0) {
	    std::cout << "could not join thread" << std::endl;
	    exit(1);
	}
	ops += p[t].ops;
	bytesLoaded += p[t].bytesLoaded;
	bytesStored += p[t].bytesStored;
    }

    auto end_time = chrono::high_resolution_clock::now();
    long long nanos = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time).count();
    
    double elapsed = (double) nanos / 1.0E9;

    //    cout << "elapsed time in nanoseconds = " << nanos << ", in seconds = " << elapsed << endl;

    double flops = (ops / elapsed) / 1000000000.0;
    double efficiency = (flops / maxFlops) * 100.0;
    cout << "correlate took " << elapsed << " s, max Gflops = " << maxFlops << ", achieved " << flops << " Gflops, " << efficiency << " % efficiency" << endl;
    double gbsLoad = (double) (bytesLoaded / (1024.0 * 1024.0 * 1024.0)) / elapsed;
    double mbsStore = (double) (bytesStored / (1024.0 * 1024.0)) / elapsed;

    cout << "throughput: " << gbsLoad << " GB/s load, " << mbsStore << " MB/s store" << endl;

    std::cout.precision(15);
    std::cout << visibilities[VISIBILITIES_INDEX(BASELINE(0, 0), 7, 0, 0, 0)] << ", " 
	      << visibilities[VISIBILITIES_INDEX(BASELINE(0, 0), 7, 0, 0, 1)]
	      << std::endl;

    std::cout << visibilities[VISIBILITIES_INDEX(BASELINE(0, 4), 7, 0, 0, 0)]  << ", " 
	      << visibilities[VISIBILITIES_INDEX(BASELINE(0, 4), 7, 0, 0, 1)] 
	      << std::endl;

    std::cout << visibilities[VISIBILITIES_INDEX(BASELINE(0, 4), 7, 0, 1, 0)]  << ", " 
	      << visibilities[VISIBILITIES_INDEX(BASELINE(0, 4), 7, 0, 1, 1)] 
	      << std::endl;

#if 0
    complex<float> *checkVis = new complex<float>[visibilitiesSize / sizeof(complex<float>)];
    correlateOnHost(samples[0], checkVis);
    std::cout << checkVis[BASELINE(0, 0)][7][0][0] << std::endl;
    std::cout << checkVis[BASELINE(0, 4)][7][0][0] << std::endl;
    std::cout << checkVis[BASELINE(0, 4)][7][0][1] << std::endl;

    for (unsigned channel = 0; channel < nrChannels; channel ++)
	for (unsigned baseline = 0; baseline < nrBaselines; baseline ++)
	    for (unsigned pol0 = 0; pol0 < 2; pol0 ++)
		for (unsigned pol1 = 0; pol1 < 2; pol1 ++) {
		    size_t index = VISIBILITIES_INDEX(baseline, channel, pol0, pol1);

		    if (visibilities[index] != checkVis[index])
			std::cout << "channel = " << channel << ", baseline = " << baseline << ", pol = " << pol0 << '/' << pol1 << ": " << visibilities[index] << " != " << checkVis[index] << std::endl;
		}

    delete [] checkVis;
#endif

    delete[] samples;
    delete[] visibilities;

    return 0;
}
