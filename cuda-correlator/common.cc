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
#include <iostream>
#include <cassert>
#include <pmmintrin.h> // sse3

#include "gpu_correlator.h"
#include "common.h"

// #include <complex> does not work in -deviceemu mode
// so we include our own simple version
#include "gpu_complex.h"

#if ! USE_MULTITHREADED_MEMORY_CONVERSION
// assume 2 polarizations
void convertBufferLayout(const complex<float>* src, complex<float>* dst)
{
#if TIME_CONVERSION
    timer conversionTimer("conversion");
    conversionTimer.start();
#endif

#if 0 // debug version
    for(unsigned channel=0; channel<nrChannels; channel++) {
	for(unsigned time=0; time<nrTimes; time++) {
	    for(unsigned station=0; station<nrStations; station++) {
		dst[ALTERNATE_SAMPLE_INDEX(station, channel, time, 0)] = src[NORMAL_SAMPLE_INDEX(station, channel, time, 0)];
		dst[ALTERNATE_SAMPLE_INDEX(station, channel, time, 1)] = src[NORMAL_SAMPLE_INDEX(station, channel, time, 1)];
	    }
	}
    }
#else
    unsigned destIndex = 0;
    for(unsigned channel=0; channel<nrChannels; channel++) {
	for(unsigned time=0; time<nrTimes; time++) {
	    unsigned srcIndex = NORMAL_SAMPLE_INDEX(0, channel, time, 0);
	    for(unsigned station=0; station<nrStations; station++) {
		__m128 sample = _mm_load_ps((const float*)(src + srcIndex));
		_mm_store_ps((float*)(dst + destIndex), sample);

		destIndex += 2;
		srcIndex += nrTimesWidth * 2;
	    }
	}
	destIndex += (nrTimesWidth - nrTimes) * nrStations * 2;
    }
#endif

#if TIME_CONVERSION
    conversionTimer.stop();
    std::cout << "normal to alternate conversion took: " << conversionTimer << std::endl;
#endif
}

#else // USE_MULTITHREADED_MEMORY_CONVERSION

unsigned nrThreads = MULTITHREADED_MEMORY_CONVERSION_THREADS;

class params {
public:
    const complex<float>* src;
    complex<float>* dst;
    unsigned myThreadID;
};

void* convertBufferLayoutThread(void* data)
{
    params* p = (params*) data;

    for(unsigned channel=p->myThreadID; channel<nrChannels; channel+=nrThreads) {
	unsigned destIndex = ALTERNATE_SAMPLE_INDEX(0, channel, 0, 0);
	for(unsigned time=0; time<nrTimes; time++) {
	    unsigned srcIndex = NORMAL_SAMPLE_INDEX(0, channel, time, 0);
	    for(unsigned station=0; station<nrStations; station++) {
		__m128 sample = _mm_load_ps((const float*)(p->src + srcIndex));
		_mm_store_ps((float*)(p->dst + destIndex), sample);

		destIndex += 2;
		srcIndex += nrTimesWidth * 2;
	    }
	}
    }

    return 0;
}

void convertBufferLayout(const complex<float>* src, complex<float>* dst)
{
    pthread_t threads[nrThreads];
    params p[nrThreads];

#if TIME_CONVERSION
    timer conversionTimer("conversion");
    conversionTimer.start();
#endif

    for(unsigned i=0; i<nrThreads; i++) {
	p[i].src = src;
	p[i].dst = dst;
	p[i].myThreadID = i;

	if (pthread_create(&threads[i], 0, convertBufferLayoutThread, &p[i]) != 0) {
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

#if TIME_CONVERSION
    conversionTimer.stop();
    std::cout << "normal to alternate conversion took: " << conversionTimer << std::endl;
#endif
}
#endif //USE_MULTITHREADED_MEMORY_CONVERSION

void printVisibility(unsigned channel, unsigned cell, unsigned sXStart, unsigned sYStart, 
		     unsigned statX, unsigned statY, unsigned baseline, unsigned pol0, unsigned pol1,
		     complex<float> *hostVisibilities, complex<float> *checkVis, size_t index)
{
    std::cout.width(3);
    std::cout << "chan ";
    std::cout.width(3);
    std::cout << channel;
    std::cout << ", cellNr ";
    std::cout.width(5);
    std::cout << cell;
    std::cout << ", XStart ";
    std::cout.width(2);
    std::cout << sXStart;
    std::cout << ", YStart ";
    std::cout.width(2);
    std::cout << sYStart;
    std::cout << ", sX ";
    std::cout.width(2);
    std::cout << statX;
    std::cout << ", sY ";
    std::cout.width(2);
    std::cout << statY;
    std::cout << ", bl ";
    std::cout.width(5);
    std::cout << baseline;
    std::cout << ", pol " << pol0 << '/' << pol1 << ": (";
    std::cout.precision(10);
    std::cout.width(10);
    std::cout << hostVisibilities[index].real;
    std::cout << ", ";
    std::cout.precision(10);
    std::cout.width(10);
    std::cout << hostVisibilities[index].imag;
    std::cout << ") should be (";
    std::cout.precision(10);
    std::cout.width(10);
    std::cout << checkVis[index].real;
    std::cout << ", ";
    std::cout.precision(10);
    std::cout.width(10);
    std::cout << checkVis[index].imag;
    std::cout << ")";
    if (hostVisibilities[index] != checkVis[index]) {
	std::cout << " !!";
    } else {
	std::cout << " OK";
    }
    std::cout << std::endl;
}

#if CALCULATE_MISSING
void checkResult(size_t visibilitiesSize, complex<float> *hostSamples,
		 complex<float> *hostVisibilities, unsigned nrBaselines, unsigned nrCells)
{
    std::cout << "checking result..." << std::endl;

    complex<float> *checkVis = new complex<float>[visibilitiesSize / sizeof(complex<float>)];
    correlateOnHost(hostSamples, checkVis);

    bool errorFound = false;
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned stat1 = 0; stat1 < nrStations; stat1 ++) {
	    for (unsigned stat0 = 0; stat0 <= stat1; stat0 ++) {
		for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
		    for (unsigned pol1 = 0; pol1 < 2; pol1 ++) { 
			unsigned baseline = BASELINE(stat0, stat1);
			size_t index = VISIBILITIES_INDEX(baseline, channel, pol0, pol1);
			if (hostVisibilities[index] != checkVis[index]) {
			    errorFound = true;
			}
			if(printResult) {
			    printVisibility(channel, 0, 0, 0, stat0, stat1, baseline, pol0, pol1, hostVisibilities, checkVis, index);
			}
		    }
		}
	    }
	}
    }

    if(errorFound) {
	std::cout << "ERRORS in result." << std::endl;
    } else {
	std::cout << "result OK." << std::endl;
    }

    delete [] checkVis;
}
#else
void checkResult(size_t visibilitiesSize, complex<float> *hostSamples,
		 complex<float> *hostVisibilities, unsigned nrBaselines, unsigned nrCells)
{
    std::cout << "checking result..." << std::endl;

    complex<float> *checkVis = new complex<float>[visibilitiesSize / sizeof(complex<float>)];
    correlateOnHost(hostSamples, checkVis);

    bool errorFound = false;
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for(unsigned cell = 0; cell < nrCells; cell++) {
	    int sXStart = cellToStatXHost[cell];
	    int sYStart = cellToStatYHost[cell];

	    for(unsigned h = 0; h<cellHeight; h++) {
		for(unsigned w = 0; w<cellWidth; w++) {
		    int statX = sXStart + w;
		    int statY = sYStart + h;
		    unsigned baseline = BASELINE(statX, statY);

		    for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
			for (unsigned pol1 = 0; pol1 < 2; pol1 ++) {
			    size_t index = VISIBILITIES_INDEX(baseline, channel, pol0, pol1);
			    if (hostVisibilities[index] != checkVis[index]) {
				errorFound = true;
			    }
			    if(printResult) {
				printVisibility(channel, cell, sXStart, sYStart, statX, statY, baseline, pol0, pol1, hostVisibilities, checkVis, index);

			    }
			}
		    }
		}
	    }
	}
    }

    if(errorFound) {
	std::cout << "ERRORS in result." << std::endl;
    } else {
	std::cout << "result OK." << std::endl;
    }

    delete [] checkVis;
}
#endif // CALCULATE_MISSING


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
