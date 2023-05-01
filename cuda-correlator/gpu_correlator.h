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
// TODO implement DO_STORES and use TEXTURE_CACHE with other sizes!
// TODO alternate + spilling
#ifndef _GPU_CORRELATOR_H
#define _GPU_CORRELATOR_H

// #include <complex> does not work in -deviceemu mode
// so we include our own simple version
#include "gpu_complex.h"

// global options
#define CALCULATE_MISSING                       0
#define DO_IO                                   1
#define USE_STREAMING                           0

// memory options
#define USE_ALTERNATE_MEMORY_LAYOUT             0
#define DO_MEMORY_CONVERSION                    0
#define CONVERT_ON_GPU                          0
#define USE_MULTITHREADED_MEMORY_CONVERSION     1
#define MULTITHREADED_MEMORY_CONVERSION_THREADS 1
#define TIME_CONVERSION                         1

#define USE_WRITE_COMBINING                     1

// kernel options
#define USE_PREFETCHING                         0
#define USE_TEXTURE_CACHE                       1
#define DO_LOADS                                1
#define DO_STORES                               1
#define SYNC_THREADS                            1

// used for calculating theoretical maximum FLOPS
#define NR_PROCESSORS_PER_MULTIPROCESSOR 8 // according to manual, 1.0 - 1.3 capability

//#define MAX_DEVICE_MEM_BANDWIDTH 141.7 // GTX280
#define MAX_DEVICE_MEM_BANDWIDTH 102 // Tesla 1060C

#define MAX_CELLS 4096     // used for fixed-sized cellToStat tables
#define MAX_THREADS 512    // used for spilling version
#define MAX_STATIONS 64   // used for spilling version
#define MAX_CHANNELS 256  // used for spilling version

#define BASELINE(STATION_1, STATION_2)			\
    ((STATION_2) * ((STATION_2) + 1) / 2 + (STATION_1))

// Alternate layout is:
// channel  time station polarization complex
#define ALTERNATE_SAMPLE_INDEX(STATION, CHANNEL, TIME, POLARIZATION)		\
    ((((CHANNEL) * nrTimesWidth + (TIME)) * nrStations + (STATION)) * 2 + (POLARIZATION))

// Normal layout is:
// channel station time polarization complex
#define NORMAL_SAMPLE_INDEX(STATION, CHANNEL, TIME, POLARIZATION)		\
    ((((CHANNEL) * nrStations + (STATION)) * nrTimesWidth + (TIME)) * 2 + (POLARIZATION))

#if USE_ALTERNATE_MEMORY_LAYOUT
#define SAMPLE_INDEX ALTERNATE_SAMPLE_INDEX
#else
#define SAMPLE_INDEX NORMAL_SAMPLE_INDEX
#endif

#define VISIBILITIES_INDEX(BASELINE, CHANNEL, POLARIZATION_1, POLARIZATION_2) \
    ((((BASELINE) * nrChannels + (CHANNEL)) * 2 + (POLARIZATION_1)) * 2 + (POLARIZATION_2))

#define LOOP_COUNT(NR_CELLS, NR_THREADS) (((NR_CELLS) + (NR_THREADS) - 1) / (NR_THREADS) * (NR_THREADS))

#if USE_PREFETCHING
#if USE_ALTERNATE_MEMORY_LAYOUT
#define PREFETCH()						     \
	unsigned prefetchIndex = SAMPLE_INDEX(0, channel, 0, 0) / 2; \
	for (unsigned i = myThread; i < nrStations; i += nrThreads) { \
	    float4 sample = tex1Dfetch(tex0, prefetchIndex + i); \
	}
#else // !USE_ALTERNATE_MEMORY_LAYOUT
#define PREFETCH()						     \
	unsigned prefetchIndex = SAMPLE_INDEX(0, channel, 0, 0) / 2; \
	for (unsigned i = myThread; i < nrStations; i += nrThreads) { \
	    float4 sample = tex1Dfetch(tex0, prefetchIndex + i * nrTimesWidth); \
	}
#endif // USE_ALTERNATE_MEMORY_LAYOUT
#else // !USE_PREFETCHING
#define PREFETCH()
#endif // USE_PREFETCHING

#define CORRELATOR_1x1           1
#define CORRELATOR_1x1_SHAREDMEM 2
#define CORRELATOR_1x2           3
#define CORRELATOR_1x3           4
#define CORRELATOR_1x4           5
#define CORRELATOR_1x6           6
#define CORRELATOR_2x2           7
#define CORRELATOR_2x2_NEW       8
#define CORRELATOR_3x2           9
#define CORRELATOR_3x2_SPILL     10
#define CORRELATOR_3x3           11
#define CORRELATOR_3x3_SHAREDMEM 12
#define CORRELATOR_4x3           13

extern unsigned nrStations;
extern unsigned nrTimes, nrTimesWidth;
extern unsigned nrChannels;
extern unsigned nrPolarizations;

extern unsigned char missedStatXHost[4*MAX_STATIONS], missedStatYHost[4*MAX_STATIONS];

extern void correlateOnHost(const complex<float> *samples, complex<float> *visibilities);
extern void correlateMissingOnHost(const unsigned nrMissed, const complex<float> *samples, complex<float> *visibilities);
unsigned long long calcNrOps();
unsigned fillCellToStatTable();

#endif // _GPU_CORRELATOR_H
