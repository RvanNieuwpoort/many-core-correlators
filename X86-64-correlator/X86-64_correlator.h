#ifndef _X86_64_CORRELATOR_H
#define _X86_64_CORRELATOR_H

#include "../common/common.h"

#if defined(__AVX512F__)
#define VECTOR_WIDTH_IN_FLOATS 16 // 512 bit
#elif defined(__AVX2__)
#define VECTOR_WIDTH_IN_FLOATS 8 // 256 bit
#elif defined(__AVX__)
#define VECTOR_WIDTH_IN_FLOATS 8 // 156 bit
#elif defined(__SSE3__)
#define VECTOR_WIDTH_IN_FLOATS 4 // 128 bit
#else
#error unsupported
#endif

#define ALIGNMENT 256 // in bytes

#define MAX_CELLS 2048

#define CORRELATOR_1X1_SSE3      5
#define CORRELATOR_1X1_TIME_SSE3 6
#define CORRELATOR_2X2_SSE3      7
#define CORRELATOR_2X2_TIME_SSE3 8
#define CORRELATOR_3X2_TIME_SSE3 9
#define CORRELATOR_1X1_TIME_AVX512 10


extern unsigned long long cpuCorrelator_1x1_sse3(const float* __restrict__ samples, float* __restrict__ visibilities,
						 unsigned nrTimes, const unsigned nrTimesWidth,
						 const unsigned nrStations, const unsigned nrChannels,
						 unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long cpuCorrelator_1x1_time_sse3(const float* __restrict__ samples, float* __restrict__ visibilities,
						      const unsigned nrTimes, const unsigned nrTimesWidth,
						      const unsigned nrStations, const unsigned nrChannels,
						      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long cpuCorrelator_2x2_sse3(const float* __restrict__ samples, float* __restrict__ visibilities,
						 const unsigned nrTimes, const unsigned nrTimesWidth,
						 const unsigned nrStations, const unsigned nrChannels,
						 unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long cpuCorrelator_2x2_time_sse3(const float* __restrict__ samples, float* __restrict__ visibilities, 
						      const unsigned nrTimes, const unsigned nrTimesWidth, 
						      const unsigned nrStations, const unsigned nrChannels,
						      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long cpuCorrelator_3x2_time_sse3(const float* __restrict__ samples, float* __restrict__ visibilities, 
						      const unsigned nrTimes, const unsigned nrTimesWidth, 
						      const unsigned nrStations, const unsigned nrChannels,
						      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long cpuCorrelator_1x1_time_avx512(const float* __restrict__ samples, float* __restrict__ visibilities,
						      const unsigned nrTimes, const unsigned nrTimesWidth,
						      const unsigned nrStations, const unsigned nrChannels,
						      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

#endif // _X86_64_CORRELATOR_H
