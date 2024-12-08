#include "../common/common.h"
#include "X86-64_correlator.h"

#include <complex>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <chrono>
#include <emmintrin.h>

using namespace std;

#define VECTOR_WIDTH_IN_FLOATS 4

static const bool verbose = true;
static const bool validateResults = true;
static const unsigned nrStations = 64;
static const unsigned nrBaselines = NR_BASELINES(nrStations);
static const unsigned nrTimes = 768, nrTimesWidth = 768; // 770
static const unsigned nrChannels = 128;
static const unsigned iter = 1;
static const unsigned nrThreads = 32;


void print_m128(__m128 a, const char* name)
{
    float result[VECTOR_WIDTH_IN_FLOATS];
    for(int i=0; i<VECTOR_WIDTH_IN_FLOATS; i++) {
	result[i] = -1.0f;
    }

    _mm_store_ps(result, a);

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
  __m128 a = _mm_set_ps1(0.2E-9f);
  __m128 b = _mm_set_ps1(0.3E-9f);
  __m128 c = _mm_set_ps1(0.4E-9f);
  __m128 d = _mm_set_ps1(0.5E-9f);
  __m128 e = _mm_set_ps1(0.6E-9f);
  __m128 f = _mm_set_ps1(0.7E-9f);
  __m128 g = _mm_set_ps1(0.8E-9f);
  __m128 h = _mm_set_ps1(0.9E-9f);
  __m128 i = _mm_set_ps1(0.10E-9f);
  __m128 j = _mm_set_ps1(0.11E-9f);
  __m128 k = _mm_set_ps1(0.12E-9f);
  __m128 l = _mm_set_ps1(0.13E-9f);
  __m128 m = _mm_set_ps1(0.14E-9f);
  __m128 n = _mm_set_ps1(0.15E-9f);
  __m128 o = _mm_set_ps1(0.16E-9f);
  __m128 p = _mm_set_ps1(0.17E-9f);

  volatile __m128 result;

  for (unsigned long long x = 0; x < 1000000000L; x++) {
    a += a * b;
    b += b * c;
    c += c * d;
    d += d * e;
    e += e * f;
    f += f * g;
    g += g * h;
    h += h * i;
    i += i * j;
    j += j * k;
    k += k * l;
    l += l * m;
    m += m * n;
    n += n * p;
    o += o * p;
    p += p * a;
  }

  result = a + b + c + d + e + f + g + h + i + j + k + l + m + n + o + p;
  *(__m128*)data = result;

    // Print results just to make sure the compiler doesn't optimize everything away.
  print_m128(a, "a");
  print_m128(b, "b");
  print_m128(c, "c");
  print_m128(d, "d");
  print_m128(e, "e");
  print_m128(f, "f");
  print_m128(g, "g");
  print_m128(h, "h");
  print_m128(i, "i");
  print_m128(j, "j");
  print_m128(k, "k");
  print_m128(l, "l");
  print_m128(m, "m");
  print_m128(n, "n");
  print_m128(o, "o");
  print_m128(p, "p");
  
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
    case CORRELATOR_1X1_SSE3:
        cout << "1x1-SSE3";
	break;
    case CORRELATOR_1X1_TIME_SSE3:
	cout << "1x1-time-SSE3";
	break;
    case CORRELATOR_2X2_SSE3:
	cout << "2x2-SSE3";
	break;
    case CORRELATOR_2X2_TIME_SSE3:
	cout << "2x2-time-SSE3";
	break;
    case CORRELATOR_3X2_TIME_SSE3:
	cout << "3x3-time-SSE3";
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
	case CORRELATOR_1X1_SSE3:
	    p->ops = cpuCorrelator_1x1_sse3(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_1X1_TIME_SSE3:
	    p->ops = cpuCorrelator_1x1_time_sse3(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_2X2_SSE3:
	    p->ops = cpuCorrelator_2x2_sse3(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_2X2_TIME_SSE3:
	    p->ops = cpuCorrelator_2x2_time_sse3(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	    break;
	case CORRELATOR_3X2_TIME_SSE3:
	    p->ops = cpuCorrelator_3x2_time_sse3(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
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
    cout << "timer resolution is " << chrono::high_resolution_clock::period::den << " ticks per second" << endl;
    cout << "Configuration: " << nrStations << " stations, " << nrBaselines <<
	    " baselines, " << nrChannels << " channels, " << nrTimes << " time samples" << endl;

    const unsigned arraySize = nrStations*nrChannels*nrTimesWidth*NR_POLARIZATIONS*2;
    const unsigned visArraySize = nrBaselines*nrChannels*NR_POLARIZATIONS*NR_POLARIZATIONS*2;

    cout << "sample array size = " << ((arraySize * nrThreads * sizeof(float))/(1024*1024)) << " mbytes, "
	 << "vis array size = " << ((visArraySize * nrThreads * sizeof(float))/(1024*1024)) << " mbytes" << endl;

    __m128 dummyData;
    double maxGflops = computeMaxGflops(nrThreads, calcMaxFlops, &dummyData);
    cout << "peak flops with " << nrThreads << " threads is: " << maxGflops << " gflops" << std::endl;

    float* samples = new float[nrThreads*nrStations*nrChannels*nrTimesWidth*NR_POLARIZATIONS*2];
    float* visibilities= new float[nrThreads*nrBaselines*nrChannels*NR_POLARIZATIONS*NR_POLARIZATIONS*2];
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

    spawnCorrelatorThreads(CORRELATOR_1X1_SSE3, runCorrelator, samples, arraySize,
			   visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			   nrThreads, iter, maxGflops, verbose, validateResults);

    spawnCorrelatorThreads(CORRELATOR_1X1_TIME_SSE3, runCorrelator, samples, arraySize,
			   visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			   nrThreads, iter, maxGflops, verbose, validateResults);

    spawnCorrelatorThreads(CORRELATOR_2X2_SSE3, runCorrelator, samples, arraySize,
			   visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			   nrThreads, iter, maxGflops, verbose, validateResults);

    spawnCorrelatorThreads(CORRELATOR_2X2_TIME_SSE3, runCorrelator, samples, arraySize,
			   visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			   nrThreads, iter, maxGflops, verbose, validateResults);

    spawnCorrelatorThreads(CORRELATOR_3X2_TIME_SSE3, runCorrelator, samples, arraySize,
			   visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			   nrThreads, iter, maxGflops, verbose, validateResults);

    endCommon();

    delete[] samples;
    delete[] visibilities;

    return 0;
}
