#include "../common/common.h"
#include "X86-64_correlator.h"

#include <iostream>
#include <chrono>
#include <immintrin.h>

using namespace std;

static const bool verbose = true;
static const bool validateResults = true;
static const unsigned nrStations = 256;
static const unsigned nrBaselines = NR_BASELINES(nrStations);
static const unsigned nrTimes = 768, nrTimesWidth = 768; // 770
static const unsigned nrChannels = 400;
static const unsigned iter = 10000;
static const unsigned nrThreads = 64;

void printVectorType()
{
#if defined(__AVX512F__)
    cout << "AVX512"; // 512 bit
#elif defined(__AVX2__)
    cout << "AVX2"; // 256 bit
#elif defined(__AVX__)
    cout << "AVX";  // 256 bit
#elif defined(__SSE3__)
    cout << "SSE3"; // 128 bit
#else
    cout << "NONE";
#endif
}


#if defined(__AVX512F__)
void* calcMaxFlops(void* data)
{
    // There are 32 registers. we use 16 varables here
    __m512 a = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 b = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 c = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 d = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 e = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 f = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 g = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 h = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 i = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 j = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 k = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 l = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 m = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 n = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 o = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);
    __m512 p = _mm512_set_ps(0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f, 0.2E-6f);

    // For AVX512, we can also use the minus. That's great for multiplying complex numbers
    for (unsigned long long x = 0; x < 1E9L; x++) {
	a += a * a;
	b += b * c;
	c += c * d;
	d -= d * e;
	e += e * f;
	f += f * g;
	g += g * h;
	h -= h * i;
	i += i * j;
	j += j * k;
	k += k * l;
	l -= l * m;
	m += m * n;
	n += n * o;
	o += o * p;
	p -= p * a;
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
    // the volatile is to make sure the compiler doesn't optimiza everything away
    volatile __m512 result = a + b + c + d + e + f + g + h + i + j + k + l + m + n + o + p; 
#pragma GCC diagnostic pop
    
    FlopTestParams* params = (FlopTestParams*) data;
    params->gFlops = VECTOR_WIDTH_IN_FLOATS * 16 * 2; // 16 lines * FMAs (2 ops) vector width
    
    return 0;
}
#elif defined(__SSE3__)
void print_m128(const __m128 a, const char* name)
{
    float result[VECTOR_WIDTH_IN_FLOATS];
    for(int i=0; i<VECTOR_WIDTH_IN_FLOATS; i++) {
	result[i] = -1.0f;
    }

    _mm_store_ps(result, a);

    cout << "float vec " << name << " = { ";
    for(int i=0; i<VECTOR_WIDTH_IN_FLOATS; i++) {
	cout << result[i] << " ";
    }
    cout << " }" << endl;
}

void* calcMaxFlops(void* data)
{
    // we use 15 varables here, with 16 we have one register to short (noted when inspecting assembly)
    __m128 a = _mm_set_ps1(0.2E-6f);
    __m128 b = _mm_set_ps1(0.3E-6f);
    __m128 c = _mm_set_ps1(0.4E-6f);
    __m128 d = _mm_set_ps1(0.5E-6f);
    __m128 e = _mm_set_ps1(0.6E-6f);
    __m128 f = _mm_set_ps1(0.7E-6f);
    __m128 g = _mm_set_ps1(0.8E-6f);
    __m128 h = _mm_set_ps1(0.9E-6f);
    __m128 i = _mm_set_ps1(0.10E-6f);
    __m128 j = _mm_set_ps1(0.11E-6f);
    __m128 k = _mm_set_ps1(0.12E-6f);
    __m128 l = _mm_set_ps1(0.13E-6f);
    __m128 m = _mm_set_ps1(0.14E-6f);
    __m128 n = _mm_set_ps1(0.15E-6f);
    __m128 o = _mm_set_ps1(0.16E-6f);
  
    for (unsigned long long x = 0; x < 1E9L; x++) {
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
	n += n * o;
	o += o * a;
    }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
    // the volatile is to make sure the compiler doesnoptimiza everything away
    volatile __m128 result = a + b + c + d + e + f + g + h + i + j + k + l + m + n + o; 
#pragma GCC diagnostic pop

    FlopTestParams* params = (FlopTestParams*) data;
    params->gFlops = VECTOR_WIDTH_IN_FLOATS * 15 * 2; // 15 lines * FMAs (2 ops) vector width
    
    return 0;
}
#else
#error unsupported
#endif // SSE3

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
	cout << "3x2-time-SSE3";
	break;
    case CORRELATOR_1X1_TIME_AVX512:
	cout << "1x1-time-AVX512";
	break;
    default:
	cout << "illegal correlator" << endl;
	exit(66);
    }
}

void* runCorrelator(void* data)
{
    CorrelatorParams* p = (CorrelatorParams*) data;

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
#if defined(__AVX512F__)
    case CORRELATOR_1X1_TIME_AVX512:
	p->ops = cpuCorrelator_1x1_time_avx512(p->samples, p->visibilities, nrTimes, nrTimesWidth, nrStations, nrChannels, &p->bytesLoaded, &p->bytesStored);
	break;
#endif // defined(__AVX512F__)
    default:
	cout << "illegal correlator" << endl;
	exit(66);
    }
    
#if PRINT_RESULT
    printResult(visibilities, nrChannels, nrBaselines);
#endif

    return 0;
}

int main()
{
    if(nrStations > MAX_NR_STATIONS) {
	cerr << "Cannot use more than " << MAX_NR_STATIONS << " stations." << endl;
    }

    cout << "using ";
    printVectorType();
    cout << " vector extensions" << endl;

    cout << "timer resolution is " << chrono::high_resolution_clock::period::den << " ticks per second" << endl;
    cout << "Configuration: " << nrStations << " stations, " << nrBaselines <<
	" baselines, " << nrChannels << " channels, " << nrTimes << " time samples, " << iter << " iterations" << endl;

    const size_t arraySize = nrStations*nrChannels*nrTimesWidth*NR_POLARIZATIONS*2;
    const size_t visArraySize = nrBaselines*nrChannels*NR_POLARIZATIONS*NR_POLARIZATIONS*2;

    cout << "sample array size = " << arraySize * nrThreads << " elements = " << (((double)arraySize * nrThreads * sizeof(float))/(1024.0*1024.0*1024.0)) << " gbytes, "
	 << "vis array size = " << visArraySize * nrThreads << " elements = " << (((double)visArraySize * nrThreads * sizeof(float))/(1024.0*1024.0)) << " mbytes" << endl;
   
//    double maxGflops = computeMaxGflops(nrThreads, calcMaxFlops, 0);
    double maxGflops = (2800.0 * 32 * 8 * 2)/1000.0; // For DAS-6: 2800 Mhz, 32 cores, 8-wide vectors, x2 for FMA
    cout << "peak flops with " << nrThreads << " threads is: " << maxGflops << " gflops" << std::endl;

    float* samples = new (align_val_t{ALIGNMENT}) float[nrThreads*arraySize];
    float* visibilities= new (align_val_t{ALIGNMENT}) float[nrThreads*visArraySize];

    for(size_t i=0; i<iter; i++) {
	initSamples(samples, nrThreads, nrTimes, nrTimesWidth, nrStations, nrChannels, arraySize);

	spawnCorrelatorThreads(CORRELATOR_1X1, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);
	
	spawnCorrelatorThreads(CORRELATOR_2X1, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);
	
	spawnCorrelatorThreads(CORRELATOR_2X2, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);

#if defined(__SSE3__)
	spawnCorrelatorThreads(CORRELATOR_1X1_SSE3, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);
	
	spawnCorrelatorThreads(CORRELATOR_1X1_TIME_SSE3, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);
	
	spawnCorrelatorThreads(CORRELATOR_2X2_SSE3, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);
	
	spawnCorrelatorThreads(CORRELATOR_2X2_TIME_SSE3, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);
	
	spawnCorrelatorThreads(CORRELATOR_3X2_TIME_SSE3, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);
#endif //  defined(__SSE3__)
#if defined(__AVX512F__)
	spawnCorrelatorThreads(CORRELATOR_1X1_TIME_AVX512, runCorrelator, samples, arraySize,
			       visibilities, visArraySize, nrTimes, nrStations, nrChannels,
			       nrThreads, maxGflops, verbose, validateResults);
#endif // defined(__AVX512F__)
    }
    
    endCommon();

    delete[] samples;
    delete[] visibilities;

    return 0;
}
