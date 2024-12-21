#ifndef _COMMON_CORRELATOR_H
#define _COMMON_CORRELATOR_H

#define CORRELATOR_REFERENCE     1
#define CORRELATOR_1X1           2
#define CORRELATOR_2X1           3
#define CORRELATOR_2X2           4

#define MAX_NR_STATIONS 64
#define MAX_NR_CELLS (MAX_NR_STATIONS * MAX_NR_STATIONS / 2)
#define NR_POLARIZATIONS 2

#define BASELINE(STATION_1, STATION_2) ((STATION_2) * ((STATION_2) + 1) / 2 + (STATION_1))
#define NR_BASELINES(nrStations) (nrStations * (nrStations + 1) / 2)

#define SAMPLE_INDEX(STATION, CHANNEL, TIME, POLARIZATION, REAL_OR_IMAG)		\
    (((((CHANNEL) * nrStations + (STATION)) * nrTimesWidth + (TIME)) * 2 + (POLARIZATION))*2 + (REAL_OR_IMAG))

#define VISIBILITIES_INDEX(BASELINE, CHANNEL, POLARIZATION_1, POLARIZATION_2, REAL_OR_IMAG) \
    (((((BASELINE) * nrChannels + (CHANNEL)) * 2 + (POLARIZATION_1)) * 2 + (POLARIZATION_2)) * 2 + (REAL_OR_IMAG))


extern void endCommon(void);
extern void initSamples(float* __restrict__ samples,
			const unsigned nrThreads, const unsigned nrTimes, const unsigned nrTimesWidth,
			const unsigned nrStations, const unsigned nrChannels, const unsigned arraySize);

extern unsigned long long referenceCorrelator(const float* __restrict__ samples, float* __restrict__ visibilities, 
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned calcNrCells(const unsigned w, const unsigned h, const unsigned nrStations);

extern void printCorrelatorType(int correlatorType);
extern void printResult(const float* visibilities, const unsigned nrChannels, const unsigned nrBaselines);

extern double computeMaxGflops(const unsigned nrThreads, void *(*runMaxFlopsTest) (void *), void* data);

extern void spawnCorrelatorThreads(int correlatorType, void *(*runCorrelator) (void *),
				   const float* __restrict__ samples, const unsigned arraySize,
				   float* __restrict__ visibilities, const unsigned visArraySize,
				   const unsigned nrTimes, const unsigned nrStations, const unsigned nrChannels,
				   const unsigned nrThreads, const unsigned iter, double maxFlops,
				   const bool verbose, const bool validateResults);

extern unsigned long long computeMissedBaselines(const float* __restrict__ samples, float* __restrict__ visibilities, 
					  const bool* __restrict__ missedBaselines,
					  const unsigned nrTimes, const unsigned nrTimesWidth,
					  const unsigned nrStations, const unsigned nrChannels,
					  unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern void checkResult(const float* __restrict__ samples, void *(*runCorrelator) (void *),
		 const float* __restrict__ visibilities, 
		 const unsigned nrThreads,
		 const unsigned nrTimes, const unsigned nrStations,
		 const unsigned nrChannels, const unsigned arraySize, const unsigned visArraySize);

extern unsigned long long cpuCorrelator_1x1(const float* __restrict__ samples, float* __restrict__ visibilities,
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long cpuCorrelator_2x1(const float* __restrict__ samples, float* __restrict__ visibilities,
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long cpuCorrelator_2x2(const float* __restrict__ samples, float* __restrict__ visibilities,
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

class CorrelatorParams {
public:
    int correlatorType;
    bool verbose;
    unsigned iter;
    const float* __restrict__ samples; 
    float* visibilities;
    unsigned long long ops, bytesLoaded, bytesStored;
};

class FlopTestParams {
public:
    unsigned long long gFlops;
};

#endif // _COMMON_CORRELATOR_H
