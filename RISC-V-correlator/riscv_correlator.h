#ifndef _RISCV_CORRELATOR_H
#define _RISCV_CORRELATOR_H

#define CORRELATOR_REFERENCE     1
#define CORRELATOR_1X1           2
#define CORRELATOR_2X2           3

#define MAX_NR_SATIONS 64
#define MAX_CELLS (MAX_NR_SATIONS * MAX_NR_SATIONS / 2)

#define BASELINE(STATION_1, STATION_2)			\
    ((STATION_2) * ((STATION_2) + 1) / 2 + (STATION_1))

#define SAMPLE_INDEX(STATION, CHANNEL, TIME, POLARIZATION, REAL_OR_IMAG)		\
    (((((CHANNEL) * nrStations + (STATION)) * nrTimesWidth + (TIME)) * 2 + (POLARIZATION))*2 + (REAL_OR_IMAG))

#define VISIBILITIES_INDEX(BASELINE, CHANNEL, POLARIZATION_1, POLARIZATION_2, REAL_OR_IMAG) \
    (((((BASELINE) * nrChannels + (CHANNEL)) * 2 + (POLARIZATION_1)) * 2 + (POLARIZATION_2)) * 2 + (REAL_OR_IMAG))


extern unsigned long long referenceCorrelator(const float* __restrict__ samples, float* __restrict__ visibilities, 
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long riscvCorrelator_1x1(const float* __restrict__ samples, float* __restrict__ visibilities,
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long riscvCorrelator_2x2(const float* __restrict__ samples, float* __restrict__ visibilities,
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

#endif // _RISCV_CORRELATOR_H
