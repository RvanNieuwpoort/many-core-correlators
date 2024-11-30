#ifndef _RISCV_CORRELATOR_H
#define _RISCV_CORRELATOR_H

#define MAX_CELLS 2048

#define CORRELATOR_REFERENCE     1
#define CORRELATOR_1X1           2


#define BASELINE(STATION_1, STATION_2)			\
    ((STATION_2) * ((STATION_2) + 1) / 2 + (STATION_1))

#define SAMPLE_INDEX(STATION, CHANNEL, TIME, POLARIZATION, REAL_OR_IMAG)		\
    (((((CHANNEL) * nrStations + (STATION)) * nrTimesWidth + (TIME)) * 2 + (POLARIZATION))*2 + (REAL_OR_IMAG))

#define VISIBILITIES_INDEX(BASELINE, CHANNEL, POLARIZATION_1, POLARIZATION_2, REAL_OR_IMAG) \
    (((((BASELINE) * nrChannels + (CHANNEL)) * 2 + (POLARIZATION_1)) * 2 + (POLARIZATION_2)) * 2 + (REAL_OR_IMAG))

#define LOOP_COUNT(NR_CELLS, NR_THREADS) ((NR_CELLS) + (NR_THREADS) - 1) / (NR_THREADS) * (NR_THREADS);

extern const unsigned nrStations;
extern const unsigned nrTimes, nrTimesWidth;
extern const unsigned nrChannels;
extern const unsigned nrPolarizations;

extern unsigned long long referenceCorrelator(float* samples, float* visibilities, 
					      unsigned nrTimes, unsigned nrTimesWidth, unsigned nrStations, unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long riscvCorrelator_1x1(float* samples, float* visibilities,
					    unsigned nrTimes, unsigned nrTimesWidth, unsigned nrStations, unsigned nrChannels,
					    unsigned long long* bytesLoaded, unsigned long long* bytesStored);


#endif // _RISCV_CORRELATOR_H
