#ifndef _RISCV_CORRELATOR_H
#define _RISCV_CORRELATOR_H

#define CORRELATOR_REFERENCE     1
#define CORRELATOR_1X1           2
#define CORRELATOR_2X1           3
#define CORRELATOR_2X2           4

extern unsigned long long riscvCorrelator_1x1(const float* __restrict__ samples, float* __restrict__ visibilities,
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long riscvCorrelator_2x2(const float* __restrict__ samples, float* __restrict__ visibilities,
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

extern unsigned long long riscvCorrelator_2x1(const float* __restrict__ samples, float* __restrict__ visibilities,
					      const unsigned nrTimes, const unsigned nrTimesWidth,
					      const unsigned nrStations, const unsigned nrChannels,
					      unsigned long long* bytesLoaded, unsigned long long* bytesStored);

#endif // _RISCV_CORRELATOR_H
