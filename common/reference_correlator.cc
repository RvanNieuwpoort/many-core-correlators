#include <cstdlib>
#include <iostream>

#include "../common/common.h"

using namespace std;

unsigned long long referenceCorrelator(const float* __restrict__ samples, float* __restrict__ visibilities, 
				       const unsigned nrTimes, const unsigned nrTimesWidth, 
				       const unsigned nrStations, const unsigned nrChannels,
				       unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned nrBaselines = NR_BASELINES(nrStations);

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned stat1 = 0; stat1 < nrStations; stat1 ++) {
	    for (unsigned stat0 = 0; stat0 <= stat1; stat0 ++) {
		for (unsigned time = 0; time < nrTimes; time ++) {
		    for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
			for (unsigned pol1 = 0; pol1 < 2; pol1 ++) { 
			    unsigned baseline = BASELINE(stat0, stat1);
			    size_t stat0_index_real = SAMPLE_INDEX(stat0, channel, time, pol0, 0);
			    size_t stat0_index_imag = SAMPLE_INDEX(stat0, channel, time, pol0, 1);
			    size_t stat1_index_real = SAMPLE_INDEX(stat1, channel, time, pol1, 0);
			    size_t stat1_index_imag = SAMPLE_INDEX(stat1, channel, time, pol1, 1);

			    float stat0_real = samples[stat0_index_real];
			    float stat0_imag = samples[stat0_index_imag];
			    float stat1_real = samples[stat1_index_real];
			    float stat1_imag = samples[stat1_index_imag];

			    float res_real = stat0_real *  stat1_real + stat0_imag * stat1_imag;
			    float res_imag = stat0_imag *  stat1_real - stat0_real * stat1_imag;

			    size_t vis_index_real = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 0);
			    size_t vis_index_imag = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 1);

			    visibilities[vis_index_real] += res_real;
			    visibilities[vis_index_imag] += res_imag;
			}
		    }
		}
	    }
	}
    }

    unsigned long long ops = nrChannels * nrBaselines * nrTimes * 16L * 2L;

    *bytesLoaded = nrChannels * nrBaselines * nrTimes * 2L * 2L * 6L * sizeof(float); // 2 samples r + i, vis r + i
    *bytesStored = nrChannels * nrBaselines * nrTimes * 2L * 2L * 2L * sizeof(float); // 2 vis r + i

    return ops;
}
