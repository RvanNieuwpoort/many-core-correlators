#include <cstdlib>
#include <iostream>
#include <string.h>

#include "gpu_correlator.h"
#include "gpu_complex.h"


void correlateOnHost(const complex<float> *samples, complex<float> *visibilities)
{
    unsigned nrBaselines = nrStations * (nrStations + 1) / 2;
    memset(visibilities, 0, nrBaselines * nrChannels * 2 * 2 * sizeof(complex<float>));

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned stat1 = 0; stat1 < nrStations; stat1 ++) {
	    for (unsigned stat0 = 0; stat0 <= stat1; stat0 ++) {
		for (unsigned time = 0; time < nrTimes; time ++) {
		    for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
			for (unsigned pol1 = 0; pol1 < 2; pol1 ++) { 
			    visibilities[VISIBILITIES_INDEX(BASELINE(stat0, stat1), channel, pol0, pol1)]
				+=  samples[SAMPLE_INDEX(stat0, channel, time, pol0)] 
				* ~(samples[SAMPLE_INDEX(stat1, channel, time, pol1)]);
			}
		    }
		}
	    }
	}
    }
}

void correlateMissingOnHost(const unsigned nrMissed, const complex<float> *samples, complex<float> *visibilities)
{
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for(unsigned missed = 0; missed<nrMissed; missed++) {
	    unsigned stat0 = missedStatXHost[missed];
	    unsigned stat1 = missedStatXHost[missed];

	    for (unsigned time = 0; time < nrTimes; time ++) {
		for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
		    for (unsigned pol1 = 0; pol1 < 2; pol1 ++) { 
			visibilities[VISIBILITIES_INDEX(BASELINE(stat0, stat1), channel, pol0, pol1)]
			    +=  samples[SAMPLE_INDEX(stat0, channel, time, pol0)] 
			    * ~(samples[SAMPLE_INDEX(stat1, channel, time, pol1)]);
		    }
		}
	    }
	}
    }
}
