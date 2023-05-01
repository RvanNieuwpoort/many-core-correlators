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
