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

#include "cpu_correlator.h"

unsigned long long referenceCorrelator(float* samples, float* visibilities, 
				       unsigned nrTimes, unsigned nrTimesWidth, 
				       unsigned nrStations, unsigned nrChannels,
				       unsigned long long* bytesLoaded, unsigned long long* bytesStored)
{
    unsigned nrBaselines = nrStations * (nrStations + 1) / 2;

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned stat1 = 0; stat1 < nrStations; stat1 ++) {
	    for (unsigned stat0 = 0; stat0 <= stat1; stat0 ++) {
		for (unsigned time = 0; time < nrTimes; time ++) {
		    for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
			for (unsigned pol1 = 0; pol1 < 2; pol1 ++) { 
			    unsigned baseline = BASELINE(stat0, stat1);
			    unsigned vis_index_real = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 0);
			    unsigned vis_index_imag = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 1);

			    unsigned stat0_index_real = SAMPLE_INDEX(stat0, channel, time, pol0, 0);
			    unsigned stat0_index_imag = SAMPLE_INDEX(stat0, channel, time, pol0, 1);
			    unsigned stat1_index_real = SAMPLE_INDEX(stat1, channel, time, pol1, 0);
			    unsigned stat1_index_imag = SAMPLE_INDEX(stat1, channel, time, pol1, 1);

			    float stat0_real = samples[stat0_index_real];
			    float stat0_imag = samples[stat0_index_imag];
			    float stat1_real = samples[stat1_index_real];
			    float stat1_imag = samples[stat1_index_imag];

			    float res_real = stat0_real * -stat1_real - stat0_imag * -stat1_imag;
			    float res_imag = stat0_real * -stat1_imag + stat0_imag * stat1_real;

			    visibilities[vis_index_real] += res_real;
			    visibilities[vis_index_imag] += res_imag;
			}
		    }
		}
	    }
	}
    }

    *bytesLoaded = nrChannels * nrBaselines * nrTimes * 2L * 2L * 6 * sizeof(float); // 2 samples r + i, vis r + i
    *bytesStored = nrChannels * nrBaselines * nrTimes * 2L * 2L * 2 * sizeof(float); // 2 vis r + i

    return nrChannels * nrBaselines * nrTimes * 2L * 2L * 8;
}
