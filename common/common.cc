#include "common.h"
#include <iostream>


unsigned calcNrCells(const unsigned w, const unsigned h, const unsigned nrStations) 
{
    unsigned nrCells = 0;
    for (unsigned stationY = h; stationY < nrStations - h + 1; stationY += h) {
	for (unsigned stationX = 0; stationX + w <= stationY; stationX += w) {
	    nrCells++;
	}
    }
    return nrCells;
}

void printResult(const float* visibilities, const unsigned nrChannels, const unsigned nrBaselines)
{
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned baseline = 0; baseline < nrBaselines; baseline ++) {
	    for (unsigned pol0 = 0; pol0 < NR_POLARIZATIONS; pol0 ++) {
		for (unsigned pol1 = 0; pol1 < NR_POLARIZATIONS; pol1 ++) {
		    std::cout.precision(15);
		    std::cout << visibilities[VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 0)] << ", " 
			      << visibilities[VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 1)]
			      << std::endl;
		    
		}
	    }
	}
    }
}
