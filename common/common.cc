#include "common.h"
#include <cstdio>

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
