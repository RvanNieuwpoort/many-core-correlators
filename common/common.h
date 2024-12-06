#ifndef _COMMON_CORRELATOR_H
#define _COMMON_CORRELATOR_H

#define CORRELATOR_REFERENCE     1
#define CORRELATOR_1X1           2
#define CORRELATOR_2X1           3
#define CORRELATOR_2X2           4

#define MAX_NR_STATIONS 64
#define MAX_NR_CELLS (MAX_NR_STATIONS * MAX_NR_STATIONS / 2)

#define BASELINE(STATION_1, STATION_2) ((STATION_2) * ((STATION_2) + 1) / 2 + (STATION_1))
#define NR_BASELINES(nrStations) (nrStations * (nrStations + 1) / 2)

extern unsigned calcNrCells(const unsigned w, const unsigned h, const unsigned nrStations);


#endif // _COMMON_CORRELATOR_H
