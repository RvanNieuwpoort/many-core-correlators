#ifndef _COMMON_H
#define _COMMON_H

extern unsigned fillCellToStatTable(unsigned w, unsigned h, unsigned nrStations, char* cellToStatX, char* cellToStatY);
extern unsigned calcNrCells(unsigned w, unsigned h, unsigned nrStations);
extern unsigned power (unsigned base, unsigned n);
extern unsigned calcNrOutputsPerCell(unsigned w, unsigned h, unsigned nrPolarizations);

extern void checkResult(size_t visibilitiesSize, complex<float> *hostSamples,
			complex<float> *hostVisibilities, unsigned nrBaselines, unsigned nrCells);

extern void convertBufferLayout(const complex<float> *src, complex<float>* dst);

extern unsigned cellWidth;
extern unsigned cellHeight;
extern bool printResult;
extern unsigned char cellToStatXHost[MAX_CELLS], cellToStatYHost[MAX_CELLS];

extern unsigned nrStations;
extern unsigned nrTimes;
extern unsigned nrTimesWidth;
extern unsigned nrChannels;
extern unsigned nrPolarizations;


#endif // _COMMON_H
