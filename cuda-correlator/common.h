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
