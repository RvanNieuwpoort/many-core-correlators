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
#ifndef _CORRELATOR_H
#define _CORRELATOR_H

extern void correlate_1x1(float* baselineToStat1Host, float* baselineToStat2Host, float* hostSamples, float* hostVisReal, float* hostVisImag);
extern void correlate_1x1_vec(float* baselineToStat1Host, float* baselineToStat2Host, float* hostSamples, float* hostVisReal, float* hostVisImag);

extern void correlate_2x2_vec(float* cellToStatXHost, float* cellToStatYHost, float* hostSamples, float* hostVisReal, float* hostVisImag);

extern void startLoadTimer();
extern void stopLoadTimer();

extern void startCorrelateTimer();
extern void stopCorrelateTimer();

extern void startStoreTimer();
extern void stopStoreTimer();

#endif // _CORRELATOR_H
