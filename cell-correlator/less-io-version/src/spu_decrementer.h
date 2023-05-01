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
#ifndef SPU_DECREMENTER_H
#define SPU_DECREMENTER_H

#include <spu_intrinsics.h>

static inline void decrementer_init (void) {
    spu_writech(SPU_WrDec, -1);
}

static inline unsigned int decrementer_get (void) {
    return spu_readch(SPU_RdDec);
}

float decrementer_seconds (unsigned int t);

float decrementer_msecs (unsigned int t);

float decrementer_cycles (unsigned int t);

void decrementer_sleep (unsigned int t);

#endif // SPU_DECREMENTER_H
