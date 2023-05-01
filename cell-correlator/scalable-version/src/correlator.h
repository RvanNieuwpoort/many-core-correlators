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
#ifndef CORRELATOR_H_
#define CORRELATOR_H_

#define CELL_1X1 1
#define CELL_2X2 2
#define CELL_3X2 3
#define CELL_3X3 4
#define CELL_4X3 5

#define CELL_SIZE CELL_4X3

#if CELL_SIZE == CELL_1X1
#define CELL_WIDTH  1
#define CELL_HEIGHT 1
#elif CELL_SIZE == CELL_2X2
#define CELL_WIDTH  2
#define CELL_HEIGHT 2
#elif CELL_SIZE == CELL_3X2
#define CELL_WIDTH  3
#define CELL_HEIGHT 2
#elif CELL_SIZE == CELL_3X3
#define CELL_WIDTH  3
#define CELL_HEIGHT 3
#elif CELL_SIZE == CELL_4X3
#define CELL_WIDTH  4
#define CELL_HEIGHT 3
#else
#error illegal cell size
#endif

// ranges are for LOFAR
#define NR_STATIONS 64 // in the range 40 -- 77

#define NR_CHANNELS 120 // Always 256

#define NR_POLARISATIONS 2 // X and Y

#define MINOR_INTEGRATION_TIME 64 // was 64

#define MAJOR_INTEGRATION_TIME 768 // samples 768 is nice multiple of 16, represents one second

#define ITERATIONS 100

#define NR_BASELINES (NR_STATIONS * (NR_STATIONS + 1) / 2 + NR_STATIONS)

#define COMPLEX_SIZE 2

#define REAL 0
#define IMG 1

#define NR_SAMPLE_BUFFERS	2
#define NR_CORRELATION_BUFFERS	3

#define VISIBILITY_SIZE (NR_POLARISATIONS * NR_POLARISATIONS * COMPLEX_SIZE * sizeof(float))


#define MAX_BASELINES_IN_STRIP (NR_STATIONS * CELL_HEIGHT)
#define CELL_STORAGE_SIZE (CELL_WIDTH * CELL_HEIGHT * VISIBILITY_SIZE)

#define BEGIN_STATION  (NR_STATIONS % 2 ? 1 : 2) // first station to start with in the triangle

#define NR_SPUS 6 // max 6 on the ps3

#define SPU_FREQ 3.192
#define TOTAL_HOST_MEM_BANDWIDTH 25.8 // GB/s

#define FLOPS_PER_CYCLE 8 // 4 madds

#define BASELINE(station1, station2) ((station2) * ((station2) + 1) / 2 + (station1))

typedef float ppu_samples_type[NR_CHANNELS][NR_STATIONS][MAJOR_INTEGRATION_TIME][NR_POLARISATIONS][COMPLEX_SIZE] __attribute__ ((aligned(128)));
typedef float ppu_visibilities_type[NR_CHANNELS][NR_BASELINES][NR_POLARISATIONS][NR_POLARISATIONS][COMPLEX_SIZE] __attribute__ ((aligned(128)));

typedef struct spu_arguments {
    unsigned	spu_id;
    unsigned	ppu_samples; /* ppu_samples_type * */
    unsigned	ppu_visibilities; /* ppu_visibilities_type */
    char	pad[116];
} spu_arguments_type;

typedef struct spu_dma_list_elt {
  unsigned int size;
  unsigned int ea_low;
} spu_dma_list_elt_type;

#endif /*CORRELATOR_H_*/

