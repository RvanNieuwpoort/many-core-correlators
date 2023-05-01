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
// this code assumes the arrays are always multiples of 16 byes
/*
     visibilities are stored at the following indices (called baselines):

S1 | 0  1  2  3  4  5  6  7
     AA                       0
     .. AA                    1
     03 04 AA                 2
     06 07 .. AA              3
     10 11 12 13 AA           4  |
     15 16 17 18 .. AA        5  | this is a strip we keep in memory: one contiguous row from 10 .. 18
     21 22 23 24 25 26 AA     6
     28 29 30 31 32 33 .. AA  7
                              S2
*/

#include "correlator.h"
#include "spu_decrementer.h"

#include <stdio.h>
#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#define DO_DMA 1

static float samples[NR_SAMPLE_BUFFERS][NR_STATIONS][MINOR_INTEGRATION_TIME][NR_POLARISATIONS][COMPLEX_SIZE] __attribute__ ((aligned(128)));
static float visibilities[NR_CORRELATION_BUFFERS][MAX_BASELINES_IN_STRIP][NR_POLARISATIONS][NR_POLARISATIONS][COMPLEX_SIZE] __attribute__ ((aligned(128)));

static struct spu_dma_list_elt samples_dma_list[NR_STATIONS] __attribute__ ((aligned(128)));

static struct spu_arguments spu_arguments __attribute__ ((aligned(128)));

static unsigned first_channel, last_channel;

extern void oneByOne(vector float samples1[MINOR_INTEGRATION_TIME],
              vector float samples2[MINOR_INTEGRATION_TIME],
              vector float* visibilities1);

extern void twoByTwo(vector float samples1[MINOR_INTEGRATION_TIME],
              vector float samples2[MINOR_INTEGRATION_TIME],
              vector float samples3[MINOR_INTEGRATION_TIME],
              vector float samples4[MINOR_INTEGRATION_TIME],
              vector float* visibilities1,
              vector float* visibilities2,
              vector float* visibilities3,
              vector float* visibilities4);

extern void threeByTwo(
              vector float samples0[MINOR_INTEGRATION_TIME], // x + 0 
              vector float samples1[MINOR_INTEGRATION_TIME], // x + 1 
              vector float samples2[MINOR_INTEGRATION_TIME], // x + 2 
              vector float samples3[MINOR_INTEGRATION_TIME], // y + 0 
              vector float samples4[MINOR_INTEGRATION_TIME], // y + 1
              vector float* vis03,
              vector float* vis13,
              vector float* vis23,
              vector float* vis04,
              vector float* vis14,
              vector float* vis24);

extern void threeByThree(
              vector float samples0[MINOR_INTEGRATION_TIME], // x + 0 
              vector float samples1[MINOR_INTEGRATION_TIME], // x + 1 
              vector float samples2[MINOR_INTEGRATION_TIME], // x + 2 
              vector float samples3[MINOR_INTEGRATION_TIME], // y + 0 
              vector float samples4[MINOR_INTEGRATION_TIME], // y + 1
              vector float samples5[MINOR_INTEGRATION_TIME], // y + 2
              vector float* vis03,
              vector float* vis13,
              vector float* vis23,
              vector float* vis04,
              vector float* vis14,
              vector float* vis24,
              vector float* vis05,
              vector float* vis15,
              vector float* vis25);

void fourByThree(
              vector float samples0[MINOR_INTEGRATION_TIME], // x + 0 
              vector float samples1[MINOR_INTEGRATION_TIME], // x + 1 
              vector float samples2[MINOR_INTEGRATION_TIME], // x + 2 
              vector float samples3[MINOR_INTEGRATION_TIME], // x + 3 
              vector float samples4[MINOR_INTEGRATION_TIME], // y + 0 
              vector float samples5[MINOR_INTEGRATION_TIME], // y + 1
              vector float samples6[MINOR_INTEGRATION_TIME], // y + 2
              vector float* vis04,
              vector float* vis14,
              vector float* vis24,
              vector float* vis34,
              vector float* vis05,
              vector float* vis15,
              vector float* vis25,
              vector float* vis35,
              vector float* vis06,
              vector float* vis16,
              vector float* vis26,
              vector float* vis36);

static unsigned long long calc_ops() {
	unsigned long long ops = 0;
	for (int stationY = BEGIN_STATION; stationY < NR_STATIONS; stationY += CELL_HEIGHT) {    	
          for (int stationX = 0; stationX + CELL_WIDTH <= stationY; stationX += CELL_WIDTH) {
       	    ops += 8L * CELL_WIDTH * CELL_HEIGHT * NR_POLARISATIONS * NR_POLARISATIONS;
          }
	}
	
	return ops * MAJOR_INTEGRATION_TIME * (last_channel - first_channel);
}

static unsigned long long calc_loads() { // in bytes
	unsigned long long loads = 0;
	for (int stationY = BEGIN_STATION; stationY < NR_STATIONS; stationY += CELL_HEIGHT) {    	
          for (int stationX = 0; stationX + CELL_WIDTH <= stationY; stationX += CELL_WIDTH) {
       	    loads += (CELL_WIDTH + CELL_HEIGHT) * sizeof(vector float);
          }
	}
	
	return loads * MAJOR_INTEGRATION_TIME * (last_channel - first_channel);
}

void print_vector(char *var, vector float val) {
	printf("Vector %s is: {%f, %f, %f, %f}\n", var, spu_extract(val, 0),
	 spu_extract(val, 1), spu_extract(val, 2), spu_extract(val, 3));
}

static void barrier(void) {
  spu_write_out_mbox(0);
  spu_read_in_mbox();
}

static void init(unsigned long long argp) {
  mfc_get(&spu_arguments, (unsigned) argp, sizeof(spu_arguments), 0, 0, 0);
  mfc_write_tag_mask(1 << 0);
  mfc_read_tag_status_all();

  first_channel =  spu_arguments.spu_id      * NR_CHANNELS / NR_SPUS;
  last_channel  = (spu_arguments.spu_id + 1) * NR_CHANNELS / NR_SPUS;

  for(int i=0; i<NR_STATIONS; i++) {
    samples_dma_list[i].size = sizeof(samples[0][0]);
  }

  if(spu_arguments.spu_id == 0) {
    printf("SPU sample dma size = %ld bytes\n", sizeof(samples[0][0]));
    printf("SPU in buffers = %ld KB @ %p, out buffers = %ld B @ %p\n", sizeof(samples) / 1024, samples, sizeof(visibilities), visibilities);
  }

  printf("I am spu %d, calculating channels %3d - %3d\n", spu_arguments.spu_id, first_channel, last_channel);
}

static void printSamples(float samples[NR_STATIONS][MINOR_INTEGRATION_TIME][NR_POLARISATIONS][COMPLEX_SIZE])
{
	for(unsigned station = 0; station < NR_STATIONS; station++) {
		for(unsigned time=0; time<MINOR_INTEGRATION_TIME; time++) {
			printf("station %d, time %d: x(%f %f) y(%f %f)\n",
				station, time,
				samples[station][time][0][0],
				samples[station][time][0][1],
				samples[station][time][1][0],
				samples[station][time][1][1]);
		}
	}
}

static void printVisibilities(int buf)
{
		for(int width=0; width < MAX_BASELINES_IN_STRIP; width++) {
			fprintf(stderr, "vis: xx(%f %f) xy(%f %f) yx(%f %f) yy(%f %f)\n",
				visibilities[buf][width][0][0][0], visibilities[buf][width][0][0][1], 
				visibilities[buf][width][0][1][0], visibilities[buf][width][0][1][1], 
				visibilities[buf][width][1][0][0], visibilities[buf][width][1][0][1], 
				visibilities[buf][width][1][1][0], visibilities[buf][width][1][1][1]
			);
		}
}

// samples are non-contiguous in mem of the ppu...
static inline void load_samples(int channel, int time, int buffer) {
#if DO_DMA
  for(int stat=0; stat<NR_STATIONS; stat++) {
    samples_dma_list[stat].ea_low = (unsigned) (*((ppu_samples_type*)spu_arguments.ppu_samples))[channel][stat][time];
//	    printf("SPE sample for station %u, time %u = %x %u, size = %u\n", stat, time, samples_dma_list[stat].ea_low, samples_dma_list[stat].ea_low, samples_dma_list[stat].size);
  }

  mfc_getl(samples[buffer], 0, samples_dma_list, sizeof(struct spu_dma_list_elt) * NR_STATIONS, buffer, 0, 0);
#endif
}

static inline void wait_for_dma_samples(int buffer) {
#if DO_DMA
  mfc_write_tag_mask(1 << buffer);
  mfc_read_tag_status_all();
#endif
}

/* mfc_get (local store space, main mem address, bytes to DMA, tag associated to DMA (031), optional, optional) */
static inline void load_visibilities(unsigned visibilities_ptr, int size, int buffer) {
#if DO_DMA
    mfc_get(visibilities[buffer], visibilities_ptr, size, buffer + NR_SAMPLE_BUFFERS, 0, 0);
#endif
}

static inline void wait_for_dma_visibilities1(int buffer) {
#if DO_DMA
  mfc_write_tag_mask(1 << (buffer + NR_SAMPLE_BUFFERS));
  mfc_read_tag_status_all();
#endif
}

static inline void wait_for_dma_visibilities2(int buffer1, int buffer2) {
#if DO_DMA
  mfc_write_tag_mask(1 << (buffer1 + NR_SAMPLE_BUFFERS) | 1 << (buffer2 + NR_SAMPLE_BUFFERS));
  mfc_read_tag_status_all();
#endif
}

static inline void store_visibilities(unsigned visibilities_ptr, int size, int buffer) {
#if DO_DMA
  mfc_put(visibilities[buffer], visibilities_ptr, size, buffer + NR_SAMPLE_BUFFERS, 0, 0);
#endif
}

static void correlateChannel(int samplesBuffer, unsigned channelvisibilities) {
  unsigned buf = 0;
  unsigned visibilitiesSize = 0;
  unsigned visStripStart = BASELINE(0, BEGIN_STATION);
  unsigned visibilities_ptr = channelvisibilities + visStripStart * VISIBILITY_SIZE;

  load_visibilities(visibilities_ptr, visibilitiesSize, 0);

  for (unsigned stationY = BEGIN_STATION; stationY < NR_STATIONS; stationY += CELL_HEIGHT) {
    unsigned next_buf = buf + 1 == 3 ? 0 : buf + 1;
    unsigned next_stationY = stationY + CELL_HEIGHT;
    unsigned next_visStripStart = BASELINE(0, next_stationY);
    unsigned next_ptr = channelvisibilities + next_visStripStart * VISIBILITY_SIZE;

    unsigned nrCellsInNextStrip = (next_stationY) / CELL_WIDTH;
    unsigned lastBaselineInNextStrip = BASELINE(0, next_stationY + (CELL_HEIGHT-1)) + CELL_WIDTH * nrCellsInNextStrip - 1;
    unsigned next_size = (lastBaselineInNextStrip - next_visStripStart + 1) * VISIBILITY_SIZE;

    // prefetch the next buffer
    if(stationY + CELL_HEIGHT < NR_STATIONS) {
      wait_for_dma_visibilities2(buf, next_buf);
      load_visibilities(next_ptr, next_size, next_buf);
    } else { // last one
      wait_for_dma_visibilities1(buf);
    }

    for (unsigned stationX = 0; stationX + CELL_WIDTH <= stationY; stationX += CELL_WIDTH) {
      unsigned visStrip1 = BASELINE(stationX, stationY) - visStripStart;

#if CELL_SIZE == CELL_1X1
      oneByOne(
               (vector float*) samples[samplesBuffer][stationX], 
	       (vector float*) samples[samplesBuffer][stationY], 
	       (vector float*) visibilities[buf][visStrip1]);
#elif CELL_SIZE == CELL_2X2
      unsigned visStrip2 = BASELINE(stationX, stationY+1) - visStripStart;
      twoByTwo(
               (vector float*) samples[samplesBuffer][stationX+0], 
	       (vector float*) samples[samplesBuffer][stationX+1],
	       (vector float*) samples[samplesBuffer][stationY+0], 
	       (vector float*) samples[samplesBuffer][stationY+1], 
	       (vector float*) visibilities[buf][visStrip1+0], 
	       (vector float*) visibilities[buf][visStrip1+1], 
	       (vector float*) visibilities[buf][visStrip2+0], 
	       (vector float*) visibilities[buf][visStrip2+1]); 
#elif CELL_SIZE == CELL_3X2
      unsigned visStrip2 = BASELINE(stationX, stationY+1) - visStripStart;
      threeByTwo(
               (vector float*) samples[samplesBuffer][stationX+0], 
	       (vector float*) samples[samplesBuffer][stationX+1],
	       (vector float*) samples[samplesBuffer][stationX+2],
	       (vector float*) samples[samplesBuffer][stationY+0], 
	       (vector float*) samples[samplesBuffer][stationY+1], 
	       (vector float*) visibilities[buf][visStrip1+0], 
	       (vector float*) visibilities[buf][visStrip1+1], 
	       (vector float*) visibilities[buf][visStrip1+2], 
	       (vector float*) visibilities[buf][visStrip2+0], 
	       (vector float*) visibilities[buf][visStrip2+1], 
	       (vector float*) visibilities[buf][visStrip2+2]);
#elif CELL_SIZE == CELL_3X3
      unsigned visStrip2 = BASELINE(stationX, stationY+1) - visStripStart;
      unsigned visStrip3 = BASELINE(stationX, stationY+2) - visStripStart;
      threeByThree(
               (vector float*) samples[samplesBuffer][stationX+0], 
	       (vector float*) samples[samplesBuffer][stationX+1],
	       (vector float*) samples[samplesBuffer][stationX+2],
	       (vector float*) samples[samplesBuffer][stationY+0], 
	       (vector float*) samples[samplesBuffer][stationY+1], 
	       (vector float*) samples[samplesBuffer][stationY+2], 
	       (vector float*) visibilities[buf][visStrip1+0], 
	       (vector float*) visibilities[buf][visStrip1+1], 
	       (vector float*) visibilities[buf][visStrip1+2], 
	       (vector float*) visibilities[buf][visStrip2+0], 
	       (vector float*) visibilities[buf][visStrip2+1], 
	       (vector float*) visibilities[buf][visStrip2+2],
	       (vector float*) visibilities[buf][visStrip3+0], 
	       (vector float*) visibilities[buf][visStrip3+1], 
	       (vector float*) visibilities[buf][visStrip3+2]);
#elif CELL_SIZE == CELL_4X3
      unsigned visStrip2 = BASELINE(stationX, stationY+1) - visStripStart;
      unsigned visStrip3 = BASELINE(stationX, stationY+2) - visStripStart;
      fourByThree(
               (vector float*) samples[samplesBuffer][stationX+0], 
	       (vector float*) samples[samplesBuffer][stationX+1],
	       (vector float*) samples[samplesBuffer][stationX+2],
	       (vector float*) samples[samplesBuffer][stationX+3],
	       (vector float*) samples[samplesBuffer][stationY+0], 
	       (vector float*) samples[samplesBuffer][stationY+1], 
	       (vector float*) samples[samplesBuffer][stationY+2], 
	       (vector float*) visibilities[buf][visStrip1+0], 
	       (vector float*) visibilities[buf][visStrip1+1], 
	       (vector float*) visibilities[buf][visStrip1+2], 
	       (vector float*) visibilities[buf][visStrip1+3], 
	       (vector float*) visibilities[buf][visStrip2+0], 
	       (vector float*) visibilities[buf][visStrip2+1], 
	       (vector float*) visibilities[buf][visStrip2+2],
	       (vector float*) visibilities[buf][visStrip2+3],
	       (vector float*) visibilities[buf][visStrip3+0], 
	       (vector float*) visibilities[buf][visStrip3+1], 
	       (vector float*) visibilities[buf][visStrip3+2],
	       (vector float*) visibilities[buf][visStrip3+3]);
#else
#error illegal cell size
#endif
    }

//    printVisibilities(buf);
    
    // and store the strip again.
    store_visibilities(visibilities_ptr, visibilitiesSize, buf);

    buf = next_buf;
    visibilities_ptr = next_ptr;
    visStripStart = next_visStripStart;
    visibilitiesSize = next_size;
  }

  // wait for the last store
  buf = buf == 0 ? 3 : buf - 1;
  wait_for_dma_visibilities1(buf);
}

static void correlate() {
  int buffer = 0;
  unsigned channelvisibilities = spu_arguments.ppu_visibilities + first_channel * NR_BASELINES * VISIBILITY_SIZE;
  unsigned int timeIterations = MAJOR_INTEGRATION_TIME / MINOR_INTEGRATION_TIME;

  for (unsigned channel = first_channel; channel < last_channel; channel++) {
      // load first samples
      load_samples(first_channel, 0, buffer);
      
      for (unsigned t = 0; t < timeIterations; t++) {
	// prefetch samples for the next time slot (except for the last iteration)
	if(t < timeIterations - 1) {
	  load_samples(channel, (t+1) * MINOR_INTEGRATION_TIME, !buffer);
	}

	wait_for_dma_samples(buffer);

//	printSamples(samples[buffer]);

	correlateChannel(buffer, channelvisibilities);

	buffer = !buffer;
      }

      channelvisibilities += NR_BASELINES * VISIBILITY_SIZE;
  }
}

int main(unsigned long long id, unsigned long long argp, unsigned long long envp) {
  id = id; // fix compiler warning
  envp = envp; // fix compiler warning

  init(argp);
  decrementer_init();

  //  barrier();

  unsigned int start = decrementer_get();
  for(int i=0; i<ITERATIONS; i++) {
    correlate();
  }
  unsigned int end = decrementer_get();

  //  barrier();

  double seconds = decrementer_seconds(start - end);
  unsigned long long ops = calc_ops() * ITERATIONS;
  double gflops = (ops / seconds) / 1000000000.0;
  double efficiency = (gflops / (SPU_FREQ * FLOPS_PER_CYCLE)) * 100.0;
  unsigned long long loads = calc_loads() * ITERATIONS;
  double gb = (double)loads / (1024.0*1024.0*1024.0);
  double throughput = gb / seconds;

  printf("SPU SECONDS = %f, ops = %lld, %6.3f gflops, efficiency is %6.3f %%, loads = %6.3f GB, %6.3f GB/s, efficiency is %6.3f %%\n", 
	seconds, ops, gflops, efficiency,
	gb, throughput, throughput / TOTAL_HOST_MEM_BANDWIDTH * 100.0);

  return 0;
}

