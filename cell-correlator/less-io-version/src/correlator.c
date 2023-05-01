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
#include "correlator.h"
#include "timer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <libspe2.h>
#include <pthread.h>
#include <complex.h>

typedef struct ppu_thread_data {
  unsigned int thread_id;
  spe_context_ptr_t spe_context;
  pthread_t pthread;
  unsigned int entry;
  spu_arguments_type* argp;
} ppu_thread_data_type;

ppu_thread_data_type threads[NR_SPUS];

ppu_samples_type	samples __attribute__ ((aligned(128)));
ppu_visibilities_type	visibilities __attribute__ ((aligned(128)));
ppu_zeros_type		zeros __attribute__ ((aligned(128)));
ppu_dma_lists_type	dma_lists __attribute__ ((aligned(128)));

spu_arguments_type	spu_arguments[NR_SPUS] __attribute__ ((aligned(128)));

static void barrier() {
  unsigned int dummy;

  for (unsigned spe_id = 0; spe_id < NR_SPUS; spe_id ++) {
    while (spe_out_mbox_status(threads[spe_id].spe_context) == 0)
      ;
    spe_out_mbox_read(threads[spe_id].spe_context, &dummy, 1);
  }

  for (unsigned spe_id = 0; spe_id < NR_SPUS; spe_id ++) {
    while (spe_in_mbox_status(threads[spe_id].spe_context) == 0)
      ;

    spe_in_mbox_write(threads[spe_id].spe_context, &dummy, 1, SPE_MBOX_ANY_NONBLOCKING);
  }
}

static void init()
{
 // init samples
 for (int station = 0; station < NR_STATIONS; station++)
   for (int channel = 0; channel < NR_CHANNELS; channel++)
     for (int t = 0; t < MAJOR_INTEGRATION_TIME; t++)
       for (int pol = 0; pol < NR_POLARISATIONS; pol ++)
	 * (complex float *) (samples[channel][station][t][pol]) = round(32 * drand48() - 16) + I * round(32 * drand48() - 16);

 memset(visibilities, 0, sizeof visibilities);
 memset(zeros, 0, sizeof zeros);

 for (unsigned channel = 0; channel < NR_CHANNELS; channel ++)
   for (unsigned time = 0; time < MAJOR_INTEGRATION_TIME / MINOR_INTEGRATION_TIME; time ++)
     for (unsigned stat = 0; stat < NR_STATIONS; stat ++) {
       dma_lists[channel][time][stat].ea_low = (unsigned) samples[channel][stat][time * MINOR_INTEGRATION_TIME];
       dma_lists[channel][time][stat].size   = sizeof(vector float) * MINOR_INTEGRATION_TIME;
     }
 
 int bufSize = sizeof(samples) / 1024;
 int outBufSize = sizeof(visibilities) / 1024;
 printf("HOST in buffer = %d KB @ %p, out buffer = %d KB @ %p\n", bufSize, samples, outBufSize, visibilities); 
 printf("host sample ptr = %p %u\n", samples, (int)samples);
}

void *ppu_pthread_function(void *arg)
{
  ppu_thread_data_type* datap = (ppu_thread_data_type*) arg;

  if (spe_context_run(datap->spe_context, &datap->entry, 0, datap->argp, NULL, NULL) < 0) {
    perror ("Failed running context\n");
    exit (1);
  }

  return 0;
}

void start_threads()
{
  extern spe_program_handle_t spu_correlator;

  for(int i = 0; i < NR_SPUS; i ++) {
    spu_arguments[i].spu_id = i;
    spu_arguments[i].ppu_samples = (unsigned) &samples;
    spu_arguments[i].ppu_visibilities = (unsigned) &visibilities;
    spu_arguments[i].ppu_zeros = (unsigned) &zeros;
    spu_arguments[i].ppu_dma_lists = (unsigned) &dma_lists;
    threads[i].thread_id = i;
    threads[i].argp = &spu_arguments[i];
    threads[i].entry = SPE_DEFAULT_ENTRY; // just run main()

    // use SPE_NOSCHED to make sure threads don't migrate
    if ((threads[i].spe_context = spe_context_create (0, 0)) == 0) {
      perror ("Failed creating context");
      exit (1);
    }

    /* Load program into context */
    if (spe_program_load(threads[i].spe_context, &spu_correlator) != 0) {
      perror ("Failed loading program");
      exit (1);
    }

    /* Create thread for each SPE context */
    if (pthread_create(&threads[i].pthread, 0, &ppu_pthread_function, &threads[i])) {
      perror ("Failed creating thread");
      exit (1);
    }
  }
}

static void end(void)
{
  /* Wait for SPU-thread to complete execution. */
  for (int i = 0; i < NR_SPUS; i ++) {
    if (pthread_join(threads[i].pthread, 0)) {
      perror("Failed pthread_join");
      exit (1);
    }
  }
}


static void check_visibilities()
{
  printf("checking visibilities ...\n");
  for (unsigned channel = 0; channel < NR_CHANNELS; channel ++)
    for (unsigned stat2 = 0, baseline = 0; stat2 < NR_STATIONS; stat2 ++)
      for (unsigned stat1 = 0; stat1 <= stat2; stat1 ++, baseline ++)
	for (unsigned pol1 = 0; pol1 < NR_POLARISATIONS; pol1 ++)
	  for (unsigned pol2 = 0; pol2 < NR_POLARISATIONS; pol2 ++) {
	    complex float sum = 0;

	    for (unsigned time = 0; time < MAJOR_INTEGRATION_TIME; time ++)
	      sum += * (complex float *) samples[channel][stat1][time][pol1] *
		conj(* (complex float *) samples[channel][stat2][time][pol2]);

	    if (creal(sum) != visibilities[channel][baseline][0][pol2][pol1] ||
	        cimag(sum) != visibilities[channel][baseline][1][pol2][pol1])
	      printf("ch: %u bl: %u%c-%u%c (%u) (%.0f, %.0f) - (%.0f, %.0f)\n",
		     channel,
		     stat1, 'X' + pol1, stat2, 'X' + pol2, baseline,
		     creal(sum), cimag(sum),
		     visibilities[channel][baseline][0][pol2][pol1],
		     visibilities[channel][baseline][1][pol2][pol1]);
	  }
  printf("checking completed ...\n");
}


int main(int argc, char** argv)
{
  struct timer timer;
  timer_init(&timer);
  
  init();
  
  start_threads();
  //  barrier();
  timer_start(&timer);
  //  barrier();
  //  timer_stop(&timer);
  
  end();
  check_visibilities();

  timer_stop(&timer);

  //double seconds = timer_time_in_seconds(&timer);
  return 0;
}
