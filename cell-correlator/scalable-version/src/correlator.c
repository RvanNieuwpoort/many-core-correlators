#include "correlator.h"
#include "timer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <libspe2.h>
#include <pthread.h>

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

static void init() {
 // init samples
 for (int station = 0; station < NR_STATIONS; station++) {
   for (int channel = 0; channel < NR_CHANNELS; channel++) {
     for (int t = 0; t < MAJOR_INTEGRATION_TIME; t++) {
	 samples[channel][station][t][0][REAL] = 1.0f;
	 samples[channel][station][t][0][IMG]  = 2.0f;
	 samples[channel][station][t][1][REAL] = 3.0f;
	 samples[channel][station][t][1][IMG]  = 4.0f;
     }
   }
 }

 // init visibilities (result)
 for (int channel = 0; channel < NR_CHANNELS; channel++) {
   for (int bl = 0; bl < NR_BASELINES; bl++) {
     for (int polarisation1 = 0; polarisation1 < NR_POLARISATIONS; polarisation1++) {
       for (int polarisation2 = 0; polarisation2 < NR_POLARISATIONS; polarisation2++) {
	 visibilities[channel][bl][polarisation1][polarisation2][REAL] = 0.0f;
	 visibilities[channel][bl][polarisation1][polarisation2][IMG]  = 0.0f;
       }
     }
   }
 }
 
 int bufSize = sizeof(samples) / 1024;
 int outBufSize = sizeof(visibilities) / 1024;
 printf("HOST in buffer = %d KB @ %p, out buffer = %d KB @ %p\n", bufSize, samples, outBufSize, visibilities); 
 printf("host sample ptr = %p %u\n", samples, (int)samples);
}

void* ppu_pthread_function(void *arg) {
  ppu_thread_data_type* datap = (ppu_thread_data_type*) arg;

  if (spe_context_run(datap->spe_context, &datap->entry, 0, datap->argp, NULL, NULL) < 0) {
    perror ("Failed running context\n");
    exit (1);
  }

  pthread_exit(NULL);
}

void start_threads()
{
  extern spe_program_handle_t spu_correlator;

  for(int i=0; i<NR_SPUS; i++) {
    spu_arguments[i].spu_id = i;
    spu_arguments[i].ppu_samples = (unsigned) &samples;
    spu_arguments[i].ppu_visibilities = (unsigned) &visibilities;
    threads[i].thread_id = i;
    threads[i].argp = &spu_arguments[i];
    threads[i].entry = SPE_DEFAULT_ENTRY; // just run main()

    // use SPE_NOSCHED to make sure threads don't migrate
    if ((threads[i].spe_context = spe_context_create (SPE_NOSCHED, NULL)) == NULL) {
      perror ("Failed creating context");
      exit (1);
    }

    /* Load program into context */
    if (spe_program_load(threads[i].spe_context, &spu_correlator) != 0) {
      perror ("Failed loading program");
      exit (1);
    }

    /* Create thread for each SPE context */
    if (pthread_create(&threads[i].pthread, NULL, &ppu_pthread_function, &threads[i])) {
      perror ("Failed creating thread");
      exit (1);
    }
  }
}

static void end(void) {
  /* Wait for SPU-thread to complete execution. */
  for (int i=0; i<NR_SPUS; i++) {
    if (pthread_join(threads[i].pthread, NULL)) {
      perror("Failed pthread_join");
      exit (1);
    }
  }
}

static unsigned long long calc_ops() {
    unsigned long long ops = 0;

    for (int station2 = BEGIN_STATION; station2 < NR_STATIONS; station2 += CELL_HEIGHT) {
        for (int station1 = 0; station1 + CELL_WIDTH <= station2; station1 += CELL_WIDTH) {
       	    ops += 8L * CELL_HEIGHT * CELL_WIDTH * NR_POLARISATIONS * NR_POLARISATIONS;
        }
    }
	
    return ops * MAJOR_INTEGRATION_TIME * NR_CHANNELS;
}

static unsigned long long calc_loads() { // in bytes
	unsigned long long loads = 0;
	for (int stationY = BEGIN_STATION; stationY < NR_STATIONS; stationY += CELL_HEIGHT) {    	
          for (int stationX = 0; stationX + CELL_WIDTH <= stationY; stationX += CELL_WIDTH) {
       	    loads += (CELL_WIDTH + CELL_HEIGHT) * sizeof(vector float);
          }
	}
	
	return loads * MAJOR_INTEGRATION_TIME * NR_CHANNELS;
}

static unsigned calcNrCells() {
	unsigned nrCells = 0;
    for (int station2 = BEGIN_STATION; station2 < NR_STATIONS; station2 += CELL_HEIGHT) {
        for (int station1 = 0; station1 + CELL_WIDTH <= station2; station1 += CELL_WIDTH) {
       	    nrCells++;
        }
    }
	
    return nrCells;
}


void printResult()
{
	for(unsigned channel = 0; channel < NR_CHANNELS; channel++) {
		for(unsigned baseline = 0; baseline < NR_BASELINES; baseline++) {
			printf("channel %u, baseline %u: (%f %f) (%f %f) (%f %f) (%f %f)\n",
				channel, baseline, 
				visibilities[channel][baseline][0][0][0],
				visibilities[channel][baseline][0][0][1],
				visibilities[channel][baseline][0][1][0],
				visibilities[channel][baseline][0][1][1],
				visibilities[channel][baseline][1][0][0],
				visibilities[channel][baseline][1][0][1],
				visibilities[channel][baseline][1][1][0],
				visibilities[channel][baseline][1][1][1]);
		}
	}


}

int main(int argc, char** argv) {
  struct timer timer;
  timer_init(&timer);
  
  init();
  
  unsigned nrCells = calcNrCells();
  printf("using %d x %d correlator, nr cells = %u\n", CELL_WIDTH, CELL_HEIGHT, nrCells);

  start_threads();
  //  barrier();
  timer_start(&timer);
  //  barrier();
  //  timer_stop(&timer);
  
  end();

  timer_stop(&timer);

 // printResult();
  
  unsigned long long ops = calc_ops() * ITERATIONS;
  double seconds = timer_time_in_seconds(&timer);
  double gflops = (ops / seconds) / 1000000000.0;
  double efficiency = (gflops / (SPU_FREQ * FLOPS_PER_CYCLE * NR_SPUS)) * 100;
  
  unsigned long long loads = calc_loads() * ITERATIONS;
  double gb = (double)loads / (1024.0*1024.0*1024.0);
  double throughput = gb / seconds;

  printf("performed %lld operations in %6.6f seconds, %6.3f GFLOPS, efficiency is %6.3f %%, loaded %6.3f GB, %6.3f GB/s, efficiency is %6.3f\n", ops, seconds, gflops, efficiency, gb, throughput, (throughput / TOTAL_HOST_MEM_BANDWIDTH ) * 100.0);
  
  return 0;
}
