#ifndef CORRELATOR_H_
#define CORRELATOR_H_


// ranges are for LOFAR
#define NR_STATIONS 64 // in the range 40 -- 77

#define NR_CHANNELS 120 // Always 256

#define NR_POLARISATIONS 2 // X and Y

#define MINOR_INTEGRATION_TIME 64 // was 64

#define MAJOR_INTEGRATION_TIME 768 // samples 768 is nice multiple of 16, represents one second

#define NR_BASELINES (NR_STATIONS * (NR_STATIONS + 1) / 2 + NR_STATIONS)
// #define NR_4X4_BLOCKS ((NR_STATIONS / 4) * ((NR_STATIONS / 4) - 1) / 2)

#define COMPLEX_SIZE 2

#define REAL 0
#define IMAG 1

#define DOUBLE_BUFFERING	2


#define NR_SPUS 6 // max 6 on the ps3

#define likely(p)	__builtin_expect((p), 1)
#define unlikely(p)	__builtin_expect((p), 0)

typedef struct spu_dma_list_elt {
  unsigned int size;
  unsigned int ea_low;
} spu_dma_list_elt_type;

typedef float  ppu_samples_type[NR_CHANNELS][NR_STATIONS][MAJOR_INTEGRATION_TIME][NR_POLARISATIONS][COMPLEX_SIZE] __attribute__ ((aligned(128)));
typedef float  ppu_visibilities_type[NR_CHANNELS][NR_BASELINES][COMPLEX_SIZE][NR_POLARISATIONS][NR_POLARISATIONS] __attribute__ ((aligned(128)));
typedef char   ppu_zeros_type[16384] __attribute__ ((aligned(128)));
typedef struct spu_dma_list_elt ppu_dma_lists_type[NR_CHANNELS][MAJOR_INTEGRATION_TIME / MINOR_INTEGRATION_TIME][NR_STATIONS] __attribute__ ((aligned(128)));

typedef struct spu_arguments {
    unsigned	spu_id;
    unsigned	ppu_samples; /* ppu_samples_type * */
    unsigned	ppu_visibilities; /* ppu_visibilities_type */
    unsigned	ppu_zeros;
    unsigned	ppu_dma_lists;
    char	pad[108];
} spu_arguments_type;

#endif /*CORRELATOR_H_*/

