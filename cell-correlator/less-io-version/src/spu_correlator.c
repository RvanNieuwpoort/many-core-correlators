#include "correlator.h"
#include "spu_decrementer.h"

#include <stdio.h>
#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#define DO_DMA

vector float samples[DOUBLE_BUFFERING][NR_STATIONS][MINOR_INTEGRATION_TIME] __attribute__ ((aligned(128)));
vector float visibilities[NR_BASELINES][COMPLEX_SIZE] __attribute__ ((aligned(128)));

#if defined DO_DMA
struct spu_dma_list_elt samples_dma_list[DOUBLE_BUFFERING][NR_STATIONS] __attribute__ ((aligned(128)));
struct spu_dma_list_elt visibilities_dma_list[NR_CHANNELS / NR_SPUS][(NR_BASELINES * 32 + 16383) / 16384] __attribute__ ((aligned(128))); // FIXME: assert(NR_CHANNELS % NR_SPUS == 0)
struct spu_dma_list_elt zero_dma_list[(NR_BASELINES * 32 + 16383) / 16384] __attribute__ ((aligned(128)));
#endif

#if 0
unsigned stat0s[NR_4X4_BLOCKS], statAs[NR_4X4_BLOCKS];
unsigned baselines[NR_4X4_BLOCKS];
#endif

static struct spu_arguments spu_arguments __attribute__ ((aligned(128)));

static unsigned first_channel, last_channel;


static inline void print_vector_float(char *var, vector float val)
{
  printf("Vector %s is: {%f, %f, %f, %f}\n",
	 var,
	 spu_extract(val, 0),
	 spu_extract(val, 1),
	 spu_extract(val, 2),
	 spu_extract(val, 3));
}


static inline void print_vector_unsigned(char *var, vector unsigned val)
{
  printf("Vector %s is: {0x%08X, 0x%08X, 0x%08X, 0x%08X}\n",
	 var,
	 spu_extract(val, 0),
	 spu_extract(val, 1),
	 spu_extract(val, 2),
	 spu_extract(val, 3));
}


static inline void barrier(void)
{
  spu_write_out_mbox(0);
  spu_read_in_mbox();
}


#if defined DO_DMA

static inline void load_samples(unsigned channel, unsigned time, unsigned buffer)
{
  // first load DMA list, then DMA all samples
  unsigned tag = buffer;

  mfc_get(samples_dma_list[buffer], (unsigned) (*((ppu_dma_lists_type *) spu_arguments.ppu_dma_lists))[channel][time / MINOR_INTEGRATION_TIME], sizeof samples_dma_list[0], tag, 0, 0);
  mfc_getlf(samples[buffer], 0, samples_dma_list[buffer], sizeof samples_dma_list[0], tag, 0, 0);
}


static inline void store_visibilities(unsigned channel)
{
  // first store all visibilities, then clear memory
  unsigned tag = 2;

  mfc_putl(visibilities, 0, visibilities_dma_list[channel - first_channel], sizeof visibilities_dma_list[0], tag, 0, 0);
  mfc_getlf(visibilities, 0, zero_dma_list, sizeof zero_dma_list, tag, 0, 0);
}


static inline void wait_for_dma_samples(unsigned buffer)
{
  unsigned tag = buffer;

  mfc_write_tag_mask(1 << tag);
  mfc_read_tag_status_all();
}


static inline void wait_for_dma_visibilities()
{
  unsigned tag = 2;

  mfc_write_tag_mask(1 << tag);
  mfc_read_tag_status_all();
}

#endif


static const vector unsigned char permute_0202 = (const vector unsigned char) {
  0x00, 0x01, 0x02, 0x03,
  0x08, 0x09, 0x0A, 0x0B,
  0x00, 0x01, 0x02, 0x03,
  0x08, 0x09, 0x0A, 0x0B,
};

static const vector unsigned char permute_0022 = (const vector unsigned char) {
  0x00, 0x01, 0x02, 0x03,
  0x00, 0x01, 0x02, 0x03,
  0x08, 0x09, 0x0A, 0x0B,
  0x08, 0x09, 0x0A, 0x0B,
};

static const vector unsigned char permute_1313 = (const vector unsigned char) {
  0x04, 0x05, 0x06, 0x07,
  0x0C, 0x0D, 0x0E, 0x0F,
  0x04, 0x05, 0x06, 0x07,
  0x0C, 0x0D, 0x0E, 0x0F,
};

static const vector unsigned char permute_1133 = (const vector unsigned char) {
  0x04, 0x05, 0x06, 0x07,
  0x04, 0x05, 0x06, 0x07,
  0x0C, 0x0D, 0x0E, 0x0F,
  0x0C, 0x0D, 0x0E, 0x0F,
};


static void correlate_all_4x4_blocks(
  vector float samples[NR_STATIONS][MINOR_INTEGRATION_TIME],
  vector float visibilities[NR_BASELINES][COMPLEX_SIZE]
)
{
  vector float vis_A0_r, vis_A1_r, vis_A2_r, vis_A3_r;
  vector float vis_B0_r, vis_B1_r, vis_B2_r, vis_B3_r;
  vector float vis_C0_r, vis_C1_r, vis_C2_r, vis_C3_r;
  vector float vis_D0_r, vis_D1_r, vis_D2_r, vis_D3_r;
  vector float vis_A0_i, vis_A1_i, vis_A2_i, vis_A3_i;
  vector float vis_B0_i, vis_B1_i, vis_B2_i, vis_B3_i;
  vector float vis_C0_i, vis_C1_i, vis_C2_i, vis_C3_i;
  vector float vis_D0_i, vis_D1_i, vis_D2_i, vis_D3_i;

  vector float sample_0x0y0x0y_r, sample_1x1y1x1y_r, sample_2x2y2x2y_r, sample_3x3y3x3y_r, sample_AxAxAyAy_r, sample_BxBxByBy_r, sample_CxCxCyCy_r, sample_DxDxDyDy_r;
  vector float sample_0x0y0x0y_i, sample_1x1y1x1y_i, sample_2x2y2x2y_i, sample_3x3y3x3y_i, sample_AxAxAyAy_i, sample_BxBxByBy_i, sample_CxCxCyCy_i, sample_DxDxDyDy_i;

  vector float sample_0, sample_1, sample_2, sample_3;
  vector float sample_A, sample_B, sample_C, sample_D;

  //for (unsigned _4x4_block = 0; _4x4_block < NR_4X4_BLOCKS; _4x4_block ++) 
    //unsigned stat0 = stat0s[_4x4_block];
    //unsigned statA = statAs[_4x4_block];
    //unsigned bl    = baselines[_4x4_block];
    //printf("statA = %u, stat0 = %u, bl = %u\n", statA, stat0, bl);

  for (unsigned statA = 4; statA < NR_STATIONS; statA += 4) {
    unsigned bl = statA * (statA + 1) / 2;

    sample_A = samples[statA + 0][0];
    sample_B = samples[statA + 1][0];
    sample_C = samples[statA + 2][0];
    sample_D = samples[statA + 3][0];
    sample_0 = samples[0][0];
    sample_1 = samples[1][0];
    sample_2 = samples[2][0];
    sample_3 = samples[3][0];

    for (unsigned stat0 = 0; stat0 < statA; stat0 += 4, bl += 4) {
      vis_A0_r = visibilities[bl + 0 * statA + 0 + 0][REAL];
      vis_A0_i = visibilities[bl + 0 * statA + 0 + 0][IMAG];
      vis_A1_r = visibilities[bl + 0 * statA + 0 + 1][REAL];
      vis_A1_i = visibilities[bl + 0 * statA + 0 + 1][IMAG];
      vis_A2_r = visibilities[bl + 0 * statA + 0 + 2][REAL];
      vis_A2_i = visibilities[bl + 0 * statA + 0 + 2][IMAG];
      vis_A3_r = visibilities[bl + 0 * statA + 0 + 3][REAL];
      vis_A3_i = visibilities[bl + 0 * statA + 0 + 3][IMAG];
      vis_B0_r = visibilities[bl + 1 * statA + 1 + 0][REAL];
      vis_B0_i = visibilities[bl + 1 * statA + 1 + 0][IMAG];
      vis_B1_r = visibilities[bl + 1 * statA + 1 + 1][REAL];
      vis_B1_i = visibilities[bl + 1 * statA + 1 + 1][IMAG];
      vis_B2_r = visibilities[bl + 1 * statA + 1 + 2][REAL];
      vis_B2_i = visibilities[bl + 1 * statA + 1 + 2][IMAG];
      vis_B3_r = visibilities[bl + 1 * statA + 1 + 3][REAL];
      vis_B3_i = visibilities[bl + 1 * statA + 1 + 3][IMAG];
      vis_C0_r = visibilities[bl + 2 * statA + 3 + 0][REAL];
      vis_C0_i = visibilities[bl + 2 * statA + 3 + 0][IMAG];
      vis_C1_r = visibilities[bl + 2 * statA + 3 + 1][REAL];
      vis_C1_i = visibilities[bl + 2 * statA + 3 + 1][IMAG];
      vis_C2_r = visibilities[bl + 2 * statA + 3 + 2][REAL];
      vis_C2_i = visibilities[bl + 2 * statA + 3 + 2][IMAG];
      vis_C3_r = visibilities[bl + 2 * statA + 3 + 3][REAL];
      vis_C3_i = visibilities[bl + 2 * statA + 3 + 3][IMAG];
      vis_D0_r = visibilities[bl + 3 * statA + 6 + 0][REAL];
      vis_D0_i = visibilities[bl + 3 * statA + 6 + 0][IMAG];
      vis_D1_r = visibilities[bl + 3 * statA + 6 + 1][REAL];
      vis_D1_i = visibilities[bl + 3 * statA + 6 + 1][IMAG];
      vis_D2_r = visibilities[bl + 3 * statA + 6 + 2][REAL];
      vis_D2_i = visibilities[bl + 3 * statA + 6 + 2][IMAG];
      vis_D3_r = visibilities[bl + 3 * statA + 6 + 3][REAL];
      vis_D3_i = visibilities[bl + 3 * statA + 6 + 3][IMAG];

      sample_0x0y0x0y_r = spu_shuffle(sample_0, sample_0, permute_0202);
      sample_1x1y1x1y_r = spu_shuffle(sample_1, sample_1, permute_0202);
      sample_2x2y2x2y_r = spu_shuffle(sample_2, sample_2, permute_0202);
      sample_3x3y3x3y_r = spu_shuffle(sample_3, sample_3, permute_0202);
      sample_AxAxAyAy_r = spu_shuffle(sample_A, sample_A, permute_0022);
      sample_BxBxByBy_r = spu_shuffle(sample_B, sample_B, permute_0022);
      sample_CxCxCyCy_r = spu_shuffle(sample_C, sample_C, permute_0022);
      sample_DxDxDyDy_r = spu_shuffle(sample_D, sample_D, permute_0022);
      sample_0x0y0x0y_i = spu_shuffle(sample_0, sample_0, permute_1313);
      sample_1x1y1x1y_i = spu_shuffle(sample_1, sample_1, permute_1313);
      sample_2x2y2x2y_i = spu_shuffle(sample_2, sample_2, permute_1313);
      sample_3x3y3x3y_i = spu_shuffle(sample_3, sample_3, permute_1313);
      sample_AxAxAyAy_i = spu_shuffle(sample_A, sample_A, permute_1133);
      sample_BxBxByBy_i = spu_shuffle(sample_B, sample_B, permute_1133);
      sample_CxCxCyCy_i = spu_shuffle(sample_C, sample_C, permute_1133);
      sample_DxDxDyDy_i = spu_shuffle(sample_D, sample_D, permute_1133);

      vis_A0_r = spu_madd(sample_AxAxAyAy_r, sample_0x0y0x0y_r, vis_A0_r);
      vis_A1_r = spu_madd(sample_AxAxAyAy_r, sample_1x1y1x1y_r, vis_A1_r);
      vis_A2_r = spu_madd(sample_AxAxAyAy_r, sample_2x2y2x2y_r, vis_A2_r);
      vis_A3_r = spu_madd(sample_AxAxAyAy_r, sample_3x3y3x3y_r, vis_A3_r);
      vis_B0_r = spu_madd(sample_BxBxByBy_r, sample_0x0y0x0y_r, vis_B0_r);
      vis_B1_r = spu_madd(sample_BxBxByBy_r, sample_1x1y1x1y_r, vis_B1_r);
      vis_B2_r = spu_madd(sample_BxBxByBy_r, sample_2x2y2x2y_r, vis_B2_r);
      vis_B3_r = spu_madd(sample_BxBxByBy_r, sample_3x3y3x3y_r, vis_B3_r);
      vis_C0_r = spu_madd(sample_CxCxCyCy_r, sample_0x0y0x0y_r, vis_C0_r);
      vis_C1_r = spu_madd(sample_CxCxCyCy_r, sample_1x1y1x1y_r, vis_C1_r);
      vis_C2_r = spu_madd(sample_CxCxCyCy_r, sample_2x2y2x2y_r, vis_C2_r);
      vis_C3_r = spu_madd(sample_CxCxCyCy_r, sample_3x3y3x3y_r, vis_C3_r);
      vis_D0_r = spu_madd(sample_DxDxDyDy_r, sample_0x0y0x0y_r, vis_D0_r);
      vis_D1_r = spu_madd(sample_DxDxDyDy_r, sample_1x1y1x1y_r, vis_D1_r);
      vis_D2_r = spu_madd(sample_DxDxDyDy_r, sample_2x2y2x2y_r, vis_D2_r);
      vis_D3_r = spu_madd(sample_DxDxDyDy_r, sample_3x3y3x3y_r, vis_D3_r);

      for (int t = 1; t < MINOR_INTEGRATION_TIME; t += 1) {
	vis_A0_i = spu_madd(sample_AxAxAyAy_r, sample_0x0y0x0y_i, vis_A0_i);
	vis_A1_i = spu_madd(sample_AxAxAyAy_r, sample_1x1y1x1y_i, vis_A1_i);
	vis_A2_i = spu_madd(sample_AxAxAyAy_r, sample_2x2y2x2y_i, vis_A2_i);
	vis_A3_i = spu_madd(sample_AxAxAyAy_r, sample_3x3y3x3y_i, vis_A3_i);
	vis_B0_i = spu_madd(sample_BxBxByBy_r, sample_0x0y0x0y_i, vis_B0_i);
	vis_B1_i = spu_madd(sample_BxBxByBy_r, sample_1x1y1x1y_i, vis_B1_i);
	vis_B2_i = spu_madd(sample_BxBxByBy_r, sample_2x2y2x2y_i, vis_B2_i);
	vis_B3_i = spu_madd(sample_BxBxByBy_r, sample_3x3y3x3y_i, vis_B3_i);
	vis_C0_i = spu_madd(sample_CxCxCyCy_r, sample_0x0y0x0y_i, vis_C0_i);
	vis_C1_i = spu_madd(sample_CxCxCyCy_r, sample_1x1y1x1y_i, vis_C1_i);
	vis_C2_i = spu_madd(sample_CxCxCyCy_r, sample_2x2y2x2y_i, vis_C2_i);
	vis_C3_i = spu_madd(sample_CxCxCyCy_r, sample_3x3y3x3y_i, vis_C3_i);
	vis_D0_i = spu_madd(sample_DxDxDyDy_r, sample_0x0y0x0y_i, vis_D0_i);
	vis_D1_i = spu_madd(sample_DxDxDyDy_r, sample_1x1y1x1y_i, vis_D1_i);
	vis_D2_i = spu_madd(sample_DxDxDyDy_r, sample_2x2y2x2y_i, vis_D2_i);
	vis_D3_i = spu_madd(sample_DxDxDyDy_r, sample_3x3y3x3y_i, vis_D3_i);

	sample_0 = samples[stat0 + 0][t];
	sample_1 = samples[stat0 + 1][t];
	sample_2 = samples[stat0 + 2][t];
	sample_3 = samples[stat0 + 3][t];
	sample_A = samples[statA + 0][t];
	sample_B = samples[statA + 1][t];
	sample_C = samples[statA + 2][t];
	sample_D = samples[statA + 3][t];

	vis_A0_r = spu_madd(sample_AxAxAyAy_i, sample_0x0y0x0y_i, vis_A0_r);
	vis_A1_r = spu_madd(sample_AxAxAyAy_i, sample_1x1y1x1y_i, vis_A1_r);
	vis_A2_r = spu_madd(sample_AxAxAyAy_i, sample_2x2y2x2y_i, vis_A2_r);
	vis_A3_r = spu_madd(sample_AxAxAyAy_i, sample_3x3y3x3y_i, vis_A3_r);
	vis_B0_r = spu_madd(sample_BxBxByBy_i, sample_0x0y0x0y_i, vis_B0_r);
	vis_B1_r = spu_madd(sample_BxBxByBy_i, sample_1x1y1x1y_i, vis_B1_r);
	vis_B2_r = spu_madd(sample_BxBxByBy_i, sample_2x2y2x2y_i, vis_B2_r);
	vis_B3_r = spu_madd(sample_BxBxByBy_i, sample_3x3y3x3y_i, vis_B3_r);
	vis_C0_r = spu_madd(sample_CxCxCyCy_i, sample_0x0y0x0y_i, vis_C0_r);
	vis_C1_r = spu_madd(sample_CxCxCyCy_i, sample_1x1y1x1y_i, vis_C1_r);
	vis_C2_r = spu_madd(sample_CxCxCyCy_i, sample_2x2y2x2y_i, vis_C2_r);
	vis_C3_r = spu_madd(sample_CxCxCyCy_i, sample_3x3y3x3y_i, vis_C3_r);
	vis_D0_r = spu_madd(sample_DxDxDyDy_i, sample_0x0y0x0y_i, vis_D0_r);
	vis_D1_r = spu_madd(sample_DxDxDyDy_i, sample_1x1y1x1y_i, vis_D1_r);
	vis_D2_r = spu_madd(sample_DxDxDyDy_i, sample_2x2y2x2y_i, vis_D2_r);
	vis_D3_r = spu_madd(sample_DxDxDyDy_i, sample_3x3y3x3y_i, vis_D3_r);

	vis_A0_i = spu_nmsub(sample_AxAxAyAy_i, sample_0x0y0x0y_r, vis_A0_i);
	vis_A1_i = spu_nmsub(sample_AxAxAyAy_i, sample_1x1y1x1y_r, vis_A1_i);
	vis_A2_i = spu_nmsub(sample_AxAxAyAy_i, sample_2x2y2x2y_r, vis_A2_i);
	vis_A3_i = spu_nmsub(sample_AxAxAyAy_i, sample_3x3y3x3y_r, vis_A3_i);
	vis_B0_i = spu_nmsub(sample_BxBxByBy_i, sample_0x0y0x0y_r, vis_B0_i);
	vis_B1_i = spu_nmsub(sample_BxBxByBy_i, sample_1x1y1x1y_r, vis_B1_i);
	vis_B2_i = spu_nmsub(sample_BxBxByBy_i, sample_2x2y2x2y_r, vis_B2_i);
	vis_B3_i = spu_nmsub(sample_BxBxByBy_i, sample_3x3y3x3y_r, vis_B3_i);
	vis_C0_i = spu_nmsub(sample_CxCxCyCy_i, sample_0x0y0x0y_r, vis_C0_i);
	vis_C1_i = spu_nmsub(sample_CxCxCyCy_i, sample_1x1y1x1y_r, vis_C1_i);
	vis_C2_i = spu_nmsub(sample_CxCxCyCy_i, sample_2x2y2x2y_r, vis_C2_i);
	vis_C3_i = spu_nmsub(sample_CxCxCyCy_i, sample_3x3y3x3y_r, vis_C3_i);
	vis_D0_i = spu_nmsub(sample_DxDxDyDy_i, sample_0x0y0x0y_r, vis_D0_i);
	vis_D1_i = spu_nmsub(sample_DxDxDyDy_i, sample_1x1y1x1y_r, vis_D1_i);
	vis_D2_i = spu_nmsub(sample_DxDxDyDy_i, sample_2x2y2x2y_r, vis_D2_i);
	vis_D3_i = spu_nmsub(sample_DxDxDyDy_i, sample_3x3y3x3y_r, vis_D3_i);

	sample_0x0y0x0y_r = spu_shuffle(sample_0, sample_0, permute_0202);
	sample_1x1y1x1y_r = spu_shuffle(sample_1, sample_1, permute_0202);
	sample_2x2y2x2y_r = spu_shuffle(sample_2, sample_2, permute_0202);
	sample_3x3y3x3y_r = spu_shuffle(sample_3, sample_3, permute_0202);
	sample_AxAxAyAy_r = spu_shuffle(sample_A, sample_A, permute_0022);
	sample_BxBxByBy_r = spu_shuffle(sample_B, sample_B, permute_0022);
	sample_CxCxCyCy_r = spu_shuffle(sample_C, sample_C, permute_0022);
	sample_DxDxDyDy_r = spu_shuffle(sample_D, sample_D, permute_0022);
	sample_0x0y0x0y_i = spu_shuffle(sample_0, sample_0, permute_1313);
	sample_1x1y1x1y_i = spu_shuffle(sample_1, sample_1, permute_1313);
	sample_2x2y2x2y_i = spu_shuffle(sample_2, sample_2, permute_1313);
	sample_3x3y3x3y_i = spu_shuffle(sample_3, sample_3, permute_1313);
	sample_AxAxAyAy_i = spu_shuffle(sample_A, sample_A, permute_1133);
	sample_BxBxByBy_i = spu_shuffle(sample_B, sample_B, permute_1133);
	sample_CxCxCyCy_i = spu_shuffle(sample_C, sample_C, permute_1133);
	sample_DxDxDyDy_i = spu_shuffle(sample_D, sample_D, permute_1133);

	vis_A0_r = spu_madd(sample_AxAxAyAy_r, sample_0x0y0x0y_r, vis_A0_r);
	vis_A1_r = spu_madd(sample_AxAxAyAy_r, sample_1x1y1x1y_r, vis_A1_r);
	vis_A2_r = spu_madd(sample_AxAxAyAy_r, sample_2x2y2x2y_r, vis_A2_r);
	vis_A3_r = spu_madd(sample_AxAxAyAy_r, sample_3x3y3x3y_r, vis_A3_r);
	vis_B0_r = spu_madd(sample_BxBxByBy_r, sample_0x0y0x0y_r, vis_B0_r);
	vis_B1_r = spu_madd(sample_BxBxByBy_r, sample_1x1y1x1y_r, vis_B1_r);
	vis_B2_r = spu_madd(sample_BxBxByBy_r, sample_2x2y2x2y_r, vis_B2_r);
	vis_B3_r = spu_madd(sample_BxBxByBy_r, sample_3x3y3x3y_r, vis_B3_r);
	vis_C0_r = spu_madd(sample_CxCxCyCy_r, sample_0x0y0x0y_r, vis_C0_r);
	vis_C1_r = spu_madd(sample_CxCxCyCy_r, sample_1x1y1x1y_r, vis_C1_r);
	vis_C2_r = spu_madd(sample_CxCxCyCy_r, sample_2x2y2x2y_r, vis_C2_r);
	vis_C3_r = spu_madd(sample_CxCxCyCy_r, sample_3x3y3x3y_r, vis_C3_r);
	vis_D0_r = spu_madd(sample_DxDxDyDy_r, sample_0x0y0x0y_r, vis_D0_r);
	vis_D1_r = spu_madd(sample_DxDxDyDy_r, sample_1x1y1x1y_r, vis_D1_r);
	vis_D2_r = spu_madd(sample_DxDxDyDy_r, sample_2x2y2x2y_r, vis_D2_r);
	vis_D3_r = spu_madd(sample_DxDxDyDy_r, sample_3x3y3x3y_r, vis_D3_r);
      }

      sample_A = samples[statA + 0][0];
      sample_B = samples[statA + 1][0];
      sample_C = samples[statA + 2][0];
      sample_D = samples[statA + 3][0];
      sample_0 = samples[stat0 + 4][0];
      sample_1 = samples[stat0 + 5][0];
      sample_2 = samples[stat0 + 6][0];
      sample_3 = samples[stat0 + 7][0];

      vis_A0_i = spu_madd(sample_AxAxAyAy_r, sample_0x0y0x0y_i, vis_A0_i);
      vis_A1_i = spu_madd(sample_AxAxAyAy_r, sample_1x1y1x1y_i, vis_A1_i);
      vis_A2_i = spu_madd(sample_AxAxAyAy_r, sample_2x2y2x2y_i, vis_A2_i);
      vis_A3_i = spu_madd(sample_AxAxAyAy_r, sample_3x3y3x3y_i, vis_A3_i);
      vis_B0_i = spu_madd(sample_BxBxByBy_r, sample_0x0y0x0y_i, vis_B0_i);
      vis_B1_i = spu_madd(sample_BxBxByBy_r, sample_1x1y1x1y_i, vis_B1_i);
      vis_B2_i = spu_madd(sample_BxBxByBy_r, sample_2x2y2x2y_i, vis_B2_i);
      vis_B3_i = spu_madd(sample_BxBxByBy_r, sample_3x3y3x3y_i, vis_B3_i);
      vis_C0_i = spu_madd(sample_CxCxCyCy_r, sample_0x0y0x0y_i, vis_C0_i);
      vis_C1_i = spu_madd(sample_CxCxCyCy_r, sample_1x1y1x1y_i, vis_C1_i);
      vis_C2_i = spu_madd(sample_CxCxCyCy_r, sample_2x2y2x2y_i, vis_C2_i);
      vis_C3_i = spu_madd(sample_CxCxCyCy_r, sample_3x3y3x3y_i, vis_C3_i);
      vis_D0_i = spu_madd(sample_DxDxDyDy_r, sample_0x0y0x0y_i, vis_D0_i);
      vis_D1_i = spu_madd(sample_DxDxDyDy_r, sample_1x1y1x1y_i, vis_D1_i);
      vis_D2_i = spu_madd(sample_DxDxDyDy_r, sample_2x2y2x2y_i, vis_D2_i);
      vis_D3_i = spu_madd(sample_DxDxDyDy_r, sample_3x3y3x3y_i, vis_D3_i);

      vis_A0_r = spu_madd(sample_AxAxAyAy_i, sample_0x0y0x0y_i, vis_A0_r);
      vis_A1_r = spu_madd(sample_AxAxAyAy_i, sample_1x1y1x1y_i, vis_A1_r);
      vis_A2_r = spu_madd(sample_AxAxAyAy_i, sample_2x2y2x2y_i, vis_A2_r);
      vis_A3_r = spu_madd(sample_AxAxAyAy_i, sample_3x3y3x3y_i, vis_A3_r);
      vis_B0_r = spu_madd(sample_BxBxByBy_i, sample_0x0y0x0y_i, vis_B0_r);
      vis_B1_r = spu_madd(sample_BxBxByBy_i, sample_1x1y1x1y_i, vis_B1_r);
      vis_B2_r = spu_madd(sample_BxBxByBy_i, sample_2x2y2x2y_i, vis_B2_r);
      vis_B3_r = spu_madd(sample_BxBxByBy_i, sample_3x3y3x3y_i, vis_B3_r);
      vis_C0_r = spu_madd(sample_CxCxCyCy_i, sample_0x0y0x0y_i, vis_C0_r);
      vis_C1_r = spu_madd(sample_CxCxCyCy_i, sample_1x1y1x1y_i, vis_C1_r);
      vis_C2_r = spu_madd(sample_CxCxCyCy_i, sample_2x2y2x2y_i, vis_C2_r);
      vis_C3_r = spu_madd(sample_CxCxCyCy_i, sample_3x3y3x3y_i, vis_C3_r);
      vis_D0_r = spu_madd(sample_DxDxDyDy_i, sample_0x0y0x0y_i, vis_D0_r);
      vis_D1_r = spu_madd(sample_DxDxDyDy_i, sample_1x1y1x1y_i, vis_D1_r);
      vis_D2_r = spu_madd(sample_DxDxDyDy_i, sample_2x2y2x2y_i, vis_D2_r);
      vis_D3_r = spu_madd(sample_DxDxDyDy_i, sample_3x3y3x3y_i, vis_D3_r);

      vis_A0_i = spu_nmsub(sample_AxAxAyAy_i, sample_0x0y0x0y_r, vis_A0_i);
      vis_A1_i = spu_nmsub(sample_AxAxAyAy_i, sample_1x1y1x1y_r, vis_A1_i);
      vis_A2_i = spu_nmsub(sample_AxAxAyAy_i, sample_2x2y2x2y_r, vis_A2_i);
      vis_A3_i = spu_nmsub(sample_AxAxAyAy_i, sample_3x3y3x3y_r, vis_A3_i);
      vis_B0_i = spu_nmsub(sample_BxBxByBy_i, sample_0x0y0x0y_r, vis_B0_i);
      vis_B1_i = spu_nmsub(sample_BxBxByBy_i, sample_1x1y1x1y_r, vis_B1_i);
      vis_B2_i = spu_nmsub(sample_BxBxByBy_i, sample_2x2y2x2y_r, vis_B2_i);
      vis_B3_i = spu_nmsub(sample_BxBxByBy_i, sample_3x3y3x3y_r, vis_B3_i);
      vis_C0_i = spu_nmsub(sample_CxCxCyCy_i, sample_0x0y0x0y_r, vis_C0_i);
      vis_C1_i = spu_nmsub(sample_CxCxCyCy_i, sample_1x1y1x1y_r, vis_C1_i);
      vis_C2_i = spu_nmsub(sample_CxCxCyCy_i, sample_2x2y2x2y_r, vis_C2_i);
      vis_C3_i = spu_nmsub(sample_CxCxCyCy_i, sample_3x3y3x3y_r, vis_C3_i);
      vis_D0_i = spu_nmsub(sample_DxDxDyDy_i, sample_0x0y0x0y_r, vis_D0_i);
      vis_D1_i = spu_nmsub(sample_DxDxDyDy_i, sample_1x1y1x1y_r, vis_D1_i);
      vis_D2_i = spu_nmsub(sample_DxDxDyDy_i, sample_2x2y2x2y_r, vis_D2_i);
      vis_D3_i = spu_nmsub(sample_DxDxDyDy_i, sample_3x3y3x3y_r, vis_D3_i);

      visibilities[bl + 0 * statA + 0 + 0][REAL] = vis_A0_r;
      visibilities[bl + 0 * statA + 0 + 1][REAL] = vis_A1_r;
      visibilities[bl + 0 * statA + 0 + 2][REAL] = vis_A2_r;
      visibilities[bl + 0 * statA + 0 + 3][REAL] = vis_A3_r;
      visibilities[bl + 1 * statA + 1 + 0][REAL] = vis_B0_r;
      visibilities[bl + 1 * statA + 1 + 1][REAL] = vis_B1_r;
      visibilities[bl + 1 * statA + 1 + 2][REAL] = vis_B2_r;
      visibilities[bl + 1 * statA + 1 + 3][REAL] = vis_B3_r;
      visibilities[bl + 2 * statA + 3 + 0][REAL] = vis_C0_r;
      visibilities[bl + 2 * statA + 3 + 1][REAL] = vis_C1_r;
      visibilities[bl + 2 * statA + 3 + 2][REAL] = vis_C2_r;
      visibilities[bl + 2 * statA + 3 + 3][REAL] = vis_C3_r;
      visibilities[bl + 3 * statA + 6 + 0][REAL] = vis_D0_r;
      visibilities[bl + 3 * statA + 6 + 1][REAL] = vis_D1_r;
      visibilities[bl + 3 * statA + 6 + 2][REAL] = vis_D2_r;
      visibilities[bl + 3 * statA + 6 + 3][REAL] = vis_D3_r;
      visibilities[bl + 0 * statA + 0 + 0][IMAG] = vis_A0_i;
      visibilities[bl + 0 * statA + 0 + 1][IMAG] = vis_A1_i;
      visibilities[bl + 0 * statA + 0 + 2][IMAG] = vis_A2_i;
      visibilities[bl + 0 * statA + 0 + 3][IMAG] = vis_A3_i;
      visibilities[bl + 1 * statA + 1 + 0][IMAG] = vis_B0_i;
      visibilities[bl + 1 * statA + 1 + 1][IMAG] = vis_B1_i;
      visibilities[bl + 1 * statA + 1 + 2][IMAG] = vis_B2_i;
      visibilities[bl + 1 * statA + 1 + 3][IMAG] = vis_B3_i;
      visibilities[bl + 2 * statA + 3 + 0][IMAG] = vis_C0_i;
      visibilities[bl + 2 * statA + 3 + 1][IMAG] = vis_C1_i;
      visibilities[bl + 2 * statA + 3 + 2][IMAG] = vis_C2_i;
      visibilities[bl + 2 * statA + 3 + 3][IMAG] = vis_C3_i;
      visibilities[bl + 3 * statA + 6 + 0][IMAG] = vis_D0_i;
      visibilities[bl + 3 * statA + 6 + 1][IMAG] = vis_D1_i;
      visibilities[bl + 3 * statA + 6 + 2][IMAG] = vis_D2_i;
      visibilities[bl + 3 * statA + 6 + 3][IMAG] = vis_D3_i;
    }
  }
}


static void correlate_all_4x4_diagonals(
  vector float samples[NR_STATIONS][MINOR_INTEGRATION_TIME],
  vector float visibilities[NR_BASELINES][COMPLEX_SIZE]
)
{
  vector float vis_00_r;
  vector float vis_10_r, vis_11_r;
  vector float vis_20_r, vis_21_r, vis_22_r;
  vector float vis_30_r, vis_31_r, vis_32_r, vis_33_r;
  vector float vis_00_i;
  vector float vis_10_i, vis_11_i;
  vector float vis_20_i, vis_21_i, vis_22_i;
  vector float vis_30_i, vis_31_i, vis_32_i, vis_33_i;

  vector float sample_0x0y0x0y_r, sample_1x1y1x1y_r, sample_2x2y2x2y_r, sample_3x3y3x3y_r, sample_0x0x0y0y_r, sample_1x1x1y1y_r, sample_2x2x2y2y_r, sample_3x3x3y3y_r;
  vector float sample_0x0y0x0y_i, sample_1x1y1x1y_i, sample_2x2y2x2y_i, sample_3x3y3x3y_i, sample_0x0x0y0y_i, sample_1x1x1y1y_i, sample_2x2x2y2y_i, sample_3x3x3y3y_i;

  vector float sample_0, sample_1, sample_2, sample_3;

  //for (unsigned _4x4_block = 0; _4x4_block < NR_4X4_1LO2KS; _4x4_block ++) {
    //unsigned stat0 = stat0s[_4x4_block];
    //unsigned stat0 = stat0s[_4x4_block];
    //unsigned bl    = baselines[_4x4_block];
    //printf("stat0 = %u, stat0 = %u, bl = %u\n", stat0, stat0, bl);

  for (unsigned stat0 = 0; stat0 < NR_STATIONS; stat0 += 4) {
    unsigned bl = stat0 * (stat0 + 3) / 2;

    sample_0 = samples[stat0 + 0][0];
    sample_1 = samples[stat0 + 1][0];
    sample_2 = samples[stat0 + 2][0];
    sample_3 = samples[stat0 + 3][0];

    vis_00_r = visibilities[bl + 0 * stat0 + 0][REAL];
    vis_00_i = visibilities[bl + 0 * stat0 + 0][IMAG];
    vis_10_r = visibilities[bl + 1 * stat0 + 1][REAL];
    vis_10_i = visibilities[bl + 1 * stat0 + 1][IMAG];
    vis_11_r = visibilities[bl + 1 * stat0 + 2][REAL];
    vis_11_i = visibilities[bl + 1 * stat0 + 2][IMAG];
    vis_20_r = visibilities[bl + 2 * stat0 + 3][REAL];
    vis_20_i = visibilities[bl + 2 * stat0 + 3][IMAG];
    vis_21_r = visibilities[bl + 2 * stat0 + 4][REAL];
    vis_21_i = visibilities[bl + 2 * stat0 + 4][IMAG];
    vis_22_r = visibilities[bl + 2 * stat0 + 5][REAL];
    vis_22_i = visibilities[bl + 2 * stat0 + 5][IMAG];
    vis_30_r = visibilities[bl + 3 * stat0 + 6][REAL];
    vis_30_i = visibilities[bl + 3 * stat0 + 6][IMAG];
    vis_31_r = visibilities[bl + 3 * stat0 + 7][REAL];
    vis_31_i = visibilities[bl + 3 * stat0 + 7][IMAG];
    vis_32_r = visibilities[bl + 3 * stat0 + 8][REAL];
    vis_32_i = visibilities[bl + 3 * stat0 + 8][IMAG];
    vis_33_r = visibilities[bl + 3 * stat0 + 9][REAL];
    vis_33_i = visibilities[bl + 3 * stat0 + 9][IMAG];

    sample_0x0y0x0y_r = spu_shuffle(sample_0, sample_0, permute_0202);
    sample_1x1y1x1y_r = spu_shuffle(sample_1, sample_1, permute_0202);
    sample_2x2y2x2y_r = spu_shuffle(sample_2, sample_2, permute_0202);
    sample_3x3y3x3y_r = spu_shuffle(sample_3, sample_3, permute_0202);
    sample_0x0x0y0y_r = spu_shuffle(sample_0, sample_0, permute_0022);
    sample_1x1x1y1y_r = spu_shuffle(sample_1, sample_1, permute_0022);
    sample_2x2x2y2y_r = spu_shuffle(sample_2, sample_2, permute_0022);
    sample_3x3x3y3y_r = spu_shuffle(sample_3, sample_3, permute_0022);
    sample_0x0y0x0y_i = spu_shuffle(sample_0, sample_0, permute_1313);
    sample_1x1y1x1y_i = spu_shuffle(sample_1, sample_1, permute_1313);
    sample_2x2y2x2y_i = spu_shuffle(sample_2, sample_2, permute_1313);
    sample_3x3y3x3y_i = spu_shuffle(sample_3, sample_3, permute_1313);
    sample_0x0x0y0y_i = spu_shuffle(sample_0, sample_0, permute_1133);
    sample_1x1x1y1y_i = spu_shuffle(sample_1, sample_1, permute_1133);
    sample_2x2x2y2y_i = spu_shuffle(sample_2, sample_2, permute_1133);
    sample_3x3x3y3y_i = spu_shuffle(sample_3, sample_3, permute_1133);

    vis_00_r = spu_madd(sample_0x0x0y0y_r, sample_0x0y0x0y_r, vis_00_r);
    vis_10_r = spu_madd(sample_1x1x1y1y_r, sample_0x0y0x0y_r, vis_10_r);
    vis_11_r = spu_madd(sample_1x1x1y1y_r, sample_1x1y1x1y_r, vis_11_r);
    vis_20_r = spu_madd(sample_2x2x2y2y_r, sample_0x0y0x0y_r, vis_20_r);
    vis_21_r = spu_madd(sample_2x2x2y2y_r, sample_1x1y1x1y_r, vis_21_r);
    vis_22_r = spu_madd(sample_2x2x2y2y_r, sample_2x2y2x2y_r, vis_22_r);
    vis_30_r = spu_madd(sample_3x3x3y3y_r, sample_0x0y0x0y_r, vis_30_r);
    vis_31_r = spu_madd(sample_3x3x3y3y_r, sample_1x1y1x1y_r, vis_31_r);
    vis_32_r = spu_madd(sample_3x3x3y3y_r, sample_2x2y2x2y_r, vis_32_r);
    vis_33_r = spu_madd(sample_3x3x3y3y_r, sample_3x3y3x3y_r, vis_33_r);

    for (int t = 1; t < MINOR_INTEGRATION_TIME; t += 1) {
      sample_0 = samples[stat0 + 0][t];
      sample_1 = samples[stat0 + 1][t];
      sample_2 = samples[stat0 + 2][t];
      sample_3 = samples[stat0 + 3][t];

      vis_00_i = spu_madd(sample_0x0x0y0y_r, sample_0x0y0x0y_i, vis_00_i);
      vis_10_i = spu_madd(sample_1x1x1y1y_r, sample_0x0y0x0y_i, vis_10_i);
      vis_11_i = spu_madd(sample_1x1x1y1y_r, sample_1x1y1x1y_i, vis_11_i);
      vis_20_i = spu_madd(sample_2x2x2y2y_r, sample_0x0y0x0y_i, vis_20_i);
      vis_21_i = spu_madd(sample_2x2x2y2y_r, sample_1x1y1x1y_i, vis_21_i);
      vis_22_i = spu_madd(sample_2x2x2y2y_r, sample_2x2y2x2y_i, vis_22_i);
      vis_30_i = spu_madd(sample_3x3x3y3y_r, sample_0x0y0x0y_i, vis_30_i);
      vis_31_i = spu_madd(sample_3x3x3y3y_r, sample_1x1y1x1y_i, vis_31_i);
      vis_32_i = spu_madd(sample_3x3x3y3y_r, sample_2x2y2x2y_i, vis_32_i);
      vis_33_i = spu_madd(sample_3x3x3y3y_r, sample_3x3y3x3y_i, vis_33_i);

      vis_00_r = spu_madd(sample_0x0x0y0y_i, sample_0x0y0x0y_i, vis_00_r);
      vis_10_r = spu_madd(sample_1x1x1y1y_i, sample_0x0y0x0y_i, vis_10_r);
      vis_11_r = spu_madd(sample_1x1x1y1y_i, sample_1x1y1x1y_i, vis_11_r);
      vis_20_r = spu_madd(sample_2x2x2y2y_i, sample_0x0y0x0y_i, vis_20_r);
      vis_21_r = spu_madd(sample_2x2x2y2y_i, sample_1x1y1x1y_i, vis_21_r);
      vis_22_r = spu_madd(sample_2x2x2y2y_i, sample_2x2y2x2y_i, vis_22_r);
      vis_30_r = spu_madd(sample_3x3x3y3y_i, sample_0x0y0x0y_i, vis_30_r);
      vis_31_r = spu_madd(sample_3x3x3y3y_i, sample_1x1y1x1y_i, vis_31_r);
      vis_32_r = spu_madd(sample_3x3x3y3y_i, sample_2x2y2x2y_i, vis_32_r);
      vis_33_r = spu_madd(sample_3x3x3y3y_i, sample_3x3y3x3y_i, vis_33_r);

      vis_00_i = spu_nmsub(sample_0x0x0y0y_i, sample_0x0y0x0y_r, vis_00_i);
      vis_10_i = spu_nmsub(sample_1x1x1y1y_i, sample_0x0y0x0y_r, vis_10_i);
      vis_11_i = spu_nmsub(sample_1x1x1y1y_i, sample_1x1y1x1y_r, vis_11_i);
      vis_20_i = spu_nmsub(sample_2x2x2y2y_i, sample_0x0y0x0y_r, vis_20_i);
      vis_21_i = spu_nmsub(sample_2x2x2y2y_i, sample_1x1y1x1y_r, vis_21_i);
      vis_22_i = spu_nmsub(sample_2x2x2y2y_i, sample_2x2y2x2y_r, vis_22_i);
      vis_30_i = spu_nmsub(sample_3x3x3y3y_i, sample_0x0y0x0y_r, vis_30_i);
      vis_31_i = spu_nmsub(sample_3x3x3y3y_i, sample_1x1y1x1y_r, vis_31_i);
      vis_32_i = spu_nmsub(sample_3x3x3y3y_i, sample_2x2y2x2y_r, vis_32_i);
      vis_33_i = spu_nmsub(sample_3x3x3y3y_i, sample_3x3y3x3y_r, vis_33_i);

      sample_0x0y0x0y_r = spu_shuffle(sample_0, sample_0, permute_0202);
      sample_1x1y1x1y_r = spu_shuffle(sample_1, sample_1, permute_0202);
      sample_2x2y2x2y_r = spu_shuffle(sample_2, sample_2, permute_0202);
      sample_3x3y3x3y_r = spu_shuffle(sample_3, sample_3, permute_0202);
      sample_0x0x0y0y_r = spu_shuffle(sample_0, sample_0, permute_0022);
      sample_1x1x1y1y_r = spu_shuffle(sample_1, sample_1, permute_0022);
      sample_2x2x2y2y_r = spu_shuffle(sample_2, sample_2, permute_0022);
      sample_3x3x3y3y_r = spu_shuffle(sample_3, sample_3, permute_0022);
      sample_0x0y0x0y_i = spu_shuffle(sample_0, sample_0, permute_1313);
      sample_1x1y1x1y_i = spu_shuffle(sample_1, sample_1, permute_1313);
      sample_2x2y2x2y_i = spu_shuffle(sample_2, sample_2, permute_1313);
      sample_3x3y3x3y_i = spu_shuffle(sample_3, sample_3, permute_1313);
      sample_0x0x0y0y_i = spu_shuffle(sample_0, sample_0, permute_1133);
      sample_1x1x1y1y_i = spu_shuffle(sample_1, sample_1, permute_1133);
      sample_2x2x2y2y_i = spu_shuffle(sample_2, sample_2, permute_1133);
      sample_3x3x3y3y_i = spu_shuffle(sample_3, sample_3, permute_1133);

      vis_00_r = spu_madd(sample_0x0x0y0y_r, sample_0x0y0x0y_r, vis_00_r);
      vis_10_r = spu_madd(sample_1x1x1y1y_r, sample_0x0y0x0y_r, vis_10_r);
      vis_11_r = spu_madd(sample_1x1x1y1y_r, sample_1x1y1x1y_r, vis_11_r);
      vis_20_r = spu_madd(sample_2x2x2y2y_r, sample_0x0y0x0y_r, vis_20_r);
      vis_21_r = spu_madd(sample_2x2x2y2y_r, sample_1x1y1x1y_r, vis_21_r);
      vis_22_r = spu_madd(sample_2x2x2y2y_r, sample_2x2y2x2y_r, vis_22_r);
      vis_30_r = spu_madd(sample_3x3x3y3y_r, sample_0x0y0x0y_r, vis_30_r);
      vis_31_r = spu_madd(sample_3x3x3y3y_r, sample_1x1y1x1y_r, vis_31_r);
      vis_32_r = spu_madd(sample_3x3x3y3y_r, sample_2x2y2x2y_r, vis_32_r);
      vis_33_r = spu_madd(sample_3x3x3y3y_r, sample_3x3y3x3y_r, vis_33_r);
    }

    vis_00_i = spu_madd(sample_0x0x0y0y_r, sample_0x0y0x0y_i, vis_00_i);
    vis_10_i = spu_madd(sample_1x1x1y1y_r, sample_0x0y0x0y_i, vis_10_i);
    vis_11_i = spu_madd(sample_1x1x1y1y_r, sample_1x1y1x1y_i, vis_11_i);
    vis_20_i = spu_madd(sample_2x2x2y2y_r, sample_0x0y0x0y_i, vis_20_i);
    vis_21_i = spu_madd(sample_2x2x2y2y_r, sample_1x1y1x1y_i, vis_21_i);
    vis_22_i = spu_madd(sample_2x2x2y2y_r, sample_2x2y2x2y_i, vis_22_i);
    vis_30_i = spu_madd(sample_3x3x3y3y_r, sample_0x0y0x0y_i, vis_30_i);
    vis_31_i = spu_madd(sample_3x3x3y3y_r, sample_1x1y1x1y_i, vis_31_i);
    vis_32_i = spu_madd(sample_3x3x3y3y_r, sample_2x2y2x2y_i, vis_32_i);
    vis_33_i = spu_madd(sample_3x3x3y3y_r, sample_3x3y3x3y_i, vis_33_i);

    vis_00_r = spu_madd(sample_0x0x0y0y_i, sample_0x0y0x0y_i, vis_00_r);
    vis_10_r = spu_madd(sample_1x1x1y1y_i, sample_0x0y0x0y_i, vis_10_r);
    vis_11_r = spu_madd(sample_1x1x1y1y_i, sample_1x1y1x1y_i, vis_11_r);
    vis_20_r = spu_madd(sample_2x2x2y2y_i, sample_0x0y0x0y_i, vis_20_r);
    vis_21_r = spu_madd(sample_2x2x2y2y_i, sample_1x1y1x1y_i, vis_21_r);
    vis_22_r = spu_madd(sample_2x2x2y2y_i, sample_2x2y2x2y_i, vis_22_r);
    vis_30_r = spu_madd(sample_3x3x3y3y_i, sample_0x0y0x0y_i, vis_30_r);
    vis_31_r = spu_madd(sample_3x3x3y3y_i, sample_1x1y1x1y_i, vis_31_r);
    vis_32_r = spu_madd(sample_3x3x3y3y_i, sample_2x2y2x2y_i, vis_32_r);
    vis_33_r = spu_madd(sample_3x3x3y3y_i, sample_3x3y3x3y_i, vis_33_r);

    vis_00_i = spu_nmsub(sample_0x0x0y0y_i, sample_0x0y0x0y_r, vis_00_i);
    vis_10_i = spu_nmsub(sample_1x1x1y1y_i, sample_0x0y0x0y_r, vis_10_i);
    vis_11_i = spu_nmsub(sample_1x1x1y1y_i, sample_1x1y1x1y_r, vis_11_i);
    vis_20_i = spu_nmsub(sample_2x2x2y2y_i, sample_0x0y0x0y_r, vis_20_i);
    vis_21_i = spu_nmsub(sample_2x2x2y2y_i, sample_1x1y1x1y_r, vis_21_i);
    vis_22_i = spu_nmsub(sample_2x2x2y2y_i, sample_2x2y2x2y_r, vis_22_i);
    vis_30_i = spu_nmsub(sample_3x3x3y3y_i, sample_0x0y0x0y_r, vis_30_i);
    vis_31_i = spu_nmsub(sample_3x3x3y3y_i, sample_1x1y1x1y_r, vis_31_i);
    vis_32_i = spu_nmsub(sample_3x3x3y3y_i, sample_2x2y2x2y_r, vis_32_i);
    vis_33_i = spu_nmsub(sample_3x3x3y3y_i, sample_3x3y3x3y_r, vis_33_i);

    visibilities[bl + 0 * stat0 + 0][REAL] = vis_00_r;
    visibilities[bl + 1 * stat0 + 1][REAL] = vis_10_r;
    visibilities[bl + 1 * stat0 + 2][REAL] = vis_11_r;
    visibilities[bl + 2 * stat0 + 3][REAL] = vis_20_r;
    visibilities[bl + 2 * stat0 + 4][REAL] = vis_21_r;
    visibilities[bl + 2 * stat0 + 5][REAL] = vis_22_r;
    visibilities[bl + 3 * stat0 + 6][REAL] = vis_30_r;
    visibilities[bl + 3 * stat0 + 7][REAL] = vis_31_r;
    visibilities[bl + 3 * stat0 + 8][REAL] = vis_32_r;
    visibilities[bl + 3 * stat0 + 9][REAL] = vis_33_r;
    visibilities[bl + 0 * stat0 + 0][IMAG] = vis_00_i;
    visibilities[bl + 1 * stat0 + 1][IMAG] = vis_10_i;
    visibilities[bl + 1 * stat0 + 2][IMAG] = vis_11_i;
    visibilities[bl + 2 * stat0 + 3][IMAG] = vis_20_i;
    visibilities[bl + 2 * stat0 + 4][IMAG] = vis_21_i;
    visibilities[bl + 2 * stat0 + 5][IMAG] = vis_22_i;
    visibilities[bl + 3 * stat0 + 6][IMAG] = vis_30_i;
    visibilities[bl + 3 * stat0 + 7][IMAG] = vis_31_i;
    visibilities[bl + 3 * stat0 + 8][IMAG] = vis_32_i;
    visibilities[bl + 3 * stat0 + 9][IMAG] = vis_33_i;
  }
}


static void init(unsigned long long argp)
{
  mfc_get(&spu_arguments, (unsigned) argp, sizeof spu_arguments, 0, 0, 0);
  mfc_write_tag_mask(1 << 0);
  mfc_read_tag_status_all();

  first_channel =  spu_arguments.spu_id      * NR_CHANNELS / NR_SPUS;
  last_channel  = (spu_arguments.spu_id + 1) * NR_CHANNELS / NR_SPUS;

#if defined DO_DMA
  for (unsigned channel = first_channel; channel < NR_CHANNELS; channel ++) {
    unsigned i;

    for (i = 0; i < sizeof visibilities / 16384; i ++) {
      visibilities_dma_list[channel - first_channel][i].ea_low = (unsigned) (*((ppu_visibilities_type *) spu_arguments.ppu_visibilities))[channel] + 16384 * i;
      visibilities_dma_list[channel - first_channel][i].size   = 16384;
    }

    if (sizeof visibilities % 16384 != 0) {
      visibilities_dma_list[channel - first_channel][i].ea_low = (unsigned) (*((ppu_visibilities_type *) spu_arguments.ppu_visibilities))[channel] + 16384 * i;
      visibilities_dma_list[channel - first_channel][i].size   = sizeof visibilities % 16384;
    }
  }

  unsigned i;

  for (i = 0; i < sizeof visibilities / 16384; i ++) {
    zero_dma_list[i].ea_low = (unsigned) (*((ppu_zeros_type *) spu_arguments.ppu_zeros));
    zero_dma_list[i].size   = 16384;
  }

  if (sizeof visibilities % 16384 != 0) {
    zero_dma_list[i].ea_low = (unsigned) (*((ppu_zeros_type *) spu_arguments.ppu_zeros));
    zero_dma_list[i].size   = sizeof visibilities % 16384;
  }
#endif

  printf("I am spu %d, calculating channels %3d - %3d\n", spu_arguments.spu_id, first_channel, last_channel);

#if 0
  for (unsigned statA = 4, _4x4_block = 0; statA < NR_STATIONS; statA += 4) {
    for (unsigned stat0 = 0; stat0 < statA; stat0 += 4, ++ _4x4_block) {
      stat0s[_4x4_block]    = stat0;
      statAs[_4x4_block]    = statA;
      baselines[_4x4_block] = statA * (statA + 1) / 2 + stat0;
    }
  }
#endif
}


static inline void memzero(vector float *ptr, unsigned size)
{
  for (unsigned i = 0; i < size / 16; i += 8) {
    ptr[i + 0] = (vector float) { 0, 0, 0, 0 };
    ptr[i + 1] = (vector float) { 0, 0, 0, 0 };
    ptr[i + 2] = (vector float) { 0, 0, 0, 0 };
    ptr[i + 3] = (vector float) { 0, 0, 0, 0 };
    ptr[i + 4] = (vector float) { 0, 0, 0, 0 };
    ptr[i + 5] = (vector float) { 0, 0, 0, 0 };
    ptr[i + 6] = (vector float) { 0, 0, 0, 0 };
    ptr[i + 7] = (vector float) { 0, 0, 0, 0 };
  }
}


int main(unsigned long long id, unsigned long long argp, unsigned long long envp)
{
  id = id; // fix compiler warning
  envp = envp; // fix compiler warning

  init(argp);
  decrementer_init();

  //  barrier();

  unsigned samples_buffer = 0;

  unsigned int start = decrementer_get();

  for (unsigned i = 0; i < 10; i ++) {
#if defined DO_DMA
    load_samples(first_channel, 0, samples_buffer);
#endif

    for (unsigned channel = first_channel; likely(channel < last_channel); channel ++) {
      //memzero(visibilities, sizeof visibilities);

      for (unsigned major_time = 0; major_time < MAJOR_INTEGRATION_TIME - MINOR_INTEGRATION_TIME; major_time += MINOR_INTEGRATION_TIME) {
#if defined DO_DMA
	load_samples(channel, major_time + MINOR_INTEGRATION_TIME, samples_buffer ^ 1);
	wait_for_dma_samples(samples_buffer);
#endif

	correlate_all_4x4_diagonals(samples[samples_buffer], visibilities);
	correlate_all_4x4_blocks(samples[samples_buffer], visibilities);
	samples_buffer ^= 1;
      }

#if defined DO_DMA
      if (likely(channel != last_channel - 1))
	load_samples(channel + 1, 0, samples_buffer ^ 1);

      wait_for_dma_samples(samples_buffer);
#endif

      correlate_all_4x4_diagonals(samples[samples_buffer], visibilities);
      correlate_all_4x4_blocks(samples[samples_buffer], visibilities);
      samples_buffer ^= 1;

#if defined DO_DMA
      store_visibilities(channel);
      wait_for_dma_visibilities();
#endif
    }
  }

  unsigned int end = decrementer_get();

  //  barrier();

  double seconds = decrementer_seconds(start - end);

  printf("SPU %u: time = %.12f s\n", spu_arguments.spu_id, seconds);

  return 0;
}

