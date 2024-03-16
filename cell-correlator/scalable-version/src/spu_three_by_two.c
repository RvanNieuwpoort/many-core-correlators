#include "correlator.h"
#include <spu_intrinsics.h>

static vector unsigned char global_real_pattern = (vector unsigned char) {0, 1, 2, 3, 0, 1, 2, 3, 8, 9, 10, 11, 8, 9, 10, 11}; // float 0, 0, 2, 2. i.e., real values
static vector unsigned char global_img_pattern = (vector unsigned char) {4, 5, 6, 7, 4, 5, 6, 7, 12, 13, 14, 15, 12, 13, 14, 15 }; // float 1, 1, 3, 3. i.e., img values

static vector unsigned char global_real_pattern_2 = (vector unsigned char) {0, 1, 2, 3, 8, 9, 10, 11, 0, 1, 2, 3, 8, 9, 10, 11}; // float 0, 2, 0, 2. i.e., real values
static vector unsigned char global_img_pattern_2 =  (vector unsigned char) {4, 5, 6, 7, 12, 13, 14, 15, 4, 5, 6, 7, 12, 13, 14, 15 }; // float 1, 3, 1, 3. i.e., img values

#define innerLoop(t) \
  s0_r = spu_shuffle(samples0[t], samples0[t], real_pattern);  \
  s0_i = spu_shuffle(samples0[t], samples0[t], img_pattern);   \
  s1_r = spu_shuffle(samples1[t], samples1[t], real_pattern);  \
  s1_i = spu_shuffle(samples1[t], samples1[t], img_pattern);   \
  s2_r = spu_shuffle(samples2[t], samples2[t], real_pattern); \
  s2_i = spu_shuffle(samples2[t], samples2[t], img_pattern);	  \
  s3_r = spu_shuffle(samples3[t], samples3[t], real_pattern_2);	  \
  s3_i = spu_shuffle(samples3[t], samples3[t], img_pattern_2);	  \
  s4_r = spu_shuffle(samples4[t], samples4[t], real_pattern_2);	  \
  s4_i = spu_shuffle(samples4[t], samples4[t], img_pattern_2);	  \
 \
  vr03 = spu_madd(s0_i, s3_i, spu_madd (s0_r, s3_r, vr03));   \
  vi03 = spu_madd(s0_i, s3_r, spu_nmsub(s0_r, s3_i, vi03));   \
  vr13 = spu_madd(s1_i, s3_i, spu_madd (s1_r, s3_r, vr13));   \
  vi13 = spu_madd(s1_i, s3_r, spu_nmsub(s1_r, s3_i, vi13));   \
  vr23 = spu_madd(s2_i, s3_i, spu_madd (s2_r, s3_r, vr23));   \
  vi23 = spu_madd(s2_i, s3_r, spu_nmsub(s2_r, s3_i, vi23));   \
  vr04 = spu_madd(s0_i, s4_i, spu_madd (s0_r, s4_r, vr04));   \
  vi04 = spu_madd(s0_i, s4_r, spu_nmsub(s0_r, s4_i, vi04));   \
  vr14 = spu_madd(s1_i, s4_i, spu_madd (s1_r, s4_r, vr14));   \
  vi14 = spu_madd(s1_i, s4_r, spu_nmsub(s1_r, s4_i, vi14));   \
  vr24 = spu_madd(s2_i, s4_i, spu_madd (s2_r, s4_r, vr24));   \
  vi24 = spu_madd(s2_i, s4_r, spu_nmsub(s2_r, s4_i, vi24));   \


void threeByTwo(
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
              vector float* vis24) {

  vector unsigned char real_pattern = global_real_pattern;
  vector unsigned char img_pattern =  global_img_pattern;
  vector unsigned char real_pattern_2 = global_real_pattern_2;
  vector unsigned char img_pattern_2 =  global_img_pattern_2;

  vector float vr03 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi03 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr13 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi13 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr23 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi23 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr04 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi04 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr14 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi14 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr24 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi24 = {0.0f, 0.0f, 0.0f, 0.0f};

  vector float s0_r;
  vector float s0_i;
  vector float s1_r;
  vector float s1_i;
  vector float s2_r;
  vector float s2_i;
  vector float s3_r;
  vector float s3_i;
  vector float s4_r;
  vector float s4_i;

  for(int t=0; t<64; t+=16) {
    innerLoop(t+0);
    innerLoop(t+1);
    innerLoop(t+2);
    innerLoop(t+3);
    innerLoop(t+4);
    innerLoop(t+5);
    innerLoop(t+6);
    innerLoop(t+7);
    innerLoop(t+8);
    innerLoop(t+9);
    innerLoop(t+10);
    innerLoop(t+11);
    innerLoop(t+12);
    innerLoop(t+13);
    innerLoop(t+14);
    innerLoop(t+15);
  }
// the compiler generates horrible code if we inline everything :-(
#if 0
  innerLoop(0);
  innerLoop(1);
  innerLoop(2);
  innerLoop(3);
  innerLoop(4);
  innerLoop(5);
  innerLoop(6);
  innerLoop(7);
  innerLoop(8);
  innerLoop(9);
  innerLoop(10);
  innerLoop(11);
  innerLoop(12);
  innerLoop(13);
  innerLoop(14);
  innerLoop(15);
  innerLoop(16);
  innerLoop(17);
  innerLoop(18);
  innerLoop(19);
  innerLoop(20);
  innerLoop(21);
  innerLoop(22);
  innerLoop(23);
  innerLoop(24);
  innerLoop(25);
  innerLoop(26);
  innerLoop(27);
  innerLoop(28);
  innerLoop(29);
  innerLoop(30);
  innerLoop(31);
  innerLoop(32);
  innerLoop(33);
  innerLoop(34);
  innerLoop(35);
  innerLoop(36);
  innerLoop(37);
  innerLoop(38);
  innerLoop(39);
  innerLoop(40);
  innerLoop(41);
  innerLoop(42);
  innerLoop(43);
  innerLoop(44);
  innerLoop(45);
  innerLoop(46);
  innerLoop(47);
  innerLoop(48);
  innerLoop(49);
  innerLoop(50);
  innerLoop(51);
  innerLoop(52);
  innerLoop(53);
  innerLoop(54);
  innerLoop(55);
  innerLoop(56);
  innerLoop(57);
  innerLoop(58);
  innerLoop(59);
  innerLoop(60);
  innerLoop(61);
  innerLoop(62);
  innerLoop(63);
#endif
#if 1
  *(vis03+0) = vr03;
  *(vis03+1) = vi03;
  *(vis13+0) = vr13;
  *(vis13+1) = vi13;
  *(vis23+0) = vr23;
  *(vis23+1) = vi23;
  *(vis04+0) = vr04;
  *(vis04+1) = vi04;
  *(vis14+0) = vr14;
  *(vis14+1) = vi14;
  *(vis24+0) = vr24;
  *(vis24+1) = vi24;
#endif
}

