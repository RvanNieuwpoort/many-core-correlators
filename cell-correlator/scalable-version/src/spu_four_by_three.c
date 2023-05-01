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
  s3_r = spu_shuffle(samples3[t], samples3[t], real_pattern); \
  s3_i = spu_shuffle(samples3[t], samples3[t], img_pattern);	  \
  s4_r = spu_shuffle(samples4[t], samples4[t], real_pattern_2);	  \
  s4_i = spu_shuffle(samples4[t], samples4[t], img_pattern_2);	  \
  s5_r = spu_shuffle(samples5[t], samples5[t], real_pattern_2);	  \
  s5_i = spu_shuffle(samples5[t], samples5[t], img_pattern_2);	  \
  s6_r = spu_shuffle(samples6[t], samples6[t], real_pattern_2);	  \
  s6_i = spu_shuffle(samples6[t], samples6[t], img_pattern_2);	  \
 \
  vr04 = spu_madd(s0_i, s4_i, spu_madd (s0_r, s4_r, vr04));   \
  vi04 = spu_madd(s0_i, s4_r, spu_nmsub(s0_r, s4_i, vi04));   \
  vr14 = spu_madd(s1_i, s4_i, spu_madd (s1_r, s4_r, vr14));   \
  vi14 = spu_madd(s1_i, s4_r, spu_nmsub(s1_r, s4_i, vi14));   \
  vr24 = spu_madd(s2_i, s4_i, spu_madd (s2_r, s4_r, vr24));   \
  vi24 = spu_madd(s2_i, s4_r, spu_nmsub(s2_r, s4_i, vi24));   \
  vr34 = spu_madd(s3_i, s4_i, spu_madd (s2_r, s4_r, vr24));   \
  vi34 = spu_madd(s3_i, s4_r, spu_nmsub(s2_r, s4_i, vi24));   \
  vr05 = spu_madd(s0_i, s5_i, spu_madd (s0_r, s5_r, vr05));   \
  vi05 = spu_madd(s0_i, s5_r, spu_nmsub(s0_r, s5_i, vi05));   \
  vr15 = spu_madd(s1_i, s5_i, spu_madd (s1_r, s5_r, vr15));   \
  vi15 = spu_madd(s1_i, s5_r, spu_nmsub(s1_r, s5_i, vi15));   \
  vr25 = spu_madd(s2_i, s5_i, spu_madd (s2_r, s5_r, vr25));   \
  vi25 = spu_madd(s2_i, s5_r, spu_nmsub(s2_r, s5_i, vi25));   \
  vr35 = spu_madd(s3_i, s5_i, spu_madd (s3_r, s5_r, vr35));   \
  vi35 = spu_madd(s3_i, s5_r, spu_nmsub(s3_r, s5_i, vi35));   \
  vr06 = spu_madd(s0_i, s6_i, spu_madd (s0_r, s6_r, vr06));   \
  vi06 = spu_madd(s0_i, s6_r, spu_nmsub(s0_r, s6_i, vi06));   \
  vr16 = spu_madd(s1_i, s6_i, spu_madd (s1_r, s6_r, vr16));   \
  vi16 = spu_madd(s1_i, s6_r, spu_nmsub(s1_r, s6_i, vi16));   \
  vr26 = spu_madd(s2_i, s6_i, spu_madd (s2_r, s6_r, vr26));   \
  vi26 = spu_madd(s2_i, s6_r, spu_nmsub(s2_r, s6_i, vi26));   \
  vr36 = spu_madd(s3_i, s6_i, spu_madd (s3_r, s6_r, vr36));   \
  vi36 = spu_madd(s3_i, s6_r, spu_nmsub(s3_r, s6_i, vi36));   \


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
              vector float* vis36) {

  vector unsigned char real_pattern = global_real_pattern;
  vector unsigned char img_pattern =  global_img_pattern;
  vector unsigned char real_pattern_2 = global_real_pattern_2;
  vector unsigned char img_pattern_2 =  global_img_pattern_2;

  vector float vr04 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi04 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr14 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi14 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr24 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi24 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr34 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi34 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr05 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi05 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr15 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi15 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr25 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi25 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr35 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi35 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr06 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi06 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr16 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi16 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr26 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi26 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vr36 = {0.0f, 0.0f, 0.0f, 0.0f};
  vector float vi36 = {0.0f, 0.0f, 0.0f, 0.0f};

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
  vector float s5_r;
  vector float s5_i;
  vector float s6_r;
  vector float s6_i;

  for(int t=0; t<64; t+=8) {
    innerLoop(t+0);
    innerLoop(t+1);
    innerLoop(t+2);
    innerLoop(t+3);
    innerLoop(t+4);
    innerLoop(t+5);
    innerLoop(t+6);
    innerLoop(t+7);
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
  *(vis04+0) = vr04;
  *(vis04+1) = vi04;
  *(vis14+0) = vr14;
  *(vis14+1) = vi14;
  *(vis24+0) = vr24;
  *(vis24+1) = vi24;
  *(vis34+0) = vr34;
  *(vis34+1) = vi34;
  *(vis05+0) = vr05;
  *(vis05+1) = vi05;
  *(vis15+0) = vr15;
  *(vis15+1) = vi15;
  *(vis25+0) = vr25;
  *(vis25+1) = vi25;
  *(vis35+0) = vr35;
  *(vis35+1) = vi35;
  *(vis06+0) = vr06;
  *(vis06+1) = vi06;
  *(vis16+0) = vr16;
  *(vis16+1) = vi16;
  *(vis26+0) = vr26;
  *(vis26+1) = vi26;
  *(vis36+0) = vr36;
  *(vis36+1) = vi36;
#endif
}

