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
#ifndef _CORRELATOR_3x3_MEM_H_
#define _CORRELATOR_3x3_MEM_H_

const CALchar* Kernel_3x3_mem =
"il_ps_2_0\n"
"dcl_input_position_interp(linear_noperspective) vWinCoord0.xy__\n" // position in domain (output buffer)
"dcl_resource_id(0)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // samples
"dcl_resource_id(1)_type(1d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // cellToStatX
"dcl_resource_id(2)_type(1d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // cellToStatY
"dcl_literal l0, 0x00000000, 0x00000000, 0x00000000, 0x00000000\n"
"dcl_literal l1, 1.0, 2.0, 3.0, 4.0\n"
"dcl_literal l18, 18.0, 18.0, 18.0, 18.0\n"
"dcl_cb cb0[1]\n"

// r0 -- r9 and r30 -- r39 are the computed visibilities
"mov r0, l0\n"  // reals03
"mov r1, l0\n"  // imags03
"mov r2, l0\n"  // reals13
"mov r3, l0\n"  // imags13
"mov r4, l0\n"  // reals23
"mov r5, l0\n"  // imags23
"mov r6, l0\n"  // reals04
"mov r7, l0\n"  // imags04
"mov r30, l0\n" // reals14
"mov r31, l0\n" // imags14
"mov r32, l0\n" // reals24
"mov r33, l0\n" // imags24
"mov r34, l0\n" // reals05
"mov r35, l0\n" // imags05
"mov r36, l0\n" // reals15
"mov r37, l0\n" // imags15
"mov r38, l0\n" // reals25
"mov r39, l0\n" // imags25

                             // cb0[0].z = nrChannels
                             // cb0[0].w = nrCells
"mov r8.x, vWinCoord0.x\n"   // r8.x = cell
"mov r8.y, vWinCoord0.y\n"   // r8.y = channel
"mov r9.x, cb0[0].x\n"       // r9.x = nrStations (float)
"mov r9.y, cb0[0].y\n"       // r9.y = nrTimes (float)
"mov r9.zw, l0\n"            // r9.z = time loopcounter (0 .. nrTimes)

"ftoi r8.y, r8.y"
"itof r8.y, r8.y\n" // r8.y = real channel value (vWinCoord is in the middle of two samples, e.g. 1.5)
"ftoi r8.x, r8.x"
"itof r8.x, r8.x\n" // r8.x = real cell value

"sample_resource(1)_sampler(0)   r8.z, r8.x\n" // r8.z = stat0 (X)
"sample_resource(2)_sampler(1)   r8.w, r8.x\n" // r8.w = stat2 (Y)


// indices should be (channel * nrStations) + station, time, 0, 0
"mul r10.y, r8.y, r9.x\n"
"add r10.y, r10.y, r8.z\n"
"mov r10.xzw, l0\n" // r10 = index0: {time, (channel * nrStations) + station0+0, 0, 0}

"add r11.y, r10.y, r11.1\n"
"mov r11.xzw, l0\n" // r11 = index0: {time, (channel * nrStations) + station0+1, 0, 0}

"add r12.y, r11.y, r12.1\n"
"mov r12.xzw, l0\n" // r12 = index0: {time, (channel * nrStations) + station0+2, 0, 0}

"mul r13.y, r8.y, r9.x\n"
"add r13.y, r13.y, r8.w\n"
"mov r13.xzw, l0\n" // r13 = index0: {time, (channel * nrStations) + station1+0, 0, 0}

"add r14.y, r13.y, r14.1\n"
"mov r14.xzw, l0\n" // r14 = index0: {time, (channel * nrStations) + station1+1, 0, 0}

"add r15.y, r14.y, r15.1\n"
"mov r15.xzw, l0\n" // r15 = index0: {time, (channel * nrStations) + station1+2, 0, 0}


"whileloop\n"
"    ge r9.w, r9.z, r9.y\n"  // r9.w = temp used for comparison
"    break_logicalnz r9.w\n"

//   update time index
"    mov r10.x, r9.z\n"
"    mov r11.x, r9.z\n"
"    mov r12.x, r9.z\n"
"    mov r13.x, r9.z\n"
"    mov r14.x, r9.z\n"
"    mov r15.x, r9.z\n"

//   load the six samples
#if DO_LOADS
"    sample_resource(0)_sampler(0) r20, r10\n" // sample 0
"    sample_resource(0)_sampler(1) r21, r11\n" // sample 1
"    sample_resource(0)_sampler(2) r22, r12\n" // sample 2
"    sample_resource(0)_sampler(3) r23, r13\n" // sample 3
"    sample_resource(0)_sampler(4) r24, r14\n" // sample 4
"    sample_resource(0)_sampler(5) r25, r15\n" // sample 5
#else
"    mov r20, r10\n"
"    mov r21, r11\n"
"    mov r22, r12\n"
"    mov r23, r13\n"
"    mov r24, r14\n"
"    mov r25, r15\n"
#endif

// calculate:
// 0 x 3 -> r0, r1
// 1 x 3 -> r2, r3
// 2 x 3 -> r4, r5

// 0 x 4 -> r6, r7
// 1 x 4 -> r30, r31
// 2 x 4 -> r32, r33

// 0 x 5 -> r34, r35
// 1 x 5 -> r36, r37
// 2 x 5 -> r38, r39

#if DO_COMPUTATION
"    mad r0,  r20.xxzz, r23.xzxz, r0\n" //    reals03 += sample0.xxzz * sample3.xzxz
"    mad r0,  r20.yyww, r23.ywyw, r0\n" //    reals03 += sample0.yyww * sample3.ywyw
"    mad r1,  r20.yyww, r23.xzxz, r1\n" //    imags03 += sample0.yyww * sample3.xzxz
"    mul r50, r20.xxzz, r23.ywyw\n"     //    imags03 -= sample0.xxzz * sample3.ywyw
"    sub r1,  r1,       r50\n"

"    mad r2,  r21.xxzz, r23.xzxz, r2\n" //    reals13 += sample1.xxzz * sample3.xzxz
"    mad r2,  r21.yyww, r23.ywyw, r2\n" //    reals13 += sample1.yyww * sample3.ywyw
"    mad r3,  r21.yyww, r23.xzxz, r3\n" //    imags13 += sample1.yyww * sample3.xzxz
"    mul r50, r21.xxzz, r23.ywyw\n"     //    imags13 -= sample1.xxzz * sample3.ywyw
"    sub r3,  r3,       r50\n"

"    mad r4,  r22.xxzz, r23.xzxz, r4\n" //    reals23 += sample2.xxzz * sample3.xzxz
"    mad r4,  r22.yyww, r23.ywyw, r4\n" //    reals23 += sample2.yyww * sample3.ywyw
"    mad r5,  r22.yyww, r23.xzxz, r5\n" //    imags23 += sample2.yyww * sample3.xzxz
"    mul r50, r22.xxzz, r23.ywyw\n"     //    imags23 -= sample2.xxzz * sample3.ywyw
"    sub r5,  r5,       r50\n"


"    mad r6,  r20.xxzz, r24.xzxz, r6\n"   //    reals04 += sample0.xxzz * sample4.xzxz
"    mad r6,  r20.yyww, r24.ywyw, r6\n"   //    reals04 += sample0.yyww * sample4.ywyw
"    mad r7,  r20.yyww, r24.xzxz, r7\n"   //    imags04 += sample0.yyww * sample4.xzxz
"    mul r50, r20.xxzz, r24.ywyw\n"       //    imags04 -= sample0.xxzz * sample4.ywyw
"    sub r7,  r7,       r50\n"

"    mad r30,  r21.xxzz, r24.xzxz, r30\n" //    reals14 += sample1.xxzz * sample4.xzxz
"    mad r30,  r21.yyww, r24.ywyw, r30\n" //    reals14 += sample1.yyww * sample4.ywyw
"    mad r31,  r21.yyww, r24.xzxz, r31\n" //    imags14 += sample1.yyww * sample4.xzxz
"    mul r50,  r21.xxzz, r24.ywyw\n"      //    imags14 -= sample1.xxzz * sample4.ywyw
"    sub r31,  r31,      r50\n"

"    mad r32,  r22.xxzz, r24.xzxz, r32\n" //    reals24 += sample2.xxzz * sample4.xzxz
"    mad r32,  r22.yyww, r24.ywyw, r32\n" //    reals24 += sample2.yyww * sample4.ywyw
"    mad r33,  r22.yyww, r24.xzxz, r33\n" //    imags24 += sample2.yyww * sample4.xzxz
"    mul r50,  r22.xxzz, r24.ywyw\n"      //    imags24 -= sample2.xxzz * sample4.ywyw
"    sub r33,  r33,      r50\n"


"    mad r34,  r20.xxzz, r25.xzxz, r34\n" //    reals05 += sample0.xxzz * sample5.xzxz
"    mad r34,  r20.yyww, r25.ywyw, r34\n" //    reals05 += sample0.yyww * sample5.ywyw
"    mad r35,  r20.yyww, r25.xzxz, r35\n" //    imags05 += sample0.yyww * sample5.xzxz
"    mul r50,  r20.xxzz, r25.ywyw\n"      //    imags05 -= sample0.xxzz * sample5.ywyw
"    sub r35,  r35,       r50\n"

"    mad r36,  r21.xxzz, r25.xzxz, r36\n" //    reals15 += sample1.xxzz * sample5.xzxz
"    mad r36,  r21.yyww, r25.ywyw, r36\n" //    reals15 += sample1.yyww * sample5.ywyw
"    mad r37,  r21.yyww, r25.xzxz, r37\n" //    imags15 += sample1.yyww * sample5.xzxz
"    mul r50,  r21.xxzz, r25.ywyw\n"      //    imags15 -= sample1.xxzz * sample5.ywyw
"    sub r37,  r37,      r50\n"

"    mad r38,  r22.xxzz, r25.xzxz, r38\n" //    reals25 += sample2.xxzz * sample5.xzxz
"    mad r38,  r22.yyww, r25.ywyw, r38\n" //    reals25 += sample2.yyww * sample5.ywyw
"    mad r39,  r22.yyww, r25.xzxz, r39\n" //    imags25 += sample2.yyww * sample5.xzxz
"    mul r50,  r22.xxzz, r25.ywyw\n"      //    imags25 -= sample2.xxzz * sample5.ywyw
"    sub r39,  r39,      r50\n"
#else
"    mov r0, r20\n"
"    mov r1, r21\n"
"    mov r2, r22\n"
"    mov r3, r23\n"
"    mov r4, r20\n"
"    mov r5, r21\n"
"    mov r6, r22\n"
"    mov r7, r23\n"
"    mov r30, r20\n"
"    mov r31, r21\n"
"    mov r32, r22\n"
"    mov r33, r23\n"
"    mov r34, r20\n"
"    mov r35, r21\n"
"    mov r36, r22\n"
"    mov r37, r23\n"
"    mov r38, r23\n"
"    mov r39, r23\n"
#endif // DO_COMPUTATION

"    add r9.z, r9.z, r9.1\n" // increase loop counter
"endloop\n"

// store

// pos = (nrCells * channel + myCell) * 18 (18 output values)
"mul r40.x, cb0[0].w, r8.y\n"
"add r40.x, r40.x, r8.x\n"
"mul r40.x, r40.x, l18.x\n"
"ftoi r40.x, r40.x\n"

"mov g[r40.x +  0], r0\n"
"mov g[r40.x +  1], r1\n"
"mov g[r40.x +  2], r2\n"
"mov g[r40.x +  3], r3\n"
"mov g[r40.x +  4], r4\n"
"mov g[r40.x +  5], r5\n"
"mov g[r40.x +  6], r6\n"
"mov g[r40.x +  7], r7\n"
"mov g[r40.x +  8], r30\n"
"mov g[r40.x +  9], r31\n"
"mov g[r40.x + 10], r32\n"
"mov g[r40.x + 11], r33\n"
"mov g[r40.x + 12], r34\n"
"mov g[r40.x + 13], r35\n"
"mov g[r40.x + 14], r36\n"
"mov g[r40.x + 15], r37\n"
"mov g[r40.x + 16], r38\n"
"mov g[r40.x + 17], r39\n"

"end\n";


#endif // _CORRELATOR_3x3_MEM_H_
