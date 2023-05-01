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
#ifndef _CORRELATOR_1x1_H_
#define _CORRELATOR_1x1_H_

const CALchar* Kernel_1x1 =
"il_ps_2_0\n"
"dcl_input_position_interp(linear_noperspective) vWinCoord0.xy__\n" // position in domain (output buffer)
"dcl_output_generic o0\n"
"dcl_output_generic o1\n"
"dcl_resource_id(0)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // samples
"dcl_resource_id(1)_type(1d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // cellToStatX
"dcl_resource_id(2)_type(1d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // cellToStatY
"dcl_literal l0, 0x00000000, 0x00000000, 0x00000000, 0x00000000\n"
"dcl_literal l1, 1.0, 2.0, 3.0, 4.0\n"
"dcl_cb cb0[1]\n"

// o0 -- o7 are the computed visibilities
"mov o0, l0\n" // reals01
"mov o1, l0\n" // imags01

"mov r8.x, vWinCoord0.x\n"   // r8.x = cell
"mov r8.y00, vWinCoord0.y\n" // r8.y = channel
"mov r9.x, cb0[0].x\n"       // r9.x = nrStations (float)
"mov r9.y, cb0[0].y\n"       // r9.y = nrTimes (float)
"mov r9.zw, l0\n"            // r9.z = time loopcounter (0 .. nrTimes)

"sample_resource(1)_sampler(0)   r8.z, r8.x\n" // r8.z = stat0 (X)
"sample_resource(2)_sampler(1)   r8.w, r8.x\n" // r8.w = stat1 (Y)

// indices should be (channel * nrStations) + station, time, 0, 0
"mul r10.y, r8.y, r9.x\n"
"add r10.y, r10.y, r8.z\n"
"mov r10.xzw, l0\n" // r10 = index0: {time, (channel * nrStations) + station0, 0, 0}

"mul r12.y, r8.y, r9.x\n"
"add r12.y, r12.y, r8.w\n"
"mov r12.xzw, l0\n" // r10 = index0: {time, (channel * nrStations) + station1, 0, 0}

"whileloop\n"
"    ge r9.w, r9.z, r9.y\n"  // r9.w = temp used for comparison
"    break_logicalnz r9.w\n"

//   update time index
"    mov r10.x, r9.z\n"
"    mov r12.x, r9.z\n"

//   load the four samples
#if DO_LOADS
"    sample_resource(0)_sampler(0) r20, r10\n" // sample 0
"    sample_resource(0)_sampler(2) r22, r12\n" // sample 2
#else
"    mov r20, r10\n"
"    mov r22, r12\n"
#endif

#if DO_COMPUTATION
"    mad o0,  r20.xxzz, r22.xzxz, o0\n" //    reals02 += sample0.xxzz * sample2.xzxz
"    mad o0,  r20.yyww, r22.ywyw, o0\n" //    reals02 += sample0.yyww * sample2.ywyw
"    mad o1,  r20.yyww, r22.xzxz, o1\n" //    imags02 += sample0.yyww * sample2.xzxz
"    mul r30, r20.xxzz, r22.ywyw\n"     //    imags02 -= sample0.xxzz * sample2.ywyw
"    sub o1,  o1,       r30\n"
#else
"    mov o0, r20\n"
"    mov o1, r21\n"
#endif // DO_COMPUTATION

"    add r9.z, r9.z, r9.1\n" // increase loop counter
"endloop\n"

"end\n";


#endif // _CORRELATOR_1x1_H_
