#ifndef _CORRELATOR_1x1_MEM_H_
#define _CORRELATOR_1x1_MEM_H_

// split the number of channels over several buffers.
// the number of channels is easier than the number of stations,
// all operations here work on a single channel, not on a single station.

const CALchar* Kernel_1x1_mem =
"il_ps_2_0\n"
"dcl_input_position_interp(linear_noperspective) vWinCoord0.xy__\n" // position in domain (output buffer)
"dcl_resource_id(0)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // samples
"dcl_resource_id(1)_type(1d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // cellToStatX
"dcl_resource_id(2)_type(1d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // cellToStatY
"dcl_literal l0, 0x00000000, 0x00000000, 0x00000000, 0x00000000\n"
"dcl_literal l1, 1.0, 2.0, 3.0, 4.0\n"
"dcl_literal l8, 8.0, 8.0, 8.0, 8.0\n"
"dcl_cb cb0[1]\n"

// r0 -- r7 are the computed visibilities
"mov r0, l0\n" // reals02
"mov r1, l0\n" // imags02

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
"mad r10.y, r8.y, r9.x, r8.z\n"
"mov r10.xzw, l0\n" // r10 = index0: {time, (channel * nrStations) + station0+0, 0, 0}

"mad r12.y, r8.y, r9.x, r8.w\n"
"mov r12.xzw, l0\n" // r12 = index0: {time, (channel * nrStations) + station1+0, 0, 0}

"whileloop\n"
"    ge r9.w, r9.z, r9.y\n"  // r9.w = temp used for comparison
"    break_logicalnz r9.w\n"

//   update time index
"    mov r10.x, r9.z\n"
"    mov r12.x, r9.z\n"

//   load the two samples
#if DO_LOADS
"    sample_resource(0)_sampler(0) r20, r10\n" // sample 0
"    sample_resource(0)_sampler(2) r22, r12\n" // sample 2
#else
"    mov r20, r10\n"
"    mov r22, r12\n"
#endif

#if DO_COMPUTATION
"    mad r0,  r20.xxzz, r22.xzxz, r0\n"  //    reals02 += sample0.xxzz * sample2.xzxz
"    mad r0,  r20.yyww, r22.ywyw, r0\n"  //    reals02 += sample0.yyww * sample2.ywyw
"    mad r1,  r20.yyww, r22.xzxz, r1\n"  //    imags02 += sample0.yyww * sample2.xzxz
"    mul r30, r20.xxzz, r22.ywyw\n"      //    imags02 -= sample0.xxzz * sample2.ywyw
"    sub r1,  r1,       r30\n"
#else
"    mov r0, r20\n"
"    mov r1, r22\n"
#endif // DO_COMPUTATION

"    add r9.z, r9.z, r9.1\n" // increase loop counter
"endloop\n"

// store

// pos = (nrCells * channel + myCell) * 8 (8 output registers)
"mad r40.x, cb0[0].w, r8.y, r8.x\n"
"mul r40.x, r40.x, l8.x\n"
"ftoi r40.x, r40.x\n"

#if 1
"mov g[r40.x + 0], r0\n"
"mov g[r40.x + 1], r1\n"
#endif

"end\n";


#endif // _CORRELATOR_1x1_MEM_H_
