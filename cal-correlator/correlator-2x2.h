#ifndef _CORRELATOR_2x2_H_
#define _CORRELATOR_2x2_H_

const CALchar* Kernel_2x2 =
"il_ps_2_0\n"
"dcl_input_position_interp(linear_noperspective) vWinCoord0.xy__\n" // position in domain (output buffer)
"dcl_output_generic o0\n"
"dcl_output_generic o1\n"
"dcl_output_generic o2\n"
"dcl_output_generic o3\n"
"dcl_output_generic o4\n"
"dcl_output_generic o5\n"
"dcl_output_generic o6\n"
"dcl_output_generic o7\n"
"dcl_resource_id(0)_type(2d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // samples
"dcl_resource_id(1)_type(1d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // cellToStatX
"dcl_resource_id(2)_type(1d,unnorm)_fmtx(float)_fmty(float)_fmtz(float)_fmtw(float)\n" // cellToStatY
"dcl_literal l0, 0x00000000, 0x00000000, 0x00000000, 0x00000000\n"
"dcl_literal l1, 1.0, 2.0, 3.0, 4.0\n"
"dcl_cb cb0[1]\n"

// o0 -- o7 are the computed visibilities
"mov o0, l0\n" // reals02
"mov o1, l0\n" // imags02
"mov o2, l0\n" // reals12
"mov o3, l0\n" // imags12
"mov o4, l0\n" // reals03
"mov o5, l0\n" // imags03
"mov o6, l0\n" // reals13
"mov o7, l0\n" // imags13

"mov r8.x, vWinCoord0.x\n"   // r8.x = cell
"mov r8.y00, vWinCoord0.y\n" // r8.y = channel
"mov r9.x, cb0[0].x\n"       // r9.x = nrStations (float)
"mov r9.y, cb0[0].y\n"       // r9.y = nrTimes (float)
"mov r9.zw, l0\n"            // r9.z = time loopcounter (0 .. nrTimes)

"sample_resource(1)_sampler(0)   r8.z, r8.x\n" // r8.z = stat0 (X)
"sample_resource(2)_sampler(1)   r8.w, r8.x\n" // r8.w = stat2 (Y)

// indices should be (channel * nrStations) + station, time, 0, 0
"mul r10.y, r8.y, r9.x\n"
"add r10.y, r10.y, r8.z\n"
"mov r10.xzw, l0\n" // r10 = index0: {time, (channel * nrStations) + station0+0, 0, 0}

"mov r11.y, r10.y\n"
"add r11.y, r11.y, r11.1\n"
"mov r11.xzw, l0\n" // r10 = index0: {time, (channel * nrStations) + station0+1, 0, 0}

"mul r12.y, r8.y, r9.x\n"
"add r12.y, r12.y, r8.w\n"
"mov r12.xzw, l0\n" // r10 = index0: {time, (channel * nrStations) + station1+0, 0, 0}

"mov r13.y, r13.y\n"
"add r13.y, r13.y, r13.1\n"
"mov r13.xzw, l0\n" // r10 = index0: {time, (channel * nrStations) + station1+1, 0, 0}

"whileloop\n"
"    ge r9.w, r9.z, r9.y\n"  // r9.w = temp used for comparison
"    break_logicalnz r9.w\n"

//   update time index
"    mov r10.x, r9.z\n"
"    mov r11.x, r9.z\n"
"    mov r12.x, r9.z\n"
"    mov r13.x, r9.z\n"

//   load the four samples
#if DO_LOADS
"    sample_resource(0)_sampler(0) r20, r10\n" // sample 0
"    sample_resource(0)_sampler(1) r21, r11\n" // sample 1
"    sample_resource(0)_sampler(2) r22, r12\n" // sample 2
"    sample_resource(0)_sampler(3) r23, r13\n" // sample 3
#else
"    mov r20, r10\n"
"    mov r21, r11\n"
"    mov r22, r12\n"
"    mov r23, r13\n"
#endif

#if DO_COMPUTATION
"    mad o0,  r20.xxzz, r22.xzxz, o0\n" //    reals02 += sample0.xxzz * sample2.xzxz
"    mad o0,  r20.yyww, r22.ywyw, o0\n" //    reals02 += sample0.yyww * sample2.ywyw
"    mad o1,  r20.yyww, r22.xzxz, o1\n" //    imags02 += sample0.yyww * sample2.xzxz
"    mul r30, r20.xxzz, r22.ywyw\n"     //    imags02 -= sample0.xxzz * sample2.ywyw
"    sub o1,  o1,       r30\n"

    // the load of sample 1 can be done here

"    mad o2,  r21.xxzz, r22.xzxz, o2\n" //    reals12 += sample1.xxzz * sample2.xzxz
"    mad o2,  r21.yyww, r22.ywyw, o2\n" //    reals12 += sample1.yyww * sample2.ywyw
"    mad o3,  r21.yyww, r22.xzxz, o3\n" //    imags12 += sample1.yyww * sample2.xzxz
"    mul r30, r21.xxzz, r22.ywyw\n"     //    imags12 -= sample1.xxzz * sample2.ywyw
"    sub o3,  o3,       r30\n"

    // load of sample3 can be done here, reusing sample2's registers

"    mad o4,  r20.xxzz, r23.xzxz, o4\n" //    reals03 += sample0.xxzz * sample3.xzxz
"    mad o4,  r20.yyww, r23.ywyw, o4\n" //    reals03 += sample0.yyww * sample3.ywyw
"    mad o5,  r20.yyww, r23.xzxz, o5\n" //    imags03 += sample0.yyww * sample3.xzxz
"    mul r30, r20.xxzz, r23.ywyw\n"     //    imags03 -= sample0.xxzz * sample3.ywyw
"    sub o5,  o5,       r30\n"

"    mad o6,  r21.xxzz, r23.xzxz, o6\n" //    reals13 += sample1.xxzz * sample3.xzxz
"    mad o6,  r21.yyww, r23.ywyw, o6\n" //    reals13 += sample1.yyww * sample3.ywyw
"    mad o7,  r21.yyww, r23.xzxz, o7\n" //    imags13 += sample1.yyww * sample3.xzxz
"    mul r30, r21.xxzz, r23.ywyw\n"     //    imags13 -= sample1.xxzz * sample3.ywyw
"    sub o7,  o7,       r30\n"
#else
"    mov o0, r20\n"
"    mov o1, r21\n"
"    mov o2, r22\n"
"    mov o3, r23\n"
"    mov o4, r20\n"
"    mov o5, r21\n"
"    mov o6, r22\n"
"    mov o7, r23\n"
#endif // DO_COMPUTATION

"    add r9.z, r9.z, r9.1\n" // increase loop counter
"endloop\n"

"end\n";


#endif // _CORRELATOR_2x2_H_
