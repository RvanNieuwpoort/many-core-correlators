__kernel void correlate_2x2_standard(__global float4 *samples,
				     __global float4 *visibilities,
				     __global const unsigned *cellToStatX, __global const unsigned *cellToStatY,
				     const unsigned nrStations, const unsigned nrTimes, const unsigned nrTimesWidth, 
				     const unsigned nrChannels, const unsigned nrCells)
{
  const unsigned int lid = get_local_id(0);
  const unsigned int lsize = get_local_size(0);
  const unsigned int gid = get_group_id(0);
  const unsigned int gsize = get_num_groups(0);

  for (unsigned channel = gid; channel < nrChannels; channel += gsize) {       
    for (unsigned cell = lid; cell < nrCells; cell += lsize) {
      
      unsigned stat0 = cellToStatX[cell];
      unsigned stat1 = stat0+1;
      unsigned stat2 = cellToStatY[cell];
      unsigned stat3 = stat2+1;

      unsigned index0 = (channel * nrStations + stat0) * nrTimesWidth;
      unsigned index1 = (channel * nrStations + stat1) * nrTimesWidth;
      unsigned index2 = (channel * nrStations + stat2) * nrTimesWidth;
      unsigned index3 = (channel * nrStations + stat3) * nrTimesWidth;

      float4 real02 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      float4 imag02 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      float4 real12 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      float4 imag12 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      float4 real03 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      float4 imag03 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      float4 real13 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      float4 imag13 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

      for (unsigned time = 0; time < nrTimes; time ++) {
	float4 sample0 = samples[index0 + time];
	float4 sample1 = samples[index1 + time];
	float4 sample2 = samples[index2 + time];
	float4 sample3 = samples[index3 + time];

	real02 += sample0.xxzz * sample2.xzxz;
	real02 += sample0.yyww * sample2.ywyw;
	imag02 += sample0.yyww * sample2.xzxz;
	imag02 -= sample0.xxzz * sample2.ywyw;

	real12 += sample1.xxzz * sample2.xzxz;
	real12 += sample1.yyww * sample2.ywyw;
	imag12 += sample1.yyww * sample2.xzxz;
	imag12 -= sample1.xxzz * sample2.ywyw;

	real03 += sample0.xxzz * sample3.xzxz;
	real03 += sample0.yyww * sample3.ywyw;
	imag03 += sample0.yyww * sample3.xzxz;
	imag03 -= sample0.xxzz * sample3.ywyw;

	real13 += sample1.xxzz * sample3.xzxz;
	real13 += sample1.yyww * sample3.ywyw;
	imag13 += sample1.yyww * sample3.xzxz;
	imag13 -= sample1.xxzz * sample3.ywyw;
      }

      unsigned baseline = stat2 * (stat2 + 1) / 2 + stat0;
      unsigned visibilityIndex = (baseline * nrChannels + channel) * 2;
      __global float *dst = (__global float*) (visibilities + visibilityIndex);
      dst[0] = real02.x; 
      dst[1] = imag02.x;
      dst[2] = real02.y;
      dst[3] = imag02.y;
      dst[4] = real02.z;
      dst[5] = imag02.z;
      dst[6] = real02.w;
      dst[7] = imag02.w;

      baseline = stat2 * (stat2 + 1) / 2 + stat1;
      visibilityIndex = (baseline * nrChannels + channel) * 2;
      dst = (__global float*) (visibilities + visibilityIndex);
      dst[0] = real12.x; 
      dst[1] = imag12.x;
      dst[2] = real12.y;
      dst[3] = imag12.y;
      dst[4] = real12.z;
      dst[5] = imag12.z;
      dst[6] = real12.w;
      dst[7] = imag12.w;

      baseline = stat3 * (stat3 + 1) / 2 + stat0;
      visibilityIndex = (baseline * nrChannels + channel) * 2;
      dst = (__global float*) (visibilities + visibilityIndex);
      dst[0] = real03.x; 
      dst[1] = imag03.x;
      dst[2] = real03.y;
      dst[3] = imag03.y;
      dst[4] = real03.z;
      dst[5] = imag03.z;
      dst[6] = real03.w;
      dst[7] = imag03.w;

      baseline = stat3 * (stat3 + 1) / 2 + stat1;
      visibilityIndex = (baseline * nrChannels + channel) * 2;
      dst = (__global float*) (visibilities + visibilityIndex);
      dst[0] = real13.x; 
      dst[1] = imag13.x;
      dst[2] = real13.y;
      dst[3] = imag13.y;
      dst[4] = real13.z;
      dst[5] = imag13.z;
      dst[6] = real13.w;
      dst[7] = imag13.w;
    }
  }
}
