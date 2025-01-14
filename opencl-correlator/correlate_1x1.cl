__kernel void correlate_1x1(__read_only image3d_t samples,
			    sampler_t sampler,		
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
      unsigned stat1 = cellToStatY[cell];
      float4 real = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
      float4 imag = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

      for (unsigned time = 0; time < nrTimes; time ++) {
	float4 sample0 = read_imagef(samples, sampler, (float4)(channel, stat0, time*4, 1.0f));
	float4 sample1 = read_imagef(samples, sampler, (float4)(channel, stat1, time*4, 1.0f));

	real += sample0.xxzz * sample1.xzxz;
	real += sample0.yyww * sample1.ywyw;
	imag += sample0.yyww * sample1.xzxz;
	imag -= sample0.xxzz * sample1.ywyw;
      }

      unsigned baseline = stat1 * (stat1 + 1) / 2 + stat0;
      unsigned visibilityIndex = (baseline * nrChannels + channel) * 2;
      __global float *dst = (__global float*) (visibilities + visibilityIndex);
      dst[0] = real.x; 
      dst[1] = imag.x;
      dst[2] = real.y;
      dst[3] = imag.y;
      dst[4] = real.z;
      dst[5] = imag.z;
      dst[6] = real.w;
      dst[7] = imag.w;
    }
  }
}
