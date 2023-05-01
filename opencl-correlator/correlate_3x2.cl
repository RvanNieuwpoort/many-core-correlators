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
__kernel void correlate_3x2(__global const float4 *samples, __global float4 *visibilities, 
			    __global const unsigned *cellToStatX, __global const unsigned *cellToStatY,
			    const unsigned nrStations, const unsigned nrTimes, const unsigned nrTimesWidth, const unsigned nrChannels, const unsigned nrCells)
{
  const unsigned int lid = get_local_id(0);
  const unsigned int lsize = get_local_size(0);
  const unsigned int gid = get_group_id(0);
  const unsigned int gsize = get_num_groups(0);

  for (unsigned channel = gid; channel < nrChannels; channel += gsize) {       
    for (unsigned cell = lid; cell < nrCells; cell += lsize) {
      
      unsigned stat0 = cellToStatX[cell];
      unsigned stat1 = cellToStatY[cell];
      
      unsigned index0 = (channel * nrStations + stat0) * nrTimesWidth;
      unsigned index1 = (channel * nrStations + stat1) * nrTimesWidth;
      float xxr = 0, xxi = 0, xyr = 0, xyi = 0, yxr = 0, yxi = 0, yyr = 0, yyi = 0;
      
      for (unsigned time = 0; time < nrTimes; time ++) {
	float4 sample0 = samples[index0 + time];
	float4 sample1 = samples[index1 + time];

	xxr += sample0.x * sample1.x;
	xxi += sample0.y * sample1.x;
	xyr += sample0.x * sample1.z;
	xyi += sample0.y * sample1.z;
	yxr += sample0.z * sample1.x;
	yxi += sample0.w * sample1.x;
	yyr += sample0.z * sample1.z;
	yyi += sample0.w * sample1.z;

	xxr += sample0.y * sample1.y;
	xxi -= sample0.x * sample1.y;
	xyr += sample0.y * sample1.w;
	xyi -= sample0.x * sample1.w;
	yxr += sample0.w * sample1.y;
	yxi -= sample0.z * sample1.y;
	yyr += sample0.w * sample1.w;
	yyi -= sample0.z * sample1.w;
      }

      unsigned baseline = stat1 * (stat1 + 1) / 2 + stat0;
      unsigned visibilityIndex = (baseline * nrChannels + channel) * 2;
      __global float *dst = (__global float*) (visibilities + visibilityIndex);
      dst[0] = xxr;
      dst[1] = xxi;
      dst[2] = xyr;
      dst[3] = xyi;
      dst[4] = yxr;
      dst[5] = yxi;
      dst[6] = yyr;
      dst[7] = yyi;
    }
  }
}
