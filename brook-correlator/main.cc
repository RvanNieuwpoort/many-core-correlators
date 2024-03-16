#include "timer.h"
#include "brook_correlator.h"

#include <stdlib.h>
#include <string.h>

unsigned int nrStations;
unsigned int nrTimes;
unsigned int nrChannels;
unsigned int nrBaselines;
unsigned int nrCells;

timer loadTimer("load");
timer correlateTimer("correlate");
timer storeTimer("store");

void startLoadTimer()
{
    loadTimer.start();
}

void stopLoadTimer()
{
    loadTimer.stop();
}

void startCorrelateTimer()
{
    correlateTimer.start();
}

void stopCorrelateTimer()
{
    correlateTimer.stop();
}

void startStoreTimer()
{
    storeTimer.start();
}

void stopStoreTimer()
{
    storeTimer.stop();
}



void fillBaselineToStatTable(unsigned nrStations, float* baselineToStat1Host, float* baselineToStat2Host)
{
    unsigned baseline;
    unsigned stat1, stat2;

    for (stat2 = baseline = 0; stat2 < nrStations; stat2 ++) {
	for (stat1 = 0; stat1 <= stat2; stat1 ++, baseline ++) {
	    baselineToStat1Host[baseline] = stat1;
	    baselineToStat2Host[baseline] = stat2;
	}
    }
}

void printBaselineToStatTable(unsigned nrBaselines, float* baselineToStat1Host, float* baselineToStat2Host)
{
  unsigned baseline;
  for(baseline = 0; baseline < nrBaselines; baseline++) {
    printf("baseline %d, stat1 = %d, stat2 = %d\n", baseline, (int)baselineToStat1Host[baseline], (int)baselineToStat2Host[baseline]);
  }
}

static unsigned fillCellToStatTable(float* cellToStatX, float* cellToStatY)
{
    unsigned nrCells, stat0, stat2;

    for (stat2 = nrStations % 2 ? 1 : 2, nrCells = 0; stat2 < nrStations; stat2 += 2) {
	for (stat0 = 0; stat0 + 2 <= stat2; stat0 += 2, nrCells ++) {
	    cellToStatX[nrCells] = stat0;
	    cellToStatY[nrCells] = stat2;
	}
    }

    return nrCells;
}

static unsigned calcNrCells()
{
    unsigned nrCells = 0;

    for (unsigned stat2 = nrStations % 2 ? 1 : 2; stat2 < nrStations; stat2 += 2) {
	for (unsigned stat0 = 0; stat0 + 2 <= stat2; stat0 += 2, nrCells ++);
    }

    return nrCells;
}


void printResult(int error, unsigned nrBaselines, unsigned nrChannels, unsigned nrPolarizations, float* visReal, float* visImag)
{
  unsigned baseline, channel, pol;
  unsigned index = 0;

  for(baseline = 0; baseline < nrBaselines; baseline++) {
    for(channel = 0; channel < nrChannels; channel++) {
      for(pol = 0; pol < nrPolarizations * nrPolarizations; pol++) {
	  if(error) {
	      fprintf(stderr, "baseline %d, channel %d, pol %d (%f %f)\n", 
		      baseline, channel, pol, visReal[index], visImag[index]);
	  } else {
	      fprintf(stdout, "baseline %d, channel %d, pol %d (%f %f)\n", 
		      baseline, channel, pol, visReal[index], visImag[index]);
	  }
	index++;
      }
    }
  }
}

void initSamples(unsigned nrSamples, float* samples)
{
  unsigned sample;
  float val = 0.0;
  for(sample = 0; sample < nrSamples*2; sample += 4) {
      samples[sample+0] = val++; // r X
      samples[sample+1] = val++; // i X
      samples[sample+2] = val++; // r Y
      samples[sample+3] = val++; // i Y
  }
}

unsigned visibilityPos(unsigned baseline, unsigned channel, unsigned nrChannels)
{
  return (baseline * nrChannels + channel) * 4;
}

unsigned samplePos(unsigned station, unsigned channel, unsigned time, unsigned nrStations, unsigned nrTimes) {
  return ((channel * nrStations + station) * nrTimes + time) * 2 * 2;
}

void hostCorrelator(unsigned nrBaselines, unsigned nrStations, unsigned nrChannels, unsigned nrTimes, 
		    float* baselineToStat1, float* baselineToStat2, float* samples, float* visReal, float* visImag) 
{
  unsigned baseline, channel, time, visPos;

  for(baseline = 0; baseline < nrBaselines; baseline++) {
    for(channel = 0; channel < nrChannels; channel++) {
      unsigned stat1 = (int) baselineToStat1[baseline];
      unsigned stat2 = (int) baselineToStat2[baseline];

      float xxr = 0.0f, xxi = 0.0f, xyr = 0.0f, xyi = 0.0f, yxr = 0.0f, yxi = 0.0f, yyr = 0.0f, yyi = 0.0f;

      for(time = 0; time < nrTimes; time++) {
	// get samples
	unsigned sample1Pos = samplePos(stat1, channel, time, nrStations, nrTimes);
	float sample1_xr = samples[sample1Pos + 0];
	float sample1_xi = samples[sample1Pos + 1];
	float sample1_yr = samples[sample1Pos + 2];
	float sample1_yi = samples[sample1Pos + 3];

	unsigned sample2Pos = samplePos(stat2, channel, time, nrStations, nrTimes);
	float sample2_xr = samples[sample2Pos + 0];
	float sample2_xi = samples[sample2Pos + 1];
	float sample2_yr = samples[sample2Pos + 2];
	float sample2_yi = samples[sample2Pos + 3];

	xxr += sample1_xr * sample2_xr;
	xxr += sample1_xi * sample2_xi;
	xxi += sample1_xi * sample2_xr;
	xxi -= sample1_xr * sample2_xi;
	xyr += sample1_xr * sample2_yr;
	xyr += sample1_xi * sample2_yi;
	xyi += sample1_xi * sample2_yr;
	xyi -= sample1_xr * sample2_yi;
	yxr += sample1_yr * sample2_xr;
	yxr += sample1_yi * sample2_xi;
	yxi += sample1_yi * sample2_xr;
	yxi -= sample1_yr * sample2_xi;
	yyr += sample1_yr * sample2_yr;
	yyr += sample1_yi * sample2_yi;
	yyi += sample1_yi * sample2_yr;
	yyi -= sample1_yr * sample2_yi;
      }

      // and store
      visPos = visibilityPos(baseline, channel, nrChannels);

      visReal[visPos + 0] = xxr;
      visReal[visPos + 1] = xyr;
      visReal[visPos + 2] = yxr;
      visReal[visPos + 3] = yyr;
      visImag[visPos + 0] = xxi;
      visImag[visPos + 1] = xyi;
      visImag[visPos + 2] = yxi;
      visImag[visPos + 3] = yyi;
    }
  }
}

int main(int argc, char** argv)
{
  nrStations      = 2;
  nrTimes	       = 768;
  nrChannels      = 2;
  nrBaselines    = nrStations * (nrStations + 1) / 2;

  unsigned int nrPolarizations = 2;
  unsigned int nrSamples      = nrStations * nrChannels * nrTimes * nrPolarizations;
  unsigned int nrVisibilities = nrBaselines * nrChannels * nrPolarizations * nrPolarizations;
  nrCells = calcNrCells();

  std::cout << "running with " << nrStations << " stations, " 
	    << nrBaselines << " baselines, " << nrCells << " 2x2 cells, " << nrChannels << " channels" << std::endl;
  std::cout << "sample buf = " << (sizeof(float) * nrSamples * 2.0) / (1024.0 * 1024.0) << " MB" << std::endl;
  std::cout << "visibilities buf = " << (nrVisibilities * sizeof(float) * 2.0) / (1024.0 * 1024.0) << " MB" << std::endl;
  std::cout << "nr visibilities = " << nrVisibilities << " (" << nrBaselines << " x " << nrChannels << ") " << std::endl;

  float* hostSamples = (float*) malloc(sizeof *hostSamples * nrSamples * 2);
  float* hostVisReal = (float*) malloc(nrVisibilities * sizeof *hostVisReal);
  float* hostVisImag = (float*) malloc(nrVisibilities * sizeof *hostVisImag);
  float* baselineToStat1Host = (float*) malloc(nrBaselines * sizeof *baselineToStat1Host);
  float* baselineToStat2Host = (float*) malloc(nrBaselines * sizeof *baselineToStat2Host);
  float* cellToStatXHost = (float*) malloc(sizeof(float) * nrCells);
  float* cellToStatYHost = (float*) malloc(sizeof(float) * nrCells);

  unsigned long long ops;
  double elapsed;
  double flops;
  double loadTime;
  double storeTime;
  double totalTime;
  double totalFlops;


  fillBaselineToStatTable(nrStations, baselineToStat1Host, baselineToStat2Host);
//  printBaselineToStatTable(nrBaselines, baselineToStat1Host, baselineToStat2Host);

  fillCellToStatTable(cellToStatXHost, cellToStatYHost);

  initSamples(nrSamples, hostSamples);
  memset(hostVisReal, 0, nrVisibilities * sizeof *hostVisReal);
  memset(hostVisImag, 0, nrVisibilities * sizeof *hostVisImag);

  correlate_1x1(baselineToStat1Host, baselineToStat2Host, hostSamples, hostVisReal, hostVisImag);

  fprintf(stdout, "gpu result:\n");
  printResult(0, nrBaselines, nrChannels, nrPolarizations, hostVisReal, hostVisImag);

  memset(hostVisReal, 0, nrVisibilities * sizeof *hostVisReal);
  memset(hostVisImag, 0, nrVisibilities * sizeof *hostVisImag);

  hostCorrelator(nrBaselines, nrStations, nrChannels, nrTimes, baselineToStat1Host, baselineToStat2Host, hostSamples, hostVisReal, hostVisImag);

  fprintf(stderr, "cpu result:\n");
  printResult(1, nrBaselines, nrChannels, nrPolarizations, hostVisReal, hostVisImag);

  // ops only valid for 1x1
  ops = (unsigned long long) nrChannels * nrBaselines * nrTimes * 16L * 2L;
  elapsed = correlateTimer.getTimeInSeconds();
  flops = (ops / elapsed) / 1000000000.0;  
  loadTime = loadTimer.getTimeInSeconds();
  storeTime = storeTimer.getTimeInSeconds();
  totalTime = elapsed + loadTime + storeTime;
  totalFlops = (ops / totalTime) / 1000000000.0;  

  std::cout << "ops = " << ops << std::endl;
  std::cout << "correlate_1x1 took " << elapsed << " s, achieved " << flops << " Gflops" 
	    << ", load time = " << loadTime 
	    << ", store time = " << storeTime 
	    << ", total FLOPS including mem I/O = " << totalFlops << std::endl;
  loadTimer.reset(); correlateTimer.reset(); storeTimer.reset();
  memset(hostVisReal, 0, nrVisibilities * sizeof *hostVisReal);
  memset(hostVisImag, 0, nrVisibilities * sizeof *hostVisImag);

  correlate_1x1_vec(baselineToStat1Host, baselineToStat2Host, hostSamples, hostVisReal, hostVisImag);

  ops = nrChannels * nrBaselines * nrTimes * 16L * 2L;
  elapsed = correlateTimer.getTimeInSeconds();
  flops = (ops / elapsed) / 1000000000.0;  
  loadTime = loadTimer.getTimeInSeconds();
  storeTime = storeTimer.getTimeInSeconds();
  totalTime = elapsed + loadTime + storeTimer.getTimeInSeconds();
  totalFlops = (ops / totalTime) / 1000000000.0;  
  std::cout << "ops = " << ops << std::endl;
  std::cout << "correlate_1x1_vec took " << elapsed << " s, achieved " << flops << " Gflops" 
	    << ", load time = " << loadTime 
	    << ", store time = " << storeTime 
	    << ", total FLOPS including mem I/O = " << totalFlops << std::endl;
  loadTimer.reset(); correlateTimer.reset(); storeTimer.reset();
  memset(hostVisReal, 0, nrVisibilities * sizeof *hostVisReal);
  memset(hostVisImag, 0, nrVisibilities * sizeof *hostVisImag);

  correlate_2x2_vec(cellToStatXHost, cellToStatYHost, hostSamples, hostVisReal, hostVisImag);

  ops =  nrChannels * nrCells * nrTimes * 16L * 4L * 2L;
  elapsed = correlateTimer.getTimeInSeconds();
  flops = (ops / elapsed) / 1000000000.0;  
  loadTime = loadTimer.getTimeInSeconds();
  storeTime = storeTimer.getTimeInSeconds();
  totalTime = elapsed + loadTime + storeTimer.getTimeInSeconds();
  totalFlops = (ops / totalTime) / 1000000000.0;  
  std::cout << "ops = " << ops << std::endl;
  std::cout << "correlate_2x2_vec took " << elapsed << " s, achieved " << flops << " Gflops" 
	    << ", load time = " << loadTime 
	    << ", store time = " << storeTime 
	    << ", total FLOPS including mem I/O = " << totalFlops << std::endl;
  loadTimer.reset(); correlateTimer.reset(); storeTimer.reset();

//  fprintf(stdout, "gpu result:\n");
//  printResult(0, nrBaselines, nrChannels, nrPolarizations, hostVisReal, hostVisImag);

//  hostCorrelator(nrBaselines, nrStations, nrChannels, nrTimes, baselineToStat1Host, baselineToStat2Host, hostSamples, hostVisReal, hostVisImag);

//  fprintf(stderr, "cpu result:\n");
//  printResult(1, nrBaselines, nrChannels, nrPolarizations, hostVisReal, hostVisImag);

  free(hostSamples);
  free(hostVisReal);
  free(hostVisImag);
  free(baselineToStat1Host);
  free(baselineToStat2Host);
}
