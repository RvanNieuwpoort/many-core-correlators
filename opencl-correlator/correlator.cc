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
#include <CL/cl.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include "gpu_complex.h"
#include "correlator.h"
#include "timer.h"

using namespace std;

const bool useGPU = true;
const bool printResult = false;
const bool verify = false;

const unsigned nrStations	 = 64;
const unsigned nrTimes	 = 768, nrTimesWidth = 768; // 770
const unsigned nrChannels	 = 200;
const unsigned nrPolarizations = 2;
const unsigned iter = 100;
const unsigned cellWidth = 2;
const unsigned cellHeight = 2;

unsigned* cellToStatX;
unsigned* cellToStatY;


void correlateOnHost(const complex<float> *samples, complex<float> *visibilities)
{
    unsigned nrBaselines = nrStations * (nrStations + 1) / 2;
    memset(visibilities, 0, nrBaselines * nrChannels * 2 * 2 * sizeof(complex<float>));

    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for (unsigned stat1 = 0; stat1 < nrStations; stat1 ++) {
	    for (unsigned stat0 = 0; stat0 <= stat1; stat0 ++) {
		for (unsigned time = 0; time < nrTimes; time ++) {
		    for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
			for (unsigned pol1 = 0; pol1 < 2; pol1 ++) { 
			    size_t index = VISIBILITIES_INDEX(BASELINE(stat0, stat1), channel, pol0, pol1, 0)/2;
			    visibilities[index]
				+=  samples[SAMPLE_INDEX(stat0, channel, time, pol0, 0)/2] 
				* ~(samples[SAMPLE_INDEX(stat1, channel, time, pol1, 0)/2]);
			}
		    }
		}
	    }
	}
    }
}

string* openKernel(const char* fileName)
{
  ifstream myfile(fileName);
  string line;
  string* all = new string;
  if (myfile.is_open())
  {
    while (!myfile.eof() )
    {
      getline (myfile,line);
      all->append(line);
    }
    myfile.close();
    return all;
  } else {
    return 0;
  }
}

unsigned calcNrCells(const unsigned w, const unsigned h, const unsigned nrStations) 
{
    unsigned nrCells = 0;
    for(int statY = nrStations - h; statY >= 0; statY -= h) {
      for(int statX = 0; statX + (int)w - 1 <= statY; statX += w) {
	    nrCells++;
	}
    }

    return nrCells;
}

unsigned fillCellToStatTable(const unsigned w, const unsigned h, const unsigned nrStations, 
			     cl_mem deviceCellToStatX, cl_mem deviceCellToStatY, cl_command_queue commands) 
{
  unsigned nrCells = calcNrCells(w, h, nrStations);
  cellToStatX = new unsigned[nrCells];
  cellToStatY = new unsigned[nrCells];

  unsigned index = 0;
  for(int statY = nrStations - h; statY >= 0; statY -= h) {
    for(int statX = 0; statX + (int)w - 1 <= statY; statX += w) {
      cellToStatX[index] = statX;
      cellToStatY[index] = statY;
      index++;
    }
  }

  cl_int err = clEnqueueWriteBuffer(commands, deviceCellToStatX, CL_TRUE, 0, nrCells * sizeof(unsigned), (void *)cellToStatX, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    printf("Error: Failed to write to cellToStat table on device!\n");
    return EXIT_FAILURE;
  }

  err = clEnqueueWriteBuffer(commands, deviceCellToStatY, CL_TRUE, 0, nrCells * sizeof(unsigned), (void *)cellToStatY, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    printf("Error: Failed to write to cellToStat table on device!\n");
    return EXIT_FAILURE;
  }

  return nrCells;
}

void printVisibility(unsigned channel, unsigned cell, unsigned sXStart, unsigned sYStart, 
		     unsigned statX, unsigned statY, unsigned baseline, unsigned pol0, unsigned pol1,
		     complex<float> *hostVisibilities, complex<float> *checkVis, size_t index)
{
    cout.width(3);
    cout << "chan ";
    cout.width(3);
    cout << channel;
    cout << ", cellNr ";
    cout.width(5);
    cout << cell;
    cout << ", XStart ";
    cout.width(2);
    cout << sXStart;
    cout << ", YStart ";
    cout.width(2);
    cout << sYStart;
    cout << ", sX ";
    cout.width(2);
    cout << statX;
    cout << ", sY ";
    cout.width(2);
    cout << statY;
    cout << ", bl ";
    cout.width(5);
    cout << baseline;
    cout << ", pol " << pol0 << '/' << pol1 << ": (";
    cout.precision(10);
    cout.width(10);
    cout << hostVisibilities[index].real;
    cout << ", ";
    cout.precision(10);
    cout.width(10);
    cout << hostVisibilities[index].imag;
    cout << ") should be (";
    cout.precision(10);
    cout.width(10);
    cout << checkVis[index].real;
    cout << ", ";
    cout.precision(10);
    cout.width(10);
    cout << checkVis[index].imag;
    cout << ")";
    if (hostVisibilities[index] != checkVis[index]) {
	cout << " !!";
    } else {
	cout << " OK";
    }
    cout << endl;
}

void checkResult(size_t visibilitiesSize, complex<float> *hostSamples,
		 complex<float> *hostVisibilities, unsigned nrBaselines, unsigned nrCells)
{
    cout << "checking result..." << endl;
    unsigned cellWidth = 1;
    unsigned cellHeight = 1;

    complex<float> *checkVis = new complex<float>[visibilitiesSize / 2];
    correlateOnHost(hostSamples, checkVis);

    bool errorFound = false;
    for (unsigned channel = 0; channel < nrChannels; channel ++) {
	for(unsigned cell = 0; cell < nrCells; cell++) {
	    int sXStart = cellToStatX[cell];
	    int sYStart = cellToStatY[cell];

	    for(unsigned h = 0; h<cellHeight; h++) {
		for(unsigned w = 0; w<cellWidth; w++) {
		    int statX = sXStart + w;
		    int statY = sYStart + h;
		    unsigned baseline = BASELINE(statX, statY);

		    for (unsigned pol0 = 0; pol0 < 2; pol0 ++) {
			for (unsigned pol1 = 0; pol1 < 2; pol1 ++) {
			    size_t index = VISIBILITIES_INDEX(baseline, channel, pol0, pol1, 0)/2;
			    if (hostVisibilities[index] != checkVis[index]) {
				errorFound = true;
			    }
			    if(printResult) {
				printVisibility(channel, cell, sXStart, sYStart, statX, statY, baseline, pol0, pol1, hostVisibilities, checkVis, index);

			    }
			}
		    }
		}
	    }
	}
    }

    if(errorFound) {
	cout << "ERRORS in result." << endl;
    } else {
	cout << "result OK." << endl;
    }

    delete [] checkVis;
}

cl_mem allocateDeviceSamples(cl_context context, cl_sampler* sampler, size_t samplesSize)
{
  cl_int err;
  cl_mem deviceSamples;

#if ! USE_TEXTURE
  *sampler = NULL;
  // Create the samples buffer on the device
  deviceSamples = clCreateBuffer(context, CL_MEM_READ_WRITE, samplesSize * sizeof(float), NULL, &err);
  if (err != CL_SUCCESS || !deviceSamples) {
    printf("Error: Failed to allocate sample data buffer on device!\n");
    exit(EXIT_FAILURE);
  }
#else
  cerr << "samples image size = " << nrChannels << " x " << nrStations << " x " << (nrTimesWidth*nrPolarizations*2) << endl;
  // create a 3d image containing the samples in device memory
  cl_image_format sample_image_format;
  sample_image_format.image_channel_order = CL_R;
  sample_image_format.image_channel_data_type = CL_FLOAT;
  deviceSamples = clCreateImage3D (context, CL_MEM_READ_ONLY, &sample_image_format,
					       nrChannels, nrStations, nrTimesWidth*nrPolarizations*2,
					       0 /* image row pitch */, 0 /* image slice pitch */,
					       NULL /* host ptr */, &err);
  if(err != CL_SUCCESS) {
    printf("Error: Failed to create 3d image for the samples: %d!\n", err);
    exit(EXIT_FAILURE);
  }

  *sampler = clCreateSampler(context, CL_FALSE, CL_ADDRESS_NONE, CL_FILTER_NEAREST, &err);
  if (err != CL_SUCCESS) {
    printf("Error: Failed to create sampler: %d!\n", err);
    exit(EXIT_FAILURE);
  }
#endif // USE_TEXTURE

  return deviceSamples;
}

void asyncCopySamples(cl_command_queue commands, cl_mem deviceSamples, size_t samplesSize, void* samples, cl_event* event)
{
  cl_int err;

#if !USE_TEXTURE
  // transfer the samples to the device
  err = clEnqueueWriteBuffer(commands, deviceSamples, CL_TRUE, 0, samplesSize * sizeof(float), samples, 0, NULL, event);
#else
  size_t origin[3] = {0, 0, 0};
  size_t region[3] = {nrChannels, nrStations, nrTimesWidth*nrPolarizations*2};
  err = clEnqueueWriteImage(commands, deviceSamples, CL_FALSE /* !blocking_write */,
			    origin, region, 0 /* input_row_pitch */,
			    0 /* input_slice_pitch */,
			    samples,
			    0 /* num_events_in_wait_list */,
			    NULL /* cl_event *event_wait_list */,
			    event);
#endif // USE_TEXTURE

    if (err != CL_SUCCESS) {
	printf("Error: Failed to write to sample data buffer to device: %d!\n", err);
	exit(EXIT_FAILURE);
    }
}

void loadKernel(cl_context context, cl_device_id device_id, cl_program* program, cl_kernel* kernel, const char* suffix)
{
  cl_int            err;

  ostringstream os;
  os << "correlate_" << cellWidth << "x" << cellHeight << suffix;
  string kernelName = os.str();

  // Load the compute program from disk into a cstring buffer
  string fileName(kernelName);
  fileName.append(".cl");
  printf("Loading program '%s'...\n", fileName.data());
  string* source = openKernel(fileName.data());
  if(!source) {
    printf("Error: Failed to load compute program from file!\n");
    exit(EXIT_FAILURE);
  }

  cerr << "creating program" << endl;

  // Create the compute program from the source buffer
  const char* sourceString = source->data();
  *program = clCreateProgramWithSource(context, 1, (const char **) &sourceString, NULL, &err);
  if (!*program || err != CL_SUCCESS) {
    printf("Error: Failed to create compute program!\n");
    exit(EXIT_FAILURE);
  }
  delete source;

  cerr << "building kernel" << endl;

  // Build the program executable
  err = clBuildProgram(*program, 0, NULL, "-cl-mad-enable -cl-strict-aliasing -cl-fast-relaxed-math", NULL, NULL);
  if (err != CL_SUCCESS) {
    size_t len;
    char buffer[2048];
    printf("Error: Failed to build program executable: %d!\n", err);
    clGetProgramBuildInfo(*program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
    printf("errors: %s\n", buffer);
    exit(EXIT_FAILURE);
  }
  
  // Create the compute kernel from within the program
  *kernel = clCreateKernel(*program, kernelName.data(), &err);
  if (!*kernel || err != CL_SUCCESS) {
      printf("Error: Failed to create compute kernel %s!: %d\n", kernelName.data(), err);
    exit(EXIT_FAILURE);
  }

  cerr << "kernel built successfully" << endl;


  size_t len;
  size_t binSize;

  err = clGetProgramInfo(*program, CL_PROGRAM_BINARY_SIZES, sizeof(binSize), &binSize, &len);
  if (err != CL_SUCCESS) {
      printf("Error: Could not get program info: %d\n", err);
      exit(EXIT_FAILURE);
  }

  cerr << " kernel binary is " << binSize << " bytes" << endl;
#if 0
  char* buffer = new char[binSize*2];

  err = clGetProgramInfo(*program, CL_PROGRAM_BINARIES, binSize, &buffer, &len);
  if (err != CL_SUCCESS) {
      printf("Error: Could not get program info: %d\n", err);
      exit(EXIT_FAILURE);
  }

  cout << buffer << endl;

  delete[] buffer;
#endif
}


int main(int argc, char **argv)
{
  cl_int            err;
  cl_device_id      device_id;
  cl_command_queue  commands;
  cl_context        context;
  cl_mem            deviceSamples;
  cl_mem            deviceVisibilities;
  cl_mem            deviceCellToStatX;
  cl_mem            deviceCellToStatY;
  cl_program        program;
  cl_kernel         kernel;
  size_t            local_work_size = 128;
  size_t            global_work_size = 10 * local_work_size;

  const unsigned nrCells = calcNrCells(cellWidth, cellHeight, nrStations);
  const unsigned nrBaselines = nrStations * (nrStations + 1) / 2;    
  const size_t samplesSize = nrChannels*nrStations*nrTimesWidth*nrPolarizations*2; // in floats
  const size_t visibilitiesSize = nrBaselines*nrChannels*nrPolarizations*nrPolarizations*2; // in floats

  cout << "nrCells = " << nrCells << ", nrBaselines = " << nrBaselines 
       << ", samplesSize = " << ((samplesSize*sizeof(float)) / (1024*1024))
       << " MB, visibilitiesSize = " << ((visibilitiesSize*sizeof(float)) / (1024*1024)) << " MB" << endl;

  cout << "work size (global x local) = " << (global_work_size/local_work_size) << " x " << local_work_size << endl;

  cellToStatX = new unsigned[nrCells];
  cellToStatY = new unsigned[nrCells];

  float* samples = new float[samplesSize];
  float* visibilities= new float[visibilitiesSize];
  memset(visibilities, 0, visibilitiesSize * sizeof(float));
  
  // initialize the samples so some pseudo-random values
  for (unsigned channel = 0; channel < nrChannels; channel ++) {
    for (unsigned stat = 0; stat < nrStations; stat ++) {
      for (unsigned time = 0; time < nrTimes; time ++) {
	samples[SAMPLE_INDEX(stat, channel, time, 0, 0)] = time % 8;
	samples[SAMPLE_INDEX(stat, channel, time, 0, 1)] = stat;
	samples[SAMPLE_INDEX(stat, channel, time, 1, 0)] = time % 8 + 1;
	samples[SAMPLE_INDEX(stat, channel, time, 1, 1)] = channel%50;
      }
    }
  }

  // Connect to a compute device
  if(useGPU) {
	  std::cerr << "using GPU" << std::endl;
	  err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);
  } else {
	  std::cerr << "using CPU" << std::endl;
	  err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
  }
  if(err != CL_SUCCESS) {
    printf("Error: Failed to create a device!\n");
    return EXIT_FAILURE;
  }

  // Create a compute context
  context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &err);
  if(!context) {
    printf("Error: Failed to create a compute context!\n");
    return EXIT_FAILURE;
  }

  // Create a command commands
  commands = clCreateCommandQueue(context, device_id, 0, &err);
  if(!commands) {
    printf("Error: Failed to create a command!\n");
    return EXIT_FAILURE;
  }

#if !USE_TEXTURE
  loadKernel(context, device_id, &program, &kernel, "_standard");
#else
  loadKernel(context, device_id, &program, &kernel, "");
#endif // USE_TEXTURE

  cl_sampler sampler;
  deviceSamples = allocateDeviceSamples(context, &sampler, samplesSize);

  // Create the visibilities buffer on the device
  deviceVisibilities = clCreateBuffer(context, CL_MEM_READ_WRITE, visibilitiesSize * sizeof(float), NULL, NULL);
  if (!deviceVisibilities) {
    printf("Error: Failed to allocate visibilities buffer on device!\n");
    return EXIT_FAILURE;
  }

  // Create the cellToStatX buffer on the device
  deviceCellToStatX = clCreateBuffer(context, CL_MEM_READ_WRITE, nrCells * sizeof(unsigned), NULL, NULL);
  if (!deviceCellToStatX) {
    printf("Error: Failed to allocate cellToStatX buffer on device!\n");
    return EXIT_FAILURE;
  }

  // Create the cellToStatY buffer on the device
  deviceCellToStatY = clCreateBuffer(context, CL_MEM_READ_WRITE, nrCells * sizeof(unsigned), NULL, NULL);
  if (!deviceCellToStatY) {
    printf("Error: Failed to allocate cellToStatY buffer on device!\n");
    return EXIT_FAILURE;
  }

  cerr << "all buffers allocated successfully" << endl;

  fillCellToStatTable(cellWidth, cellHeight, nrStations, deviceCellToStatX, deviceCellToStatY, commands);

  cerr << "cellToStat tables initialized successfully" << endl;

  unsigned paramNr=0;
  err =  clSetKernelArg(kernel, paramNr++, sizeof(cl_mem), &deviceSamples);
#if USE_TEXTURE
  err |= clSetKernelArg(kernel, paramNr++, sizeof(cl_sampler), &sampler);
#endif // USE_TEXTURE
  err |= clSetKernelArg(kernel, paramNr++, sizeof(cl_mem), &deviceVisibilities);
  err |= clSetKernelArg(kernel, paramNr++, sizeof(cl_mem), &deviceCellToStatX);
  err |= clSetKernelArg(kernel, paramNr++, sizeof(cl_mem), &deviceCellToStatY);
  err |= clSetKernelArg(kernel, paramNr++, sizeof(unsigned), &nrStations);
  err |= clSetKernelArg(kernel, paramNr++, sizeof(unsigned), &nrTimes);
  err |= clSetKernelArg(kernel, paramNr++, sizeof(unsigned), &nrTimesWidth);
  err |= clSetKernelArg(kernel, paramNr++, sizeof(unsigned), &nrChannels);
  err |= clSetKernelArg(kernel, paramNr++, sizeof(unsigned), &nrCells);
  if (err != CL_SUCCESS) {
    printf("Error: Failed to set kernel arguments: %d!\n", err);
    return EXIT_FAILURE;
  }

  cl_event copyEvent;
  asyncCopySamples(commands, deviceSamples, samplesSize, samples, &copyEvent);
  clWaitForEvents(1, &copyEvent);

  cerr << "running kernel " << iter << " time(s)" << endl;

  clFinish(commands);

  timer kernelTimer("kernel");
  kernelTimer.start();

  for(unsigned i=0; i<iter; i++) {
      // And run the kernel
      err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);
      if (err != CL_SUCCESS) {
	  printf("Error: Failed to execute kernel: %d!\n", err);
	  return EXIT_FAILURE;
      }
      clFinish(commands);
  }

  kernelTimer.stop();
  cout << kernelTimer;

  unsigned long long ops =  (unsigned long long) iter * nrChannels * nrCells * cellWidth * cellHeight
      * nrTimes * 4L * 2L * 4L;

  double perf = (double) ops / kernelTimer.getTimeInSeconds();
  perf /= 1e9; // gigaflops

  cout << "ops = " << ops << endl;

  cout << "achieved " << perf << " gflops" << endl; 

  // Copy back the visibilities
  err = clEnqueueReadBuffer(commands, deviceVisibilities, CL_TRUE, 0, visibilitiesSize * sizeof(float), visibilities, 0, NULL, NULL);
  if (err) {
    printf("Error: Failed to read back visibilities from the device!\n");
    return EXIT_FAILURE;
  }

  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseMemObject(deviceSamples);
  clReleaseMemObject(deviceVisibilities);
  clReleaseCommandQueue(commands);
  clReleaseContext(context);
 
#if USE_TEXTURE
 clReleaseSampler(sampler);
#endif

 if(verify) {
  checkResult(visibilitiesSize, (complex<float>*)samples, (complex<float>*)visibilities, nrBaselines, nrCells);
 }

  delete[] samples;
  delete[] visibilities;
  delete[] cellToStatX;
  delete[] cellToStatY;

  return 0;
}
