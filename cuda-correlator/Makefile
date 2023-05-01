NVCC		= nvcc
CC		= g++
#CU_FLAGS	= -O3 -g -deviceemu --host-compilation c++
CU_FLAGS	= -O3 -g --host-compilation c++ --ptxas-options=-v -LIST:source=on --maxrregcount 128 -arch compute_13 -code sm_13
CC_FLAGS	= -O3 -msse3 -Wall -march=nocona -m64 -mfpmath=sse

CU_SOURCES	= main.cu
CC_SOURCES	= timer.cc host_correlator.cc common.cc

CU_OBJECTS	= $(CU_SOURCES:%.cu=%.o)
CU_PTX		= $(CU_SOURCES:%.cu=%.ptx)
CC_OBJECTS	= $(CC_SOURCES:%.cc=%.o)

%.o:		%.cu
		$(NVCC) $(CU_FLAGS) -c $< -o $@

%.o:		%.cc
		$(CC) $(CC_FLAGS) -c $< -o $@

%.ptx:		%.cu
		$(NVCC) $(CU_FLAGS) --ptx $< -o $@

correlator:	$(CU_OBJECTS) $(CC_OBJECTS)
		$(NVCC) $^ -o $@

ptx:		$(CU_PTX) 

clean:
		rm -f *.o main.ptx correlator *~
