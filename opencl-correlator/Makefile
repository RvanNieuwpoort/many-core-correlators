CC		= g++-4.3

#AMD, 64 bit GPU
CC_FLAGS	= -I$(ATISTREAMSDKROOT)/include -O3 -msse3 -Wall -march=nocona -m64 -mfpmath=sse -g
CC_LINK_FLAGS   = -m64 -laticalcl -laticalrt -L$(ATISTREAMSDKROOT)/lib/x86_64 -lOpenCL

# AMD, 64 bit CPU
#CC_FLAGS	= -I$(ATISTREAMSDKROOT)/include -O3 -msse3 -Wall -march=nocona -m64 -mfpmath=sse -g
#CC_LINK_FLAGS   = -m64 -lpthread -L$(ATISTREAMSDKROOT)/lib/x86_64 -lOpenCL

# AMD, 32 bit CPU
#CC_FLAGS	= -I$(ATISTREAMSDKROOT)/include -O3 -msse3 -Wall -mfpmath=sse -g
#CC_LINK_FLAGS   = -lpthread -L$(ATISTREAMSDKROOT)/lib/x86 -lOpenCL

# NVIDIA, 64 bit GPU
#CC_FLAGS	= -I$(NVIDIA_OPENCL)/OpenCL/common/inc -O3 -msse3 -Wall -march=nocona -m64 -mfpmath=sse -g
##CC_FLAGS	= -I$(NVIDIA_OPENCL)/OpenCL/common/inc -Wall -m64 -g
#CC_LINK_FLAGS   = -m64 -lpthread -L$(NVIDIA_OPENCL)/OpenCL/common/lib/Linux64/ -lOpenCL


CC_SOURCES	= correlator.cc timer.cc

CC_OBJECTS	= $(CC_SOURCES:%.cc=%.o)

CC_AS		= $(CC_SOURCES:%.cc=%.s)

%.o:		%.cc
		$(CC) $(CC_FLAGS) -c $< -o $@

%.s:		%.cc
		$(CC) -S $(CC_FLAGS) -c $< -o $@

correlator:	$(CC_OBJECTS)
		$(CC) $(CC_LINK_FLAGS) $^ -o $@

s:		$(CC_AS)
#		$(CC) $(CC_FLAGS) -S cpu_correlator_1x1.cc -o cpu_correlator_1x1.s
#		$(CC) $(CC_FLAGS) -S cpu_correlator_1x1_sse3.cc -o cpu_correlator_1x1_sse3.s

clean:
		rm -f *.o correlator *~ core *.ll
