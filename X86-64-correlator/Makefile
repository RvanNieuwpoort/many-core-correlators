CC		= g++
ICC		= /opt/intel/cc/10.1.012/bin/icc
CC_FLAGS	= -O3 -msse3 -Wall -march=nocona -m64 -mfpmath=sse -g
#CC_FLAGS	= -O3 -Wall 
ICC_FLAGS	= -O3 -msse3 -Wall -march=nocona
CC_LINK_FLAGS   = -m64 -lpthread

CC_SOURCES	= main.cc cpu_correlator_1x1.cc cpu_correlator_1x1_sse3.cc cpu_correlator_1x1_sse3.cc cpu_correlator_1x1_time_sse3.cc reference_correlator.cc cpu_correlator_2x2_sse3.cc cpu_correlator_2x2_time_sse3.cc cpu_correlator_3x2_time_sse3.cc

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
		rm -f *.o correlator *~
