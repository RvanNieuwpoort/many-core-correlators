CC		= g++
#CC		= clang++
#CC_FLAGS	= -std=c++17 -O3 -msse3 -Wall -march=nocona -m64 -mfpmath=sse
#CC_FLAGS	= -std=c++17 -O3 -mavx512f -Wall -march=nocona -m64 -mfpmath=sse
CC_FLAGS	= -std=c++20 -O3 -mavx512f -Wall -march=native -m64 -mfpmath=sse
CC_LINK_FLAGS   = -m64 -lpthread

CC_SOURCES	= ../common/common.cc ../common/reference_correlator.cc \
		  ../common/cpu_correlator_1x1.cc ../common/cpu_correlator_2x1.cc \
		  ../common/cpu_correlator_2x2.cc main.cc \
		  cpu_correlator_1x1_sse3.cc cpu_correlator_1x1_time_sse3.cc \
		  cpu_correlator_2x2_sse3.cc cpu_correlator_2x2_time_sse3.cc \
		  cpu_correlator_3x2_time_sse3.cc \
		  cpu_correlator_1x1_time_avx512.cc 

CC_OBJECTS	= $(CC_SOURCES:%.cc=%.o)

CC_AS		= $(CC_SOURCES:%.cc=%.s)

%.o:		%.cc
		$(CC) $(CC_FLAGS) -c $< -o $(notdir $@)

%.s:		%.cc
		$(CC) -S $(CC_FLAGS) -c $< -o $@

correlator:	$(CC_OBJECTS)
		$(CC) $(CC_LINK_FLAGS) *.o -o $@

s:		$(CC_AS)
#		$(CC) $(CC_FLAGS) -S cpu_correlator_1x1.cc -o cpu_correlator_1x1.s
#		$(CC) $(CC_FLAGS) -S cpu_correlator_1x1_sse3.cc -o cpu_correlator_1x1_sse3.s

clean:
		rm -f *.o correlator *~
