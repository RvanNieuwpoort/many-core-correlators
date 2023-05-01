BR_INCLUDE	= /usr/local/atibrook/sdk/include
BR_LIB		= /usr/local/atibrook/sdk/lib
BRCC		= /usr/local/atibrook/sdk/bin/brcc
CC		= g++
BRCC_FLAGS	= -p cpu
CC_FLAGS	= -O2 -Wall -Wfloat-equal -Wpointer-arith  -g3 -ffor-scope

BR_SOURCES	= correlator.br
CC_SOURCES	= main.cc timer.cc

BR_CPP_FILES	= $(BR_SOURCES:%.br=%.cpp)
BR_OBJECTS	= $(BR_SOURCES:%.br=%.o)
CC_OBJECTS	= $(CC_SOURCES:%.cc=%.o)

%.cpp:		%.br
		$(BRCC) $(BRCC_FLAGS) $<

%.o:		%.cc
		$(CC) $(CC_FLAGS) -c $< -o $@

%.o:		%.cpp
		$(CC) $(CC_FLAGS) -I $(BR_INCLUDE) -c $< -o $@

all:		$(CC_OBJECTS) $(BR_OBJECTS) 
		$(CC) $^ -L $(BR_LIB) -lbrook -o correlator

clean:
		rm -f *.o correlator *~ $(BR_CPP_FILES) *_gpu.hpp *.il *.hlsl correlator.h correlator_gpu.h

