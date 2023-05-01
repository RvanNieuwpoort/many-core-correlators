CAL_INCLUDE	= /usr/local/ati-stream-sdk-v2.0-beta4-lnx64/include
CAL_LIB		= /usr/local/ati-stream-sdk-v2.0-beta4-lnx64/lib/x86_64
CC		= g++
CC_FLAGS	= -g -Wall -g3

CC_SOURCES	= correlator.cc timer.cc common.cc

CC_OBJECTS	= $(CC_SOURCES:%.cc=%.o)

%.o:		%.cc
		$(CC) $(CC_FLAGS) -I $(CAL_INCLUDE) -c $< -o $@

correlator:	$(CC_OBJECTS)
		$(CC) $^ -L $(CAL_LIB) -l aticalcl -l aticalrt -o $@

clean:
		rm -f *.o correlator core *~
