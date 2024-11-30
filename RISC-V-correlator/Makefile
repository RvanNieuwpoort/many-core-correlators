CC		= g++
CC_FLAGS	= -O3 -x c++ -march=rv64gcv_zba_zbb_zbs -Wall 

#CC		= clang
#CC_FLAGS	= -O3 -x c++ --target=riscv64 -march=rv64gcv_zba_zbb_zbs -Wall

CC_LINK_FLAGS   = -lpthread

CC_SOURCES	= main.cc riscv_correlator_1x1.cc reference_correlator.cc 

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

clean:
		rm -f *.o correlator *~
