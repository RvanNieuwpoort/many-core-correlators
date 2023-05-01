#SPU_INCLUDE_DIR =	/usr/local/cellsim/sysroot/usr/include
SPU_INCLUDE_DIR =	/opt/cell/sdk/usr/spu/include
SPU_LIB_DIR =		.
SPU_LD_FLAGS =		
SPU_CC =		spu-gcc
SPU_CFLAGS =		-std=gnu99 -W -Wall -Winline -g -O3 -funroll-loops \
                        -fcall-used-80 -fcall-used-81 -fcall-used-82\
                        -fcall-used-83 -fcall-used-84 -fcall-used-85\
                        -fcall-used-86 -fcall-used-87 -fcall-used-88\
                        -fcall-used-89 -fcall-used-90 -fcall-used-91\
                        -fcall-used-92 -fcall-used-93 -fcall-used-94\
                        -fcall-used-95 -fcall-used-96 -fcall-used-97\
                        -fcall-used-98 -fcall-used-99 -fcall-used-100\
                        -fcall-used-101 -fcall-used-102 -fcall-used-103\
                        -fcall-used-104 -fcall-used-105 -fcall-used-106\
                        -fcall-used-107 -fcall-used-108 -fcall-used-109\
                        -fcall-used-110 -fcall-used-111 -fcall-used-112\
                        -fcall-used-113 -fcall-used-114 -fcall-used-115\
                        -fcall-used-116 -fcall-used-117 -fcall-used-118\
                        -fcall-used-119 -fcall-used-120 -fcall-used-121\
                        -fcall-used-122 -fcall-used-123 -fcall-used-123\
                        -fcall-used-124 -fcall-used-125 -fcall-used-126

PPU_CC =		ppu32-gcc
PPU_CFLAGS =		-Wall -std=gnu99 -Winline -m32  -mabi=altivec -maltivec -g -O3
PPU_INCLUDE_DIR =	/opt/cell/sdk/usr/include
PPU_LIB_DIR =		/opt/cell/sdk/usr/lib
#PPU_INCLUDE_DIR =	/usr/local/cellsim/sysroot/usr/include
#PPU_LIB_DIR =		/usr/local/cellsim/sysroot/usr/lib
PPU_LD_FLAGS =		-lspe2 -lpthread -m32 -Wl,-m,elf32ppc


tmp/spu_%.o:	src/spu_%.S
		$(SPU_CC) $(SPU_CFLAGS) -I$(SPU_INCLUDE_DIR) -c $< -o $@

tmp/spu_%.o:	src/spu_%.c
		$(SPU_CC) $(SPU_CFLAGS) -I$(SPU_INCLUDE_DIR) -c $< -o $@

tmp/%.o:	src/%.c
		$(PPU_CC) $(PPU_CFLAGS) -c $< -o $@

tmp/%.obj:	tmp/%
		embedspu -m32 $* $< $@

all:		init bin/correlator

bin/correlator:	tmp/correlator.o tmp/timer.o tmp/spu_correlator.obj 
		$(PPU_CC) $(PPU_CFLAGS) -L$(PPU_LIB_DIR) -R$(PPU_LIB_DIR) $^ $(PPU_LD_FLAGS) -o $@
#		rm -rf tmp

tmp/spu_correlator:	tmp/spu_correlator.o tmp/spu_one_by_one.o tmp/spu_two_by_two.o tmp/spu_three_by_two.o tmp/spu_three_by_three.o tmp/spu_four_by_three.o tmp/spu_decrementer.o
		$(SPU_CC) $(SPU_CFLAGS) -Wl,-N     tmp/spu_correlator.o tmp/spu_one_by_one.o tmp/spu_two_by_two.o tmp/spu_three_by_two.o tmp/spu_three_by_three.o tmp/spu_four_by_three.o tmp/spu_decrementer.o -o $@

init:		clean
		mkdir tmp; mkdir bin

clean::
		$(RM) -rf tmp bin

