#ifndef SPU_DECREMENTER_H
#define SPU_DECREMENTER_H

#include <spu_intrinsics.h>

static inline void decrementer_init (void) {
    spu_writech(SPU_WrDec, -1);
}

static inline unsigned int decrementer_get (void) {
    return spu_readch(SPU_RdDec);
}

float decrementer_seconds (unsigned int t);

float decrementer_msecs (unsigned int t);

float decrementer_cycles (unsigned int t);

void decrementer_sleep (unsigned int t);

#endif // SPU_DECREMENTER_H
