#include "spu_decrementer.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static float timebase;
static float GHz;

static void get_MHz(void) {
    FILE *cpuinfo;
    char buffer [256], *ptr;

    if ((cpuinfo = fopen("/proc/cpuinfo", "r")) != 0) {
      while (fgets(buffer, 256, cpuinfo) != 0) {
        if (strncmp("timebase", buffer, 8) == 0 && (ptr = strchr(buffer, ':')) != 0) {
            timebase = atof(ptr + 2);
        }
        if (strncmp("clock", buffer, 5) == 0 && (ptr = strchr(buffer, ':')) != 0) {
            GHz = atof(ptr + 2);
        }
      }
        fclose(cpuinfo);
    }

//    printf("SPU timebase = %f, GHz = %f\n", timebase, GHz);
}


float decrementer_seconds (unsigned int t) {
  if(timebase == 0) {
    get_MHz();
  }
    return 1.0 * t / timebase;
}

float decrementer_msecs (unsigned int t) {
    return decrementer_seconds(t) * 1000.0;
}

float decrementer_cycles (unsigned int t) {
  return t * (GHz * 1000000 / timebase);
}

void decrementer_sleep (unsigned int t) {
    unsigned int start = decrementer_get(), target;
    if (start < t) {
        printf("SPU decrementer sleep: Requested: %u, remaining: %u!\n",
               t, start);
        target = 0;
    } else {
        target = start - t;
    }
    while (decrementer_get() > target);
}
