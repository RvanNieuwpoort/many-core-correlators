/*
Copyright (C) 2009 Rob van Nieuwpoort & John Romein
Astron
P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, nieuwpoort@astron.nl

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
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
