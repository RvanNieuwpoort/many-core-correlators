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
#include "timer.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


unsigned long divisor = 0UL;
float	      Hz, KHz, MHz, GHz;


static void get_MHz(void) {
#if defined __linux__
    FILE *cpuinfo;
    char buffer [256], *ptr;

    if ((cpuinfo = fopen("/proc/cpuinfo", "r")) != 0) {
      while (fgets(buffer, 256, cpuinfo) != 0) {
#if defined __PPC__
        if (strncmp("timebase", buffer, 8) == 0 && (ptr = strchr(buffer, ':')) != 0) {
            Hz = atof(ptr + 2);

            KHz = Hz / 1000.0;
            MHz = Hz / 1000000.0;
            GHz = Hz / 1000000000.0;
            break;
        }
#else
        if (strncmp("cpu MHz", buffer, 7) == 0 && (ptr = strchr(buffer, ':')) != 0) {
            MHz = atof(ptr + 2);
            GHz = MHz / 1000.0;
            KHz = MHz * 1000.0;
            Hz  = MHz * 1000000.0;
            break;
        }
#endif
      }
        fclose(cpuinfo);
    }
#endif
}

static void print_time(const char *name, unsigned long long time) {
    static const char *prefixes [] = { "n", "u", "m", " ", "k", 0
                                     };
    const char	      **prefix	   = prefixes;
    float	      dtime	   = time / GHz;

    while (dtime >= 999.5 && prefix[1] != 0)
        dtime /= 1000, prefix ++;

    printf("%s = %4.3g %ss", name, dtime, *prefix);
}

double timer_time_in_seconds(struct timer *timer) {
    if (MHz == 0)
        get_MHz();

    if (timer->count <= 0) {
        printf("not used\n");
        return 0.0;
    }

    return (timer->total_time_low / GHz) / 1000000000.0;
}

void timer_print(struct timer *timer, const char *timer_name) {
    if (MHz == 0)
        get_MHz();

    printf("%-25s: ", timer_name);

    if (timer->count > 0) {
        unsigned long long avg = timer->total_time / timer->count;

        print_time("avg", avg);
        print_time(", total", timer->total_time);
        printf(", count = %9llu\n", timer->count);
    } else
        printf("not used\n");
}
