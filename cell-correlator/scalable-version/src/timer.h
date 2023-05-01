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
#ifndef timer_h
#define timer_h

#if defined __cplusplus
extern "C" {
#endif

#if defined __ia64__ && defined __INTEL_COMPILER
#include <ia64regs.h>
#endif


struct timer
{
	
union {
      long long	   total_time;
      struct {
#if defined __PPC__
	int	   total_time_high, total_time_low;
#else
	int	   total_time_low, total_time_high;
#endif
      };
};

#if defined __i386__ && defined __INTEL_COMPILER && defined _OPENMP
    union {
      unsigned long long count;
      struct {
	int count_low, count_high;
      };
    };
#else
    unsigned long long count;
#endif
	
	
};
	
extern float Hz, KHz, MHz, GHz;

__inline__ static void timer_init(struct timer *timer)
{
    timer->total_time = 0;
    timer->count      = 0;
}

__inline__ static void timer_start(struct timer *timer)
  {
#if defined __x86_64__ && defined __INTEL_COMPILER && defined _OPENMP
    asm volatile
    (
	"rdtsc\n\t"
	"shlq $32,%%rdx\n\t"
	"leaq (%%rax,%%rdx),%%rax\n\t"
	"lock;subq %%rax,%0"
    :
	"+m" (timer->total_time)
    :
    :
	"rax", "rdx"
    );
#elif defined __i386__ && defined __INTEL_COMPILER && defined _OPENMP
    asm volatile
    (
	"rdtsc\n\t"
	"lock;subl %%eax,%0\n\t"
	"lock;sbbl %%edx,%1"
    :
	"+m" (timer->total_time_low), "+m" (timer->total_time_high)
    :
    :
	"eax", "edx"
    );
#elif (defined __i386__ || defined __x86_64__) && (defined __GNUC__ || defined __INTEL_COMPILER)
    asm volatile
    (
	"rdtsc\n\t"
	"subl %%eax, %0\n\t"
	"sbbl %%edx, %1"
    :
	"+m" (timer->total_time_low), "+m" (timer->total_time_high)
    :
    :
	"eax", "edx"
    );
#elif (defined __i386__ || defined __x86_64__) && defined __PATHSCALE__
    unsigned eax, edx;

    asm volatile ("rdtsc" : "=a" (eax), "=d" (edx));

    timer->total_time -= ((unsigned long long) edx << 32) + eax;
#elif defined __ia64__ && defined __INTEL_COMPILER
    timer->total_time -= __getReg(_IA64_REG_AR_ITC);
#elif defined __ia64__ && defined __GNUC__
    long long time;
    asm volatile ("mov %0=ar.itc" : "=r" (time));
    timer->total_time -= time;
#elif defined __PPC__ && (defined __GNUC__ || defined __xlC__)
    int high, low, retry;

    asm
    (
	"0:\n\t"
	"mfspr %0,269\n\t"
	"mfspr %1,268\n\t"
	"mfspr %2,269\n\t"
	"cmpw %2,%0\n\t"
	"bne 0b\n\t"
	"subfc %3,%1,%3\n\t"
	"subfe %4,%0,%4"
    :
	"=r" (high), "=r" (low), "=r" (retry),
	"=r" (timer->total_time_low), "=r" (timer->total_time_high)
    :
	"3" (timer->total_time_low), "4" (timer->total_time_high)
    );
#endif
  }


  __inline__ static void timer_stop(struct timer* timer)
  {
#if defined __x86_64__ && defined __INTEL_COMPILER && defined _OPENMP
    asm volatile
    (
	"rdtsc\n\t"
	"shlq $32,%%rdx\n\t"
	"leaq (%%rax,%%rdx),%%rax\n\t"
	"lock;addq %%rax,%0"
    :
	"+m" (timer->total_time)
    :
    :
	"rax", "rdx"
    );
#elif defined __i386__ && defined __INTEL_COMPILER && defined _OPENMP
    asm volatile
    (
	"rdtsc\n\t"
	"lock;addl %%eax, %0\n\t"
	"lock;adcl %%edx, %1"
    :
	"+m" (timer->total_time_low), "+m" (timer->total_time_high)
    :
    :
	"eax", "edx"
    );
#elif (defined __i386__ || defined __x86_64__) && (defined __GNUC__ || defined __INTEL_COMPILER)
    asm volatile
    (
	"rdtsc\n\t"
	"addl %%eax, %0\n\t"
	"adcl %%edx, %1"
    :
	"+m" (timer->total_time_low), "+m" (timer->total_time_high)
    :
    :
	"eax", "edx"
    );
#elif (defined __i386__ || defined __x86_64__) && defined __PATHSCALE__
    unsigned eax, edx;

    asm volatile ("rdtsc\n\t" : "=a" (eax), "=d" (edx));
    timer->total_time += ((unsigned long long) edx << 32) + eax;
#elif defined __ia64__ && defined __INTEL_COMPILER
    timer->total_time += __getReg(_IA64_REG_AR_ITC);
#elif defined __ia64__ && defined __GNUC__
    long long time;
    asm volatile ("mov %0=ar.itc" : "=r" (time));
    timer->total_time += time;
#elif defined __PPC__ && (defined __GNUC__ || defined __xlC__)
    int high, low, retry;

    asm
    (
	"0:\n\t"
	"mfspr %0,269\n\t"
	"mfspr %1,268\n\t"
	"mfspr %2,269\n\t"
	"cmpw %2,%0\n\t"
	"bne 0b\n\t"
	"addc %3,%3,%1\n\t"
	"adde %4,%4,%0"
    :
	"=r" (high), "=r" (low), "=r" (retry),
	"=r" (timer->total_time_low), "=r" (timer->total_time_high)
    :
	"3" (timer->total_time_low), "4" (timer->total_time_high)
    );
#endif

#if defined __x86_64__ && defined __INTEL_COMPILER && defined _OPENMP
    asm volatile ("lock;addq $1,%0" : "+m" (timer->count));
#elif defined __i386__ && defined __INTEL_COMPILER && defined _OPENMP
    asm volatile
    (
	"lock;addl $1,%0\n\t"
	"lock;adcl $0,%1"
    :
	"+m" (timer->count_low), "+m" (timer->count_high)
    );
#else
    ++ timer->count;
#endif
  }

extern void timer_print(struct timer *, const char *name);
extern double timer_time_in_seconds(struct timer *timer);

#if defined __cplusplus
}
#endif

#endif
