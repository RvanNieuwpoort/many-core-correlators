These are example codes for a software correlator. 
The package includes source code for:

- Intel X86 compatible general purpose platform
- The Sony/Toshiba/IBM Cell/B.E. processor (used in the PlayStation3)
- NVIDIA graphics processors (GPUs)
- ATI graphics processors (GPUs)
- An OpenCL version that is portable to all these platforms.

The correlator is a signal processing application, that multiplies all
pairs of complex input signals. At Astron, we use the correlator for
the LOFAR radio telescope. See www.lofar.org or www.astron.nl. This
telescope consists of tens of thousands of small antennas which are
combined in software to form one large "virtual dish". The correlator
is the most time consuming step in the real-time pipeline. Its computational
complexity is quadratic in the number of receivers. For more
information, see the included paper. If you want to cite it, the
reference is:

Rob V. van Nieuwpoort and John W. Romein:
Correlating Radio Astronomy Signals with Many-Core Hardware
Springer International Journal of Parallel Programming, Special Issue on NY-2009 
Intl. Conf. on Supercomputing, Volume 39, Number 1, 88-114, DOI: 10.1007/s10766-010-0144-3, 2011.

Cheers,

Rob van Nieuwpoort
www.vannieuwpoort.com
