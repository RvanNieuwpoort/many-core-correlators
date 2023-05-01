Dear reader,

This is the example source code for a software correlator. 
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
is the most time consuming step in this process. Its computational
complexity is quadratic in the number of receivers. For more
information, see the included paper. If you want to cite it, the
reference is:

Rob V. van Nieuwpoort and John W. Romein:
Using Many-Core Hardware to Correlate Radio Astronomy Signals,
Proceedings of the ACM International Conference on Supercomputing (ICS'09), 
pp. 440-449, June 8-12, 2009, Yorktown Heights, New York, USA. 

The production correlator for LOFAR currently runs on an IBM
BlueGene/P supercomputer. If you would like to have the source code
for that version, please send me an email.  For questions or more
information, feel free to contact me.  Enjoy!

Cheers,

Rob van Nieuwpoort
nieuwpoort@astron.nl or rob@cs.vu.nl or rob@vannieuwpoort.com
http://www.astron.nl/~nieuwpoort or http://www.cs.vu.nl/~rob or
http://www.vannieuwpoort.com
