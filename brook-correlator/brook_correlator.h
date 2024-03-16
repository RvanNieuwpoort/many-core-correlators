#ifndef _CORRELATOR_H
#define _CORRELATOR_H

extern void correlate_1x1(float* baselineToStat1Host, float* baselineToStat2Host, float* hostSamples, float* hostVisReal, float* hostVisImag);
extern void correlate_1x1_vec(float* baselineToStat1Host, float* baselineToStat2Host, float* hostSamples, float* hostVisReal, float* hostVisImag);

extern void correlate_2x2_vec(float* cellToStatXHost, float* cellToStatYHost, float* hostSamples, float* hostVisReal, float* hostVisImag);

extern void startLoadTimer();
extern void stopLoadTimer();

extern void startCorrelateTimer();
extern void stopCorrelateTimer();

extern void startStoreTimer();
extern void stopStoreTimer();

#endif // _CORRELATOR_H
