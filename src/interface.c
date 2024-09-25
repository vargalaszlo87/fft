#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stddef.h>

#include "../include/fft.h"

// interfaces
void makeWindow(complex double* signal, fftSetup* setup, int windowType) {
    void (*window[])(complex double*, int) = {Hamming, Hanning, Blackman, Bartlett};
    windowType = (windowType > enumSize - 1 || windowType < 0) ? 0 : windowType;
    window[windowType](signal, setup->samplingSize);
}

void makeFFT(complex double* signal, fftSetup* setup) {
    fft(signal, setup->samplingSize);
}

double *getAmplitude(complex double* signal, fftSetup *s) {
    return results(signal, s->samplingPeriod, s->samplingSize, 1);
}

double *getFrequency(complex double* signal, fftSetup *s) {
    return results(signal, s->samplingPeriod, s->samplingSize, 0);
}