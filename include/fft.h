#ifndef FFT_H
#define FFT_H

#ifndef M_PI
#define M_PI 3.1415926538
#endif

typedef struct fft_ {
    double samplingPeriod;
    int samplingSize;
} fftSetup;

enum {
    HAMMING,
    HANNING,
    BLACKMAN,
    BARTLETT,
    enumSize
};

#endif