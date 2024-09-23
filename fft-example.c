#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stddef.h>

// make a sine-wave
#define SAMPLE	4096
#define FREQUENCY	5e2
#define FREQUENCY2	3.5e2
#define AMPLITUDE	2
#define dt	1e-4

// define the parameters of fft
#define FS		1/dt

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

// windows
void Hamming(complex double* signal, int size) {
	int i = -1;
    while (++i < size)
        signal[i] *= 0.54 - 0.46 * cos(2.0 * M_PI * i / (size - 1));
}

void Hanning(complex double* signal, int size) {
    int i = -1;
    while (++i < size)
        signal[i] *= 0.5 * (1.0 - cos(2.0 * M_PI * i / (size - 1)));
}

void Blackman(complex double* signal, int size) {
    int i = -1;
    while (++i < size)
        signal[i] *= 0.42 - 0.5 * cos(2.0 * M_PI * i / (size - 1)) + 0.08 * cos(4.0 * M_PI * i / (size - 1));
}

void Bartlett(complex double* signal, int size) {
    int i = -1;
    while (++i < size)
        signal[i] *= 1.0 - fabs((i - (size - 1) / 2.0) / ((size - 1) / 2.0));
}

void fft(complex double* x, int size) {
    if (size <= 1)
        return;
	complex double* even = (complex double *)calloc(size / 2, sizeof(complex double));
    complex double* odd = (complex double *)calloc(size / 2, sizeof(complex double));

    int i = -1;
    while (++i < size / 2 )
        even[i] = x[i * 2],
        odd[i] = x[i * 2 + 1];

    fft(even, size / 2);
    fft(odd, size / 2);

    i = -1;
    while (++i < size / 2 ) {
        complex double t = cexp(-2.0 * I * M_PI * i / size) * odd[i];
        x[i] = even[i] + t;
        x[i + size / 2] = even[i] - t;
    }

    free(even);
    free(odd);
}

double* results(complex double* signal, double samplingPeriod, int size, int type) {
    double* out = (double*)calloc(size / 2, sizeof(double));
    if (!out) {
        fprintf(stderr, "Failed memory allocation.\n");
        exit(1);
    }
    double samplingFrequency = 1 / samplingPeriod;
    for (int k = 0; k < SAMPLE / 2; k++) {
        *(out + k) = (type) ? sqrt(creal(signal[k]) * creal(signal[k]) + cimag(signal[k]) * cimag(signal[k])) : (double)k * samplingFrequency / size;
    }
    return out;
}

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

int main(int argc, char* argv[]) {

// make signal
	complex double signal[SAMPLE];
	for (int i = 0 ; i < SAMPLE; i++) {
		signal[i] = AMPLITUDE * cexp(I * 2 * M_PI * FREQUENCY * dt * i) + AMPLITUDE * 0.2 * cexp(I * 2 * M_PI * FREQUENCY2 * dt * i);
	}

// use it

    // setup (sampling period, sample size)
    fftSetup set = {1e-4, 4096};

    // making
	makeWindow(signal, &set, BLACKMAN);
	makeFFT(signal, &set);

	// results
	double* amplitude = getAmplitude(signal, &set);
	double* frequency = getFrequency(signal, &set);

    // show
    printf("Frequency (Hz) | Amplitude\n");
    for (int k = 0; k < SAMPLE / 2; k++) {  // first N/2
        printf("%f Hz | %f\n", frequency[k], amplitude[k]);
    }


    // destroy
    free(&signal);
    free(&set);


	return 0;
}
