#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stddef.h>

#include "../include/fft.h"

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
    for (int k = 0; k < size / 2; k++) {
        *(out + k) = (type) ? sqrt(creal(signal[k]) * creal(signal[k]) + cimag(signal[k]) * cimag(signal[k])) : (double)k * samplingFrequency / size;
    }
    return out;
}