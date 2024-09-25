#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stddef.h>

#include "../include/fft.h"
#include "../src/fft.c"  

#define EPSILON 1e-6

void test_makeWindow() {
    fftSetup setup = { .samplingPeriod = 1e-6, .samplingSize = 4 };
    complex double signal[] = {1.0, 2.0, 3.0, 4.0};
    complex double expectedHamming[] = {
        1.0 * (0.54 - 0.46 * cos(0)),
        2.0 * (0.54 - 0.46 * cos(2 * M_PI / 3)),
        3.0 * (0.54 - 0.46 * cos(4 * M_PI / 3)),
        4.0 * (0.54 - 0.46 * cos(2 * M_PI))
    };

    makeWindow(signal, &setup, 0);  
    for (int i = 0; i < setup.samplingSize; i++) {
        assert(fabs(creal(signal[i]) - creal(expectedHamming[i])) < EPSILON);
    }

    printf("makeWindow Hamming test passed.\n");

}

void test_makeFFT() {
    fftSetup setup = { .samplingPeriod = 1e-6, .samplingSize = 4 };
    complex double signal[] = {1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I};
    complex double expected[] = {1.0 + 0.0*I, 1.0 + 0.0*I, 1.0 + 0.0*I, 1.0 + 0.0*I};

    makeFFT(signal, &setup);  // FFT végrehajtása
    for (int i = 0; i < setup.samplingSize; i++) {
        assert(fabs(creal(signal[i]) - creal(expected[i])) < EPSILON);
        assert(fabs(cimag(signal[i]) - cimag(expected[i])) < EPSILON);
    }

    printf("makeFFT test passed.\n");
}

void test_getAmplitude() {
    fftSetup setup = { .samplingPeriod = 1e-6, .samplingSize = 4 };
    complex double signal[] = {1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I};
    
    double *amplitudes = getAmplitude(signal, &setup);  
    double expected[] = {1.0, 1.0, 1.0, 1.0};

    for (int i = 0; i < setup.samplingSize; i++) {
        assert(fabs(amplitudes[i] - expected[i]) < EPSILON);
    }

    free(amplitudes); 
    printf("getAmplitude test passed.\n");
}

void test_getFrequency() {
    fftSetup setup = { .samplingPeriod = 1e-6, .samplingSize = 4 };
    complex double signal[] = {1.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I, 0.0 + 0.0*I};
    
    double *frequencies = getFrequency(signal, &setup);  // Frekvencia számítás
    double expected[] = {0.0, 1.0 / (setup.samplingPeriod * setup.samplingSize), 
                         2.0 / (setup.samplingPeriod * setup.samplingSize), 
                         3.0 / (setup.samplingPeriod * setup.samplingSize)};

    for (int i = 0; i < setup.samplingSize; i++) {
        assert(fabs(frequencies[i] - expected[i]) < EPSILON);
    }

    free(frequencies);  
    printf("getFrequency test passed.\n");
}

int main() {
    test_makeWindow();
    test_makeFFT();
    test_getAmplitude();
    test_getFrequency();

    printf("All tests passed.\n");
    return 0;
}