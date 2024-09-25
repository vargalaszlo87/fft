#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "../src/fft.c"  

#define EPSILON 1e-6

void test_hamming_window() {
    int size = 4;
    complex double signal[] = {1.0, 2.0, 3.0, 4.0};
    complex double expected[] = {0.54, 1.08, 1.62, 2.16}; 

    Hamming(signal, size);

    for (int i = 0; i < size; i++) {
        assert(fabs(creal(signal[i]) - creal(expected[i])) < EPSILON);
    }

    printf("Hamming window test passed.\n");
}

void test_hanning_window() {
    int size = 4;
    complex double signal[] = {1.0, 2.0, 3.0, 4.0};
    complex double expected[] = {0.0, 1.5, 3.0, 0.0}; 

    Hanning(signal, size);

    for (int i = 0; i < size; i++) {
        assert(fabs(creal(signal[i]) - creal(expected[i])) < EPSILON);
    }

    printf("Hanning window test passed.\n");
}

void test_blackman_window() {
    int size = 4;
    complex double signal[] = {1.0, 2.0, 3.0, 4.0};
    complex double expected[] = {0.42, 1.26, 1.62, 2.16}; 

    Blackman(signal, size);

    for (int i = 0; i < size; i++) {
        assert(fabs(creal(signal[i]) - creal(expected[i])) < EPSILON);
    }

    printf("Blackman window test passed.\n");
}

int compare_complex(complex double a, complex double b, double tolerance) {
    return (cabs(a - b) < tolerance);
}

void test_fft() {
    complex double x[] = {1.0 + 0.0*I, 1.0 + 0.0*I};
    int size = 2;

    complex double expected[] = {2.0 + 0.0*I, 0.0 + 0.0*I};

    fft(x, size);

    double tolerance = 1e-6;
    for (int i = 0; i < size; i++) {
        assert(compare_complex(x[i], expected[i], tolerance));
    }
}


int main() {
    test_hamming_window();
    test_hanning_window();
    test_blackman_window();
    test_fft();

    printf("All tests passed.\n");
    return 0;
}

