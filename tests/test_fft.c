#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stddef.h>
#include <assert.h>

#include "../include/fft.h"
#include "../src/fft.c"  

#define EPSILON 1e-6

void test_hamming_window() {
    int size = 4;
    complex double signal[] = {1.0, 2.0, 3.0, 4.0};
    complex double expected[] = {
        1.0 * (0.54 - 0.46 * cos(0)),               // w(0)
        2.0 * (0.54 - 0.46 * cos(2 * M_PI / 3)),    // w(1)
        3.0 * (0.54 - 0.46 * cos(4 * M_PI / 3)),    // w(2)
        4.0 * (0.54 - 0.46 * cos(2 * M_PI))         // w(3)
    };

    Hamming(signal, size);

    for (int i = 0; i < size; i++) {
        assert(fabs(creal(signal[i]) - creal(expected[i])) < EPSILON);
    }

    printf("Hamming window test passed.\n");
}

void test_hanning_window() {
    int size = 4;
    complex double signal[] = {1.0, 2.0, 3.0, 4.0};
    complex double expected[] = {
        1.0 * (0.5 * (1 - cos(0))),               // w(0)
        2.0 * (0.5 * (1 - cos(2 * M_PI / 3))),    // w(1)
        3.0 * (0.5 * (1 - cos(4 * M_PI / 3))),    // w(2)
        4.0 * (0.5 * (1 - cos(2 * M_PI)))         // w(3)
    };

    Hanning(signal, size);

    for (int i = 0; i < size; i++) {
        assert(fabs(creal(signal[i]) - creal(expected[i])) < EPSILON);
    }

    printf("Hanning window test passed.\n");
}

void test_blackman_window() {
    int size = 4;
    complex double signal[] = {1.0, 2.0, 3.0, 4.0};
    complex double expected[] = {
        1.0 * (0.42 - 0.5 * cos(0) + 0.08 * cos(0)),                    // w(0)
        2.0 * (0.42 - 0.5 * cos(2 * M_PI / 3) + 0.08 * cos(4 * M_PI / 3)),  // w(1)
        3.0 * (0.42 - 0.5 * cos(4 * M_PI / 3) + 0.08 * cos(8 * M_PI / 3)),  // w(2)
        4.0 * (0.42 - 0.5 * cos(2 * M_PI) + 0.08 * cos(4 * M_PI))           // w(3)
    };

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

