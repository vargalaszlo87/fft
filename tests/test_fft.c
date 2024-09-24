// tests/test_fft.c
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <assert.h>
#include "../src/fft.c"

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
    test_fft();
    printf("All tests passed!\n");
    return 0;
}