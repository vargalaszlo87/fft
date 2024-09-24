#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.1415926538
#endif

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