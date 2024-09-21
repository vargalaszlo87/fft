#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


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


void windowHamming(complex double* signal, int size) {
    for (int i = 0; i < size; i++)
        signal[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / (size - 1));
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

int main(int argc, char* argv[]) {

// make signal
	complex double signal[SAMPLE];
	for (int i = 0 ; i < SAMPLE; i++) {
		signal[i] = AMPLITUDE * cexp(I * 2 * M_PI * FREQUENCY * dt * i) + AMPLITUDE * 0.2 * cexp(I * 2 * M_PI * FREQUENCY2 * dt * i);
	}	
	
// use Hamming window and make FFT
	windowHamming(signal, SAMPLE);
	fft(signal, SAMPLE);
	
// show results
    printf("Frequency (Hz) | Amplitude\n");
    for (int k = 0; k < SAMPLE / 2; k++) {  // first N/2
        double frequency = (double)k * FS / SAMPLE;
        double amplitude = sqrt(creal(signal[k]) * creal(signal[k]) + cimag(signal[k]) * cimag(signal[k]));
        printf("%f Hz | %f\n", frequency, amplitude);
    }	
	
	return 0;
}