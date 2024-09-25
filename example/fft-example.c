#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stddef.h>

#include "../include/fft.h"
#include "../src/fft.c"
#include "../src/interface.c"

// make a sine-wave
#define SAMPLE	4096
#define FREQUENCY	5e2
#define FREQUENCY2	3.5e2
#define AMPLITUDE	2
#define dt	1e-4

// define the parameters of fft
#define FS		1/dt

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
	makeWindow(signal, &set, BARTLETT);
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
    free(amplitude);
    free(frequency);
    amplitude = NULL;
    frequency = NULL;

	return 0;
}
