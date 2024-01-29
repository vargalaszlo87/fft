#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

/* language */
#define EN_FAIL_MEMORY_01 "Memory allocation failed."

#define WINDOW_SIZE 2048 
#define LAST_SUBHARMONIC 15

#ifndef M_PI
#define M_PI 3.1415926535
#endif

/* ABLAKOK */

void windowHamming(double *window, int size) {
	int i = -1;
	while (++i < size)
        window[i] = 0.54 - 0.46 * cos(2.0 * M_PI * i / (size - 1));
}

void windowHanning(double* window, int size) {
	int i = -1;
	while (++i < size)
        window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (size - 1)));
}

void windowBlackman(double* window, int size) {
	int i = -1;
	while (++i < size)
        window[i] = 0.42 - 0.5 * cos(2.0 * M_PI * i / (size - 1)) + 0.08 * cos(4.0 * M_PI * i / (size - 1));
}

void windowBartlett(double* window, int size) {
	int i = -1;
	while (++i < size)
        window[i] = 1.0 - fabs((i - (size - 1) / 2.0) / ((size - 1) / 2.0));
}

/* ALULÁTERSZTÕ SZÛRÕK */

void LPFMovingAVG(double* x, int size, int windowSize) {
    double sum = 0.0;
	int i = -1;
	while (++i < size) {
        sum += x[i];
        if (i >= windowSize) {
            x[i - windowSize / 2] = sum / windowSize;
            sum -= x[i - windowSize];
        }		
	}    
}

int maxIndex(double *a, int n) {
	double temp = a[0];
	int i = -1, j = 0;
	while (++i < n) {
		if (a[i] > temp)
			temp = a[i],
			j = i;
	}
	return j;
}

int minDifferentIndex(double *f, int N, double needle) {
	printf ("-------------------%lf--------------------\n", needle);
	int j = 0;
	double nearest = f[0];
	double minDiff	= fabs(needle - nearest);
	int i = 0;
	while (++i < N / 2) {
		double diff = fabs(needle - f[i]);
		
		if (diff < minDiff) {
			minDiff = diff;
			j = i;
		}
	
	}
	printf ("%d: %lf\n",j, f[j]);
	return j;
}
	
void fft(complex double* x, int N) {
    if (N <= 1) 
        return;

    complex double* even = (complex double *)calloc(N / 2, sizeof(complex double));
    complex double* odd = (complex double *)calloc(N / 2, sizeof(complex double));

    int i = -1;
    while (++i < N / 2 )
    	even[i] = x[i * 2],
    	odd[i] = x[i * 2 + 1];
  	
    fft(even, N / 2);
    fft(odd, N / 2);
    
	i = -1;
    while (++i < N / 2 ) {
        complex double t = cexp(-2.0 * I * M_PI * i / N) * odd[i];
        x[i] = even[i] + t;
        x[i + N / 2] = even[i] - t;
    }
    free(even);
    free(odd);
}

void fftCalcAmplitude(complex double *comp, double *amplitude, int N, int linlog) {
	int i = -1;
	while (++i < N/2)
		amplitude[i] = (linlog) ?  20.0 * log10(cabs(comp[i]) * 2.0 / N) : cabs(comp[i]) * 2.0 / N;
}

void fftCalcFrequency(double *frequency, int N, double samplerate) {
	int i = -1;
	while (++i < N/2)
		frequency[i] = (double)i * samplerate / N;	
}

void fftResultSegmentation(double *frequency, double *amplitude, int N) {
	/* tömb az alharmonikusoknak */
	double s[LAST_SUBHARMONIC/2] = {0};
	/* melyik a fõkomponens */
    int index = maxIndex(amplitude, N/2);
    	
	int i = -1;
	int j = 3;
	while (++i < sizeof(s)/sizeof(double)) {
		double temp = j * frequency[index];
		int nearest = minDifferentIndex(frequency, N, temp);
/*		s[i] = temp;*/

		printf ("%lf - %lf\n", temp, frequency[nearest]);

		j += 2;
	}
	
	
	/* melyik a fõkomponens */
/*    int index = maxIndex(amplitudes, N/2);
  */  
    

	
}

int main() {
    int N = WINDOW_SIZE;

    double* window = calloc(WINDOW_SIZE, sizeof(double));
    if (window == NULL) 
    	fprintf (stderr, "%s\n", EN_FAIL_MEMORY_01),
    	exit(EXIT_FAILURE);
    windowHanning(window, WINDOW_SIZE);

    complex double* x = (complex double*)calloc(N, sizeof(complex double));
    if (x == NULL)
    	fprintf (stderr, "%s\n", EN_FAIL_MEMORY_01),
    	exit(EXIT_FAILURE);

	/* dev: dummy dates */
    double sampleRate = 8000.0;
    double frequency = 50.0;
    double amplitude = 400.0;	
    int i = -1;
     while (++i < N )
        x[i] = amplitude * sin(2.0 * M_PI * frequency * i / sampleRate) + I * window[i] * 0.0;
    /* emd */

	/* make fft */
    fft(x, N);

	/* calc amplitudes and frequencies */
    double* amplitudes = (double*)calloc(N / 2,sizeof(double));
    if (amplitudes == NULL) 
    	fprintf (stderr, "%s\n", EN_FAIL_MEMORY_01),
    	exit(EXIT_FAILURE);	
	
    double* frequencies = (double*)calloc(N / 2,sizeof(double));
    if (frequencies == NULL)
    	fprintf (stderr, "%s\n", EN_FAIL_MEMORY_01),
    	exit(EXIT_FAILURE);
    
    fftCalcAmplitude(x, amplitudes, N, 0);
    fftCalcFrequency(frequencies, N, sampleRate);
    
   /* segmentation */
   
   fftResultSegmentation(frequencies, amplitudes, N);
    
    
    /* kiirat */

    for (i = 0; i < N / 2; i++) {
    	/*double frequency = (double)i * sampleRate / N;*/
        printf("%lf - Amplitude[%d] = %.2f\n", frequencies[i], i, amplitudes[i]); 
    }

	free(x);
    free(amplitudes);

    return 0;
}

