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

void windowGaussian(double* window, int size, double sigma) {
    int i = -1;
    double norm = 0.0;
    while (++i < size)	
        window[i] = exp(-0.5 * pow((i - (size - 1) / 2.0) / (sigma * (size - 1) / 2.0), 2)),
        norm += window[i];
	i = -1;
    while (++i < size)	
        window[i] /= norm;
}
