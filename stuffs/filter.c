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
