#include <check.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include "../src/fft.c"

// Segédfüggvény a komplex számok összehasonlítására
int compare_complex(complex double a, complex double b, double tolerance) {
    return (cabs(a - b) < tolerance);
}

// Teszt a 2 pontos FFT futtatására
START_TEST(test_fft_simple)
{
    // Egyszerű 2 pontos bemenet: [1 + 0i, 1 + 0i]
    complex double x[] = {1.0 + 0.0*I, 1.0 + 0.0*I};
    int size = 2;

    // Elvárt eredmény: [2 + 0i, 0 + 0i]
    complex double expected[] = {2.0 + 0.0*I, 0.0 + 0.0*I};

    // FFT futtatása
    fft(x, size);

    // Eredmények ellenőrzése
    double tolerance = 1e-6;
    for (int i = 0; i < size; i++) {
        ck_assert(compare_complex(x[i], expected[i], tolerance));
    }
}
END_TEST

Suite* fft_suite(void) {
    Suite *s;
    TCase *tc_core;

    s = suite_create("FFT");
    tc_core = tcase_create("Core");

    // FFT teszt hozzáadása
    tcase_add_test(tc_core, test_fft_simple);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void) {
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = fft_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? 0 : 1;
}