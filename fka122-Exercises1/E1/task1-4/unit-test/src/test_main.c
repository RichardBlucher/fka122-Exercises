#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "test_main.h"
#include "tools.h"

/* *************************************
 * The test framework forces you
 * to work with global variables which
 * is a bad pattern.
 *
 * If you are interested in digging
 * deeper into testing, we recommend
 * cmocka whcih is a better designed
 * testing framework. However, it requiers
 * a somewhat more advanced understanding
 * of C.
 * ************************************/


/* ************************************
 * Gloabl variables makes the accesible
 * int the testing scope
 * ************************************/
#define SIZE 4
#define ARRAY_SIZE_N 3
#define ARRAY_SIZE_M 4
#define ARRAY_SIZE_K 5

/* ************************************
 * Here follows the test we want
 * to run
 * ***********************************/
START_TEST(test_elementwise_addition)
{
    // assert re and correct answer is within 0.001f //
    double *v1 = malloc(sizeof(double) * SIZE);
    double *v2 = malloc(sizeof(double) * SIZE);
    double *result_add = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v1[i] = 0.1*i;
        v2[i] = 0.01*i;
    }
    elementwise_addition(result_add, v1, v2, SIZE);
    ck_assert_double_eq_tol(result_add[0], 0, 1e-6);
    ck_assert_double_eq_tol(result_add[1], 0.11, 1e-6);
    ck_assert_double_eq_tol(result_add[2], 0.22, 1e-6);
    ck_assert_double_eq_tol(result_add[3], 0.33, 1e-6);
    free(v1); v1 = NULL;
    free(v2); v2 = NULL;
    free(result_add); result_add = NULL;
}

START_TEST(test_addition_with_constant)
{
    // assert re and correct answer is within 0.001f //
    double *v = malloc(sizeof(double) * SIZE);
    double *result_add = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v[i] = 0.1*i;
    }

    addition_with_constant(result_add, v, 2.2, SIZE);
    ck_assert_double_eq_tol(result_add[0], 2.2, 1e-6);
    ck_assert_double_eq_tol(result_add[1], 2.3, 1e-6);
    ck_assert_double_eq_tol(result_add[2], 2.4, 1e-6);
    ck_assert_double_eq_tol(result_add[3], 2.5, 1e-6);
    free(v); v = NULL;
    free(result_add); result_add = NULL;
}

START_TEST(test_elementwise_multiplication)
{
    // assert re and correct answer is within 0.001f //
    double * v1 = malloc(sizeof(double) * SIZE);
    double * v2 = malloc(sizeof(double) * SIZE);
    double *result_mul = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v1[i] = 0.1*i;
        v2[i] = 0.01*i;
    }
    elementwise_multiplication(result_mul, v1, v2, SIZE);
    ck_assert_double_eq_tol(result_mul[0], 0., 1e-6);
    ck_assert_double_eq_tol(result_mul[1], 0.001, 1e-6);
    ck_assert_double_eq_tol(result_mul[2], 0.004, 1e-6);
    ck_assert_double_eq_tol(result_mul[3], 0.009, 1e-6);
    free(v1); v1 = NULL;
    free(v2); v2 = NULL;
    free(result_mul); result_mul = NULL;
}

START_TEST(test_multiplication_with_constant)
{
    double *v = malloc(sizeof(double) * SIZE);
    double *result_mul = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v[i] = 0.1*i;
    }
    multiplication_with_constant(result_mul, v, 2.2, SIZE);

    ck_assert_double_eq_tol(result_mul[0], 0., 1e-6);
    ck_assert_double_eq_tol(result_mul[1], 0.22, 1e-6);
    ck_assert_double_eq_tol(result_mul[2], 0.44, 1e-6);
    ck_assert_double_eq_tol(result_mul[3], 0.66, 1e-6);
    free(v); v = NULL;
    free(result_mul); result_mul = NULL;
}

START_TEST(test_dot_product)
{
    // assert re and correct answer is within 0.001f //
    double *v1 = malloc(sizeof(double) * SIZE);
    double *v2 = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v1[i] = 0.1*i;
        v2[i] = 0.01*i;
    }
    double result_dot_prd = dot_product(v1, v2, SIZE);
    ck_assert_double_eq_tol(result_dot_prd, 0.014, 1e-6);
    free(v1); v1 = NULL;
    free(v2); v2 = NULL;
}

void print_matrix(
                   double **matrix,
                   unsigned int n,
                   unsigned int m)
{
    printf("===============================\n");
    for(int i = 0; i < n; ++i){
        printf("[");
        for(int k = 0; k < m; ++k){
            printf("%.02f ", matrix[i][k]);
        }
        printf("]\n");
    }
    printf("===============================\n");
}

START_TEST(test_matrix_vector_multiplication)
{
    double *vector_R = (double *) malloc(sizeof(double) * ARRAY_SIZE_N);
    double **matrix_A = create_2D_array(ARRAY_SIZE_N, ARRAY_SIZE_M);
    double *vector_B = (double *) malloc(sizeof(double) * ARRAY_SIZE_M);
    for(int i = 0; i < ARRAY_SIZE_N; ++i){
        // Write some random garbage to the result vector
        vector_R[i] = 0.123 * i - 0.01425 * i * i;
        for(int j = 0; j < ARRAY_SIZE_M; ++j){
            matrix_A[i][j] = i;
        }
    }
    for(int i = 0; i < ARRAY_SIZE_M; ++i){
        vector_B[i] = 2 + i;
    }

    matrix_vector_multiplication(vector_R, matrix_A, vector_B,
                                 ARRAY_SIZE_N, ARRAY_SIZE_M);

    // assert re and correct answer is within 0.001f //

    ck_assert_double_eq_tol(vector_R[0], 0, 1e-6);
    ck_assert_double_eq_tol(vector_R[1], 14.f, 1e-6);
    ck_assert_double_eq_tol(vector_R[2], 28.f, 1e-6);

    destroy_2D_array(matrix_A, ARRAY_SIZE_N); matrix_A = NULL;
    free(vector_B); vector_B = NULL;
    free(vector_R); vector_R = NULL;
}

START_TEST(test_matrix_matrix_multiplication)
{
    double **matrix_R = create_2D_array(ARRAY_SIZE_N, ARRAY_SIZE_K);
    double **matrix_A = create_2D_array(ARRAY_SIZE_N, ARRAY_SIZE_M);
    double **matrix_B = create_2D_array(ARRAY_SIZE_M, ARRAY_SIZE_K);
    for(int i = 0; i < ARRAY_SIZE_N; ++i){
        for(int j = 0; j < ARRAY_SIZE_K; ++j){
            // Write some random garbage to the result vector
            matrix_R[i][j] = 0.123 * i - 0.01425 * i * i + 0.431 * j;
        }
    }
    for(int i = 0; i < ARRAY_SIZE_N; ++i){
        for(int j = 0; j < ARRAY_SIZE_M; ++j){
            matrix_A[i][j] = i;
        }
    }
    for(int i = 0; i < ARRAY_SIZE_M; ++i){
        for(int j = 0; j < ARRAY_SIZE_K; ++j){
            matrix_B[i][j] = j;
        }
    }

    matrix_matrix_multiplication(matrix_R, matrix_A, matrix_B,
                          ARRAY_SIZE_N, ARRAY_SIZE_M, ARRAY_SIZE_K);

    printf("Matrix A:\n");
    print_matrix(matrix_A, ARRAY_SIZE_N, ARRAY_SIZE_M);
    printf("\n");
    printf("Matrix B:\n");
    print_matrix(matrix_B, ARRAY_SIZE_M, ARRAY_SIZE_K);
    printf("\n");
    printf("Matrix R:\n");
    print_matrix(matrix_R, ARRAY_SIZE_N, ARRAY_SIZE_K);

    // assert re and correct answer is within 0.001f //
    ck_assert_double_eq_tol(matrix_R[0][0], 0, 1e-6);
    ck_assert_double_eq_tol(matrix_R[0][1], 0.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[0][2], 0.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[0][3], 0.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[0][4], 0.f, 1e-6);

    ck_assert_double_eq_tol(matrix_R[1][0], 0, 1e-6);
    ck_assert_double_eq_tol(matrix_R[1][1], 4.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[1][2], 8.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[1][3], 12.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[1][4], 16.f, 1e-6);

    ck_assert_double_eq_tol(matrix_R[2][0], 0, 1e-6);
    ck_assert_double_eq_tol(matrix_R[2][1], 8.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[2][2], 16.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[2][3], 24.f, 1e-6);
    ck_assert_double_eq_tol(matrix_R[2][4], 32.f, 1e-6);

    destroy_2D_array(matrix_A, ARRAY_SIZE_N); matrix_A = NULL;
    destroy_2D_array(matrix_B, ARRAY_SIZE_M); matrix_B = NULL;
    destroy_2D_array(matrix_R, ARRAY_SIZE_N); matrix_R = NULL;
}

START_TEST(test_vector_length)
{
    // assert re and correct answer is within 0.001f //
    double *v = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v[i] = 0.1*i;
    }
    double v_norm = vector_norm(v, SIZE);
    ck_assert_double_eq_tol(v_norm, 0.3741657, 1e-6);
    free(v); v = NULL;
}

START_TEST(test_vector_normed)
{
    // assert re and correct answer is within 0.001f //
    double *v_normed = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v_normed[i] = 0.1*i;
    }
    normalize_vector(v_normed, SIZE);
    ck_assert_double_eq_tol(v_normed[0], 0.f, 1e-6);
    ck_assert_double_eq_tol(v_normed[1], 0.26726124f, 1e-6);
    ck_assert_double_eq_tol(v_normed[2], 0.53452248f, 1e-6);
    ck_assert_double_eq_tol(v_normed[3], 0.80178373f, 1e-6);
    free(v_normed); v_normed = NULL;
}

START_TEST(test_average_and_std)
{
    // assert re and correct answer is within 0.001f //
    double *v = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v[i] = 0.01*i;
    }
    double v_average = average(v, SIZE);
    double v_std = standard_deviation(v, SIZE);
    ck_assert_double_eq_tol(v_std, 0.011180339f, 1e-6);
    ck_assert_double_eq_tol(v_average, 0.015f, 1e-6);
    free(v); v = NULL;
}

START_TEST(test_distance_between_vectors)
{
    double *v1 = malloc(sizeof(double) * SIZE);
    double *v2 = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v1[i] = 0.1*i;
        v2[i] = 0.01*i;
    }
    double v1_v2_distance = distance_between_vectors(v1, v2, SIZE);
    ck_assert_double_eq_tol(v1_v2_distance, 0.33674916, 1e-6);
    free(v1); v1 = NULL;
    free(v2); v2 = NULL;
}

START_TEST(test_cumulative_integration)
{
    double *v = malloc(sizeof(double) * SIZE);
    double *res = malloc(sizeof(double) * SIZE);
    for(int i = 0; i < SIZE; i++){
        v[i] = ((double)i / SIZE)*((double)i / SIZE);
    }
    double dx = (double)1 / SIZE;
    cumulative_integration(res, v, dx, SIZE);
    ck_assert_double_eq_tol(res[0], 0, 1e-6);
    ck_assert_double_eq_tol(res[1], 0.0078125, 1e-6);
    ck_assert_double_eq_tol(res[2], 0.046875, 1e-6);
    ck_assert_double_eq_tol(res[3], 0.1484375, 1e-6);
    free(v); free(res);
}

START_TEST(test_write_xyz)
{
    double **positions = create_2D_array(SIZE, 3);
    double **velocities = create_2D_array(SIZE, 3);
    positions[0][0] = 0.;
    positions[0][1] = 0.;
    positions[0][2] = 0.;

    positions[1][0] = 0.;
    positions[1][1] = 2.025;
    positions[1][2] = 2.025;

    positions[2][0] = 2.025;
    positions[2][1] = 2.025;
    positions[2][2] = 0.;

    positions[3][0] = 2.025;
    positions[3][1] = 0.;
    positions[3][2] = 2.025;

    memset(velocities[0], 0, 3 * SIZE * sizeof(double));

    FILE *fp = fopen("atoms.extxyz", "w");
    write_xyz(fp, "Al", positions, velocities, 4.05, SIZE);
    fclose(fp);

    double **read_positions = create_2D_array(SIZE, 3);
    double **read_velocities = create_2D_array(SIZE, 3);
    FILE *fp_read = fopen("atoms.extxyz", "r");
    char read_symbol[3];
    double read_alat;
    read_xyz(fp_read, read_symbol, read_positions, read_velocities, &read_alat);
    fclose(fp_read);

    for(int i = 0; i < SIZE; ++i){
        for(int j = 0; j < 3; ++j){
            ck_assert_double_eq_tol(read_positions[i][j], positions[i][j], 1e-6);
            ck_assert_double_eq_tol(read_velocities[i][j], velocities[i][j], 1e-6);
        }
    }

    destroy_2D_array(positions, SIZE); positions = NULL;
    destroy_2D_array(velocities, SIZE); velocities = NULL;
    destroy_2D_array(read_positions, SIZE); read_positions = NULL;
    destroy_2D_array(read_velocities, SIZE);read_velocities = NULL;
}

START_TEST(test_fourier_transform)
{
    double signal[256];
    double power_spectrum[256];
    double freq[256];
    for(int i = 0; i < sizeof(signal) / sizeof(signal[0]); ++i){
        signal[i] = sin(i);
    }
    powerspectrum(power_spectrum, signal, sizeof(signal) / sizeof(signal[0]), 1);
    fft_freq(freq, sizeof(freq) / sizeof(freq[0]), 1);
    ck_assert_double_eq_tol(freq[41], 1.006291, 1e-6);
    ck_assert_double_eq_tol(freq[215], -1.006291, 1e-6);
    ck_assert_double_eq_tol(power_spectrum[41], 50.978375, 1e-6);
    ck_assert_double_eq_tol(power_spectrum[215], 50.978375, 1e-6);
}

int
main()
{
    test_setup("testing", "core");

    // Tests
    add_test(test_elementwise_addition);
    add_test(test_addition_with_constant);
    add_test(test_elementwise_multiplication);
    add_test(test_multiplication_with_constant);
    add_test(test_dot_product);
    add_test(test_matrix_vector_multiplication);
    add_test(test_matrix_matrix_multiplication);
    add_test(test_vector_length);
    add_test(test_vector_normed);
    add_test(test_average_and_std);
    add_test(test_distance_between_vectors);
    add_test(test_cumulative_integration);
    add_test(test_write_xyz);
    add_test(test_fourier_transform);

    test_teardown();
    return 0;
}
