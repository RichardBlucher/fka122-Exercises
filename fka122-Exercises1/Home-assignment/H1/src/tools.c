#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>

#include "tools.h"

void
elementwise_addition(
                     double *res,
                     double *v1,
                     double *v2,
                     unsigned int len
                    )
{
    for (unsigned int i = 0; i < len; ++i){
        res[i] = v1[i] + v2[i];
    }
}

void
elementwise_multiplication(
                           double *res,
                           double *v1,
                           double *v2,
                           unsigned int len
                          )
{
    for (unsigned int i = 0; i < len; ++i){
        res[i] = v1[i] * v2[i];
    }
}

void
addition_with_constant(
                       double *res,
                       double *v,
                       double constant,
                       unsigned int len)
{
    for (unsigned int i = 0; i < len; ++i){
        res[i] = v[i] + constant;
    }
}


void
multiplication_with_constant(
                             double *res,
                             double *v,
                             double constant,
                             unsigned int len)
{
    for (unsigned int i = 0; i < len; ++i){
        res[i] = v[i] * constant;
    }
}

double
dot_product(
            double *v1,
            double *v2,
            unsigned int len
           )
{
    double res = 0;
    for (unsigned int i = 0; i < len; ++i){
        res += v1[i] * v2[i];
    }
    return res;
}

double **
create_2D_array(
                unsigned int row_size,
                unsigned int column_size
               )
{
    if (row_size == 0 || column_size == 0) {
        return NULL;
    }

    double **array2d = malloc(row_size * sizeof(double *));
    if (!array2d) {
        return NULL;
    }

    array2d[0] = calloc(row_size * column_size, sizeof(double));
    if (!array2d[0]) {
        free(array2d);
        return NULL;
    }

    for (unsigned int i = 1; i < row_size; i++) {
        array2d[i] = array2d[0] + i * column_size;
    }

    return array2d;
}

void
destroy_2D_array(
                 double **array,
                 unsigned int n
                )

{
    if (!array) return;
    free(array[0]);   
    free(array);      
}

void
matrix_vector_multiplication(
                             double *result,
                             double **A,
                             double *b,
                             unsigned int n,
                             unsigned int m
                            )
{
    for(unsigned int i = 0; i < n; i++){
        result[i] = 0.;
        for(unsigned int j = 0; j < m; j++){
            result[i] += A[i][j] * b[j];
        }
    }
}

void
matrix_matrix_multiplication(
                             double **result,
                             double **A,
                             double **B,
                             unsigned int n,
                             unsigned int m,
                             unsigned int k
                            )
{
    for(unsigned int i = 0; i < n; i++){
       
        for(unsigned int j = 0; j < k; j++){
            result[i][j] = 0.;
            for(unsigned int p = 0; p < m; p++){
                result[i][j] += A[i][p] * B[p][j];
            }
            
        }
    }
}

double
vector_norm(
            double *v1,
            unsigned int len
           )
{
    double sum = 0.;
    for(unsigned int i = 0; i < len; i++){
        sum += v1[i] * v1[i];
    }
    return sqrt(sum);
}


void
normalize_vector(
                 double *v1,
                 unsigned int len
                )
{
    multiplication_with_constant(v1, v1, 1.0/vector_norm(v1, len), len);
}

double
average(
        double *v1,
        unsigned int len
       )
{
    double sum = 0.;
    for(unsigned int i = 0; i < len; i++){
        sum += v1[i];
    }
    return sum/len;
}


double
standard_deviation(
                       double *v1,
                       unsigned int len
                  )
{
    double avg = average(v1, len);
    double sum = 0.;
    for(unsigned int i = 0; i < len; i++){
        sum += pow(v1[i]-avg, 2);
    }
    return sqrt(sum/len);
}

double
distance_between_vectors(
                         double *v1,
                         double *v2,
                         unsigned int len
                        )
{
    double sum = 0;
    for (unsigned int i = 0; i < len; i++){
        sum += pow(v1[i]-v2[i], 2);
    }
    double root = sqrt(sum);
    return root;
}

void
cumulative_integration(
                       double *res,
                       double *v,
                       double dx,
                       unsigned int v_len
                      )
{
    res[0] = 0.0;

    for (unsigned int i = 1; i < v_len; i++) {
        res[i] = res[i-1] + 0.5 * (v[i] + v[i-1]) * dx;
    }
}


void
write_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double alat,
          int natoms)
{
    fprintf(fp, "%i\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" ", natoms, alat, alat, alat);
    fprintf(fp, "Properties=species:S:1:pos:R:3:vel:R:3 pbc=\"T T T\"\n");
    for(int i = 0; i < natoms; ++i){
        fprintf(fp, "%s %.15f %.15f %.15f %.15f %.15f %.15f\n",
                symbol,
                positions[i][0], positions[i][1], positions[i][2],
                velocities[i][0], velocities[i][1], velocities[i][2]);
    }
}

void fft_freq(
          double *res,
              int n,
              double timestep)
{
    const double factor = 2.0 * M_PI / (n * timestep);

    for (int k = 0; k < n; ++k) {
        int kk = (k <= n/2) ? k : (k - n);  // NumPy frequency index mapping
        res[k] = kk * factor;              // Angular frequency ω = 2π f
    }
}

/* Freely given functions */
void
skip_line(FILE *fp)
{
    int c;
    while (c = fgetc(fp), c != '\n' && c != EOF);
}

void
read_xyz(
         FILE *fp,
         char *symbol,
         double **positions,
         double **velocities,
         double *alat)
{
    int natoms;
    if(fscanf(fp, "%i\nLattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" ", &natoms, alat, alat, alat) == 0){
        perror("Error");
    }
    skip_line(fp);
    for(int i = 0; i < natoms; ++i){
        fscanf(fp, "%s %lf %lf %lf ",
                symbol, &positions[i][0], &positions[i][1], &positions[i][2]);
        fscanf(fp, "%lf %lf %lf\n",
                &velocities[i][0], &velocities[i][1], &velocities[i][2]);
    }
}

void powerspectrum(
           double *res,
           double *signal,
           int n,
                   double timestep)
{
    /* Declaration of variables */
    double *complex_coefficient = malloc(sizeof(double) * 2*n); // array for the complex fft data
    double *data_cp = malloc(sizeof(double) * n);

    /*make copy of data to avoid messing with data in the transform*/
    for (int i = 0; i < n; i++) {
    data_cp[i] = signal[i];
    }

    /* Declare wavetable and workspace for fft */
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;

    /* Allocate space for wavetable and workspace for fft */
    work = gsl_fft_real_workspace_alloc(n);
    real = gsl_fft_real_wavetable_alloc(n);

    /* Do the fft*/
    gsl_fft_real_transform(data_cp, 1, n, real, work);

    /* Unpack the output into array with alternating real and imaginary part */
    gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);

    /*fill the output powspec_data with the powerspectrum */
    for (int i = 0; i < n; i++) {
    res[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1]);
    res[i] *= timestep / n;
    }

    /* Free memory of wavetable and workspace */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
    free(complex_coefficient);
    free(data_cp);
}
