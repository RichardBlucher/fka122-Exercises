

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"


#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "attempt.h"

#define PI 3.141592653589

/*
 * The following struct is defined in the header file
 * Use it to store the result

typedef struct {
    double integral;
    double error;
} result_t;
 */


/* **********************************************
 *
 * Perform Monte Carlo integration of the given
 * integral without using importance sampling.
 *
 * Parameters
 * ----------
 *  N - Number of points to sample
 *  k - GSL random number generator object
 *
 * Returns
 * -------
 *  Struct with the result
 *
 * **********************************************/
result_t MC_without_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    double sum = 0.0;
    double sum_sq = 0.0;

    for (int i = 0; i < N; i++) {
        double x = gsl_rng_uniform(k); 
        double fx = x*(1-x);              

        sum += fx;
        sum_sq += fx * fx;
    }

    double mean = sum / N;
    double variance = (sum_sq / N) - (mean * mean);
    double error = sqrt(variance / N);

    result.integral = mean;
    result.error = error;

    return result;
}

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "attempt.h"

#define PI 3.141592653589

/*
 * The following struct is defined in the header file
 * Use it to store the result

typedef struct {
    double integral;
    double error;
} result_t;
 */


/* **********************************************
 *
 * Perform Monte Carlo integration of the given
 * integral using importance sampling.
 *
 * Parameters
 * ----------
 *  N - Number of points to sample
 *  k - GSL random number generator object
 *
 * Returns
 * -------
 *  Struct with the result
 *
 * **********************************************/
result_t MC_with_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    double sum = 0.0;
    double sum_sq = 0.0;

    for (int i = 0; i < N; i++) {
        double u = gsl_rng_uniform(k); 
        double x = (2/PI) * asin(sqrt(u));
        double fx = x*(1-x);
        double px = sin(PI * x)*(PI / 2);
        double gx = fx / px;          

        sum += gx;
        sum_sq += gx * gx;
    }

    double mean = sum / N;
    double variance = (sum_sq / N) - (mean * mean);
    double error = sqrt(variance / N);

    result.integral = mean;
    result.error = error;

    return result;
}





int save_2d_csv(const char *filename, double **data,
                unsigned int rows, unsigned int cols)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen");
        return -1;
    }

    for (unsigned int i = 0; i < rows; ++i) {
        for (unsigned int j = 0; j < cols; ++j) {
            fprintf(fp, "%.10e", data[i][j]);
            if (j + 1 < cols) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}


int main()
{
    return 0;
}