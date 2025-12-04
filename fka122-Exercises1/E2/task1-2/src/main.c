

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include <gsl/gsl_rng.h>

#define PI 3.141592653589

/*
 * The following struct is defined in the header file
 * Use it to store the result
 */
typedef struct {
    double integral;
    double error;
} result_t;



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
    char filename[256];

    int max_N_sample_exp = 4;
    int N_samples = 0;

    double** results_wo_importance = create_2D_array(max_N_sample_exp, 3);
    double** results_w_importance = create_2D_array(max_N_sample_exp, 3);

    // GSL random number generator object set up 
    const gsl_rng_type * T;
    gsl_rng * r;


    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // typically the seed is set to the current time.
    //time_t seed = time(NULL);
    int seed = 42;
    gsl_rng_set(r, seed); 

    // loop for every different number of samples
    for(int i = 0; i<max_N_sample_exp; i++){
        N_samples = pow(10, i+1);

        result_t result_wo_importance = MC_without_importance_sampling(N_samples, r);
        result_t result_w_importance = MC_with_importance_sampling(N_samples, r);

        results_wo_importance[i][0] = N_samples;
        results_wo_importance[i][1] = result_wo_importance.integral;
        results_wo_importance[i][2] = result_wo_importance.error;

        results_w_importance[i][0] = N_samples;
        results_w_importance[i][1] = result_w_importance.integral;
        results_w_importance[i][2] = result_w_importance.error;




    }

    snprintf(filename, sizeof(filename), "results_wo_importance.csv");

    save_2d_csv(filename,
                results_wo_importance,
                max_N_sample_exp,
                3);


    snprintf(filename, sizeof(filename), "results_w_importance.csv");

    save_2d_csv(filename,
                results_w_importance,
                max_N_sample_exp,
                3);

    return 0;
}