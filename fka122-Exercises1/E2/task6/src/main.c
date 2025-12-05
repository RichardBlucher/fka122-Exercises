

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include <gsl/gsl_rng.h>

//#include "attempt.h"
#define PI 3.141592653589

/*
 * The following struct is defined in the header file
 * Use it to store the result
*/
typedef struct{
    double probability;
    double function_value;
    int accepted;
} result_t;


/* ***************************************
*
* Perfmore the calculation of the
* unormalized probability function, i.e.,
* the weight used in the MCMC routine
*
* Parameters
* ----------
*  x - current walker position, size = 3
*
* Returns
* -------
* The weight of the current position
*
* ***************************************/
double weight(
              double *x)
{
    double exponent = (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    double weight_value = exp(-exponent);
    return weight_value;
}

/* ***************************************
*
* Perfmore the calculation of the
* function to be sampled
*
* Parameters
* ----------
*  x - current walker position, size = 3
*
* Returns
* -------
* The function value
*
* ***************************************/
double function(
                double *x)
{
    double function_value = x[0]*x[0] + x[1]*x[1]*x[0]*x[0] + x[2]*x[2]*x[0]*x[0]*x[1]*x[1];
    return function_value;
}

/* ***************************************
*
* Perfmore the calculation of the
* function to be sampled
*
* - Use gsl_rng_uniform to displace the
*   walker.
* - The walker should be displaced
*   in the order x,y,z.
* - The draw for the acceptance condition
*   should be done after the draws
*   for the displacing of the walker
*
*
* Parameters
* ----------
*  x - current walker position, size = 3
*  delta - step size
*  k - GSL random number generator object
*
* Returns
* -------
* - The function should update the x parameter
*   to reflect if the move was accepted or
*   rejected.
* - result should contain the probability
*   of the "exiting" x parameter.
* - result should should contain the function
*   value of the "exiting" x parameter.
* - result should be 1 if the move was accepted
*   and 0 if it was rejected.
*
* ***************************************/
 result_t MCMC_step_displace_all(
                                 double *x,
                                 double delta,
                                 gsl_rng *k)
{
    result_t result; 
    double x_new[3];
    for (int i = 0; i < 3; i++){
        double r = gsl_rng_uniform(k);
        x_new[i] = x[i] + (r - 0.5) * delta;
    }
    double pm = weight(x);
    double pt = weight(x_new);

    double q = pt/pm;
    double epsilon = gsl_rng_uniform(k);
    if (epsilon < q){
        result.accepted = 1;
        x[0] = x_new[0];
        x[1] = x_new[1];
        x[2] = x_new[2];
    }
    else{
        result.accepted = 0;
    }


    result.probability = weight(x);
    result.function_value = function(x);
    //result.accepted = 0;

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
    double delta = 2.1;
    const unsigned int n_steps = 100000;
    double **positions = create_2D_array(n_steps, 3);
    double *x = malloc(sizeof(double) * 3);
    int accepted_moves = 0;
    double **function_values = create_2D_array(n_steps, 2);
    gsl_rng *k = gsl_rng_alloc(gsl_rng_mt19937);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    for (unsigned int i = 0; i < n_steps; i++){
        result_t res = MCMC_step_displace_all(x, delta, k);
        positions[i][0] = x[0];
        positions[i][1] = x[1];
        positions[i][2] = x[2];
        function_values[i][0] = res.function_value;
        function_values[i][1] = res.probability;
        accepted_moves += res.accepted;
    }
    printf("Acceptance ratio: %.5f\n", (double)accepted_moves / n_steps);
    double sum = 0.0;
    double sum_sq = 0.0;
    for(unsigned int i = 1000; i < n_steps; i++){
        sum += function_values[i][0];
        sum_sq += function_values[i][0] * function_values[i][0];
    }
    double mean = sum / n_steps;
    double variance = sum_sq / n_steps - (mean*mean);
    double error = sqrt(variance / n_steps);
    printf("Estimated integral: %.5f +/- %.5f\n", mean, error);


    //save_2d_csv("positions.csv", positions, n_steps, 3);
    //save_2d_csv("function_values.csv", function_values, n_steps, 2);

    destroy_2D_array(positions, n_steps);
    destroy_2D_array(function_values, n_steps);
    gsl_rng_free(k);
    free(x);


    return 0;
}