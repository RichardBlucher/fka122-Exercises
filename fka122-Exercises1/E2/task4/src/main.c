

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
    return 0;
}