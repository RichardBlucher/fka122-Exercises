#include "attempt.h"
#include <math.h>
#include <stdlib.h>
/*
 * Implement the functions in this file.
 * Utilize the fact that the transformation can be written as
 * a matrix multiplication. To this end, you can use the following
 * functions from C4. No need to declare them, Yata already knows about them
 *

double **create_2D_array(unsigned int n, unsigned int m);

void destroy_2D_array(double **array, unsigned int n);

void matrix_vector_multiplication(double *result, double **A, double *b, unsigned int n, unsigned int m);

 **/

/**
 * Transform to normal modes.
 *
 * Parameters
 * ----------
 *  positions - Vector of positions
 *  Q - Vector of Q-coordinates
 *  N - Number of particles (length of vectors)
 *
*/
void transform_to_normal_modes(double *positions, double *Q, const unsigned int N)
{
    for(unsigned int k = 0; k < N; k++) {
        Q[k] = 0;
        for(unsigned int i = 0; i < N; i++) {
            Q[k] +=  positions[i]*sin((M_PI*(i+1)*(k+1))/(N+1));
        }
        Q[k] *= sqrt(2.0/(N+1));
    }
    

}



/**
 * Calculate the normal mode energies
 *
 * Parameters
 * ----------
 *  energies - Vector where energies will be written
 *  positions - Vector of positions
 *  velocities - Vector of positions
 *  N - Number of particles (length of vectors)
 *
*/
void calculate_normal_mode_energies(double *energies, double *positions, double *velocities, const unsigned int N)
{
    double *Q = (double*)calloc(N, sizeof(double));
    double *P = (double*)calloc(N, sizeof(double));
    
    transform_to_normal_modes(positions, Q, N);
    transform_to_normal_modes(velocities, P, N);

    for(unsigned int k=0; k<N;k++){
        double omega = 2 * sin(((k+1)*M_PI)/(2*(N+1)));
        energies[k] = 0.5 * (P[k]*P[k] + omega*omega*Q[k]*Q[k]);
    }
}
