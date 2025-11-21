

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"




/**
 * Calculate the potential energy
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  positions - Vector of positions
 *  kappa - Spring constant
 * Returns
 * -------
 *  Potential energy
 *
*/
double calculate_potential_energy(double *positions, double kappa)
{
    const unsigned int N = 3;  // There are three particles

    double potential_energy = 0.0;

    for(unsigned int i = 0; i <N-1 ; i++) {
        // Calculate potential energy contribution from each spring
        potential_energy += 0.5 * kappa * (positions[i+1] - positions[i]) * (positions[i+1] - positions[i]);

    }

    return potential_energy;
}

/**
 * Calculate the kinetic energy
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  velocities - Vector of velocities
 *  masses - Vector of masses
 *
 * Returns
 * -------
 *  Kinetic energy
 *
*/
double calculate_kinetic_energy(double *velocities, double *masses)
{
    const unsigned int N = 3;  // There are three particles

    double kinetic_energy = 0.0;

    for(unsigned int i = 0; i < N; i++) {
        // Calculate kinetic energy contribution from each particle
        kinetic_energy += 0.5 * masses[i] * velocities[i] * velocities[i];
        
    }
    return kinetic_energy;
}




/**
 * Calculate the acceleration.
 *
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  alpha - Anharmonicity constant
 *  N - Number of atoms
 *
*/
void calculate_acceleration(double *accelerations, double *positions, double alpha, const unsigned int N)
{
    // Write to the accelerations vector
    accelerations[0] = positions[1]-2*positions[0]+alpha*((positions[1]-positions[0])*(positions[1]-positions[0])-positions[0]*positions[0]);
    accelerations[N-1] = positions[N-2]-2*positions[N-1]+alpha*(-(positions[N-1]-positions[N-2])*(positions[N-1]-positions[N-2])+positions[N-1]*positions[N-1]);

    double simple_part = 0;
    double complex_part = 0;
    for (unsigned int i = 1; i<N-1; i++){
        simple_part = positions[i+1] - 2*positions[i] + positions[i-1];
        complex_part = alpha * ((positions[i+1]-positions[i])*(positions[i+1]-positions[i])-
                                (positions[i]-positions[i-1])*(positions[i]-positions[i-1]));
        accelerations[i] = simple_part + complex_part;
    }


}

/**
 * Perform one velocity Verlet step
 *
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  alpha - Anharmonicity constant
 *  timestep - Time step
 *  N - Number of atoms
 *
*/
void velocity_verlet_one_step(double *accelerations, double *positions, double *velocities,
                              double alpha, double timestep, const unsigned int N)
{
    // Write to accelerations, positions and velocities vectors
    
    for(unsigned int i = 0; i < N; i++) {

        velocities[i] += 0.5 * accelerations[i] * timestep;
        positions[i] += velocities[i] * timestep;
       
    }
    calculate_acceleration(accelerations, positions, alpha, N);
    for(unsigned int i = 0; i < N; i++) {
        velocities[i] += 0.5 * accelerations[i] * timestep;
    }
}

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

int main()
{
    //gott make this shit work
    //constants
    //const double kappa = 99.864;
    const double timestep = 0.1;
    const unsigned int N = 32;
    const unsigned int steps = 250000;
    const double alpha = 0.1;

    //initialize positions
    double *positions = (double*)calloc(N ,sizeof(double));
    

    //initialize velocities
    double *velocities = (double*)calloc(N,sizeof(double));
    double *P = (double*)calloc(N,sizeof(double));
    P[0] = sqrt(2);
    transform_to_normal_modes(P, velocities, N);
    

    //initialize masses

    


    //initialize accelarations
    double *accelarations = calloc(N, sizeof(double));
    calculate_acceleration(accelarations, positions, alpha, N);

    

    double **energies_time = create_2D_array(steps, N);

    double *energies = (double*)calloc(N, sizeof(double));

    

    for(unsigned int i = 0; i < steps; i++){
        
        calculate_normal_mode_energies(energies, positions, velocities, N);
        
        velocity_verlet_one_step(accelarations, positions, velocities,
        alpha, timestep, N);
        
        
        
        for(unsigned int j = 0; j < N; j++){
            energies_time[i][j] = energies[j];
        }
        

        

    }
    double total_energy = 0;
    for(unsigned int k=0; k<N;k++){
        total_energy += energies[k];

    }
    printf("Total energy: %f\n", total_energy);
    printf("Initial energy: 1\n");
    

    FILE *file_ptr;
    file_ptr = fopen("energies-task9-alpha1.csv", "w");
    if (file_ptr == NULL) {
        perror("Error opening file");
        return 1;
    }
    for(unsigned int i = 0; i < steps; i++){
        fprintf(file_ptr, "%f,%f,%f,%f,%f\n",energies_time[i][0],
                energies_time[i][1], energies_time[i][2],
                energies_time[i][3], energies_time[i][4]);
    }
    fclose(file_ptr);


    /*
    double *positions_1 = malloc(steps * sizeof(double));
    double *positions_2 = malloc(steps * sizeof(double));
    double *positions_3 = malloc(steps * sizeof(double));
    for(unsigned int i = 0; i < steps; i++){
        positions_1[i] = values_time[i][0];
        positions_2[i] = values_time[i][1];
        positions_3[i] = values_time[i][2];

    }
    double *freqs_1 = malloc(steps * sizeof(double));
    double *freqs_2 = malloc(steps * sizeof(double));
    double *freqs_3 = malloc(steps * sizeof(double));
    double *freqs = malloc(steps * sizeof(double));
    powerspectrum(freqs_1, positions_1, steps, timestep);
    powerspectrum(freqs_2, positions_2, steps, timestep);
    powerspectrum(freqs_3, positions_3, steps, timestep);
    fft_freq(freqs, steps, timestep);
    

    FILE *file_ptr;
    file_ptr = fopen("power.csv", "w");
    if (file_ptr == NULL) {
        perror("Error opening file");
        return 1;
    }
    for(unsigned int i = 0; i < steps; i++){
        fprintf(file_ptr, "%f,%f,%f,%f\n",freqs[i],
                freqs_1[i], freqs_2[i], freqs_3[i]);
    }
    fclose(file_ptr);


    file_ptr = fopen("data.csv", "w");
    if (file_ptr == NULL) {
        perror("Error opening file");
        return 1;
    }
    for(unsigned int i = 0; i < steps; i++){
        fprintf(file_ptr, "%f,%f,%f,%f,%f,%f\n",
                values_time[i][0], values_time[i][1], values_time[i][2],
                values_time[i][3], values_time[i][4], values_time[i][5]);
    }
    fclose(file_ptr);

    return 0;
    */
    
}