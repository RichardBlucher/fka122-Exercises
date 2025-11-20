/**
 

* Calculate the acceleration.
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector where accelerations are to be written
 *  positions - Vector of positions
 *  masses - Vector of masses
 *  kappa - Spring constant
 *
*/
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
 * Perform one velocity Verlet step
 *
 * All vectors are assumed to be of length 3
 *
 * Parameters
 * ----------
 *  accelerations - Vector of accelerations
 *  positions - Vector of positions
 *  velocities - Vector of velocities
 *  masses - Vector of masses
 *  kappa - Spring constant
 *  timestep - Time step
 *
*/


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

int main()
{
    //gott make this shit work
    //constants
    const double kappa = 99.864;
    const double timestep = 0.0001;
    const unsigned int N = 3;
    const unsigned int steps = 500000;

    //initialize positions
    double *positions = malloc(N * sizeof(double));
    positions[0] = 0.01;
    positions[1] = 0.005;
    positions[2] = -0.005;

    //initialize velocities
    double *velocities = calloc(N,sizeof(double));

    //initialize masses

    double *masses = malloc(N*sizeof(double));
    masses[0] = 0.001658;
    masses[2] = 0.001658;
    masses[1] = 0.00124479;


    //initialize accelarations
    double *accelarations = calloc(N, sizeof(double));
    calculate_acceleration(accelarations, positions, masses, kappa);

    double pot_energy = 0;
    double kin_energy = 0;
    double tot_energy = 0;

    double **values_time = create_2D_array(steps,2* N);

    for(unsigned int i = 0; i < steps; i++){
        velocity_verlet_one_step(accelarations, positions, velocities,
        masses, kappa, timestep);
        values_time[i][0] = positions[0];
        values_time[i][1] = positions[1];
        values_time[i][2] = positions[2];
        values_time[i][3] = velocities[0];
        values_time[i][4] = velocities[1];
        values_time[i][5] = velocities[2];

        pot_energy = calculate_potential_energy(positions, kappa);
        kin_energy = calculate_kinetic_energy(velocities, masses);
        tot_energy = pot_energy + kin_energy;

        printf(" %f %f %f %f\n", i*timestep, pot_energy, kin_energy, tot_energy);

    }
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

}