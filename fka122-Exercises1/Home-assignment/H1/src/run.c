#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"

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



int
run(
    int argc,
    char *argv[]
   )
{
    // Write your code here
    // This makes it possible to test
    // 100% of you code
    const unsigned int n_atoms = 256;
    const unsigned int N = 4;
    const double lattice_constant = 4.00;

    printf("Creating fcc lattice with %i atoms\n", n_atoms);

    double **positions = create_2D_array(n_atoms, 3);

    printf("Initializing fcc lattice\n");

    init_fcc(positions, N, lattice_constant);
    printf("Before pot\n");

    double potential = 0;
    printf("efter pot\n");

    potential = get_energy_AL(positions, N * lattice_constant, n_atoms);

    printf("Potential energy: %.10f\n", potential);

    unsigned int n_lattice_consts = 150;

    double **energy_for_const_lattice = create_2D_array(n_lattice_consts, 2);

    for(unsigned int i = 0; i < n_lattice_consts; i++){
        double lattice_const = lattice_constant + 0.001 * i;
        init_fcc(positions, N, lattice_const);
        double energy = get_energy_AL(positions, N * lattice_const, n_atoms);
        energy_for_const_lattice[i][0] = lattice_const;
        energy_for_const_lattice[i][1] = energy;
        printf("Lattice const: %.4f Energy: %.10f\n", lattice_const, energy);
        
    }

    FILE *file_ptr;
    file_ptr = fopen("energies-task1.csv", "w");
    if (file_ptr == NULL) {
        perror("Error opening file");
        return 1;
    }
    for(unsigned int i = 0; i < n_lattice_consts; i++){
        fprintf(file_ptr, "%f,%f\n",energy_for_const_lattice[i][0], energy_for_const_lattice[i][1]);
    }
    fclose(file_ptr);

    destroy_2D_array(positions, n_atoms);
    destroy_2D_array(energy_for_const_lattice, n_lattice_consts);




    return 0;
}
