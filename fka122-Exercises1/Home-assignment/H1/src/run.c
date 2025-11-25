#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void velocity_verlet_one_step(double **forces, double **positions, double **velocities,double mass, double timestep, const unsigned int N, double cell_length)
{
    // Write to accelerations, positions and velocities vectors
    double inverse_mass = 1.0 / mass;
    for(unsigned int i = 0; i < N; i++) {

        velocities[i][0] += 0.5 * forces[i][0] * timestep*inverse_mass;
        velocities[i][1] += 0.5 * forces[i][1] * timestep*inverse_mass;
        velocities[i][2] += 0.5 * forces[i][2] * timestep*inverse_mass;

        positions[i][0] += velocities[i][0] * timestep;
        positions[i][1] += velocities[i][1] * timestep;
        positions[i][2] += velocities[i][2] * timestep;
       
    }
    
    //update forces/accelerations
    get_forces_AL(forces, positions, cell_length, N);

    //update velocities
    for(unsigned int i = 0; i < N; i++) {

        velocities[i][0] += 0.5 * forces[i][0] * timestep*inverse_mass;
        velocities[i][1] += 0.5 * forces[i][1] * timestep*inverse_mass;
        velocities[i][2] += 0.5 * forces[i][2] * timestep*inverse_mass;

    }


    
    
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

void task1(int n_atoms, int N)
{
    // Write your code here
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

    save_2d_csv("energies-task1.csv",
                energy_for_const_lattice,
                n_lattice_consts,
                2);

    destroy_2D_array(positions, n_atoms);
    destroy_2D_array(energy_for_const_lattice, n_lattice_consts);
}

void task2(int n_atoms,int N, double simulation_mass)
{
    // Write your code here
    double lattice_constant = 4.03;
    double total_time = 10;
    char filename[256];


    double timestep = 0;

    for(unsigned int time_step = 0; time_step < 4; time_step++){
        timestep =1/ pow(10, time_step);
        printf("Running simulation with timestep: %.10f\n", timestep);

        unsigned int n_steps = (unsigned int)(total_time / timestep);
        printf("Creating fcc lattice with %i atoms\n", n_atoms);

        double **positions = create_2D_array(n_atoms, 3);

        printf("Initializing fcc lattice\n");

        init_fcc(positions, N, lattice_constant);

        //perturb all positions by +- 0.065 from uniform random distribution
        const gsl_rng_type * T;
        gsl_rng * r;


        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);

        // typically the seed is set to the current time.
        //time_t seed = time(NULL);
        int seed = 42;
        gsl_rng_set(r, seed); 
        for(unsigned int i = 0; i < n_atoms; i++){
            for(unsigned int j = 0; j < 3; j++){
                double rand_perturb = (gsl_rng_uniform(r) * 2 - 1) * 0.065 * lattice_constant;
                positions[i][j] += rand_perturb;
            }
        }
        double **velocities = create_2D_array(n_atoms, 3);
        double **forces = create_2D_array(n_atoms, 3);
        //initialize velocities to zero
        for(unsigned int i = 0; i < n_atoms; i++){
            for(unsigned int j = 0; j < 3; j++){
                velocities[i][j] = 0.0;
            }
        }  
        //get initial forces
        get_forces_AL(forces, positions, N*lattice_constant, n_atoms);
        double **energies_time = create_2D_array(n_steps, 4);
        for(unsigned int step = 0; step < n_steps; step++){
            velocity_verlet_one_step(forces, positions, velocities, simulation_mass, timestep, n_atoms, N*lattice_constant);
            double potential = get_energy_AL(positions, N * lattice_constant, n_atoms);
            double kinetic = 0;
            for(unsigned int i = 0; i < n_atoms; i++){
                kinetic += 0.5 * simulation_mass * (velocities[i][0]*velocities[i][0] +
                                                    velocities[i][1]*velocities[i][1] +
                                                    velocities[i][2]*velocities[i][2]);  
            
            }
            energies_time[step][0] = step * timestep;
            energies_time[step][1] = potential;
            energies_time[step][2] = kinetic;
            energies_time[step][3] = potential + kinetic;

            if(step % (int)pow(10,time_step) == 0){
                printf("Step: %i Time: %.3f Potential: %.10f Kinetic: %.10f Total: %.10f\n",
                        step, step * timestep, potential, kinetic, potential + kinetic);
            }

            //printf("Step: %i Time: %.3f Potential: %.10f Kinetic: %.10f Total: %.10f\n",
                //     step, step * timestep, potential, kinetic, potential + kinetic);
            
            
        }
        //calculate temperature
        double temperature = 0.0;
        for(unsigned int step = 0; step < n_steps; step++){
            double kinetic = energies_time[step][2];
            temperature += (2.0 / 3.0) * (kinetic / n_atoms);
        }
        temperature /= n_steps;
        temperature /= 8.617333262145e-5; //convert eV to Kelvin
        printf("Average temperature: %.2f K\n", temperature);

        snprintf(filename, sizeof(filename), "output_%d.csv", time_step);

        save_2d_csv(filename,
                    energies_time,
                    n_steps,
                    4);
        gsl_rng_free(r);
        destroy_2D_array(energies_time, n_steps);
        destroy_2D_array(forces, n_atoms);
        destroy_2D_array(velocities, n_atoms);




        destroy_2D_array(positions, n_atoms);
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
    
    const double atomic_mass = 26.9815; 
    double simulation_mass = atomic_mass / 9649.0;

    //task1(n_atoms, N);
    task2(n_atoms, N, simulation_mass);




    return 0;
}
