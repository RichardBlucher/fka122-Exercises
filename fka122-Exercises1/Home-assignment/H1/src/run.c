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

double thermometer(double ** velocities, double mass, double boltzmann, int n_atoms){
    double temp = 0.0;
    for(int i = 0; i < n_atoms; i++){
        temp += velocities[i][0]*velocities[i][0];
        temp += velocities[i][1]*velocities[i][1];
        temp += velocities[i][2]*velocities[i][2];
    }
    temp /= (3*n_atoms*boltzmann);
    temp *= mass;

    return temp;
}

double barometer(double ** positions, double temp, double boltzmann, double n_atoms, double cell_lenght){
    double pressure = 0.0;
    double virial = get_virial_AL(positions, cell_lenght, n_atoms);
    pressure = n_atoms * boltzmann * temp + virial;
    pressure /= (cell_lenght*cell_lenght*cell_lenght);
    pressure *= 160.21766208; //convert to GPa
    return pressure;
}

void thermostat(double ** velocities, double target_temp, double temp, double decay_time, double timestep, int n_atoms){
    double alpha = 1 + (2 * timestep * (target_temp-temp))/(decay_time * temp);
    for(unsigned int i = 0; i< n_atoms; i++){
        velocities[i][0] *=sqrt(alpha);
        velocities[i][1] *=sqrt(alpha);
        velocities[i][2] *=sqrt(alpha);
    }
}

void barostat(double ** positions, double decay_time, double timestep, double target_pressure, double iso_comp, double pressure, int n_atoms){
    double alpha = 1 - (iso_comp * timestep * (target_pressure - pressure))/decay_time;
    for (unsigned int i = 0; i<n_atoms; i++){
        positions[i][0] *=pow(alpha, 1.0/3.0);
        positions[i][1] *=pow(alpha, 1.0/3.0);
        positions[i][2] *=pow(alpha, 1.0/3.0);
    }
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
    double total_time = 100;
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

        snprintf(filename, sizeof(filename), "2_output_%d_2.csv", time_step);

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

void task2_b(int n_atoms,int N, double simulation_mass)
{
    // Write your code here
    double lattice_constant = 4.03;
    double total_time = 100;
    char filename[256];


    double timestep = 0.01;

    for(unsigned int time_step = 0; time_step < 5; time_step++){
        timestep -= 0.002;
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

        snprintf(filename, sizeof(filename), "2b_output_%d.csv", time_step);

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

void task3(int n_atoms,int N, double simulation_mass){
    const double iso_compress = 0.01385; //GPa^-1 at 300K
    const double target_temp = 500.0 + 273.15; //Kelvin
    const double target_pressure = 0.0001; // GPa
    const double boltzmann = 8.617333262145e-5; //eV
    const double temp_delay_time = 1.0;
    const double pressure_delay_time = 1.0;

    double lattice_constant = 4.03;
    double total_time = 15;
    char filename[256];


    double timestep = 0.001;

    
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
    
    double temp = thermometer(velocities, simulation_mass, boltzmann, n_atoms);
    double pressure = barometer(positions, temp, boltzmann, n_atoms, N*lattice_constant);
    double avgDiff = N * lattice_constant;
    //get initial forces
    get_forces_AL(forces, positions, N*lattice_constant, n_atoms);
    double **energies_time = create_2D_array(n_steps, 4);
    for(unsigned int step = 0; step < n_steps; step++){


        velocity_verlet_one_step(forces, positions, velocities, simulation_mass, timestep, n_atoms, avgDiff);
        temp = thermometer(velocities, simulation_mass, boltzmann, n_atoms);
        pressure = barometer(positions, temp, boltzmann, n_atoms, avgDiff);

        //Equilibration
        thermostat(velocities, target_temp, temp, temp_delay_time, timestep, n_atoms);
        barostat(positions, pressure_delay_time, timestep, target_pressure, iso_compress, pressure, n_atoms);

        double minVals[3], maxVals[3];
    double diff[3];
    
    // Initialize min/max with the first row
    for (int col = 0; col < 3; col++) {
        minVals[col] = positions[0][col];
        maxVals[col] = positions[0][col];
    }

    // Scan remaining rows
    for (int i = 1; i < 256; i++) {
        for (int col = 0; col < 3; col++) {
            double val = positions[i][col];
            if (val < minVals[col]) minVals[col] = val;
            if (val > maxVals[col]) maxVals[col] = val;
        }
    }

    // Compute differences for each column
    for (int col = 0; col < 3; col++) {
        diff[col] = maxVals[col] - minVals[col];
    }

    // Compute the average difference
    avgDiff = (diff[0] + diff[1] + diff[2]) / 3.0;

    // Debug prints (optional)
    //printf("Diffs: %f, %f, %f\n", diff[0], diff[1], diff[2]);
    //printf("Average diff: %f\n", avgDiff);


        double potential = get_energy_AL(positions, N * lattice_constant, n_atoms);
        double kinetic = 0;
        for(unsigned int i = 0; i < n_atoms; i++){
            kinetic += 0.5 * (velocities[i][0]*velocities[i][0] + velocities[i][1]*velocities[i][1] + velocities[i][2]*velocities[i][2]) * simulation_mass;  
        
        }
        energies_time[step][0] = step * timestep;
        energies_time[step][1] = potential;
        energies_time[step][2] = kinetic;
        energies_time[step][3] = potential + kinetic;

        
        //printf("Step: %i Time: %.3f Potential: %.10f Kinetic: %.10f Total: %.10f\n", step, step * timestep, potential, kinetic, potential + kinetic);
        

        printf("Step: %i Time: %.3f Kinetic: %f Temperature: %f Pressure %f Total: %.10f\n",
                 step, step * timestep, kinetic, temp , pressure , potential + kinetic);
            
            
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

    snprintf(filename, sizeof(filename), "3_output.csv");

    save_2d_csv(filename,
                energies_time,
                n_steps,
                4);
    //gsl_rng_free(r);
    destroy_2D_array(energies_time, n_steps);
    destroy_2D_array(forces, n_atoms);
    destroy_2D_array(velocities, n_atoms);




    destroy_2D_array(positions, n_atoms);

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
    //task2(n_atoms, N, simulation_mass);
    //task2_b(n_atoms, N, simulation_mass);
    task3(n_atoms, N, simulation_mass);
    





    return 0;
}
