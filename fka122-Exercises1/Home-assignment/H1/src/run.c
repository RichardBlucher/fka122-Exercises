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
    //pressure *= 160.21766208; //convert to GPa
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

void barostat(double ** positions, double decay_time, double timestep, double target_pressure, double iso_comp, double pressure, int n_atoms, double *lattice_constant){
    double alpha = 1 - (iso_comp * timestep * (target_pressure - pressure))/(decay_time);
    *lattice_constant *= pow(alpha, 1.0/3.0);
    for (unsigned int i = 0; i<n_atoms; i++){
        positions[i][0] *=pow(alpha, 1.0/3.0);
        positions[i][1] *=pow(alpha, 1.0/3.0);
        positions[i][2] *=pow(alpha, 1.0/3.0);
    }
        
}
double MSD(double ** positions, double ** initial_positions, int n_atoms){
    double msd = 0.0;
    for(int i = 0; i < n_atoms; i++){
        /*if(i == 0){
            printf("Initial position atom 0: %f %f %f\n", initial_positions[i][0], initial_positions[i][1], initial_positions[i][2]);
            printf("Current position atom 0: %f %f %f\n", positions[i][0], positions[i][1], positions[i][2]);
        }*/
        double dx = positions[i][0] - initial_positions[i][0];
        double dy = positions[i][1] - initial_positions[i][1];
        double dz = positions[i][2] - initial_positions[i][2];
        msd += dx*dx + dy*dy + dz*dz;
    }
    msd /= n_atoms;
    return msd;
}

double velocity_correlation(double ** velocities, double** initial_velocities, int n_atoms){
    double correlation = 0;
    for(unsigned int i = 0; i<n_atoms; i++){
        double dx = velocities[i][0]*initial_velocities[i][0];
        double dy = velocities[i][1]*initial_velocities[i][1];
        double dz = velocities[i][2]*initial_velocities[i][2];
        correlation += dx+dy+dz;
    }
    return correlation;
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

    for(unsigned int time_step = 0; time_step < 4; time_step++){
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
    const double iso_compress =160.218* 0.01385; //Å^3/eV at 300K
    const double target_temp = 500.0 + 273.15; //Kelvin
    const double target_pressure = 0.0001/160.218; // eV/Å^3
    const double boltzmann = 8.617333262145e-5; //eV
    const double temp_delay_time = 15.0;
    const double pressure_delay_time = 15;

    double lattice_constant = 4.03;
    double total_time = 500+150;
    char filename[256];


    double timestep = 0.005;

    
    printf("Running simulation with timestep: %.10f\n", timestep);

    unsigned int n_steps = (unsigned int)(total_time / timestep);
    printf("Creating fcc lattice with %i atoms\n", n_atoms);

    double **positions = create_2D_array(n_atoms, 3);

    printf("Initializing fcc lattice\n");

    init_fcc(positions, N, lattice_constant);

    save_2d_csv("initial_positions.csv",
                positions,
                n_atoms,
                3);
    
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
    //get initial forces
    get_forces_AL(forces, positions, N*lattice_constant, n_atoms);
    double **simulation_values = create_2D_array(n_steps, 19);


    double **energies_time = create_2D_array(n_steps, 4);
    for(unsigned int step = 0; step < n_steps; step++){
        simulation_values[step][0] = step * timestep;
        simulation_values[step][1] = lattice_constant;
        simulation_values[step][2] = temp;
        simulation_values[step][3] = pressure;
        for(unsigned int i = 0; i < 5; i++){
            simulation_values[step][4 + i*3 + 0] = positions[i][0];
            simulation_values[step][4 + i*3 + 1] = positions[i][1];
            simulation_values[step][4 + i*3 + 2] = positions[i][2];
        }


        velocity_verlet_one_step(forces, positions, velocities, simulation_mass, timestep, n_atoms, N*lattice_constant);
        temp = thermometer(velocities, simulation_mass, boltzmann, n_atoms);
        pressure = barometer(positions, temp, boltzmann, n_atoms, N*lattice_constant);

        //Equilibration
        if(step*timestep < (unsigned int)(150)){
            thermostat(velocities, target_temp, temp, temp_delay_time, timestep, n_atoms);
            barostat(positions, pressure_delay_time, timestep, target_pressure, iso_compress, pressure, n_atoms, &lattice_constant);
        }
        //thermostat(velocities, target_temp, temp, temp_delay_time, timestep, n_atoms);
       // barostat(positions, pressure_delay_time, timestep, target_pressure, iso_compress, pressure, n_atoms, &lattice_constant);
       


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
        

        printf("Step: %i Time: %.3f Kinetic: %f Temperature: %f Pressure %f Lattice constant %f Total: %.10f\n",
                 step, step * timestep, kinetic, temp , pressure,lattice_constant , potential + kinetic);
            
            
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

    snprintf(filename, sizeof(filename), "task3_simulation_values.csv");

    save_2d_csv(filename,
                simulation_values,
                n_steps,
                19);
    //gsl_rng_free(r);
    destroy_2D_array(energies_time, n_steps);
    destroy_2D_array(forces, n_atoms);
    destroy_2D_array(velocities, n_atoms);




    destroy_2D_array(positions, n_atoms);

}

void task4(int n_atoms,int N, double simulation_mass){
    const double iso_compress =160.218* 0.01385; //Å^3/eV at 300K
    const double target_temp = 700.0 + 273.15; //Kelvin
    const double target_pressure = 0.0001/160.218; // eV/Å^3
    const double boltzmann = 8.617333262145e-5; //eV
    const double temp_delay_time = 50.0;
    const double pressure_delay_time = 0.1;

    double lattice_constant = 4.03;
    double total_time = 500;
    char filename[256];


    double timestep = 0.005;

    
    printf("Running simulation with timestep: %.10f\n", timestep);

    unsigned int n_steps = (unsigned int)(total_time / timestep);
    printf("Creating fcc lattice with %i atoms\n", n_atoms);

    double **positions = create_2D_array(n_atoms, 3);

    printf("Initializing fcc lattice\n");

    init_fcc(positions, N, lattice_constant);

    save_2d_csv("initial_positions.csv",
                positions,
                n_atoms,
                3);
    
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
    //get initial forces
    get_forces_AL(forces, positions, N*lattice_constant, n_atoms);
    double **simulation_values = create_2D_array(n_steps, 19);


    double **energies_time = create_2D_array(n_steps, 4);
    for(unsigned int step = 0; step < n_steps; step++){
        simulation_values[step][0] = step * timestep;
        simulation_values[step][1] = lattice_constant;
        simulation_values[step][2] = temp;
        simulation_values[step][3] = pressure;
        for(unsigned int i = 0; i < 5; i++){
            simulation_values[step][4 + i*3 + 0] = positions[i][0];
            simulation_values[step][4 + i*3 + 1] = positions[i][1];
            simulation_values[step][4 + i*3 + 2] = positions[i][2];
        }


        velocity_verlet_one_step(forces, positions, velocities, simulation_mass, timestep, n_atoms, N*lattice_constant);
        temp = thermometer(velocities, simulation_mass, boltzmann, n_atoms);
        pressure = barometer(positions, temp, boltzmann, n_atoms, N*lattice_constant);

        //Equilibration
        thermostat(velocities, target_temp, temp, temp_delay_time, timestep, n_atoms);
        barostat(positions, pressure_delay_time, timestep, target_pressure, iso_compress, pressure, n_atoms, &lattice_constant);
       


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
        

        printf("Step: %i Time: %.3f Kinetic: %f Temperature: %f Pressure %f Lattice constant %f Total: %.10f\n",
                 step, step * timestep, kinetic, temp , pressure,lattice_constant , potential + kinetic);
            
            
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

    snprintf(filename, sizeof(filename), "task4_simulation_values.csv");

    save_2d_csv(filename,
                simulation_values,
                n_steps,
                19);
    //gsl_rng_free(r);
    destroy_2D_array(energies_time, n_steps);
    destroy_2D_array(forces, n_atoms);
    destroy_2D_array(velocities, n_atoms);




    destroy_2D_array(positions, n_atoms);

}

void task4_melt(int n_atoms,int N, double simulation_mass){
    const double iso_compress =160.218* 0.01385; //Å^3/eV at 300K
    double target_temp = 1200.0 + 273.15; //Kelvin
    const double target_pressure = 0.0001/160.218; // eV/Å^3
    const double boltzmann = 8.617333262145e-5; //eV
    const double temp_delay_time = 15.0;
    const double pressure_delay_time = 15.0;

    double lattice_constant = 4.03;
    double total_time = 500+175;
    char filename[256];


    double timestep = 0.005;

    
    printf("Running simulation with timestep: %.10f\n", timestep);

    unsigned int n_steps = (unsigned int)(total_time / timestep);
    printf("Creating fcc lattice with %i atoms\n", n_atoms);

    double **positions = create_2D_array(n_atoms, 3);

    printf("Initializing fcc lattice\n");

    init_fcc(positions, N, lattice_constant);

    save_2d_csv("initial_positions.csv",
                positions,
                n_atoms,
                3);
    
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
    //get initial forces
    get_forces_AL(forces, positions, N*lattice_constant, n_atoms);
    double **simulation_values = create_2D_array(n_steps, 19);


    double **energies_time = create_2D_array(n_steps, 4);
    for(unsigned int step = 0; step < n_steps; step++){
        simulation_values[step][0] = step * timestep;
        simulation_values[step][1] = lattice_constant;
        simulation_values[step][2] = temp;
        simulation_values[step][3] = pressure;
        for(unsigned int i = 0; i < 5; i++){
            simulation_values[step][4 + i*3 + 0] = positions[i][0];
            simulation_values[step][4 + i*3 + 1] = positions[i][1];
            simulation_values[step][4 + i*3 + 2] = positions[i][2];
        }

        /*if((temp > 1500) && (step > 1000)){
            target_temp = 700.0 + 273.15;
            printf("melted");
        }*/
        velocity_verlet_one_step(forces, positions, velocities, simulation_mass, timestep, n_atoms, N*lattice_constant);
        temp = thermometer(velocities, simulation_mass, boltzmann, n_atoms);
        pressure = barometer(positions, temp, boltzmann, n_atoms, N*lattice_constant);

        //Equilibration
        if(step*timestep < (unsigned int)(150)){
            thermostat(velocities, target_temp, temp, temp_delay_time, timestep, n_atoms);
            barostat(positions, pressure_delay_time, timestep, target_pressure, iso_compress, pressure, n_atoms, &lattice_constant);
        }
       


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
        

        printf("Step: %i Time: %.3f Kinetic: %f Temperature: %f Pressure %f Lattice constant %f Total: %.10f\n",
                 step, step * timestep, kinetic, temp , pressure,lattice_constant , potential + kinetic);
            
            
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

    snprintf(filename, sizeof(filename), "task4_high_simulation_values.csv");

    save_2d_csv(filename,
                simulation_values,
                n_steps,
                19);
    //gsl_rng_free(r);
    destroy_2D_array(energies_time, n_steps);
    destroy_2D_array(forces, n_atoms);
    destroy_2D_array(velocities, n_atoms);




    destroy_2D_array(positions, n_atoms);

}

void task5a(int n_atoms,int N, double simulation_mass){
    const double iso_compress =160.218* 0.01385; //Å^3/eV at 300K
    const double target_temp = 500.0 + 273.15; //Kelvin
    const double target_pressure = 0.0001/160.218; // eV/Å^3
    const double boltzmann = 8.617333262145e-5; //eV
    const double temp_delay_time = 50.0;
    const double pressure_delay_time = 0.1;

    double lattice_constant = 4.03;
    double total_time = 500;
    char filename[256];


    double timestep = 0.005;

    
    printf("Running simulation with timestep: %.10f\n", timestep);

    unsigned int n_steps = (unsigned int)(total_time / timestep);
    printf("Creating fcc lattice with %i atoms\n", n_atoms);

    double **positions = create_2D_array(n_atoms, 3);

    printf("Initializing fcc lattice\n");

    init_fcc(positions, N, lattice_constant);

    double **initial_positions = create_2D_array(n_atoms, 3);

    //copy positions to initial_positions
    for(unsigned int i = 0; i < n_atoms; i++){
        for(unsigned int j = 0; j < 3; j++){
            initial_positions[i][j] = positions[i][j];
        }
    }

    save_2d_csv("initial_positions.csv",
                positions,
                n_atoms,
                3);
    double **velocities = create_2D_array(n_atoms, 3);
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
    
    
    //double **velocities = create_2D_array(n_atoms, 3);
    double **forces = create_2D_array(n_atoms, 3);
    //initialize velocities to zero
    /*for(unsigned int i = 0; i < n_atoms; i++){
        for(unsigned int j = 0; j < 3; j++){
            velocities[i][j] = 1.0;
        }
    }
    */
    
    double temp = 0;
    double pressure = 0;
    //get initial forces
    get_forces_AL(forces, positions, N*lattice_constant, n_atoms);

    int sample_rate = 20;   
    int sample_idx = 0;
    double **simulation_values = create_2D_array(n_steps/sample_rate, 2);

    double **energies_time = create_2D_array(n_steps, 4);
    for(unsigned int step = 0; step < n_steps; step++){

        if(step%sample_rate == 0){
            simulation_values[sample_idx][0] = step * timestep;
            simulation_values[sample_idx][1] = MSD(positions, initial_positions, n_atoms);
            sample_idx ++;
        }

        
        //printf("MSD at step %i: %f\n", step, simulation_values[step][1]);
        


        velocity_verlet_one_step(forces, positions, velocities, simulation_mass, timestep, n_atoms, N*lattice_constant);
        temp = thermometer(velocities, simulation_mass, boltzmann, n_atoms);
        pressure = barometer(positions, temp, boltzmann, n_atoms, N*lattice_constant);

        //Equilibration
        thermostat(velocities, target_temp, temp, temp_delay_time, timestep, n_atoms);
        barostat(positions, pressure_delay_time, timestep, target_pressure, iso_compress, pressure, n_atoms, &lattice_constant);
       

        
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
        

        //printf("Step: %i Time: %.3f Kinetic: %f Temperature: %f Pressure %f Lattice constant %f Total: %.10f\n",
                 //step, step * timestep, kinetic, temp , pressure,lattice_constant , potential + kinetic);
            
            
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

    snprintf(filename, sizeof(filename), "task5a_simulation_values.csv");

    save_2d_csv(filename,
                simulation_values,
                n_steps/sample_rate,
                2);
    //gsl_rng_free(r);
    destroy_2D_array(energies_time, n_steps);
    destroy_2D_array(forces, n_atoms);
    destroy_2D_array(velocities, n_atoms);




    destroy_2D_array(positions, n_atoms);

}

void task5b(int n_atoms,int N, double simulation_mass){
    const double iso_compress =160.218* 0.01385; //Å^3/eV at 300K
    double target_temp = 1200.0 + 273.15; //Kelvin
    const double target_temp2 = 700 +273.15; //Kelvin
    const double target_pressure = 0.0001/160.218; // eV/Å^3
    const double boltzmann = 8.617333262145e-5; //eV
    const double temp_delay_time = 50.0;
    const double pressure_delay_time = 0.1;

    double lattice_constant = 4.03;
    double total_time = 500;
    char filename[256];


    double timestep = 0.005;

    
    printf("Running simulation with timestep: %.10f\n", timestep);

    unsigned int n_steps = (unsigned int)(total_time / timestep);
    printf("Creating fcc lattice with %i atoms\n", n_atoms);

    double **positions = create_2D_array(n_atoms, 3);

    printf("Initializing fcc lattice\n");

    init_fcc(positions, N, lattice_constant);

    double **initial_positions = create_2D_array(n_atoms, 3);

    //copy positions to initial_positions
    for(unsigned int i = 0; i < n_atoms; i++){
        for(unsigned int j = 0; j < 3; j++){
            initial_positions[i][j] = positions[i][j];
        }
    }

    save_2d_csv("initial_positions.csv",
                positions,
                n_atoms,
                3);
    double **velocities = create_2D_array(n_atoms, 3);
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
    
    
    //double **velocities = create_2D_array(n_atoms, 3);
    double **forces = create_2D_array(n_atoms, 3);
    //initialize velocities to zero
    /*for(unsigned int i = 0; i < n_atoms; i++){
        for(unsigned int j = 0; j < 3; j++){
            velocities[i][j] = 1.0;
        }
    }
    */
    
    double temp = 0;
    double pressure = 0;
    //get initial forces
    get_forces_AL(forces, positions, N*lattice_constant, n_atoms);

    int sample_rate = 20;   
    int sample_idx = 0;
    double **simulation_values = create_2D_array(n_steps/sample_rate, 2);

    double **energies_time = create_2D_array(n_steps, 4);
    for(unsigned int step = 0; step < n_steps; step++){

        if(step%sample_rate == 0){
            simulation_values[sample_idx][0] = step * timestep;
            simulation_values[sample_idx][1] = MSD(positions, initial_positions, n_atoms);
            sample_idx ++;
        }

        if((temp > 1500) && (step>1000)){
            target_temp = target_temp2;
            printf("melted \n");
        }

        
        //printf("MSD at step %i: %f\n", step, simulation_values[step][1]);
        


        velocity_verlet_one_step(forces, positions, velocities, simulation_mass, timestep, n_atoms, N*lattice_constant);
        temp = thermometer(velocities, simulation_mass, boltzmann, n_atoms);
        pressure = barometer(positions, temp, boltzmann, n_atoms, N*lattice_constant);

        //Equilibration
        thermostat(velocities, target_temp, temp, temp_delay_time, timestep, n_atoms);
        barostat(positions, pressure_delay_time, timestep, target_pressure, iso_compress, pressure, n_atoms, &lattice_constant);
       

        
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
        

        //printf("Step: %i Time: %.3f Kinetic: %f Temperature: %f Pressure %f Lattice constant %f Total: %.10f\n",
                 //step, step * timestep, kinetic, temp , pressure,lattice_constant , potential + kinetic);
            
            
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

    snprintf(filename, sizeof(filename), "task5b_simulation_values.csv");

    save_2d_csv(filename,
                simulation_values,
                n_steps/sample_rate,
                2);
    //gsl_rng_free(r);
    destroy_2D_array(energies_time, n_steps);
    destroy_2D_array(forces, n_atoms);
    destroy_2D_array(velocities, n_atoms);




    destroy_2D_array(positions, n_atoms);

}


void task6(int n_atoms,int N, double simulation_mass){
    const double iso_compress =160.218* 0.01385; //Å^3/eV at 300K
    double target_temp = 1200.0 + 273.15; //Kelvin
    const double target_temp2 = 700+273.15; //Kelvin
    const double target_pressure = 0.0001/160.218; // eV/Å^3
    const double boltzmann = 8.617333262145e-5; //eV
    const double temp_delay_time = 50.0;
    const double pressure_delay_time = 0.1;

    double lattice_constant = 4.03;
    double total_time = 500;
    char filename[256];


    double timestep = 0.005;

    
    printf("Running simulation with timestep: %.10f\n", timestep);

    unsigned int n_steps = (unsigned int)(total_time / timestep);
    printf("Creating fcc lattice with %i atoms\n", n_atoms);

    double **positions = create_2D_array(n_atoms, 3);

    printf("Initializing fcc lattice\n");

    init_fcc(positions, N, lattice_constant);

    double **initial_positions = create_2D_array(n_atoms, 3);

    //copy positions to initial_positions
    for(unsigned int i = 0; i < n_atoms; i++){
        for(unsigned int j = 0; j < 3; j++){
            initial_positions[i][j] = positions[i][j];
        }
    }

    save_2d_csv("initial_positions.csv",
                positions,
                n_atoms,
                3);
    double **velocities = create_2D_array(n_atoms, 3);
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
    
    
    //double **velocities = create_2D_array(n_atoms, 3);
    double **forces = create_2D_array(n_atoms, 3);
    //initialize velocities to zero
    /*for(unsigned int i = 0; i < n_atoms; i++){
        for(unsigned int j = 0; j < 3; j++){
            velocities[i][j] = 1.0;
        }
    }
    */
    
    double temp = 0;
    double pressure = 0;
    //get initial forces
    get_forces_AL(forces, positions, N*lattice_constant, n_atoms);

    //int sample_rate = 20;   
    //int sample_idx = 0;
    double **simulation_values = create_2D_array(n_steps, 2+n_atoms*3);
    double** initial_velocities = create_2D_array(n_atoms, 3);

    double **energies_time = create_2D_array(n_steps, 4);
    for(unsigned int step = 0; step < n_steps; step++){

        

        

        
        //printf("MSD at step %i: %f\n", step, simulation_values[step][1]);
        


        velocity_verlet_one_step(forces, positions, velocities, simulation_mass, timestep, n_atoms, N*lattice_constant);
        temp = thermometer(velocities, simulation_mass, boltzmann, n_atoms);
        pressure = barometer(positions, temp, boltzmann, n_atoms, N*lattice_constant);

        //Equilibration
        thermostat(velocities, target_temp, temp, temp_delay_time, timestep, n_atoms);
        barostat(positions, pressure_delay_time, timestep, target_pressure, iso_compress, pressure, n_atoms, &lattice_constant);

        if(step == 0){
           
            for(unsigned int i = 0; i<n_atoms; i++){
                for(unsigned int j = 0; j<3; j++){
                    initial_velocities[i][j] = velocities[i][j];
                }
            }
        }
        
        simulation_values[step][0] = step * timestep;
        simulation_values[step][1] = velocity_correlation(velocities, initial_velocities, n_atoms);
        for(unsigned int i = 0; i<n_atoms;i++){
            for(unsigned int j = 0; j<3; j++){
                simulation_values[step][2 + i*3 + j] = velocities[i][j];
            }
        }
        
        

        if((temp > 1500) && (step>1000)){
            target_temp = target_temp2;
            printf("melted \n");
        }
        

       

        
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
        

        printf("Step: %i Time: %.3f Kinetic: %f Temperature: %f Pressure %f Lattice constant %f Total: %.10f\n",
                step, step * timestep, kinetic, temp , pressure,lattice_constant , potential + kinetic);
            
            
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

    snprintf(filename, sizeof(filename), "task6_simulation_values.csv");

    save_2d_csv(filename,
                simulation_values,
                n_steps,
                2+n_atoms*3);
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
    //task3(n_atoms, N, simulation_mass);
    //task4(n_atoms, N, simulation_mass);
    task4_melt(n_atoms, N, simulation_mass);
    //task5a(n_atoms, N, simulation_mass);
    //task5b(n_atoms, N, simulation_mass);
    //task6(n_atoms, N, simulation_mass);
    






    return 0;
}
