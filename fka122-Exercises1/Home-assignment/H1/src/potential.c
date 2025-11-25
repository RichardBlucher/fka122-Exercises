#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "interface.h"

void calculate(double *potential, double *virial, double **force,
           double **position, double cell_length, int nbr_atoms)
{
    /* The interface nep_wrapper is general, however, this interface is not general
     * but made for Al.
     */

    // Setup input variable //
    // box are in the following order
    // [xx yx zx xy yy zy xz yz zz]
    double box[] = {cell_length, 0., 0., 0., cell_length, 0., 0., 0., cell_length};

    // Positions are in the following order
    // [x1 ... xn, y1 ... yn, z1 ... zn]
    double _position[nbr_atoms * 3];
    for(int j = 0; j < 3; ++j){
        for(int i = 0; i < nbr_atoms; ++i){
            _position[i + j * nbr_atoms] = position[i][j];
        }
    }

    // This is a mapping to the symbol of the first atom in nep.txt //
    int type[nbr_atoms];
    memset(type, 0, sizeof(int) * nbr_atoms);


    // Setup return values //
    // setup return for virial
    double _virial[nbr_atoms * 9];
    double _potential[nbr_atoms];
    double _force[nbr_atoms * 3];

    nep_wrapper(nbr_atoms, type, box, _position, _potential, _force, _virial);

    // Handle the copying of the return values to the user //
    // Potential is per atoms
    if (potential != NULL){
        *potential = 0;
        for(int i = 0; i < nbr_atoms; ++i){
            *potential += _potential[i];
        }
    }

    // just copy the forces
    if(force != NULL){
        for(int j = 0; j < 3; ++j){
            for(int i = 0; i < nbr_atoms; ++i){
                force[i][j] = _force[j * nbr_atoms + i];
            }
        }
    }

    if(virial != NULL){
        // Virial is per atoms and order as
        // (I think...)
        // [xx1 ... xxn, yx1, ... yxn, ...]
        double _virial_[9] = {0};
        for(int j = 0; j < 9; ++j){
            for(int i = 0; i < nbr_atoms; ++i){
                _virial_[j] += _virial[j * nbr_atoms + i];
            }
        }
        // For the Al case with isotropic cell we only care
        // about one value
        *virial = _virial_[0];
    }
}

double get_virial_AL(double **position, double cell_length, int nbr_atoms)
{
    double virial = 0;
    calculate(NULL, &virial, NULL, position, cell_length, nbr_atoms);
    return virial;
}

double get_energy_AL(double **position, double cell_length, int nbr_atoms)
{
    double potential = 0;
    calculate(&potential, NULL, NULL, position, cell_length, nbr_atoms);
    return potential;
}

void get_forces_AL(double **force, double **position, double cell_length, int nbr_atoms)
{
    calculate(NULL, NULL, force, position, cell_length, nbr_atoms);
}
