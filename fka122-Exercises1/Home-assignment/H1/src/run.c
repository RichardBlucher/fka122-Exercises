#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"


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
    const double lattice_constant = 1.0;

    double **positions = create_2D_array(n_atoms, 3);

    init_fcc(positions, n_atoms, lattice_constant);

    double potential = 0;

    potential = get_energy_AL(positions, N * lattice_constant, n_atoms);

    printf("Potential energy: %.10f\n", potential);
    




    return 0;
}
