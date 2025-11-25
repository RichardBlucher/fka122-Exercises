/* compile with USE_TABLE_FOR_RADIAL_FUNCTIONS */
#include "nep.h"
#include <cstring>
#include <iostream>

#ifdef __cplusplus
extern "C"
{
#endif

/* Not something I recommend doing, but in order to avoid recreating
 * the nep object all the time we are going to make it global.
 * Not something I recommend though...
 */

NEP3 nep("nep.txt");

void nep_wrapper(int natoms, int * _type, double *_box, double *_position,
		         double *_potential, double *_force, double *_virial)
{
	std::vector<int> type(natoms);
	memcpy((void*)&type[0], (const void*)_type, natoms * sizeof(int));
	std::vector<double> box(9);
	memcpy((void*)&box[0], (const void*)_box, 9 * sizeof(double));
	std::vector<double> position(3 * natoms);
	memcpy((void*)&position[0], (const void*)_position, natoms * 3 * sizeof(double));

	std::vector<double> potential(1 * natoms);
	std::vector<double> force(3 * natoms);
	std::vector<double> virial(9 * natoms);
	nep.compute((const std::vector<int>)type,
		    (const std::vector<double>)box,
		    (const std::vector<double>)position,
		    potential, force, virial);

	memcpy(_potential, &potential[0], natoms * sizeof(double));
	memcpy(_force, &force[0], natoms * 3 * sizeof(double));
	memcpy(_virial, &virial[0], natoms * 9 * sizeof(double));
}

#ifdef __cplusplus
}
#endif
