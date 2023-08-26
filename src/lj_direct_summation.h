#ifndef __LJ_DIRECT_SUMMATION_H
#define __LJ_DIRECT_SUMMATION_H

#include "atoms.h"

// Force computation with Lennard-Jones potential (https://en.wikipedia.org/wiki/Lennard-Jones_potential). 
// Returns the potential energy of the system.
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

#endif  // __LJ_DIRECT_SUMMATION_H
