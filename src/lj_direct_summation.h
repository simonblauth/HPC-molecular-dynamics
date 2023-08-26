#ifndef __LJ_DIRECT_SUMMATION_H
#define __LJ_DIRECT_SUMMATION_H

#include "atoms.h"
#include "neighbors.h"

// Force computation with Lennard-Jones potential (https://en.wikipedia.org/wiki/Lennard-Jones_potential). 
// Returns the potential energy of the system.
double lj_direct_summation(Atoms &atoms, double epsilon, double sigma);

// Force computation with Lennard-Jones potential (https://en.wikipedia.org/wiki/Lennard-Jones_potential). 
// Returns the potential energy of the system.
double lj_direct_summation(Atoms &atoms, NeighborList &neighbor_list, double cutoff, double epsilon, double sigma);


#endif  // __LJ_DIRECT_SUMMATION_H
