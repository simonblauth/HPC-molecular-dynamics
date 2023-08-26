#ifndef __VERLET_H
#define __VERLET_H

#include "atoms.h"
#include "types.h"

// predictor step of the velocity verlet integrator (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet), implemented manually
void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep, double mass=1);
// corrector step of the velocity verlet integrator (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet), implemented manually
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass=1);
// predictor step of the velocity verlet integrator (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet), implemented with Eigen
void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep, double mass=1);
// corrector step of the velocity verlet integrator (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet), implemented with Eigen
void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double timestep, double mass=1);

// predictor step of the velocity verlet integrator (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet), implemented with Eigen
void verlet_step1(Atoms &atoms, double timestep);
// corrector step of the velocity verlet integrator (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet), implemented with Eigen
void verlet_step2(Atoms &atoms, double timestep);

#endif  // __VERLET_H
