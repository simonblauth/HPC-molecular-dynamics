#ifndef __VERLET_H
#define __VERLET_H

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep, double mass=1);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass=1);

#endif  // __VERLET_H
