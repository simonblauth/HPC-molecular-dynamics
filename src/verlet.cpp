#include "verlet.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep, double mass) {
    vx += fx * timestep / (2 * mass);
    vy += fy * timestep / (2 * mass);
    vz += fz * timestep / (2 * mass);
    x += vx * timestep;
    y += vy * timestep;
    z += vz * timestep;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass) {
    vx += fx * timestep / (2 * mass);
    vy += fy * timestep / (2 * mass);
    vz += fz * timestep / (2 * mass);
}

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep, double mass) {
    velocities += forces * timestep / (2 * mass);
    positions += velocities * timestep;
}

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double timestep, double mass) {
    velocities += forces * timestep / (2 * mass);
}

void verlet_step1(Atoms &atoms, double timestep) {
    atoms.velocities += atoms.forces * timestep / (2 * atoms.get_mass());
    atoms.positions += atoms.velocities * timestep;
}
void verlet_step2(Atoms &atoms, double timestep) {
    atoms.velocities += atoms.forces * timestep / (2 * atoms.get_mass());
}
