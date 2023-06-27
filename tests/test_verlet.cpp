#include "types.h"
#include "verlet.h"
#include <gtest/gtest.h>

TEST(VerletTest, BouncingBall) {
    double vals[9] = {0};
    double timestep = 1e-2;
    for (size_t i = 0; i < 9; i++) {
        EXPECT_EQ(vals[i], 0);
    }
    verlet_step1(vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7],
                 vals[8], timestep);
    for (size_t i = 0; i < 9; i++) {
        EXPECT_EQ(vals[i], 0);
    }
    vals[1] = 50; // set initial height
    vals[7] = -10; // set gravitational force
    int nb_steps = 100000;
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(vals[0], vals[1], vals[2], vals[3], vals[4], vals[5],
                     vals[6], vals[7], vals[8], timestep);
        // gravity stays constant
        verlet_step2(vals[3], vals[4], vals[5], vals[6], vals[7], vals[8],
                     timestep);
        if (vals[1] <= 0) {
            vals[4] *= -1; // bounce off the floor
        }
    }
    double energy_before = 500; // potential energy
    double energy_after = 10 * vals[1] + vals[4] * vals[4] / 2; // pot. + kin.
    EXPECT_NEAR(energy_after, energy_before, 1e-3);
}

TEST(VerletTest, Eigen) {
    double timestep = 1e-2;
    double g = 10;
    size_t nballs = 10;
    Positions_t positions(3, nballs);
    Velocities_t velocities(3, nballs);
    Forces_t forces(3, nballs);
    positions.fill(0);
    velocities.fill(0);
    forces.fill(0);
    forces.row(1).fill(-g); // set gravitational force
    for (size_t i = 0; i < nballs; i++) {
        positions(1, i) = i * 10; // set initial heights
    }
    int nb_steps = 100000;
    auto energies_before{g * positions.row(1)}; // m * g * h
    auto energy_before = energies_before.sum();
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(positions, velocities, forces, timestep);
        // gravity stays constant
        verlet_step2(velocities, forces, timestep);
        for (size_t i = 0; i < nballs; i++) {
            if (positions(1, i) <= 0) {
                velocities(1, i) *= -1; // bounce off the floor
            }
        }
    }
    auto energies_after{g * positions.row(1) + velocities.row(1).pow(2) / 2};
    auto energy_after = energies_after.sum();
    EXPECT_NEAR(energy_after, energy_before, 1e-3);
}
