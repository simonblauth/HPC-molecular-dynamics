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