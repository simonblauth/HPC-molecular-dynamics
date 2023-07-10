#include "lj_direct_summation.h"
#include "thermostat.h"
#include "verlet.h"
#include <gtest/gtest.h>

TEST(ThermostatTest, ThermostatOnly) {
    constexpr int nb_atoms = 10;

    Atoms atoms(nb_atoms);
    atoms.positions.setRandom();
    atoms.velocities.setRandom();
    for (size_t i = 0; i < 1000; i++) {
        berendsen_thermostat(atoms, 1.0, 0.01, 1);
    }
    EXPECT_NEAR(atoms.current_temperature(), 1.0, 0.01);
}

TEST(ThermostatTest, Simulate) {
    size_t nb_atoms = 10;

    double mass = 1;
    double sigma = 1;
    double epsilon = 1;
    double timestep = 1e-3 * std::sqrt(mass * sigma * sigma / epsilon);
    size_t max_timesteps = 10000;

    double target_temperaure = 42;
    double relaxation_time = 100 * timestep;

    Atoms atoms(nb_atoms);
    atoms.positions.setRandom();
    atoms.velocities.setRandom();
    for (size_t i = 0; i < max_timesteps; i++) {
        verlet_step1(atoms, timestep);
        lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms, timestep);
        berendsen_thermostat(atoms, target_temperaure, timestep, relaxation_time);
    }
    EXPECT_NEAR(atoms.current_temperature(), target_temperaure, 0.01);
}
