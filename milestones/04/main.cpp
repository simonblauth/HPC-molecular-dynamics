#include "atoms.h"
#include "lj_direct_summation.h"
#include "types.h"
#include "verlet.h"
#include "xyz.h"
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    fs::path p = argv[0];
    p = p.parent_path();
    auto traj_path = p;
    traj_path /= "traj.xyz";
    p /= "lj54.xyz";
    std::ofstream traj(traj_path);
    std::cout << "loaded file from: " << p << std::endl;
    auto [names, positions, velocities]{read_xyz_with_velocities(p)};
    Atoms atoms(positions, velocities);

    double mass = 1;
    double sigma = 1;
    double epsilon = 1;
    double timestep = 1e-3 * std::sqrt(mass * sigma * sigma / epsilon);
    size_t max_timesteps = 100000;

    for (size_t ts = 0; ts < max_timesteps; ts++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     mass);
        double epot = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        double ekin = atoms.kinetic_energy();
        if (ts % 100 == 0) {
            std::cout << "Total Energy: " << ekin + epot << ", ";
            std::cout << "Kinetic Energy: " << ekin << ", ";
            std::cout << "Potential Energy: " << epot << std::endl;
            write_xyz(traj, atoms);
            // TODO: plot epot in report
        }
    }
    traj.close();

    return 0;
}