#include "atoms.h"
#include "lj_direct_summation.h"
#include "thermostat.h"
#include "types.h"
#include "verlet.h"
#include "xyz.h"
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Please specify the side length of the cubic lattice (in number of atoms) to simulate." << std::endl;
        return 1;
    }
    fs::path p = argv[0];
    auto traj_path = p.parent_path();
    traj_path /= "traj.xyz";
    std::ofstream traj(traj_path);
    

    double mass = 1;
    double sigma = 1;
    double epsilon = 1;
    double timestep = 1e-3 * std::sqrt(mass * sigma * sigma / epsilon);
    size_t max_timesteps = 10000;

    double target_temperaure = 50 * 1e-5;
    double relaxation_time = 1000 * timestep;
    double lattice_distance = 1.05 * sigma;

    size_t nb_atoms_per_lattice = atoi(argv[1]);
    Atoms atoms(std::pow(nb_atoms_per_lattice, 3));
    for (size_t i = 0; i < nb_atoms_per_lattice; i++) {
        for (size_t j = 0; j < nb_atoms_per_lattice; j++) {
            for (size_t k = 0; k < nb_atoms_per_lattice; k++) {
                size_t atom_idx =
                    k + nb_atoms_per_lattice * j +
                    nb_atoms_per_lattice * nb_atoms_per_lattice * i;
                Eigen::Vector3d pos_vec = {(double) i, (double) j, (double) k};
                pos_vec *= lattice_distance;
                atoms.positions.col(atom_idx) = pos_vec;
            }
        }
    }

    for (size_t ts = 0; ts < max_timesteps; ts++) {
        if (ts % 100 == 0) {
            write_xyz(traj, atoms);
        }
        verlet_step1(atoms, timestep);
        double epot = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms, timestep);
        berendsen_thermostat(atoms, target_temperaure, timestep, relaxation_time);
        double ekin = atoms.kinetic_energy();
        if (ts % 100 == 0) {
            std::cout << "Total Energy: " << ekin + epot << ", ";
            std::cout << "Kinetic Energy: " << ekin << ", ";
            std::cout << "Potential Energy: " << epot << std::endl;
        }
    }
    // TODO: equilibrium strategy with changing relaxation time
    traj.close();

    return 0;
}