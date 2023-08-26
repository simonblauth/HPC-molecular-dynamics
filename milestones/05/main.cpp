#include "atoms.h"
#include "lj_direct_summation.h"
#include "thermostat.h"
#include "simulation_utils.h"
#include "types.h"
#include "verlet.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    fs::path filepath = argv[0];
    auto pwd = filepath.parent_path();

    // argument parsing
    argparse::ArgumentParser parser = default_parser("milestone 05");
    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }
    Writer writer(pwd, parser);

    // initialize simulation
    double sigma = parser.get<double>("--sigma");
    double epsilon = parser.get<double>("--epsilon");

    double lattice_distance = parser.get<double>("--lattice_dist") * sigma;
    size_t nb_atoms_per_lattice = parser.get<size_t>("--lattice_size");

    Atoms atoms = init_cubic_lattice(nb_atoms_per_lattice, lattice_distance);

    atoms.mass = parser.get<double>("--mass");
    double timestep = parser.get<double>("--timestep") * std::sqrt(atoms.mass * sigma * sigma / epsilon);
    size_t max_timesteps = parser.get<size_t>("--max_timesteps");

    double target_temperaure = parser.get<double>("--temperature") * 1e-5;
    double relaxation_time = parser.get<size_t>("--relaxation_time") * timestep;

    // simulate
    for (size_t ts = 0; ts < max_timesteps; ts++) {
        writer.write_traj(ts, atoms);
        verlet_step1(atoms, timestep);
        double epot = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms, timestep);
        berendsen_thermostat(atoms, target_temperaure, timestep, relaxation_time);
        double ekin = atoms.kinetic_energy();
        writer.write_energy(ts, ekin, epot);
    }
    // TODO: equilibrium strategy with changing relaxation time

    return 0;
}
