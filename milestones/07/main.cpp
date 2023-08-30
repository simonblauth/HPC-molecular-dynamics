#include "atoms.h"
#include "average.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "simulation_utils.h"
#include "thermostat.h"
#include "types.h"
#include "verlet.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

void deposit_energy(Atoms& atoms, double delta_Q) {
    double ekin = atoms.kinetic_energy();
    double lambda = std::sqrt(1 + delta_Q / ekin);
    atoms.velocities *= lambda;
}

int main(int argc, char *argv[]) {
    // argument parsing
    argparse::ArgumentParser parser = default_parser("milestone 07");
    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }

    fs::path filepath = argv[0];
    auto pwd = filepath.parent_path();
    auto input_path = parser.get<std::string>("--input");
    auto [names, positions]{read_xyz(input_path)};
    std::cout << "loaded file from: " << input_path << std::endl;

    Writer writer(pwd, parser);

    // initialize simulation
    Atoms atoms(names, positions);

    atoms.set_mass(parser.get<double>("--mass") * 103.6);
    double timestep = parser.get<double>("--timestep");
    size_t max_timesteps = parser.get<size_t>("--max_timesteps");
    size_t init_timesteps = parser.get<size_t>("--initial_relaxation");
    size_t relaxation_time = parser.get<size_t>("--relaxation_time");
    double target_temperaure = parser.get<double>("--temperature") * 1e-5;
    size_t relaxation_time_deposit = parser.get<size_t>("--relaxation_time_deposit");
    double delta_Q = parser.get<double>("--deposit_energy");
    double alpha = parser.get<double>("--smoothing");

    double cutoff = parser.get<double>("--cutoff");
    NeighborList neighbor_list(cutoff);

    // relax
    for (size_t i = 0; i < init_timesteps; i++) {
        verlet_step1(atoms, timestep);
        neighbor_list.update(atoms);
        double epot = ducastelle(atoms, neighbor_list, cutoff);
        verlet_step2(atoms, timestep);
        berendsen_thermostat(atoms, target_temperaure, timestep, relaxation_time);
    }
    

    // simulate
    bool relaxed = false;
    ExponentialAverage avg_temp(alpha, atoms.current_temperature_kelvin());
    for (size_t ts = 0; ts < max_timesteps; ts++) {
        writer.write_traj(ts, atoms);
        verlet_step1(atoms, timestep);
        neighbor_list.update(atoms);
        double epot = ducastelle(atoms, neighbor_list, cutoff);
        verlet_step2(atoms, timestep);
        double ekin = atoms.kinetic_energy();
        writer.write_stats(ts, ekin, epot, avg_temp.get());
        if (relaxed) {
            avg_temp.update(atoms.current_temperature_kelvin());
        }
        if (ts % (relaxation_time_deposit * 2) == 0) {
            deposit_energy(atoms, delta_Q);
            relaxed = false;
        } else if (ts % relaxation_time_deposit == 0) {
            relaxed = true;
        }
    }

    return 0;
}
