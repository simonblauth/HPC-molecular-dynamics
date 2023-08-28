#include "atoms.h"
#include "ducastelle.h"
#include "mpi_support.h"
#include "neighbors.h"
#include "simulation_utils.h"
#include "types.h"
#include "verlet.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    MPI::init_guard guard(&argc, &argv);
    // argument parsing
    argparse::ArgumentParser parser = default_parser("milestone 08");
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
    size_t relaxation_time = parser.get<size_t>("--relaxation_time");
    double delta_Q = parser.get<double>("--deposit_energy");

    double cutoff = parser.get<double>("--cutoff");
    NeighborList neighbor_list(cutoff);
    neighbor_list.update(atoms);

    // simulate
    bool relax = false;
    double avg_temp = 0;
    for (size_t ts = 0; ts < max_timesteps; ts++) {
        writer.write_traj(ts, atoms);
        verlet_step1(atoms, timestep);
        double epot = ducastelle(atoms, neighbor_list, cutoff);
        verlet_step2(atoms, timestep);
        double ekin = atoms.kinetic_energy();
        if (ts % relaxation_time == 0) {
            relax = !relax;
            if (relax) {
                avg_temp = 0; // compute temperature
            } else {
                // deposit energy
                double lambda = std::sqrt(1 + delta_Q / ekin);
                atoms.velocities *= lambda;
            }
        }
        if (relax) {
            // cumulative average
            size_t step = ts % relaxation_time + 1;
            double current_temp = atoms.current_temperature_kelvin();
            avg_temp += (current_temp - avg_temp) / step;
        }
        writer.write_stats(ts, ekin, epot, avg_temp);
    }

    return 0;
}
