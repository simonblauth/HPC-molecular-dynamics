#include "atoms.h"
#include "average.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "simulation_utils.h"
#include "thermostat.h"
#include "types.h"
#include "verlet.h"
#include "writer.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;


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
    Writer writer(pwd, parser);

    auto input_path = parser.get<std::string>("--input");
    auto [names, positions]{read_xyz(input_path)};

    writer.log("loaded file from: ", input_path);

    // initialize simulation
    Atoms atoms(names, positions);

    atoms.set_mass(parser.get<double>("--mass") * 103.6);
    double timestep = parser.get<double>("--timestep");
    size_t max_timesteps = parser.get<size_t>("--max_timesteps");

    size_t init_timesteps = parser.get<size_t>("--initial_relaxation");
    size_t relaxation_time = parser.get<size_t>("--relaxation_time");
    double target_temperaure = parser.get<double>("--temperature") * 1e-5;
    double relaxation_increase = parser.get<double>("--relaxation_time_factor");
    ThermostatScheduler scheduler(relaxation_increase, relaxation_time);
    Equilibrium equilibrium(relaxation_increase, relaxation_time,
                            target_temperaure, timestep, init_timesteps);

    double delta_Q = parser.get<double>("--deposit_energy");
    size_t relaxation_time_deposit = parser.get<size_t>("--relaxation_time_deposit");
    EnergyPump pump(relaxation_time_deposit, delta_Q);

    double cutoff = parser.get<double>("--cutoff");
    NeighborList neighbor_list(cutoff);

    // relax
    writer.log("Equilibriating the system...");
    for (size_t i = 0; i < init_timesteps; i++) {
        // writer.write_traj(i, atoms);
        verlet_step1(atoms, timestep);
        neighbor_list.update(atoms);
        double epot = ducastelle(atoms, neighbor_list, cutoff);
        verlet_step2(atoms, timestep);
        equilibrium.step(atoms, i, atoms.current_temperature());
    }

    // simulate
    bool relaxed = false;
    double alpha = parser.get<double>("--smoothing");
    ExponentialAverage avg_temp(alpha, atoms.current_temperature_kelvin());
    writer.log("Starting simulation");
    for (size_t ts = 0; ts < max_timesteps; ts++) {
        writer.write_traj(ts, atoms);
        verlet_step1(atoms, timestep);
        neighbor_list.update(atoms);
        double epot = ducastelle(atoms, neighbor_list, cutoff);
        verlet_step2(atoms, timestep);
        double ekin = atoms.kinetic_energy();
        writer.write_stats(ts, ekin, epot, avg_temp.get());
        if (pump.relaxed()) {
            avg_temp.update(atoms.current_temperature_kelvin());
        }
        pump.step(atoms, ts, ekin);
    }

    return 0;
}
