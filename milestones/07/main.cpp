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
    argparse::ArgumentParser parser = default_parser("milestone 07", argc, argv);

    fs::path filepath = argv[0];
    auto pwd = filepath.parent_path();
    Writer writer(pwd, parser);

    auto input_path = parser.get<std::string>("--input");
    auto [names, positions]{read_xyz(input_path)};

    writer.log("loaded file from: ", input_path);

    // initialize simulation
    Atoms atoms(names, positions);
    atoms.set_mass(parser.get<double>("--mass") * 103.6);

    SimulationParameters sim(parser);
    ThermostatScheduler scheduler(sim.relaxation_factor(), sim.relaxation_time());
    Equilibrium equilibrium(sim.relaxation_factor(), sim.relaxation_time(),
                            sim.target_temperature(), sim.timestep(), sim.init_timesteps());
    EnergyPump pump(sim.relaxation_time_deposit(), sim.delta_Q());
    NeighborList neighbor_list(sim.cutoff());

    // relax
    writer.log("Equilibriating the system...");
    for (size_t i = 0; i < sim.init_timesteps(); i++) {
        // writer.write_traj(i, atoms);
        verlet_step1(atoms, sim.timestep());
        neighbor_list.update(atoms);
        double epot = ducastelle(atoms, neighbor_list, sim.cutoff());
        verlet_step2(atoms, sim.timestep());
        equilibrium.step(atoms, i, atoms.current_temperature());
    }

    // simulate
    bool relaxed = false;
    double alpha = parser.get<double>("--smoothing");
    ExponentialAverage avg_temp(alpha, atoms.current_temperature_kelvin());
    writer.log("Starting simulation");
    for (size_t ts = 0; ts < sim.max_timesteps(); ts++) {
        writer.write_traj(ts, atoms);
        verlet_step1(atoms, sim.timestep());
        neighbor_list.update(atoms);
        double epot = ducastelle(atoms, neighbor_list, sim.cutoff());
        verlet_step2(atoms, sim.timestep());
        double ekin = atoms.kinetic_energy();
        writer.write_stats(ts, ekin, epot, avg_temp.get());
        if (pump.relaxed()) {
            avg_temp.update(atoms.current_temperature_kelvin());
        }
        pump.step(atoms, ts, ekin);
    }

    return 0;
}
