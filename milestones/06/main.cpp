#include "atoms.h"
#include "lj_direct_summation.h"
#include "neighbors.h"
#include "thermostat.h"
#include "simulation_utils.h"
#include "types.h"
#include "verlet.h"
#include "writer.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    argparse::ArgumentParser parser = default_parser("milestone 06", argc, argv);

    fs::path filepath = argv[0];
    auto pwd = filepath.parent_path();
    Writer writer(pwd, parser);

    // initialize simulation
    double sigma = parser.get<double>("--sigma");
    double epsilon = parser.get<double>("--epsilon");

    double lattice_distance = parser.get<double>("--lattice_dist") * sigma;
    size_t nb_atoms_per_lattice = parser.get<size_t>("--lattice_size");

    Atoms atoms = init_cubic_lattice(nb_atoms_per_lattice, lattice_distance);
    atoms.set_mass(parser.get<double>("--mass"));

    double timestep = parser.get<double>("--timestep") * std::sqrt(atoms.get_mass() * sigma * sigma / epsilon);
    SimulationParameters sim(parser);
    Equilibrium equilibrium(sim.relaxation_factor(), sim.relaxation_time(),
                            sim.target_temperature(), sim.timestep(), sim.max_timesteps());
    NeighborList neighbor_list;

    // simulate
    for (size_t ts = 0; ts < sim.max_timesteps(); ts++) {
        writer.write_traj(ts, atoms);
        verlet_step1(atoms, timestep);
        double epot = lj_direct_summation(atoms, neighbor_list, sim.cutoff(), epsilon, sigma);
        verlet_step2(atoms, timestep);
        double ekin = atoms.kinetic_energy();
        equilibrium.step(atoms, ts, atoms.current_temperature());
        writer.write_stats(ts, ekin, epot, atoms.current_temperature_kelvin());
    }

    return 0;
}
