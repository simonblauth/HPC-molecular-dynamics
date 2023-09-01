#include "atoms.h"
#include "average.h"
#include "domain.h"
#include "ducastelle.h"
#include "mpi_support.h"
#include "neighbors.h"
#include "simulation_utils.h"
#include "simulation_utils_mpi.h"
#include "thermostat.h"
#include "types.h"
#include "verlet.h"
#include "writer_mpi.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    MPI::init_guard guard(&argc, &argv);
    argparse::ArgumentParser parser = default_parser("milestone 09", argc, argv);

    fs::path filepath = argv[0];
    auto pwd = filepath.parent_path();
    MPIWriter writer(pwd, parser);

    auto input_path = parser.get<std::string>("--input");
    auto [names, positions]{read_xyz(input_path)};

    writer.log("loaded file from: ", input_path);

    // initialize simulation
    Atoms atoms(names, positions);
    atoms.set_mass(parser.get<double>("--mass") * 103.6);
    writer.debug("initialized atoms");
    SimulationParameters sim(parser);
    NeighborList neighbor_list(sim.cutoff());
    writer.debug("initialized neighbors");

    Equilibrium equilibrium(sim.relaxation_factor(), sim.relaxation_time(),
                            sim.target_temperature(), sim.timestep(), sim.init_timesteps());

    // domain setup
    auto domain = init_domain(atoms, parser);
    writer.debug("initialized domain");
    writer.debug("Total number of atoms: ", atoms.nb_atoms());
    domain.enable(atoms);
    writer.debug_all("Number of atoms: ", atoms.nb_atoms());
    writer.debug("enabled domain");

    // relax
    writer.log("Equilibriating the system...");
    for (size_t i = 0; i < sim.init_timesteps(); i++) {
        // writer.write_traj(i, atoms);
        verlet_step1(atoms, sim.timestep());
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * sim.cutoff());
        neighbor_list.update(atoms);
        double epot = ducastelle(atoms, neighbor_list, sim.cutoff());
        verlet_step2(atoms, sim.timestep());
        double temp_local = atoms.current_temperature(domain.nb_local());
        double temp = MPI::allreduce(temp_local, MPI_SUM, MPI_COMM_WORLD) / domain.size();
        equilibrium.step(atoms, i, temp);
    }

    // simulate
    double alpha = parser.get<double>("--smoothing");
    ExponentialAverage avg_stress(alpha);
    CumulativeAverage avg_temp(writer.get_output_interval());
    double strain = 0;
    writer.log("Starting actual simulation");
    for (size_t ts = 0; ts < sim.max_timesteps(); ts++) {
        verlet_step1(atoms, sim.timestep());
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * sim.cutoff());
        neighbor_list.update(atoms);
        double epot_local = ducastelle(atoms, neighbor_list, domain.nb_local(), sim.cutoff());
        double stress_local = compute_stress(domain, atoms);
        verlet_step2(atoms, sim.timestep());
        // double temp_scaled_local = atoms.current_temperature(domain.nb_local());
        // double temp_scaled = MPI::allreduce(temp_scaled_local, MPI_SUM, MPI_COMM_WORLD) / domain.size();
        // berendsen_thermostat(atoms, target_temperaure, timestep,
        //                      relaxation_time, temp_scaled);

        if (ts % sim.stretch_interval() == 0) {
            Eigen::Array3d new_length{
                domain.domain_length(0), domain.domain_length(1),
                domain.domain_length(2) + sim.length_increase()};
            domain.scale(atoms, new_length);
            strain += sim.length_increase();
        }

        double ekin_local = atoms.kinetic_energy(domain.nb_local());
        double temp_local = atoms.current_temperature_kelvin(domain.nb_local());

        double ekin = MPI::allreduce(ekin_local, MPI_SUM, MPI_COMM_WORLD);
        double epot = MPI::allreduce(epot_local, MPI_SUM, MPI_COMM_WORLD);

        // cumulative average over temp
        double temp = MPI::allreduce(temp_local, MPI_SUM, MPI_COMM_WORLD) / domain.size();
        avg_temp.update(temp, ts);

        // cumulative average over stress
        double stress = MPI::allreduce(stress_local, MPI_SUM, MPI_COMM_WORLD);
        stress /= (domain.domain_length(0) * domain.domain_length(1) * domain.decomposition(2));
        avg_stress.update(stress);

        if (ts % writer.get_output_interval() == 0) {
            domain.disable(atoms);
            writer.write_traj(ts, atoms);
            writer.write_stats(ts, ekin, epot, avg_temp.get(), avg_stress.get(), strain);
            domain.enable(atoms);
        }
    }

    return 0;
}
