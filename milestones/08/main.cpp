#include "atoms.h"
#include "domain.h"
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
    std::stringstream ss;
    ss << "Worker " << MPI::comm_rank(MPI_COMM_WORLD) << ": ";
    auto worker = ss.str();
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

    std::cout << worker << "loaded file from: " << input_path << std::endl;

    Writer* writer = NULL;
    if (MPI::comm_rank(MPI_COMM_WORLD) == 0) {
        writer = new Writer(pwd, parser);
    }

    size_t output_interval = parser.get<size_t>("--output_interval");

    // initialize simulation
    Atoms atoms(names, positions);
    std::cout << worker << " initialized atoms" << std::endl;

    atoms.set_mass(parser.get<double>("--mass") * 103.6);
    double timestep = parser.get<double>("--timestep");
    size_t max_timesteps = parser.get<size_t>("--max_timesteps");

    double cutoff = parser.get<double>("--cutoff");
    NeighborList neighbor_list(cutoff);
    std::cout << worker << " initialized neighbors" << std::endl;

    // domain setup
    auto max_pos = atoms.positions.rowwise().maxCoeff();
    auto min_pos = atoms.positions.rowwise().minCoeff();
    auto box_size_pre = max_pos - min_pos;
    auto offset = box_size_pre * parser.get<double>("--shift_atoms");
    auto box_size = box_size_pre + 2 * offset;
    atoms.positions.colwise() += offset;
    std::cout << worker << " computed box size: \n" << box_size << std::endl;
    

    auto ds = parser.get<std::vector<int>>("--domains");
    Domain domain(MPI_COMM_WORLD,
        {box_size(0), box_size(1), box_size(2)},
        {ds[0], ds[1], ds[2]},
        {0, 0, 0}
    );
    std::cout << worker << " initialized domain" << std::endl;

    // simulate
    domain.enable(atoms);
    std::cout << worker << " enabled domain" << std::endl;
    for (size_t ts = 0; ts < max_timesteps; ts++) {
        verlet_step1(atoms, timestep);
        domain.update_ghosts(atoms, 2 * cutoff);
        neighbor_list.update(atoms);
        double epot_local = ducastelle(atoms, neighbor_list, domain.nb_local(), cutoff);
        verlet_step2(atoms, timestep);
        domain.exchange_atoms(atoms);

        double ekin_local = atoms.kinetic_energy();

        double ekin = MPI::allreduce(ekin_local, MPI_SUM, MPI_COMM_WORLD);
        double epot = MPI::allreduce(epot_local, MPI_SUM, MPI_COMM_WORLD);

        if (ts % output_interval == 0) {
            domain.disable(atoms);
            if (MPI::comm_rank(MPI_COMM_WORLD) == 0) {
                writer->write_traj(ts, atoms);
                writer->write_stats(ts, ekin, epot, 0);
            }
            domain.enable(atoms);
        }
    }
    delete writer;

    return 0;
}
