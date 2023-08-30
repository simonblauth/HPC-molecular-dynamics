#include "atoms.h"
#include "average.h"
#include "domain.h"
#include "ducastelle.h"
#include "mpi_support.h"
#include "neighbors.h"
#include "simulation_utils.h"
#include "thermostat.h"
#include "types.h"
#include "verlet.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;


Domain init_domain(Atoms& atoms, argparse::ArgumentParser parser) {
    auto ds = parser.get<std::vector<int>>("--domains");
    auto periodic = parser.get<std::vector<int>>("--periodic");
    double shift = parser.get<double>("--shift_atoms");

    auto max_pos = atoms.positions.rowwise().maxCoeff();
    auto min_pos = atoms.positions.rowwise().minCoeff();
    // bounding box
    Eigen::Array3d bbox{max_pos(0) - min_pos(0),
                        max_pos(1) - min_pos(1),
                        max_pos(2) - min_pos(2)};
    double offset_x = bbox(0) * shift - min_pos(0);
    double offset_y = bbox(1) * shift - min_pos(1);
    double offset_z = bbox(2) * shift - min_pos(2);
    Eigen::Array3d offset{offset_x, offset_y, offset_z};
    for (size_t i = 0; i < 3; i++) {
        if (periodic[i] == 0) {
            atoms.positions.row(i) += offset(i);
        }
    }

    double box_sz_x = periodic[0] == 0 ? bbox(0) + 2 * bbox(0) * shift : max_pos(0) + min_pos(0);
    double box_sz_y = periodic[1] == 0 ? bbox(1) + 2 * bbox(1) * shift : max_pos(1) + min_pos(1);
    double box_sz_z = periodic[2] == 0 ? bbox(2) + 2 * bbox(2) * shift : max_pos(2) + min_pos(2);

    Domain domain(MPI_COMM_WORLD,
        {box_sz_x, box_sz_y, box_sz_z},
        {ds[0], ds[1], ds[2]},
        {periodic[0], periodic[1], periodic[2]}
    );
    return domain;
}

double compute_stress(const Domain& domain, const Atoms& atoms, int dim = 2) {
    double left_domain_boundary{domain.coordinate(dim) * domain.domain_length(dim) / domain.decomposition(dim)};
    auto left_mask{atoms.positions.row(dim) < left_domain_boundary};
    double stress = 0;
    for (::Eigen::Index i{0}; i < left_mask.size(); ++i) {
        if (left_mask[i]) {
            stress += atoms.forces(dim, i);
        }
    }
    return stress;
}

int main(int argc, char *argv[]) {
    MPI::init_guard guard(&argc, &argv);
    std::stringstream ss;
    ss << "Worker " << MPI::comm_rank(MPI_COMM_WORLD) << ": ";
    auto worker = ss.str();
    // argument parsing
    argparse::ArgumentParser parser = default_parser("milestone 09");
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

    double target_temperaure = parser.get<double>("--temperature") * 1e-5;
    double relaxation_time = parser.get<size_t>("--relaxation_time") * timestep;


    // domain setup
    auto domain = init_domain(atoms, parser);
    std::cout << worker << " initialized domain" << std::endl;

    // stretching
    size_t stretch_interval = parser.get<size_t>("--stretch_interval");
    double length_increase = parser.get<double>("--stretch");

    // simulate
    std::cout << worker << " natoms: " << atoms.nb_atoms() << std::endl;
    domain.enable(atoms);
    std::cout << worker << " natoms: " << atoms.nb_atoms() << std::endl;
    std::cout << worker << " enabled domain" << std::endl;
    CumulativeAverage avg_stress(output_interval);
    CumulativeAverage avg_temp(output_interval);
    for (size_t ts = 0; ts < max_timesteps; ts++) {
        verlet_step1(atoms, timestep);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * cutoff);
        neighbor_list.update(atoms);
        double epot_local = ducastelle(atoms, neighbor_list, domain.nb_local(), cutoff);
        double stress_local = compute_stress(domain, atoms);
        verlet_step2(atoms, timestep);
        double temp_scaled_local = atoms.current_temperature(domain.nb_local());
        double temp_scaled = MPI::allreduce(temp_scaled_local, MPI_SUM, MPI_COMM_WORLD) / domain.size();
        berendsen_thermostat(atoms, target_temperaure, timestep,
                             relaxation_time, temp_scaled);

        if (ts % stretch_interval == 0) {
            Eigen::Array3d new_length{
                domain.domain_length(0), domain.domain_length(1),
                domain.domain_length(2) + length_increase};
            domain.scale(atoms, new_length);
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
        avg_stress.update(stress, ts);

        if (ts % output_interval == 0) {
            domain.disable(atoms);
            if (MPI::comm_rank(MPI_COMM_WORLD) == 0) {
                writer->write_traj(ts, atoms);
                writer->write_stats(ts, ekin, epot, avg_temp.get(), avg_stress.get());
            }
            domain.enable(atoms);
            avg_stress = 0;
            avg_temp = 0;
        }
    }
    delete writer;

    return 0;
}
