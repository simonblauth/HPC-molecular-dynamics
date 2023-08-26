#include "atoms.h"
#include "lj_direct_summation.h"
#include "simulation_utils.h"
#include "types.h"
#include "verlet.h"
#include "xyz.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    fs::path filepath = argv[0];
    auto pwd = filepath.parent_path();
    auto input_path = pwd / "lj54.xyz";
    auto [names, positions, velocities]{read_xyz_with_velocities(input_path)};
    std::cout << "loaded file from: " << input_path << std::endl;
    Atoms atoms(positions, velocities);

    argparse::ArgumentParser parser = default_parser("milestone 04");
    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }
    Writer writer(pwd, parser);

    // initialize simulation
    atoms.mass = parser.get<double>("--mass");
    double sigma = parser.get<double>("--sigma");
    double epsilon = parser.get<double>("--epsilon");
    double timestep = parser.get<double>("--timestep") * std::sqrt(atoms.mass * sigma * sigma / epsilon);
    size_t max_timesteps = parser.get<size_t>("--max_timesteps");

    // simulate
    for (size_t ts = 0; ts < max_timesteps; ts++) {
        writer.write_traj(ts, atoms);
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     atoms.mass);
        double epot = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.mass);
        double ekin = atoms.kinetic_energy();
        writer.write_energy(ts, ekin, epot);
    }

    return 0;
}
