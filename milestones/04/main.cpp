#include "atoms.h"
#include "lj_direct_summation.h"
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
    auto traj_path = pwd / "traj.xyz";
    auto input_path = pwd / "lj54.xyz";
    auto csv_path = pwd / "output.csv";
    std::ofstream traj(traj_path);
    auto [names, positions, velocities]{read_xyz_with_velocities(input_path)};
    std::cout << "loaded file from: " << input_path << std::endl;
    Atoms atoms(positions, velocities);

    argparse::ArgumentParser parser("milestone 04");
    parser.add_argument("-s", "--silent")
        .help("Do not print stats.")
        .default_value(false)
        .implicit_value(true);
    parser.add_argument("--output_interval")
        .help("How many timesteps to wait between outputting stats.")
        .nargs(1)
        .default_value<size_t>(100)
        .scan<'u', size_t>();
    parser.add_argument("--csv")
        .help("If set write stats to csv. Optionally takes a path for the output file.")
        .nargs(argparse::nargs_pattern::optional);
    parser.add_argument("--max_timesteps")
        .help("The number of timesteps to run the simulation for.")
        .nargs(1)
        .default_value<size_t>(100000)
        .scan<'u', size_t>();
    parser.add_argument("-t", "--timestep")
        .help("The timestep for the simulation (multiplied by √(m ⋅ σ² / ϵ)).")
        .nargs(1)
        .default_value<double>(1e-3)
        .scan<'g', double>();
    parser.add_argument("--mass")
        .help("The number of timesteps to run the simulation for.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    parser.add_argument("--sigma")
        .help("The number of timesteps to run the simulation for.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    parser.add_argument("--epsilon")
        .help("The number of timesteps to run the simulation for.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }

    // simulation parameters
    double mass = parser.get<double>("--mass");
    double sigma = parser.get<double>("--sigma");
    double epsilon = parser.get<double>("--epsilon");
    double timestep = parser.get<double>("--timestep") * std::sqrt(mass * sigma * sigma / epsilon);
    size_t max_timesteps = parser.get<size_t>("--max_timesteps");

    // output parameters
    size_t output_interval = parser.get<size_t>("--output_interval");
    bool write_to_console = !parser.get<bool>("--silent");
    bool write_to_csv = parser.is_used("--csv");
    std::ofstream csv;
    if (write_to_csv) {
        try {
            auto p = parser.get<std::string>("--csv");
            csv_path = fs::path(p);
        } catch (const std::logic_error &e) {
            std::cout << "No csv path provided, using default path: " << csv_path << std::endl;
        }

        csv.open(csv_path);
        csv << "Timestep,Total Energy,Kinetic Energy,Potential Energy"
            << std::endl;
    }

    for (size_t ts = 0; ts < max_timesteps; ts++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     mass);
        double epot = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        double ekin = atoms.kinetic_energy();
        if (ts % output_interval == 0) {
            if (write_to_console) {
                std::cout << "Total Energy: " << ekin + epot << ", ";
                std::cout << "Kinetic Energy: " << ekin << ", ";
                std::cout << "Potential Energy: " << epot << std::endl;
            }
            if (write_to_csv) {
                csv << ts << "," << ekin + epot << "," << ekin << "," << epot << std::endl;
            }

            write_xyz(traj, atoms);
            // TODO: plot epot in report
        }
    }
    traj.close();
    csv.close();

    return 0;
}