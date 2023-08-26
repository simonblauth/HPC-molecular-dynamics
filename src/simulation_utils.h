#include "atoms.h"
#include "xyz.h"
#include <argparse/argparse.hpp>
#include <filesystem>

namespace fs = std::filesystem;

// Handles writing simulation data in various formats such as console, csv and xyz
class Writer {
private:
    std::filesystem::path csv_path;
    bool write_to_csv;
    bool write_to_console;
    size_t output_interval;
    std::ofstream csv;
    std::ofstream traj;

  public:
    Writer(fs::path pwd, argparse::ArgumentParser parser) {
        write_to_csv = parser.is_used("--csv");
        write_to_console = !parser.get<bool>("--silent");
        output_interval = parser.get<size_t>("--output_interval");
        if (write_to_csv) {
            try {
                auto p = parser.get<std::string>("--csv");
                csv_path = fs::path(p);
            } catch (const std::logic_error &e) {
                csv_path = pwd / "output.csv";
                std::cout << "No csv path provided, using default path: "
                          << csv_path << std::endl;
            }

            csv.open(csv_path);
            csv << "Timestep,Total Energy,Kinetic Energy,Potential Energy"
                << std::endl;
        }
        auto traj_path = pwd / "traj.xyz";
        traj.open(traj_path);
    }
    ~Writer() {
        csv.close();
        traj.close();
    }
    void write_energy(size_t timestep, double ekin, double epot) {
        if (timestep % output_interval == 0) {
            if (write_to_console) {
                std::cout << "Total Energy: " << ekin + epot << ", ";
                std::cout << "Kinetic Energy: " << ekin << ", ";
                std::cout << "Potential Energy: " << epot << std::endl;
            }
            if (write_to_csv) {
                csv << timestep << "," << ekin + epot << "," << ekin << "," << epot << std::endl;
            }
        }
    }
    void write_traj(size_t timestep, Atoms& atoms) {
        if (timestep % output_interval == 0) {
            write_xyz(traj, atoms);
        }
    }
};

// constructs an ArgumentParser with default arguments used for most simulations
argparse::ArgumentParser default_parser(const char* name) {
    argparse::ArgumentParser parser(name);
    // IO parameters
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
    // basic simulation parameters
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
        .help("The mass of the atoms.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    parser.add_argument("--sigma")
        .help("The σ parameter for the lennard-jones potential.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    parser.add_argument("--epsilon")
        .help("The ϵ parameter for the lennard-jones potential.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    // cubic lattice
    parser.add_argument("--lattice_size")
        .help("The side length of the cubic lattice to simulate.")
        .default_value<size_t>(3)
        .scan<'u', size_t>();
    parser.add_argument("--lattice_dist")
        .help("The distance between atoms in the cubic lattice to simulate (in multiples of σ).")
        .nargs(1)
        .default_value<double>(1.05)
        .scan<'g', double>();
    // thermostat
    parser.add_argument("--temperature")
        .help("The target temperature for the thermostat (in Kelvin).")
        .nargs(1)
        .default_value<double>(50)
        .scan<'g', double>();
    parser.add_argument("--relaxation_time")
        .help("The relaxation time for the thermostat (in Timesteps).")
        .nargs(1)
        .default_value<size_t>(1000)
        .scan<'u', size_t>();
    return parser;
}

// constructs an atom configuration in a simple cubic lattice
Atoms init_cubic_lattice(size_t nb_atoms_per_lattice, double lattice_distance) {
    Atoms atoms(std::pow(nb_atoms_per_lattice, 3));
    for (size_t i = 0; i < nb_atoms_per_lattice; i++) {
        for (size_t j = 0; j < nb_atoms_per_lattice; j++) {
            for (size_t k = 0; k < nb_atoms_per_lattice; k++) {
                size_t atom_idx =
                    k + nb_atoms_per_lattice * j +
                    nb_atoms_per_lattice * nb_atoms_per_lattice * i;
                Eigen::Vector3d pos_vec = {(double) i, (double) j, (double) k};
                pos_vec *= lattice_distance;
                atoms.positions.col(atom_idx) = pos_vec;
            }
        }
    }
    return atoms;
}
