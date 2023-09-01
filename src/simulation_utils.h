#ifndef __SIMULATION_UTILS_H
#define __SIMULATION_UTILS_H

#include "atoms.h"
#include <argparse/argparse.hpp>
#include <iostream>

// Holds the parameters of a simulation
class SimulationParameters {
  private:
    double timestep_;
    size_t max_timesteps_;
    double cutoff_;
    double target_temperature_;
    double relaxation_time_;
    double relaxation_factor_;
    size_t init_timesteps_;
    size_t stretch_interval_;
    double length_increase_;
    double delta_Q_;
    size_t relaxation_time_deposit_;

  public:
    SimulationParameters(argparse::ArgumentParser& parser) {
        timestep_ = parser.get<double>("--timestep");
        max_timesteps_ = parser.get<size_t>("--max_timesteps");
        cutoff_ = parser.get<double>("--cutoff");
        target_temperature_ = parser.get<double>("--temperature") * 1e-5;
        relaxation_time_ = parser.get<size_t>("--relaxation_time") * timestep_;
        relaxation_factor_ = parser.get<double>("--relaxation_time_factor");
        init_timesteps_ = parser.get<size_t>("--initial_relaxation");
        stretch_interval_ = parser.get<size_t>("--stretch_interval");
        length_increase_ = parser.get<double>("--stretch");
        delta_Q_ = parser.get<double>("--deposit_energy");
        relaxation_time_deposit_ = parser.get<size_t>("--relaxation_time_deposit");
    }
    ~SimulationParameters() {}
    double timestep() const { return timestep_; }
    size_t max_timesteps() const { return max_timesteps_; }
    double cutoff() const { return cutoff_; }
    double target_temperature() const { return target_temperature_; }
    double relaxation_time() const { return relaxation_time_; }
    double relaxation_factor() const { return relaxation_factor_; }
    size_t init_timesteps() const { return init_timesteps_; }
    size_t stretch_interval() const { return stretch_interval_; }
    double length_increase() const { return length_increase_; }
    double delta_Q() const { return delta_Q_; }
    size_t relaxation_time_deposit() const { return relaxation_time_deposit_; }
};


// constructs an ArgumentParser with default arguments used for most simulations
argparse::ArgumentParser default_parser(const char* name) {
    argparse::ArgumentParser parser(name);
    // IO parameters
    parser.add_argument("-s", "--silent")
        .help("Do not print stats.")
        .default_value(false)
        .implicit_value(true);
    parser.add_argument("-v", "--verbose")
        .help("Print all sorts of debug info.")
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
    parser.add_argument("--traj")
        .help("If set write trajectory to xyz. Optionally takes a path for the output file.")
        .nargs(argparse::nargs_pattern::optional);
    parser.add_argument("-i", "--input")
        .help("Takes a path for the input file.")
        .nargs(1);
    parser.add_argument("--smoothing")
        .help("Smoothing factor for exponential average.")
        .nargs(1)
        .default_value<double>(0.01)
        .scan<'g', double>();
    // basic simulation parameters
    parser.add_argument("--max_timesteps")
        .help("The number of timesteps to run the simulation for.")
        .nargs(1)
        .default_value<size_t>(100000)
        .scan<'u', size_t>();
    parser.add_argument("-t", "--timestep")
        .help("The timestep for the simulation in femtoseconds.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    parser.add_argument("--mass")
        .help("The mass of the atoms in g/mol.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    // lennard-jones
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
    // temperature
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
    parser.add_argument("--relaxation_time_factor")
        .help("Fraction of the initial relaxation time to run the thermosthat for.")
        .nargs(1)
        .default_value<double>(2.0)
        .scan<'g', double>();
    parser.add_argument("--relaxation_time_deposit")
        .help("The relaxation time for the energy deposit (in Timesteps).")
        .nargs(1)
        .default_value<size_t>(1000)
        .scan<'u', size_t>();
    parser.add_argument("--deposit_energy")
        .help("The amount of energy to deposit into the system every interval.")
        .nargs(1)
        .default_value<double>(10)
        .scan<'g', double>();
    parser.add_argument("--initial_relaxation")
        .help("Number of time steps to let the system relax initially.")
        .nargs(1)
        .default_value<size_t>(1000)
        .scan<'u', size_t>();
    // domain
    parser.add_argument("--cutoff")
        .help("The maximum distance for the neighbor search.")
        .nargs(1)
        .default_value<double>(5.0)
        .scan<'g', double>();
    parser.add_argument("--domains")
        .help("The number of domains in x, y, z direction.")
        .nargs(3)
        .default_value(std::vector<int>{1, 1, 1})
        .scan<'i', int>();
    parser.add_argument("--periodic")
        .help("The periodicity of domains in x, y, z direction, 1 means periodic, 0 not.")
        .nargs(3)
        .default_value(std::vector<int>{0, 0, 0})
        .scan<'i', int>();
    parser.add_argument("--shift_atoms")
        .help("Shift the atoms by a percentage of the box size (also increases box size) to have more space.")
        .nargs(1)
        .default_value<double>(0.1)
        .scan<'g', double>();
    parser.add_argument("--stretch")
        .help("Stretch the domain in z-direction by <stretch> angstrom.")
        .nargs(1)
        .default_value<double>(1.0)
        .scan<'g', double>();
    parser.add_argument("--stretch_interval")
        .help("Stretch the domain in z-direction every <stretch_interval> steps.")
        .nargs(1)
        .default_value<size_t>(100)
        .scan<'u', size_t>();

    return parser;
}

// constructs and initializes an ArgumentParser with default arguments used for most simulations
argparse::ArgumentParser default_parser(const char* name, int argc, const char *const *argv) {
    argparse::ArgumentParser parser = default_parser(name);
    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }
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

// TODO: Stretcher class
// TODO: stats_collector function
// TODO: move implementations to .cpp
// TODO: check default values of parser
// TODO: init atoms function

#endif // __SIMULATION_UTILS_H
