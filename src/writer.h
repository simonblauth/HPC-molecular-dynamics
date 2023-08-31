#ifndef __WRITER_H
#define __WRITER_H

#include "xyz.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

// Handles writing simulation data in various formats such as console, csv and xyz
class Writer {
  protected:
    fs::path csv_path;
    fs::path traj_path;
    bool write_to_csv;
    bool write_to_xyz;
    bool write_to_console;
    bool verbose;
    size_t output_interval;
    std::ofstream csv;
    std::ofstream traj;

  public:
    Writer() = default;
    Writer(fs::path pwd, argparse::ArgumentParser parser) {
        write_to_csv = parser.is_used("--csv");
        write_to_xyz = parser.is_used("--traj");
        write_to_console = !parser.get<bool>("--silent");
        verbose = parser.get<bool>("--verbose");
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
            csv << "Timestep,Total Energy,Kinetic Energy,Potential Energy,Temperature,Stress,Strain"
                << std::endl;
        }
        if (write_to_xyz) {
            try {
                auto p = parser.get<std::string>("--traj");
                traj_path = fs::path(p);
            } catch (const std::logic_error &e) {
                traj_path = pwd / "traj.xyz";
                std::cout << "No traj path provided, using default path: "
                          << csv_path << std::endl;
            }

            traj.open(traj_path);
        }
    }
    ~Writer() {
        csv.close();
        traj.close();
    }
    // write stats to console/csv in given intervals
    void write_stats(size_t timestep, double ekin, double epot, double temp = 0, double stress = 0, double strain = 0) {
        if (timestep % output_interval == 0) {
            if (write_to_console) {
                std::cout << "Frame: " << (int)(timestep / output_interval) << " ";
                std::cout << "Total Energy: " << ekin + epot << ", ";
                std::cout << "Kinetic Energy: " << ekin << ", ";
                std::cout << "Potential Energy: " << epot << ", ";
                std::cout << "Temperature: " << temp << ", ";
                std::cout << "Stress: " << stress << ", ";
                std::cout << "Strain: " << strain << std::endl;
            }
            if (write_to_csv) {
                csv << timestep << "," << ekin + epot << "," << ekin << "," << epot << "," << temp << "," << stress << "," << strain << std::endl;
            }
        }
    }
    // write simulation state to xyz file in given intervals
    void write_traj(size_t timestep, Atoms& atoms) {
        if (timestep % output_interval == 0) {
            write_xyz(traj, atoms);
        }
    }
    // print some info to console
    void log(std::string msg) {
        if (write_to_console) {
            std::cout << msg << std::endl;
        }
    }
    // print some info to console with an additional value
    template<typename T>
    void log(std::string msg, T value) {
        if (write_to_console) {
            std::cout << msg << value << std::endl;
        }
    }
    // print some debug info to console
    void debug(std::string msg) {
        if (verbose) {
            std::cout << msg << std::endl;
        }
    }
    // print some debug info to console with an additional value
    template<typename T>
    void debug(std::string msg, T value) {
        if (verbose) {
            std::cout << msg << value << std::endl;
        }
    }
    size_t get_output_interval() { return output_interval; }
};

#endif // __WRITER_H