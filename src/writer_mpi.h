#ifndef __WRITER_MPI_H
#define __WRITER_MPI_H

#include "mpi_support.h"
#include "writer.h"
#include "xyz.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>
#include <stddef.h>

namespace fs = std::filesystem;

// Handles writing simulation data in various formats such as console, csv and xyz. Can also deal with multiple processes.
class MPIWriter : public Writer {
private:
    int thread_id;
public:
    MPIWriter(fs::path pwd, argparse::ArgumentParser parser) : Writer(pwd, parser) {
        thread_id = MPI::comm_rank(MPI_COMM_WORLD);
    }
    ~MPIWriter() {}
    // write stats to console/csv in given intervals in main thread
    void write_stats(size_t timestep, double ekin, double epot, double temp = 0, double stress = 0, double strain = 0) {
        if (thread_id == 0) {
            Writer::write_stats(timestep, ekin, epot, temp, stress, strain);
        }
    }
    // write simulation state to xyz file in given intervals in main thread
    void write_traj(size_t timestep, Atoms& atoms) {
        if (thread_id == 0) {
            Writer::write_traj(timestep, atoms);
        }
    }
    // print some info to console in main thread
    void log(std::string msg) {
        if (thread_id == 0) {
            Writer::log(msg);
        }
    }
    // print some info to console with an additional value in main thread
    template<typename T>
    void log(std::string msg, T value) {
        if (thread_id == 0) {
            Writer::log(msg, value);
        }
    }
    // print some info to console for each worker
    void log_all(std::string msg) {
        if (write_to_console) {
            std::cout << "Worker " << thread_id << ": " << msg << std::endl;
        }
    }
    // print some info to console with an additional value for each worker
    template<typename T>
    void log_all(std::string msg, T value) {
        if (write_to_console) {
            std::cout << "Worker " << thread_id << ": "  << msg << value << std::endl;
        }
    }
    // print some debug info to console in main thread
    void debug(std::string msg) {
        if (thread_id == 0) {
            Writer::debug(msg);
        }
    }
    // print some debug info to console with an additional value in main thread
    template<typename T>
    void debug(std::string msg, T value) {
        if (thread_id == 0) {
            Writer::debug(msg, value);
        }
    }
    // print some debug info to console for each worker
    void debug_all(std::string msg) {
        if (verbose) {
            std::cout << "Worker " << thread_id << ": " << msg << std::endl;
        }
    }
    // print some debug info to console with an additional value for each worker
    template<typename T>
    void debug_all(std::string msg, T value) {
        if (verbose) {
            std::cout << "Worker " << thread_id << ": "  << msg << value << std::endl;
        }
    }

};


#endif // __WRITER_MPI_H