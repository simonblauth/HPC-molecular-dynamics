#ifndef __ATOMS_H
#define __ATOMS_H

#include <algorithm>
#include <cmath>
#include "types.h"

// Holds the current state of a simulation.
class Atoms {
  private:
    double mass_ = 1.0;
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;
    Names_t names;

    Atoms(const Atoms& other)
        : positions(3, other.nb_atoms()),
          velocities(3, other.nb_atoms()),
          forces(3, other.nb_atoms()),
          masses(other.nb_atoms()),
          names(other.nb_atoms()) {
        positions = other.positions;
        velocities = other.velocities;
        forces = other.forces;
        masses = other.masses;
        names = other.names;
    }

    Atoms(const size_t nb_atoms)
        : positions(3, nb_atoms),
          velocities(3, nb_atoms),
          forces(3, nb_atoms),
          masses(nb_atoms),
          names(nb_atoms) {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        std::fill(names.begin(), names.end(), "H");
    }

    Atoms(const Positions_t &p)
        : positions{p},
          velocities{3, p.cols()},
          forces{3, p.cols()},
          masses{p.cols()},
          names(p.cols()) {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        std::fill(names.begin(), names.end(), "H");
    }

    Atoms(const Names_t &n, Positions_t &p)
        : positions{p},
          velocities{3, p.cols()},
          forces{3, p.cols()},
          masses{p.cols()},
          names{n} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }

    Atoms(const Positions_t &p, const Velocities_t &v)
        : positions{p},
          velocities{v},
          forces{3, p.cols()},
          masses{p.cols()},
          names(p.cols()) {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setOnes();
        std::fill(names.begin(), names.end(), "H");
    }

    Atoms(const Names_t &n, const Positions_t &p, const Velocities_t &v)
        : positions{p},
          velocities{v},
          forces{3, p.cols()},
          masses{p.cols()},
          names{n} {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setOnes();
    }

    void resize(Eigen::Array3Xd &a, size_t size) {
        size_t original_size = a.cols();
        Eigen::Array3Xd tmp = a; // should do deepcopy
        a.resize(3, size); // reallocates memory
        size_t upper = std::min(original_size, size);
        a(Eigen::all, Eigen::seq(0, upper - 1)) = tmp(Eigen::all, Eigen::seq(0, upper - 1)); // copy data back
    }

    void resize(size_t size) {
        resize(positions, size);
        resize(velocities, size);
        resize(forces, size);
        masses.resize(size);
        masses.fill(mass_); // masses are the same anyways
        // no need to resize names
        // names.resize(size); // vector resize preserves data
    }

    void set_mass(double mass) {
        mass_ = mass;
        masses.fill(mass);
    }

    double get_mass() {
        return mass_;
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

    double kinetic_energy() const {
        return mass_ * velocities.cwiseAbs2().sum() / 2;
    }

    // returns temperature in 1e-5 * K, so scale with 1e5 to get K
    double current_temperature() const {
      double k_B = 8.617333262;
      return kinetic_energy() / (nb_atoms() * k_B) * 2 / 3;
    }

    // returns temperature in Kelvin
    double current_temperature_kelvin() const {
      return current_temperature() * 1e5;
    }
};

#endif // __ATOMS_H
