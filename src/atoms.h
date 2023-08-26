#ifndef __ATOMS_H
#define __ATOMS_H

#include <cmath>
#include "types.h"

// Holds the current state of a simulation.
class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    double mass = 1.0;

    Atoms(const size_t nb_atoms)
        : positions(3, nb_atoms),
          velocities(3, nb_atoms),
          forces(3, nb_atoms) {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p)
        : positions{p},
          velocities{3, p.cols()},
          forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v)
        : positions{p},
          velocities{v},
          forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

    double kinetic_energy() const {
        return mass * velocities.cwiseAbs2().sum() / 2;
    }

    double current_temperature() const {
      double k_B = 8.617333262;
      return kinetic_energy() / (nb_atoms() * k_B) * 2 / 3;
    }
};

#endif // __ATOMS_H
