#ifndef __THERMOSTAT_H
#define __THERMOSTAT_H

#include "atoms.h"

// Implementation of the https://en.wikipedia.org/wiki/Berendsen_thermostat.
void berendsen_thermostat(Atoms &atoms, double target_temperature,
                          double timestep, double relaxation_time);
void berendsen_thermostat(Atoms &atoms, double target_temperature,
                          double timestep, double relaxation_time,
                          double current_temperature);

class ThermostatScheduler {
  private:
    double factor_;
    double relaxation_time_;

  public:
    ThermostatScheduler(double factor, double relaxation_time)
        : factor_(factor),
          relaxation_time_(relaxation_time) {
    }
    double relaxation_time() {
        return relaxation_time_;
    }
    void step(size_t timestep) {
        if (timestep >= relaxation_time_) {
            relaxation_time_ *= factor_;
        }
    }
};

class Equilibrium {
  private:
    double factor_;
    double relaxation_time_;
    double target_temperature_;
    double timestep_;
    double budget_;

  public:
    Equilibrium(double factor, double relaxation_time,
                double target_temperature, double timestep, double budget)
        : factor_(factor),
          relaxation_time_(relaxation_time),
          target_temperature_(target_temperature),
          timestep_(timestep),
          budget_(budget) {}
    
    void step(Atoms& atoms, size_t timestep, double temp) {
        if (timestep < budget_ / factor_) {
            berendsen_thermostat(atoms, target_temperature_, timestep_,
                                 relaxation_time_, temp);
        }
    }
};

class EnergyPump {
  private:
    size_t interval_;
    double delta_Q_;
    double total_Q_ = 0;
    bool relaxed_ = false;

    void deposit_energy(Atoms& atoms, double ekin) {
        double lambda = std::sqrt(1 + delta_Q_ / ekin);
        atoms.velocities *= lambda;
        total_Q_ += delta_Q_;
    }

  public:
    EnergyPump(size_t interval, double delta_Q) : interval_(interval), delta_Q_(delta_Q) {}
    bool relaxed() { return relaxed_; }
    double total_Q() { return total_Q_; }

    void step(Atoms& atoms, size_t timestep, double ekin) {
        if (timestep % (interval_ * 2) == 0) {
            deposit_energy(atoms, ekin);
            relaxed_ = false;
        } else if (timestep % interval_ == 0) {
            relaxed_ = true;
        }
    }
};

#endif // __THERMOSTAT_H
