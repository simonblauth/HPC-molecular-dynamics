#include "thermostat.h"

void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep,
                          double relaxation_time) {
    double current_temperature = atoms.current_temperature();
    berendsen_thermostat(atoms, target_temperature, timestep, relaxation_time,
                         current_temperature);
}

void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep,
                          double relaxation_time, double current_temperature) {
    double lambda =
        std::sqrt(1 + (target_temperature / current_temperature - 1) *
                          timestep / relaxation_time);
    atoms.velocities *= lambda;
}
