#ifndef __THERMOSTAT_H
#define __THERMOSTAT_H

#include "atoms.h"

// Implementation of the https://en.wikipedia.org/wiki/Berendsen_thermostat.
void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep,
                          double relaxation_time);
void berendsen_thermostat(Atoms &atoms, double target_temperature,
                          double timestep, double relaxation_time,
                          double current_temperature);

#endif // __THERMOSTAT_H
