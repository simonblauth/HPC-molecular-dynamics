#ifndef __THERMOSTAT_H
#define __THERMOSTAT_H

#include "atoms.h"

// Implementation of the https://en.wikipedia.org/wiki/Berendsen_thermostat.
void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep,
                          double relaxation_time);


#endif // __THERMOSTAT_H
