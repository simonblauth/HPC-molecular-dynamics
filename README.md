# Molecular Dynamics Code for IMTEK-HPC-MD class

This repository contains code for the [HPC MD
with C++
project](https://imtek-simulation.github.io/MolecularDynamics/_project/general_remarks.html).


## Setup & Installation
First you need to make sure that `MPI` is installed on your system. If you do not have MPI, you need to set the `-DUSE_MPI=0` flag instead. This will only compile targets that do not require MPI.

Setting up the `c++` code:
```bash
cd <repository>

# Create build directory
mkdir build
cd build

# Compile with MPI
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=1 ..
make

# Run tests
make test
```
Setting up a `python` environment:
```bash
conda create -n hpc-molecular-dynamics -y
conda activate hpc-molecular-dynamics
pip install -r requirements.txt
```
The code was tested with python 3.11.3, however it should work for all python versions >=3.9

## Running simulations

## Reproducing Figures
Here I provide instructions on how to reproduce each figure from my report. For some simulations that take extremely long, the csv files from which the plots are generated are provided in this repo. However if you want to run the simulation yourself there are instructions for that as well. Note that the trajectory files are not present in this repo as they tend to grow very large for bigger simulations. 

In the provided commands, if you find expressions in `<braces>` that expression has to be provided by the user. Usually used for e.g. output file locations. All commands are assumed to be run from the root folder of this repository. Also make shure the code is compiled and the python environment is activated.

### Figure 1
```bash
python plotting/simulate_and_plot.py --mode total_energy --simulators ./build/milestones/04/04
```

### Figure 2
```bash
python plotting/simulate_and_plot.py --mode simulation_time --simulators ./build/milestones/05/05 ./build/milestones/06/06 --labels lj_naive lj_cutoff
```

### Figure 3
```bash
python plotting/simulate_and_plot.py --mode total_energy_domains --simulators ./build/milestones/08/08
```
Note that this takes a while since each simulation runs for 10000 timesteps. To get faster results, you can lower the timesteps with the `--timesteps` argument. For example the following only takes about a minute on my machine:
```bash
python plotting/simulate_and_plot.py --mode total_energy_domains --simulators ./build/milestones/08/08 --timesteps 1000
```

### Figure 4
To generate the trajectory file run the following. Make sure to specify the output files. You can then inspect the `.xyz` file in OVITO ath the timesteps mentioned in the report. The trajectory is only written every 100 timesteps, so to look at the 5000th timestep you have to look at frame 50.
```bash
./build/milestones/05/05 --lattice_size 5 --max_timesteps 20000 --timestep 1e-3 --traj <output_file.xyz> --csv <output_file.csv> --thermostat_factor 1 --relaxation_time 10 --temperature 500 --silent
```
To verify that the system remains stable without thermostat you can increase the `--thermostat_factor` to disable the thermostat after some time. If you remove the `--silent` argument you can then see the temperature development in the console.

### Figure 5
To plot the data from the provided csv files:
```bash
python plotting/plot_csv.py --csv_files milestones/07/outputs/cluster_*.csv --xdata "Total Energy" --ydata "Temperature" --zero_x --xlabel "Total Energy Increase [eV]" --ylabel "Temperature [K]" 
```
To simulate the data youself:
```bash
./build/milestones/07/07 -i .milestones/07/clusters/cluster_<cluster_size>.xyz --traj <output_file.xyz> --csv <output_file.csv> --cutoff 7 --max_timesteps 50000 --output_interval 100 --timestep 1 --mass 197 --relaxation_time 100 --deposit_energy <energy> --relaxation_time_deposit 20 --initial_relaxation 15000 --temperature 100 --smoothing 0.008 --thermostat_factor 10.0
```
Note that this only simulates one of the clusters, so you would have to run it for all five clusters. The `--deposit_energy` parameter was set to `[0.5, 0.75, 1.0, 1.5, 2.0]` for the cluster sizes `[923, 1415, 2057, 2869, 3871]` respectively. To plot the data afterwards use the first command but adjust the `--csv_files` argument according to where your output data is.

## Figure 6
To plot the data from the provided csv files:
```bash
python plotting/plot_csv.py --csv_files milestones/07/outputs/cluster_923.csv milestones/07/outputs/cluster_1415.csv milestones/07/outputs/cluster_2057.csv  milestones/07/outputs/cluster_2869.csv milestones/07/outputs/cluster_3871.csv --mode heat_capacity
```
If you want to simulate the data yourself follow the instructions from Figure 5 and adjust the `--csv_files` argument according to where your output data is.

## Figure 7
To plot the data from the provided csv files:
```bash
python plotting/plot_csv.py --csv_files milestones/07/outputs/cluster_923.csv milestones/07/outputs/cluster_1415.csv milestones/07/outputs/cluster_2057.csv  milestones/07/outputs/cluster_2869.csv milestones/07/outputs/cluster_3871.csv --mode latent_heat
```
If you want to simulate the data yourself follow the instructions from Figure 5 and adjust the `--csv_files` argument according to where your output data is.

## Figure 8
To plot the data from the provided csv files:
```bash
python plotting/plot_csv.py --csv_files milestones/07/outputs/cluster_923.csv milestones/07/outputs/cluster_1415.csv milestones/07/outputs/cluster_2057.csv  milestones/07/outputs/cluster_2869.csv milestones/07/outputs/cluster_3871.csv --mode melting_point
```
If you want to simulate the data yourself follow the instructions from Figure 5 and adjust the `--csv_files` argument according to where your output data is.

## Figure 9
To plot the data from the provided csv files:
```bash
python plotting/plot_csv.py --csv_files milestones/09/outputs/whisker_small_t0_li01.csv  milestones/09/outputs/whisker_medium_t0_li01.csv milestones/09/outputs/whisker_large_t0_li01.csv --xdata "Strain" --ydata "Stress" --xlabel "Strain [Å]" --ylabel "Stress [eV/Å^3]" --labels whisker_small whisker_medium whisker_large
```
To simulate the data yourself:
```bash
python plotting/simulate_and_plot.py --mode stress_strain --input_files ./milestones/09/whisker_small.xyz ./milestones/09/whisker_medium.xyz ./milestones/09/whisker_large.xyz --target_strains 25 25 40 --sim_only --output_dir <output_dir>
```


## Figure 10
To simulate the wire:
```bash
mpirun -n 6 ./build/milestones/09/09 -i ./milestones/09/whisker_large.xyz --traj <output_file.xyz> --csv <output_file.csv> --mass 197 --timestep 1 --max_timesteps 40000 --output_interval 100 --cutoff 7.0 --domains 1 1 6 --periodic 0 0 1 --shift_atoms 0.1 --smoothing 0.01 --stretch 0.1 --stretch_interval 100 --temperature 100 --relaxation_time 100 --initial_relaxation 10000 --thermostat_factor 10.0
```
Afterwards look at the `.xyz` file in OVITO and turn on common-neighbor-analysis.

## Figure 11
To plot the data from the provided csv files:
```bash
python plotting/plot_csv.py --csv_files milestones/09/outputs/whisker_small_t0_li01.csv milestones/09/outputs/whisker_large_t0_li01.csv milestones/09/outputs/whisker_small_t0_li005.csv milestones/09/outputs/whisker_large_t0_li005.csv --xdata "Strain" --ydata "Stress" --xlabel "Strain [Å]" --ylabel "Stress [eV/Å^3]" --labels "whisker_small, length_increase=0.1" "whisker_large, length_increase=0.1" "whisker_small, length_increase=0.05" "whisker_large, length_increase=0.05" 
```

## Figure 12
To plot the data from the provided csv files:
```bash
python plotting/plot_csv.py --csv_files milestones/09/outputs/whisker_small_t100_li005.csv  milestones/09/outputs/whisker_medium_t100_li005.csv milestones/09/outputs/whisker_large_t100_li005.csv --xdata "Strain" --ydata "Stress" --xlabel "Strain [Å]" --ylabel "Stress [eV/Å^3]" --labels whisker_small whisker_medium whisker_large
```



