import argparse
import math
import matplotlib.pyplot as plt
import os
import os.path as osp
import tempfile
import time
from plot_csv import make_plot
from typing import Optional


def total_energy(simulator: str):
    timesteps = [0.01, 0.02, 0.03, 0.04]
    labels = [f"timestep={ts}" for ts in timesteps]
    filepaths = []
    for ts in timesteps:
        tmp = tempfile.gettempdir()
        filename = f"ts-{str(ts).replace('.', '')}.csv"
        outpath = osp.join(tmp, "hpc_sims", filename)
        filepaths.append(outpath)
        os.makedirs(osp.dirname(outpath), exist_ok=True)
        os.system(f"{simulator} --max_timesteps 100 --silent --timestep {ts} --output_interval 1 --csv {outpath}")
    make_plot(filepaths, labels,
              x="Timestep", y="Total Energy",
              xlabel="Timestep", ylabel="Total Energy"
              )


def total_energy_domains(simulator: str):
    domain_decompositions = [(1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]
    ndomains = [math.prod(dd) for dd in domain_decompositions]
    labels = [f"processes={nd}, decomposition={dd}" for nd, dd in zip(ndomains, domain_decompositions)]
    filepaths = []
    for nd, dd in zip(ndomains, domain_decompositions):
        tmp = tempfile.gettempdir()
        filename = f"domains-{nd}.csv"
        outpath = osp.join(tmp, "hpc_sims", filename)
        filepaths.append(outpath)
        os.makedirs(osp.dirname(outpath), exist_ok=True)
        command = f"mpirun -n {nd} --oversubscribe {simulator} -i ./milestones/07/cluster_923.xyz --mass 197 --timestep 1 --max_timesteps 10000 --output_interval 100 --cutoff 10.0 --domains {dd[0]} {dd[1]} {dd[2]} --shift_atoms 0.1 --deposit_energy 0.0 --csv {outpath} --silent"
        scriptpath = osp.join(tmp, "hpc_sims", "doit.sh")
        with open(scriptpath, "w") as f:
            f.write(command)
        os.system(f"sh {scriptpath}")
    make_plot(filepaths, labels,
              x="Timestep", y="Total Energy",
              xlabel="Timestep", ylabel="Total Energy [eV]"
              )


def simulation_times(simulators: list[str], labels: Optional[list[str]] = None):
    if labels is None:
        labels = [f"simulator {i}" for i in range(len(simulators))]
    for label, simulator in zip(labels, simulators):
        simulation_time(simulator, label)
    # plt.semilogy()
    plt.legend()
    plt.xlabel("Number of Atoms")
    plt.ylabel("Simulation Time per Timestep [s]")
    plt.show()


def simulation_time(simulator: str, label):
    lattice_sizes = list(range(2, 11))
    nb_atoms = [math.pow(ls, 3) for ls in lattice_sizes]
    simulation_times = []
    timesteps = 10
    cutoff = 3
    for ls in lattice_sizes:
        start = time.time()
        os.system(f"{simulator} --max_timesteps {timesteps} --silent --lattice_size {ls} --cutoff {cutoff}")
        stop = time.time()
        simulation_times.append((stop - start) / timesteps)
    plt.plot(nb_atoms, simulation_times, label=label)


def temperature_over_energy(simulator: str, input_files: list[str]):
    filepaths = []
    labels = [osp.basename(path).split('.')[0] for path in filepaths]
    for infile in input_files:
        tmp = tempfile.gettempdir()
        filename = osp.basename(infile).split('.')[0] + ".csv"
        outpath = osp.join(tmp, "hpc_sims", filename)
        filepaths.append(outpath)
        os.makedirs(osp.dirname(outpath), exist_ok=True)
        # TODO: fix deposit energies
        os.system(f"{simulator} -i {infile} --cutoff 7 --max_timesteps 50000 --output_interval 100 --timestep 1 --mass 197 --relaxation_time 100 --deposit_energy 0.5 --relaxation_time_deposit 20 --csv {outpath} --initial_relaxation 15000 --temperature 100 --smoothing 0.008 --relaxation_time_increase 10.0")
    make_plot(filepaths, labels, "Total Energy", "Temperature", ylabel="Temperature [K]", xlabel="Total Energy [eV]", zero_x=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate and Plot.")
    parser.add_argument("--mode", type=str, choices=["total_energy", "total_energy_domains", "simulation_time", "temperature_over_energy"], required=True)
    parser.add_argument("--simulators", nargs='+', type=str, required=True)
    parser.add_argument("--labels", nargs='+', type=str)
    args = parser.parse_args()
    print(args.simulators)
    match args.mode:
        case "total_energy":
            total_energy(args.simulators[0])
        case "simulation_time":
            simulation_times(args.simulators, args.labels)
        case "total_energy_domains":
            total_energy_domains(args.simulators[0])
        case "temperature_over_energy":
            # TODO: input files arg
            temperature_over_energy(args.simulators[0], args.input_files)
