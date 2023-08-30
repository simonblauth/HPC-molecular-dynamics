import argparse
import math
import matplotlib.pyplot as plt
import os
import os.path as osp
import tempfile
import time
from plot_csv import plot_total_energy
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
    plot_total_energy(filepaths, labels)


def simulation_times(simulators: list[str], labels: Optional[list[str]] = None):
    if labels is None:
        labels = [f"simulator {i}" for i in range(len(simulators))]
    for label, simulator in zip(labels, simulators):
        simulation_time(simulator, label)
    # plt.semilogy()
    plt.legend()
    plt.xlabel("Number of Atoms")
    plt.ylabel("Simulation Time per Timestep in seconds")
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate and Plot.")
    parser.add_argument("--mode", type=str, choices=["total_energy", "simulation_time"], required=True)
    parser.add_argument("--simulators", nargs='+', type=str, required=True)
    parser.add_argument("--labels", nargs='+', type=str)
    args = parser.parse_args()
    print(args.simulators)
    match args.mode:
        case "total_energy":
            total_energy(args.simulators[0])
        case "simulation_time":
            simulation_times(args.simulators, args.labels)
