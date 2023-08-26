import argparse
import math
import matplotlib.pyplot as plt
import os
import os.path as osp
import tempfile
import time
from plot_energy import plot_total_energy


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


def simulation_time(simulator: str):
    lattice_sizes = list(range(2, 10))
    nb_atoms = [math.pow(ls, 3) for ls in lattice_sizes]
    simulation_times = []
    timesteps = 10
    for ls in lattice_sizes:
        start = time.time()
        os.system(f"{simulator} --max_timesteps {timesteps} --silent --lattice_size {ls}")
        stop = time.time()
        simulation_times.append((stop - start) / timesteps)
    plt.plot(nb_atoms, simulation_times)
    # plt.semilogy()
    plt.xlabel("Number of Atoms")
    plt.ylabel("Simulation Time per Timestep in seconds")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate and Plot.")
    parser.add_argument("--mode", type=str, choices=["total_energy", "simulation_time"], required=True)
    parser.add_argument("--simulator", type=str, required=True)
    args = parser.parse_args()
    match args.mode:
        case "total_energy":
            total_energy(args.simulator)
        case "simulation_time":
            simulation_time(args.simulator)
