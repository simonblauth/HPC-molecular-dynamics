import argparse
import math
import matplotlib.pyplot as plt
import os
import os.path as osp
import tempfile
import time
from plot_csv import make_plot
from typing import Optional


def total_energy(args):
    simulator = args.simulators[0]
    outdir = tempfile.gettempdir() if args.output_dir is None else args.output_dir
    timesteps = [0.01, 0.02, 0.03, 0.04]
    labels = [f"timestep={ts}" for ts in timesteps]
    filepaths = []
    for ts in timesteps:
        filename = f"ts-{str(ts).replace('.', '')}.csv"
        outpath = osp.join(outdir, "hpc_sims", filename)
        filepaths.append(outpath)
        os.makedirs(osp.dirname(outpath), exist_ok=True)
        os.system(f"{simulator} --max_timesteps 100 --silent --timestep {ts} --output_interval 1 --csv {outpath}")
    if not args.sim_only:
        make_plot(filepaths, labels,
                  x="Timestep", y="Total Energy",
                  xlabel="Timestep", ylabel="Total Energy"
                  )


def total_energy_domains(args):
    simulator = args.simulators[0]
    timesteps = args.timesteps
    outdir = tempfile.gettempdir() if args.output_dir is None else args.output_dir
    domain_decompositions = [(1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]
    ndomains = [math.prod(dd) for dd in domain_decompositions]
    labels = [f"processes={nd}, decomposition={dd}" for nd, dd in zip(ndomains, domain_decompositions)]
    input_file = osp.join("milestones", "07", "clusters", "cluster_923.xyz")
    filepaths = []
    for nd, dd in zip(ndomains, domain_decompositions):
        filename = f"domains-{nd}.csv"
        outpath = osp.join(outdir, "hpc_sims", filename)
        filepaths.append(outpath)
        os.makedirs(osp.dirname(outpath), exist_ok=True)
        command = f"mpirun -n {nd} --oversubscribe {simulator} -i {input_file} --mass 197 --timestep 1 --max_timesteps {timesteps} --output_interval 100 --cutoff 10.0 --domains {dd[0]} {dd[1]} {dd[2]} --shift_atoms 0.1 --deposit_energy 0.0 --csv {outpath} --silent"
        print(f"Simulating with {nd} workers.")
        os.system(command)
    if not args.sim_only:
        make_plot(filepaths, labels,
                  x="Timestep", y="Total Energy",
                  xlabel="Timestep", ylabel="Total Energy [eV]"
                  )


def simulation_times(args):
    simulators = args.simulators
    labels = args.labels
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


def temperature_over_energy(args):
    simulator = args.simulators[0]
    input_files = args.input_files
    outdir = tempfile.gettempdir() if args.output_dir is None else args.output_dir
    filepaths = []
    labels = [osp.basename(path).split('.')[0] for path in filepaths]
    for infile in input_files:
        basename = osp.basename(infile).split('.')[0]
        filename = basename + ".csv"
        outpath = osp.join(outdir, "hpc_sims", filename)
        filepaths.append(outpath)
        os.makedirs(osp.dirname(outpath), exist_ok=True)
        # TODO: fix deposit energies
        os.system(f"{simulator} -i {infile} --cutoff 7 --max_timesteps 50000 --output_interval 100 --timestep 1 --mass 197 --relaxation_time 100 --deposit_energy 0.5 --relaxation_time_deposit 20 --csv {outpath} --initial_relaxation 15000 --temperature 100 --smoothing 0.008 --thermostat_factor 10.0")
    if not args.sim_only:
        make_plot(filepaths, labels, "Total Energy", "Temperature", ylabel="Temperature [K]", xlabel="Total Energy [eV]", zero_x=True)


def stress_strain(args):
    target_strains = args.target_strains
    input_files = args.input_files
    outdir = tempfile.gettempdir() if args.output_dir is None else args.output_dir
    temperatures = [0, 100]
    strain_steps = [0.1, 0.05]
    strain_steps_label = ['01', '005']
    strain_interval = 100
    for target_strain, infile in zip(target_strains, input_files):
        filepaths = []
        labels = []
        for t in temperatures:
            for li, ss in zip(strain_steps_label, strain_steps):
                labels.append(f"initial_temperature={t}, length_increase={ss}")
                filename = f"{osp.basename(infile).split('.')[0]}_t{t}_li{li}"
                outpath_csv = osp.join(outdir, filename + ".csv")
                outpath_xyz = osp.join(outdir, filename + ".xyz")
                filepaths.append(outpath_csv)
                os.makedirs(osp.dirname(outpath_csv), exist_ok=True)
                timesteps = int(strain_interval * target_strain / ss)
                command = f"mpirun -n 6 ./build/milestones/09/09 -i {infile} --mass 197 --timestep 1 --max_timesteps {timesteps} --output_interval 100 --cutoff 7.0 --domains 1 1 6 --periodic 0 0 1 --shift_atoms 0.1 --smoothing 0.01 --stretch {ss} --stretch_interval {strain_interval} --temperature {t} --relaxation_time 100 --initial_relaxation 10000 --thermostat_factor 10.0 --csv {outpath_csv} --traj {outpath_xyz}"
                os.system(command)
                # print(command)
        if not args.sim_only:
            make_plot(filepaths, labels, "Strain", "Stress", "Strain [Å]", "Stress [eV/Å^3]")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate and Plot.")
    parser.add_argument("--mode", type=str, choices=["total_energy", "total_energy_domains", "simulation_time", "temperature_over_energy", "stress_strain"], required=True)
    parser.add_argument("--simulators", nargs='+', type=str)
    parser.add_argument("--labels", nargs='+', type=str)
    parser.add_argument("--input_files", nargs="+",
                        help="Path(s) to the xyz file(s) containing the data for the simulation.")
    parser.add_argument("--target_strains", nargs='+', type=int)
    parser.add_argument("--output_dir", type=str, help="Directory to write the output files. Uses temp directory by default")
    parser.add_argument("--timesteps", type=int, default=10000)
    parser.add_argument("--sim_only", default=False, action='store_true')
    args = parser.parse_args()
    # print(args.simulators)
    match args.mode:
        case "total_energy":
            total_energy(args)
        case "simulation_time":
            simulation_times(args)
        case "total_energy_domains":
            total_energy_domains(args)
        case "temperature_over_energy":
            temperature_over_energy(args)
        case "stress_strain":
            stress_strain(args)
