import os
from plot_energy import plot_total_energy

if __name__ == "__main__":
    timesteps = [0.01, 0.02, 0.03, 0.04]
    labels = [f"timestep={ts}" for ts in timesteps]
    filepaths = []
    for ts in timesteps:
        outfile = f"/tmp/hpc-sims/ts-{str(ts).replace('.', '')}.csv"
        filepaths.append(outfile)
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        filepath = os.path.dirname(os.path.abspath(__file__))
        os.system(f"{filepath}/../../build/milestones/04/04 --max_timesteps 100 --silent --timestep {ts} --output_interval 1 --csv {outfile}")
    plot_total_energy(filepaths, labels)
