import argparse
import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Total Energy over time.")
    parser.add_argument("--csv_files", nargs="+", help="Path(s) to the csv file(s) containing the data from the simulation.")
    parser.add_argument("--labels", nargs="*", help="Labels for the plots made from the csv files.")
    args = parser.parse_args()
    filepaths = args.csv_files
    labels = args.labels if args.labels is not None and len(args.labels) >= len(filepaths) else [f"label_{i}" for i in range(len(filepaths))]
    for label, filepath in zip(labels, filepaths):
        df = pd.read_csv(filepath)
        x = df["Timestep"]
        y = df["Total Energy"]
        plt.plot(x, y, label=label)
    plt.xlabel("Timestep")
    plt.ylabel("Total Energy")
    plt.legend()
    plt.show()
