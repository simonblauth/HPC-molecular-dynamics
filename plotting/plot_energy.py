import argparse
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional


def plot_from_csv(filepaths: list[str], labels: list[str], x: str, y: str):
    for label, filepath in zip(labels, filepaths):
        df = pd.read_csv(filepath)
        x_data = df[x]
        y_data = df[y]
        plt.plot(x_data, y_data, label=label)


def make_plot(filepaths: list[str], labels: list[str], x: str, y: str, xlabel: Optional[str], ylabel: Optional[str]):
    if xlabel is None:
        xlabel = x
    if ylabel is None:
        ylabel = y
    plot_from_csv(filepaths, labels, x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()


def plot_total_energy(filepaths: list[str], labels: list[str]):
    make_plot(filepaths, labels,
              x="Timestep", y="Total Energy",
              xlabel="Timestep", ylabel="Total Energy"
              )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot from csv files.")
    parser.add_argument("--csv_files", nargs="+",
                        help="Path(s) to the csv file(s) containing the data from the simulation.", required=True)
    parser.add_argument("--labels", nargs="*",
                        help="Labels for the plots made from the csv files.")
    parser.add_argument("--xlabel", type=str,
                        help="Label for the x axis")
    parser.add_argument("--ylabel", type=str,
                        help="Label for the y axis")
    parser.add_argument("--xdata", type=str,
                        help="Data field for the x axis", required=True)
    parser.add_argument("--ydata", type=str,
                        help="Data field for the y axis", required=True)
    args = parser.parse_args()
    filepaths = args.csv_files
    if args.labels is not None and len(args.labels) >= len(filepaths):
        labels = args.labels
    else:
        labels = [f"label_{i}" for i in range(len(filepaths))]
    make_plot(filepaths, labels, args.xdata,
              args.ydata, args.xlabel, args.ylabel)
