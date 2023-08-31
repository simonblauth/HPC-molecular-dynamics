import argparse
import os.path as osp
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional


def zero_shift(df: pd.DataFrame, col):
    df[col] -= df[col][0]


def slope(df: pd.DataFrame, x, y, meltpoint_lower_bound, meltpoint_upper_bound):
    dps = df.drop(df[(df[y] > meltpoint_lower_bound) &
                  (df[y] < meltpoint_upper_bound)].index)
    slopes = ((dps[y] - dps[y][0])[1:] / (dps[x] - dps[x][0]))[1:]
    return slopes.mean()


def heat_capacity(df: pd.DataFrame, meltpoint_lower_bound, meltpoint_upper_bound):
    return 1 / slope(df, "Total Energy", "Temperature", meltpoint_lower_bound, meltpoint_upper_bound)


def latent_heat(df: pd.DataFrame, meltpoint_lower_bound, meltpoint_upper_bound):
    m = slope(df, "Total Energy", "Temperature",
              meltpoint_lower_bound, meltpoint_upper_bound)
    zero_shift(df, "Total Energy")
    zero_shift(df, "Temperature")
    x = df["Total Energy"].tail(1).item()
    y = df["Temperature"].tail(1).item()
    y_0 = y - m * x
    return -y_0 / m


# checks if point (x, y) is below line with slope m and offset y_0
def isunder(x, y, m, y_0):
    return m * x + y_0 > y


def melting_point(df: pd.DataFrame, meltpoint_lower_bound, meltpoint_upper_bound):
    m = slope(df, "Total Energy", "Temperature",
              meltpoint_lower_bound, meltpoint_upper_bound)
    lh = latent_heat(df, meltpoint_lower_bound, meltpoint_upper_bound)
    start = lh / 2
    y_0 = -m * lh
    zero_shift(df, "Total Energy")
    zero_shift(df, "Temperature")
    for x, y in zip(df["Total Energy"], df["Temperature"]):
        if isunder(x, y, m, y_0):
            return y
    return df["Temperature"].tail(1).item()


def plot_heat_capacity(filepaths: list[str], labels: list[str]):
    cluster_sizes = [int(label.split('_')[1]) for label in labels]
    heat_capacities = [heat_capacity(pd.read_csv(filepath), 800, 1000) for filepath in filepaths]
    plt.plot(cluster_sizes, heat_capacities)
    plt.xlabel("Number of Atoms")
    plt.ylabel("Heat Capacity [eV/K]")
    plt.show()


def plot_latent_heat(filepaths: list[str], labels: list[str]):
    cluster_sizes = [int(label.split('_')[1]) for label in labels]
    latent_heats = [latent_heat(pd.read_csv(filepath), 800, 1000) for filepath in filepaths]
    plt.plot(cluster_sizes, latent_heats)
    plt.xlabel("Number of Atoms")
    plt.ylabel("Latent Heat [K]")
    plt.show()


def plot_melting_point(filepaths: list[str], labels: list[str]):
    cluster_sizes = [int(label.split('_')[1]) for label in labels]
    melting_points = [melting_point(pd.read_csv(filepath), 800, 1000) for filepath in filepaths]
    plt.plot(cluster_sizes, melting_points)
    plt.xlabel("Number of Atoms")
    plt.ylabel("Melting Point [K]")
    plt.show()


def plot_from_csv(filepaths: list[str], labels: list[str], x: str, y: str, zero_x=False, zero_y=False):
    for label, filepath in zip(labels, filepaths):
        df = pd.read_csv(filepath)
        x_data = df[x]
        if zero_x:
            x_data -= x_data[0]
        y_data = df[y]
        if zero_y:
            y_data -= y_data[0]
        plt.plot(x_data, y_data, label=label)


def make_plot(filepaths: list[str], labels: list[str], x: str, y: str, xlabel: Optional[str] = None, ylabel: Optional[str] = None, zero_x=False, zero_y=False):
    if xlabel is None:
        xlabel = x
    if ylabel is None:
        ylabel = y
    plot_from_csv(filepaths, labels, x, y, zero_x, zero_y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()


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
                        help="Data field for the x axis")
    parser.add_argument("--ydata", type=str,
                        help="Data field for the y axis")
    parser.add_argument("--zero_x", default=False, action='store_true',
                        help="Shift x axis to zero")
    parser.add_argument("--zero_y", default=False, action='store_true',
                        help="Shift y axis to zero")
    parser.add_argument("--mode", type=str, choices=[
                        "plain", "latent_heat", "heat_capacity", "melting_point"], default="plain")
    args = parser.parse_args()
    filepaths = args.csv_files
    if args.labels is not None and len(args.labels) >= len(filepaths):
        labels = args.labels
    else:
        labels = [osp.basename(path).split('.')[0] for path in filepaths]
    match args.mode:
        case "plain":
            make_plot(filepaths, labels, args.xdata,
                      args.ydata, args.xlabel, args.ylabel, args.zero_x, args.zero_y)
        case "latent_heat":
            plot_latent_heat(filepaths, labels)
        case "heat_capacity":
            plot_heat_capacity(filepaths, labels)
        case "melting_point":
            plot_melting_point(filepaths, labels)
