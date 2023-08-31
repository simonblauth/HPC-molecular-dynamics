#! /usr/bin/env python

from math import sqrt

import numpy as np

import ase.io as io
from ase.lattice.cubic import FaceCenteredCubic


def make_whisker(radius, size, name):
    a = FaceCenteredCubic('Au', directions=[[1, 0, -1], [0, 1, 0], [1, 0, 1]],
                      size=size)
    c = a.cell.diagonal() / 2

    dir1 = [sqrt(2), 1, 0]
    dir2 = [sqrt(2), -1, 0]
    dir3 = [0, 1, 0]

    dir1 = np.array(dir1) / np.linalg.norm(dir1)
    dir2 = np.array(dir2) / np.linalg.norm(dir2)
    dir3 = np.array(dir3) / np.linalg.norm(dir3)

    r = a.get_positions() - c
    m = np.abs(r.dot(dir1)) > radius
    m = np.logical_or(m, np.abs(r.dot(dir2)) > radius)
    m = np.logical_or(m, np.abs(r.dot(dir3)) > radius)
    del a[m]

    a.center()
    a.set_pbc([False, False, True])

    io.write(f'{name}.xyz', a)
    io.write(f'{name}.data', a, format='lammps-data')


if __name__ == "__main__":
    # radius = 15
    # size = [14, 10, 50]  # 6590 atoms

    # radius = 30
    # size = [28, 20, 100]  # 51500 atoms
    radii = [15]
    sizes = [[14, 10, 50]]
    for r, sz in zip(radii, sizes):
        make_whisker(r, sz, f"whisker_radius_{r}")
