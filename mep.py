"""
(c) 2013, 2019, 2022 Marius Retegan
License: BSD-2-Clause
Description: Create a .cube file of the molecular electrostatic
             potential (MEP) using ORCA.
Run: python mep.py basename npoints nprocs (e.g. python mep.py water 40 2).
Arguments: basename - file name without the extension;
                      this should be the same for the .gbw and .scfp.
           npoints  - number of grid points per side
                      (80 should be fine)
           nprocs - number of processors for parallel calculations
"""

#!/usr/bin/env python


def read_xyz(xyz):
    atoms = []
    x, y, z = [], [], []

    with open(xyz) as fp:
        # Skip the first two lines.
        next(fp)
        next(fp)
        for line in fp:
            data = line.split()
            atoms.append(data[0])
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))

    return atoms, np.array(x), np.array(y), np.array(z)


def read_vpot(vpot):
    v = []

    with open(vpot) as fp:
        next(fp)
        for line in fp:
            data = line.split()
            v.append(float(data[3]))

    return np.array(v)


if __name__ == "__main__":
    import os
    import shutil
    import subprocess
    import sys

    import numpy as np

    ang_to_au = 1.0 / 0.5291772083


    elements = [None,
         "H", "He",
         "Li", "Be",
         "B", "C", "N", "O", "F", "Ne",
         "Na", "Mg",
         "Al", "Si", "P", "S", "Cl", "Ar",
         "K", "Ca",
         "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
         "Ga", "Ge", "As", "Se", "Br", "Kr",
         "Rb", "Sr",
         "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
         "In", "Sn", "Sb", "Te", "I", "Xe",
         "Cs", "Ba",
         "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
         "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
         "Tl", "Pb", "Bi", "Po", "At", "Rn",
         "Fr", "Ra",
         "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
         "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub"]


    basename = sys.argv[1]
    xyz = f"{basename}.xyz"

    if not os.path.isfile(xyz):
        sys.exit("Could not find the .xyz. To quickly generate one for "
                 f"your molecule run: echo 11 | orca_plot {basename}.gbw -i.")

    atoms, x, y, z = read_xyz(xyz)

    try:
        npoints = int(sys.argv[2])
    except ValueError:
        sys.exit(f"Invalid number of points: {sys.argv[2]}")

    try:
        nprocs = int(sys.argv[3])
    except IndexError:
        nprocs = 1
    except ValueError:
        sys.exit(f"Invalid number of cpus: {sys.argv[3]}")

    natoms = len(atoms)

    extent = 7.0
    xmin = x.min() * ang_to_au - extent
    xmax = x.max() * ang_to_au + extent
    ymin = y.min() * ang_to_au - extent
    ymax = y.max() * ang_to_au + extent
    zmin = z.min() * ang_to_au - extent
    zmax = z.max() * ang_to_au + extent

    with open(f"{basename}_mep.inp", "w") as fp:
        fp.write(f"{nprocs:d}\n")
        fp.write(f"{basename}.gbw\n")
        fp.write(f"{basename}.scfp\n")
        fp.write(f"{basename}_mep.xyz\n")
        fp.write(f"{basename}_mep.out\n")

    with open(f"{basename}_mep.xyz", "w") as fp:
        fp.write(f"{npoints**3:d}\n")
        for ix in np.linspace(xmin, xmax, npoints, True):
            for iy in np.linspace(ymin, ymax, npoints, True):
                for iz in np.linspace(zmin, zmax, npoints, True):
                    fp.write(f"{ix:12.6f} {iy:12.6f} {iz:12.6f}\n")

    orca_vpot = shutil.which("orca_vpot")
    if orca_vpot is None:
        sys.exit(f"Could not find the orca_vpot executable. Make sure it is added to your PATH.")

    subprocess.check_call([f"{orca_vpot}", f"{basename}_mep.inp"])

    vpot = read_vpot(f"{basename}_mep.out")

    with open(f"{basename}_mep.cube", "w") as fp:
        fp.write("Generated with ORCA\n")
        fp.write(f"Electrostatic potential for {basename}\n")
        fp.write(f"{len(atoms):5d}{xmin:12.6f}{ymin:12.6f}{zmin:12.6f}\n")
        
        xstep = (xmax - xmin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{xstep:12.6f}{0:12.6f}{0:12.6f}\n")
        
        ystep = (ymax - ymin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{0:12.6f}{ystep:12.6f}{0:12.6f}\n")
        
        zstep = (zmax - zmin) / float(npoints - 1)
        fp.write(f"{npoints:5d}{0:12.6f}{0:12.6f}{zstep:12.6f}\n")

        for i, atom in enumerate(atoms):
            index = elements.index(atom)
            xi, yi, zi = x[i] * ang_to_au, y[i] * ang_to_au, z[i] * ang_to_au
            fp.write(f"{index:5d}{0:12.6f}{xi:12.6f}{yi:12.6f}{zi:12.6f}\n")

        m = 0
        n = 0
        vpot = np.reshape(vpot, (npoints, npoints, npoints))
        for ix in range(npoints):
            for iy in range(npoints):
                for iz in range(npoints):
                    fp.write(f"{vpot[ix][iy][iz]:14.5e}")
                    m += 1
                    n += 1
                    if (n > 5):
                        fp.write("\n")
                        n = 0
                if n != 0:
                    fp.write("\n")
                    n = 0
