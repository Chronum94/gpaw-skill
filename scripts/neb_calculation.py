#!/usr/bin/env python3
"""
Nudged Elastic Band (NEB) calculation for diffusion / reaction barriers.

Run with MPI: mpiexec -np <N_cpus> gpaw python neb_calculation.py

Requirements:
    - initial.traj and final.traj: pre-relaxed endpoint structures
    - N_IMAGES * N_CPUS_PER_IMAGE MPI processes total

Adjust N_IMAGES and the GPAW calculator settings below.
"""

from ase.io import read
from ase.neb import NEB
from ase.optimize import BFGS
from gpaw.dft import GPAW, FermiDirac


# ── configuration ────────────────────────────────────────────────────────────
N_IMAGES = 3          # number of intermediate images
ECUT     = 500        # eV
KGRID    = (4, 4, 1)  # use (k,k,k) for bulk, (k,k,1) for surfaces
XC       = 'PBE'
SMEARING = 0.1        # eV

FMAX_COARSE = 0.10    # eV/Å — first pass without climbing image
FMAX_FINE   = 0.05    # eV/Å — second pass with climbing image
# ─────────────────────────────────────────────────────────────────────────────


def build_images():
    initial = read('initial.traj')
    final   = read('final.traj')

    images = [initial.copy()] + [initial.copy() for _ in range(N_IMAGES)] + [final.copy()]
    neb = NEB(images, climb=True, parallel=True)
    neb.interpolate(method='idpp')
    return images, neb


def attach_calculators(images):
    for image in images[1:-1]:
        image.calc = GPAW(
            mode=dict(name='pw', ecut=ECUT),
            xc=XC,
            kpts=dict(size=KGRID, gamma=True),
            occupations=FermiDirac(SMEARING),
            convergence=dict(
                energy=5e-3,          # looser SCF tolerance for NEB images
                eigenstates=1e-7,
            ),
            mixer=dict(backend='msr1', beta=0.05, nmaxold=8),
            symmetry='off',
            txt='-',
        )


def run_neb():
    images, neb = build_images()
    attach_calculators(images)

    # Pass 1: without climbing image (coarse)
    opt = BFGS(neb, logfile='neb_coarse.log', trajectory='neb_coarse.traj')
    opt.run(fmax=FMAX_COARSE)

    # Pass 2: activate climbing image (fine)
    neb.climb = True
    opt = BFGS(neb, logfile='neb_fine.log', trajectory='neb_fine.traj')
    opt.run(fmax=FMAX_FINE)

    # Report energies
    energies = [img.get_potential_energy() for img in images]
    e0 = energies[0]
    print('\nNEB path energies:')
    for i, e in enumerate(energies):
        print(f'  image {i:2d}:  {e - e0:+.4f} eV')

    barrier = max(energies) - energies[0]
    print(f'\nForward barrier:  {barrier:.4f} eV')
    print(f'Reverse barrier:  {max(energies) - energies[-1]:.4f} eV')


if __name__ == '__main__':
    run_neb()
