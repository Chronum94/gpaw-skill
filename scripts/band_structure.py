#!/usr/bin/env python3
"""
Two-step band structure calculation.

Step 1: converged SCF on k-mesh  → saves Si_gs.gpw
Step 2: fixed-density band path  → saves Si_bs.gpw, bandstructure.png

Usage:
    python band_structure.py

Modify SYSTEM, PATH, NPOINTS, NBANDS below for your material.
"""

from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

# ── configuration ────────────────────────────────────────────────────────────
ELEMENT   = 'Si'
STRUCTURE = 'diamond'
A         = 5.43          # Å; set to None to use ASE default
ECUT      = 500           # eV
KGRID     = (8, 8, 8)
SMEARING  = 0.01          # eV (use 0.1 for metals)
XC        = 'PBE'

KPATH     = 'GXWKL'      # high-symmetry path label (ASE/spglib notation)
NPOINTS   = 100           # k-points along path
NBANDS    = 16            # total bands for band structure
CONV_BANDS = 8            # how many bands to converge in step 2
# ─────────────────────────────────────────────────────────────────────────────


def step1_ground_state():
    """Converged SCF on Monkhorst-Pack grid."""
    kw = dict(a=A) if A else dict()
    atoms = bulk(ELEMENT, STRUCTURE, **kw)

    calc = GPAW(
        mode=dict(name='pw', ecut=ECUT),
        xc=XC,
        kpts=dict(size=KGRID, gamma=True),
        occupations=FermiDirac(SMEARING),
        convergence=dict(energy=5e-4, density=1e-4, eigenstates=4e-8),
        txt='gs.txt',
    )
    atoms.calc = calc
    atoms.get_potential_energy()
    ef = calc.get_fermi_level()
    calc.write('gs.gpw')
    print(f'[step 1] Fermi level: {ef:.4f} eV  →  gs.gpw')
    return ef


def step2_band_structure():
    """Non-self-consistent calculation along high-symmetry path."""
    bs_calc = GPAW('gs.gpw').fixed_density(
        nbands=NBANDS,
        symmetry='off',                          # required for complete k-path
        kpts=dict(path=KPATH, npoints=NPOINTS),
        convergence=dict(bands=CONV_BANDS),
        txt='bs.txt',
    )
    bs = bs_calc.band_structure()
    bs.plot(filename='bandstructure.png', show=False, emax=10.0, emin=-6.0)
    bs_calc.write('bs.gpw')
    print(f'[step 2] Band structure saved → bs.gpw, bandstructure.png')
    return bs


if __name__ == '__main__':
    step1_ground_state()
    step2_band_structure()
