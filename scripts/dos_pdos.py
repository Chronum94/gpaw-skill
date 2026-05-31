#!/usr/bin/env python3
"""
Density of states (total + projected) from a converged GPAW calculation.

Usage:
    python dos_pdos.py [gpw_file]

Requires a converged .gpw file (from bulk_scf.py or band_structure.py).
For better PDOS resolution, use a denser k-mesh in the ground state.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from gpaw.dft import GPAW
from gpaw.dos import DOSCalculator


def compute_dos(
    gpw_file: str = 'bulk.gpw',
    emin: float = -10.0,      # eV relative to Fermi level
    emax: float =   5.0,
    npts: int   =  500,
    width: float = 0.1,       # Gaussian broadening eV; 0.0 → tetrahedron
    pdos_atoms: list = None,  # atom indices for PDOS; None = atom 0 only
    pdos_l: list = None,      # angular momenta (0=s,1=p,2=d,3=f)
    outfile: str = 'dos.png',
):
    calc = GPAW(gpw_file, txt=None)
    ef   = calc.get_fermi_level()
    nspins = calc.get_number_of_spins()

    energies = np.linspace(ef + emin, ef + emax, npts)
    dos_calc = DOSCalculator.from_calculator(calc)

    fig, ax = plt.subplots(figsize=(7, 4))

    # Total DOS (per spin)
    for spin in range(nspins):
        total = dos_calc.raw_dos(energies, spin=spin, width=width)
        label = 'Total DOS' if nspins == 1 else f'Total DOS (spin {spin})'
        sign  = 1 if spin == 0 else -1
        ax.plot(energies - ef, sign * total, label=label)

    # PDOS
    if pdos_atoms is None:
        pdos_atoms = [0]
    if pdos_l is None:
        pdos_l = [2]          # d-orbital by default

    l_names = {0: 's', 1: 'p', 2: 'd', 3: 'f'}
    atoms = calc.get_atoms()

    for a in pdos_atoms:
        symbol = atoms[a].symbol
        for l in pdos_l:
            for spin in range(nspins):
                pdos = dos_calc.raw_pdos(energies, a=a, l=l, spin=spin, width=width)
                sign = 1 if spin == 0 else -1
                lbl  = f'{symbol}({a}) {l_names[l]}'
                if nspins == 2:
                    lbl += f' ↑' if spin == 0 else ' ↓'
                ax.plot(energies - ef, sign * pdos, '--', label=lbl)

    ax.axvline(0, color='k', lw=0.8, ls='--', label='$E_F$')
    ax.set_xlabel('Energy − $E_F$ (eV)')
    ax.set_ylabel('DOS (states/eV)')
    ax.legend(fontsize=8)
    ax.set_xlim(emin, emax)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    print(f'DOS plot saved → {outfile}')


if __name__ == '__main__':
    gpw = sys.argv[1] if len(sys.argv) > 1 else 'bulk.gpw'
    compute_dos(gpw_file=gpw)
