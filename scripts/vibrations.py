#!/usr/bin/env python3
"""
Vibrational frequency calculation for a molecule or cluster.

Workflow:
  1. Relax geometry to tight fmax
  2. Compute finite-difference Hessian
  3. Report frequencies and compare to reference

Usage:
    python vibrations.py

Edit SYSTEM section below for your molecule.
"""

from math import cos, sin, pi
from ase import Atoms
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from gpaw import GPAW


# ── system configuration ─────────────────────────────────────────────────────
def build_h2o():
    d, t = 0.9575, pi / 180 * 104.51
    mol = Atoms(
        'H2O',
        positions=[(0, 0, 0), (d, 0, 0), (d * cos(t), d * sin(t), 0)],
    )
    mol.center(vacuum=4.0)
    # break cubic cell symmetry to prevent Poisson solver issues
    mol.cell[1, 1] += 0.0001
    mol.cell[2, 2] += 0.0002
    return mol

# Experimental reference (cm⁻¹) for H2O
REFERENCE = {
    'bending':     1595,
    'sym_stretch': 3657,
    'asym_stretch': 3756,
}
# ─────────────────────────────────────────────────────────────────────────────


def run_vibrations(delta: float = 0.02):
    mol = build_h2o()

    calc = GPAW(
        mode='lcao',
        basis='dzp',
        xc='PBE',
        convergence=dict(density=1e-6),
        symmetry='off',                   # mandatory for vibrational analysis
        txt='h2o.txt',
    )
    mol.calc = calc

    # Tight geometry optimization before vibrational analysis
    qn = QuasiNewton(mol, logfile='opt.log')
    qn.run(fmax=0.005)                    # 0.005 eV/Å, not 0.05
    print(f'Optimized: {mol.get_positions()}')

    # Finite-difference Hessian
    vib = Vibrations(mol, delta=delta)
    vib.run()
    vib.summary(method='frederiksen')    # Frederiksen damping for imaginary modes
    vib.write_jmol()                     # JMol trajectory for visualization

    freqs = vib.get_frequencies()
    real_freqs = freqs[freqs.imag == 0].real

    print('\nVibrational frequencies (cm⁻¹):')
    for i, f in enumerate(freqs):
        flag = '  ← imaginary' if f.imag != 0 else ''
        print(f'  mode {i:2d}:  {f.real:8.1f}{flag}')

    if REFERENCE:
        print('\nComparison with experiment:')
        for mode, ref in REFERENCE.items():
            print(f'  {mode:15s}:  expt {ref} cm⁻¹')

    return freqs


if __name__ == '__main__':
    run_vibrations()
