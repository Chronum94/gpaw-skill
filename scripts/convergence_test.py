#!/usr/bin/env python3
"""
Systematic convergence tests for plane-wave cutoff and k-point density.

Usage:
    python convergence_test.py --test ecut   # ecut convergence
    python convergence_test.py --test kpts   # k-point convergence
    python convergence_test.py --test both   # both (sequential)

Edit SYSTEM section below before running.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac


# ── system configuration ─────────────────────────────────────────────────────
ELEMENT   = 'Fe'
STRUCTURE = 'bcc'
MAGNETIC  = True
MAGMOM    = 2.2
XC        = 'PBE'
SMEARING  = 0.1       # eV

# Fixed parameters for k-point test
REF_ECUT  = 500       # eV

# Fixed parameters for ecut test
REF_KGRID = (8, 8, 8)

# Ranges to scan
ECUT_VALUES = [300, 400, 500, 600, 700]
KGRID_VALUES = [4, 6, 8, 10, 12]
# ─────────────────────────────────────────────────────────────────────────────


def make_atoms():
    atoms = bulk(ELEMENT, STRUCTURE)
    if MAGNETIC:
        atoms.set_initial_magnetic_moments([MAGMOM] * len(atoms))
    return atoms


def run_single(ecut, kgrid, label):
    atoms = make_atoms()
    calc = GPAW(
        mode=dict(name='pw', ecut=ecut),
        xc=XC,
        kpts=dict(size=kgrid, gamma=True),
        occupations=FermiDirac(SMEARING),
        spinpol=MAGNETIC,
        nbands=-10,
        mixer=dict(backend='msr1', beta=0.02, nmaxold=10, weight=100),
        convergence=dict(energy=5e-4, density=1e-5),
        txt=f'conv_{label}.txt',
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    natoms = len(atoms)
    print(f'  {label:20s}  E = {e:.6f} eV  ({e/natoms:.6f} eV/atom)')
    return e / natoms


def test_ecut():
    print(f'\n=== Ecut convergence (k-grid {REF_KGRID}) ===')
    results = {}
    for ecut in ECUT_VALUES:
        results[ecut] = run_single(ecut, REF_KGRID, f'ecut{ecut}')

    e_ref = results[ECUT_VALUES[-1]]
    print('\n  ecut (eV)   ΔE vs highest (meV/atom)')
    for ecut, e in results.items():
        print(f'  {ecut:6d}      {1000*(e - e_ref):+.2f}')

    plt.figure()
    ecuts = list(results.keys())
    energies = list(results.values())
    plt.plot(ecuts, np.array(energies) - e_ref, 'o-')
    plt.axhline(0, color='k', lw=0.5)
    plt.axhline(0.005, color='r', ls='--', label='5 meV/atom')
    plt.axhline(-0.005, color='r', ls='--')
    plt.xlabel('Plane-wave cutoff (eV)')
    plt.ylabel('ΔE (eV/atom) vs highest ecut')
    plt.legend()
    plt.savefig('convergence_ecut.png', dpi=150)
    print('Plot → convergence_ecut.png')


def test_kpts():
    print(f'\n=== k-point convergence (ecut {REF_ECUT} eV) ===')
    results = {}
    for k in KGRID_VALUES:
        results[k] = run_single(REF_ECUT, (k, k, k), f'k{k}')

    e_ref = results[KGRID_VALUES[-1]]
    print('\n  k-grid   ΔE vs densest (meV/atom)')
    for k, e in results.items():
        print(f'  {k}×{k}×{k}     {1000*(e - e_ref):+.2f}')

    plt.figure()
    ks = list(results.keys())
    energies = list(results.values())
    plt.plot(ks, np.array(energies) - e_ref, 'o-')
    plt.axhline(0, color='k', lw=0.5)
    plt.axhline(0.005, color='r', ls='--', label='5 meV/atom')
    plt.axhline(-0.005, color='r', ls='--')
    plt.xlabel('k-grid size (N×N×N)')
    plt.ylabel('ΔE (eV/atom) vs densest grid')
    plt.legend()
    plt.savefig('convergence_kpts.png', dpi=150)
    print('Plot → convergence_kpts.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--test', choices=['ecut', 'kpts', 'both'], default='both')
    args = parser.parse_args()

    if args.test in ('ecut', 'both'):
        test_ecut()
    if args.test in ('kpts', 'both'):
        test_kpts()
