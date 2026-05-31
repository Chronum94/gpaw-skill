#!/usr/bin/env python3
"""
Spin-polarized DFT+U calculation for strongly correlated systems.

Default system: antiferromagnetic NiO.

Usage:
    python magnetic_dftu.py

Edit SYSTEM section below for your material.
"""

from ase import Atoms
from gpaw import GPAW, PW, FermiDirac


# ── system configuration ─────────────────────────────────────────────────────
# NiO rock-salt antiferromagnetic (2-atom primitive cell)
A = 4.19
B = A / 2**0.5
M = 2.0          # initial magnetic moment on Ni (μB)
U_EFF = 6.0      # Hubbard U_eff (= U − J) for Ni-d in eV

def build_nio():
    return Atoms(
        'Ni2O2',
        pbc=True,
        cell=(B, B, A),
        positions=[(0,0,0), (B/2, B/2, A/2), (0,0,A/2), (B/2, B/2, 0)],
        magmoms=[M, -M, 0, 0],
    )
# ─────────────────────────────────────────────────────────────────────────────


def run_dftu(u_eff: float = U_EFF, outfile: str = 'nio_dftu'):
    atoms = build_nio()

    calc = GPAW(
        mode=dict(name='pw', ecut=600),
        xc='PBE',
        kpts=dict(size=(4, 4, 4), gamma=True),
        occupations=FermiDirac(0.05),
        spinpol=True,
        setups=dict(Ni=f':d,{u_eff}'),         # DFT+U on Ni d-orbitals
        nbands=-10,
        mixer=dict(
            backend='msr1',
            beta=0.02,
            nmaxold=10,
            weight=100,
            method='separate',
        ),
        convergence=dict(energy=5e-4, density=1e-5, eigenstates=4e-8),
        txt=f'{outfile}.txt',
    )
    atoms.calc = calc

    e = atoms.get_potential_energy()
    mm = calc.get_magnetic_moments()
    ef = calc.get_fermi_level()
    calc.write(f'{outfile}.gpw')

    print(f'U_eff = {u_eff:.1f} eV')
    print(f'Total energy: {e:.6f} eV')
    print(f'Fermi level:  {ef:.4f} eV')
    print(f'Mag. moments: {mm}')
    print(f'Total moment: {calc.get_magnetic_moment():.4f} μB')

    return e


def u_sweep(u_values=None):
    """Scan U_eff to check gap opening and moment convergence."""
    if u_values is None:
        u_values = [0, 2, 4, 6, 8]

    print('\nU_eff sweep:')
    print(f'{"U_eff":>8}  {"E (eV)":>14}  {"moment (μB)":>12}')
    for u in u_values:
        atoms = build_nio()
        calc = GPAW(
            mode=dict(name='pw', ecut=500),
            xc='PBE',
            kpts=dict(size=(4, 4, 4), gamma=True),
            occupations=FermiDirac(0.05),
            spinpol=True,
            setups=dict(Ni=f':d,{u}') if u > 0 else 'paw',
            nbands=-5,
            mixer=dict(backend='msr1', beta=0.02, nmaxold=8, weight=100),
            convergence=dict(energy=5e-4, density=1e-5),
            txt=f'nio_u{u:.0f}.txt',
        )
        atoms.calc = calc
        e = atoms.get_potential_energy()
        mm = calc.get_magnetic_moment()
        print(f'{u:>8.1f}  {e:>14.6f}  {mm:>12.4f}')


if __name__ == '__main__':
    run_dftu()
    # Uncomment to scan U values:
    # u_sweep([0, 2, 4, 6, 8])
