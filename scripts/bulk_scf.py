"""
Ground state SCF for bulk periodic systems (plane-wave mode).

Usage:
    gpaw python bulk_scf.py

Outputs:
    bulk.txt     — GPAW log file
    bulk.gpw     — restart file with converged density/wavefunctions
"""

from ase.build import bulk
from gpaw.dft import GPAW


def run_bulk_scf(
    element: str = 'Fe',
    structure: str = 'bcc',
    ecut: float = 500,
    kgrid: tuple = (8, 8, 8),
    smearing: float = 0.1,
    magnetic: bool = True,
    magmom: float = 2.2,
    xc: str = 'PBE',
    outfile: str = 'bulk',
):
    atoms = bulk(element, structure)

    if magnetic:
        atoms.set_initial_magnetic_moments([magmom] * len(atoms))

    calc = GPAW(
        mode=dict(name='pw', ecut=ecut),
        xc=xc,
        kpts=dict(size=kgrid, gamma=True),
        occupations=dict(name='fermi-dirac', width=smearing),
        spinpol=magnetic,
        nbands=-10,                           # 10 extra empty bands above Fermi
        mixer=dict(
            backend='msr1',
            beta=0.05,
            nmaxold=10,
            weight=100,
            method='fullspin' if magnetic else 'separate',
        ),
        convergence=dict(
            energy=5e-4,
            density=1e-5 if magnetic else 1e-4,
            eigenstates=4e-8,
        ),
        txt=f'{outfile}.txt',
    )
    atoms.calc = calc

    e = atoms.get_potential_energy()
    ef = calc.get_fermi_level()
    calc.write(f'{outfile}.gpw')

    print(f'Total energy:  {e:.6f} eV')
    print(f'Fermi level:   {ef:.4f} eV')

    if magnetic:
        mm = calc.get_magnetic_moments()
        print(f'Mag. moments:  {mm}')
        print(f'Total moment:  {calc.get_magnetic_moment():.3f} μB')

    return e


if __name__ == '__main__':
    run_bulk_scf()
