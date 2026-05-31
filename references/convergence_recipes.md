# GPAW Convergence Recipes

Systematic protocols for converging GPAW calculations. Run these before production calculations.

## 1. Plane-Wave Cutoff Convergence

```python
"""Test total energy vs. ecut. Run once per system type, not per structure."""
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

atoms = bulk('Fe', 'bcc')
results = {}

for ecut in [300, 400, 500, 600, 700]:
    calc = GPAW(
        mode=PW(ecut),
        xc='PBE',
        kpts=dict(size=(8, 8, 8), gamma=True),
        occupations=FermiDirac(0.1),
        txt=f'ecut_{ecut}.txt',
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    results[ecut] = e
    print(f'ecut={ecut:4d} eV  E={e:.6f} eV  ΔE={e-results[300]:.4f} eV')
```

**Decision criterion:** Choose lowest ecut where ΔE < 5 meV/atom vs. next higher value.
Typical converged values: 400–500 eV (PBE), 500–600 eV (cohesive energies/stress).

## 2. k-point Convergence

```python
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

atoms = bulk('Fe', 'bcc')
atoms.set_initial_magnetic_moments([2.2])

for k in [4, 6, 8, 10, 12]:
    calc = GPAW(
        mode=PW(500),
        xc='PBE',
        kpts=dict(size=(k, k, k), gamma=True),
        occupations=FermiDirac(0.1),
        spinpol=True,
        txt=f'kpts_{k}.txt',
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    print(f'k={k}  E={e:.6f} eV')
```

**Decision criterion:** ΔE < 5 meV/atom vs. next denser grid.

## 3. SCF Convergence Tightening

Default convergence is sufficient for total energies. Tighten for these workflows:

| Workflow | Recommended settings |
|----------|---------------------|
| Geometry optimization (forces) | `energy: 1e-5, density: 1e-5` |
| Cell relaxation (stress tensor) | `eigenstates: 1e-10` + `PW(ecut, dedecut='estimate')` |
| Vibrational modes / phonons | `density: 1e-6, symmetry='off'` |
| NEB intermediate images | `energy: 5e-3, eigenstates: 1e-7` (looser is fine) |
| Band structure (fixed_density) | `bands: N` (N = number of target bands) |
| GW / BSE input | `eigenstates: 1e-8` + `diagonalize_full_hamiltonian()` |

## 4. Smearing Width Convergence (Metals)

For final total energies, verify result is stable vs. smearing width:

```python
from gpaw import GPAW, PW, MethfesselPaxton

for width in [0.2, 0.1, 0.05, 0.02]:
    calc = GPAW(
        mode=PW(500),
        xc='PBE',
        kpts=dict(size=(12, 12, 12), gamma=True),
        occupations=MethfesselPaxton(width, order=1),
        txt=f'smear_{width}.txt',
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    e0 = calc.get_xc_difference(dict(name='PBE'))  # T→0 extrapolation
    print(f'width={width:.2f}  E={e:.6f}  E(T=0)={e0:.6f}')
```

Use `MethfesselPaxton` for metals: it gives the most accurate T→0 extrapolation.

## 5. Vacuum Convergence (Slabs / Molecules)

```python
from ase.build import fcc111
from gpaw import GPAW, PW, FermiDirac

for vac in [6, 8, 10, 12, 15]:
    slab = fcc111('Al', size=(1, 1, 4), vacuum=vac)
    calc = GPAW(
        mode=PW(500),
        xc='PBE',
        kpts=dict(size=(8, 8, 1), gamma=True),
        occupations=FermiDirac(0.1),
        txt=f'vac_{vac}.txt',
    )
    slab.calc = calc
    e = slab.get_potential_energy()
    print(f'vacuum={vac:.0f} Å  E={e:.6f} eV')
```

Typical converged vacuum: 10–12 Å for slabs, 6 Å for molecules.

## 6. LCAO Basis Convergence

```python
for basis in ['sz', 'szp', 'dz', 'dzp', 'tzdp']:
    calc = GPAW(
        mode='lcao',
        basis=basis,
        xc='PBE',
        kpts=dict(size=(4, 4, 4), gamma=True),
        txt=f'basis_{basis}.txt',
    )
    atoms.calc = calc
    e = atoms.get_potential_energy()
    print(f'basis={basis}  E={e:.6f} eV')
```

DZP is the standard production choice. TZP rarely needed except for response properties.

## 7. Grid Spacing Convergence (FD mode)

```python
for h in [0.25, 0.20, 0.18, 0.15]:
    calc = GPAW(
        mode='fd',
        h=h,
        xc='PBE',
        txt=f'h_{h:.2f}.txt',
    )
    mol.calc = calc
    e = mol.get_potential_energy()
    print(f'h={h:.2f} Å  E={e:.6f} eV')
```

h=0.20 Å is the standard; h=0.18 for tighter; h=0.15 for vibrational analysis of light molecules.

## Common Failure Checklist

Before diagnosing a convergence problem, verify:

- [ ] `magmoms` set for all magnetic atoms (not just non-zero ones)
- [ ] `spinpol=True` explicit for magnetic calculations
- [ ] `nbands` includes at least 5 empty bands above Fermi level (metals: use `nbands=-10`)
- [ ] smearing is non-zero for metals
- [ ] k-mesh not Γ-only for periodic metals
- [ ] `symmetry='off'` for calculations requiring broken symmetry
- [ ] gpts divisible by 4 (FD mode) or 8 (better Poisson convergence)
- [ ] cell not perfectly cubic for isolated molecules (shift by 0.0001 Å)
- [ ] `GPAW_SETUP_PATH` environment variable set correctly
