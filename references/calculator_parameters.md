# GPAW Calculator Parameters Reference

Complete reference for all GPAW calculator keyword arguments.

## Import Paths

```python
# New GPAW (preferred)
from gpaw.new import GPAW
from gpaw import PW, FermiDirac, MethfesselPaxton

# Old GPAW (legacy)
from gpaw import GPAW, PW, FermiDirac
```

## Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | str / dict / PW | `'fd'` | Calculation mode: `'fd'`, `'lcao'`, `PW(ecut)`, or dict |
| `xc` | str / dict | `'LDA'` | Exchange-correlation functional |
| `kpts` | sequence / dict | Γ only | k-point sampling |
| `occupations` | dict / object | Fermi-Dirac 0.1 eV | Smearing method and width |
| `convergence` | dict | See below | SCF stopping criteria |
| `mixer` | dict / object | Pulay | Density mixing settings |
| `nbands` | int | Auto | Number of bands (negative = extra on top of occupied) |
| `spinpol` | bool | Auto | Force spin-polarized calculation |
| `charge` | float | 0 | Total system charge (electrons) |
| `setups` | str / dict | `'paw'` | PAW dataset selection |
| `basis` | str / dict | `{}` | LCAO basis set specification |
| `h` | float | 0.2 | FD grid spacing (Å) |
| `gpts` | sequence | Auto | Explicit grid points (n1, n2, n3); must be divisible by 4 |
| `eigensolver` | str / dict | `'ppcg'` | Diagonalization method |
| `poissonsolver` | dict / object | FFT / FastPoisson | Poisson solver |
| `symmetry` | dict / str | `{}` | Symmetry settings |
| `parallel` | dict | `{}` | MPI parallelization |
| `maxiter` | int | 333 | Max SCF iterations |
| `random` | bool | False | Random wavefunction initialization |
| `hund` | bool | False | Hund's rule initial occupations (atoms) |
| `external` | object | None | External potential |
| `txt` | str / None | `'-'` | Log file path (`'-'` = stdout, `None` = silent) |

## Mode Specifications

```python
# Plane-wave
mode = PW(500)                                  # 500 eV cutoff
mode = PW(500, dedecut='estimate')              # Pulay stress correction
mode = dict(name="pw", ecut=500)              # dict form
mode = dict(name="pw", ecut=500, dtype="single")  # GPU single-precision

# Finite-difference
mode = 'fd'
mode = dict(name='fd', nn=3)                  # 3-point stencil

# LCAO
mode = 'lcao'
```

## k-point Specifications

```python
kpts = (4, 4, 4)                                # Monkhorst-Pack
kpts = dict(size=(4, 4, 4))                      # explicit size
kpts = dict(size=(4, 4, 4), gamma=True)       # gamma-centered
kpts = dict(density=3.5)                         # Å⁻¹, auto-size
kpts = dict(path='GXWKL', npoints=100)       # band structure path
```

## Occupations (Smearing)

```python
# Object form
from gpaw import FermiDirac, MethfesselPaxton
occupations = FermiDirac(0.1)                   # width in eV
occupations = MethfesselPaxton(0.1, order=1)    # order 1 or 2

# Dict form
occupations = dict(name='fermi-dirac', width=0.1)
occupations = dict(name='marzari-vanderbilt', width=0.1)
occupations = dict(name='methfessel-paxton', width=0.1, order=1)
occupations = dict(name='fixed-uniform')         # for open-shell atoms with hund=True
```

## Convergence Criteria

```python
convergence = dict(
    energy=5e-4,       # eV per electron (default: 5e-4)
    density=1e-4,      # electrons per electron (default: 1e-4)
    eigenstates=4e-8,  # eV² per electron (default: 4e-8)
    bands='occupied',  # 'occupied', 'all', or integer count
)
```

Tighter values for demanding calculations:
- Forces / geometry: `energy: 1e-5`
- Stress tensor / cell relaxation: `eigenstates: 1e-10`
- Vibrational modes: `density: 1e-6`
- GW/BSE input: `eigenstates: 1e-8`

## Mixer Settings

```python
mixer = dict(
    backend='msr1',    # 'msr1' (recommended) or 'pulay'
    beta=0.05,         # mixing parameter (0.02–0.25)
    nmaxold=8,         # history length (3–16)
    weight=100,        # metric weight (1 molecules, 50–200 metals)
    method='separate', # 'separate' (spin-polarized), 'sum', 'difference'
)
```

## Symmetry

```python
symmetry = dict()                        # default: point group + time reversal on
symmetry = 'off'                         # disable all symmetry
symmetry = dict(point_group=False)       # time reversal only
symmetry = dict(time_reversal=False)     # point group only
```

Always set `symmetry='off'` for: band structure paths, NEB images, vibrational analysis.

## Eigensolver

```python
eigensolver = 'ppcg'                             # default (PW + FD)
eigensolver = 'rmm-diis'
eigensolver = {'name': 'rmm-diis', 'niter': 3, 'diis-steps': 3}
eigensolver = dict(name='ppcg', niter=5)      # more iterations per step
```

Increase `niter` if eigenstates column in log never converges despite energy/density convergence.

## Poisson Solver

```python
poissonsolver = dict(name='fft')                  # periodic (default)
poissonsolver = dict(dipolelayer='xy')            # slab dipole correction
poissonsolver = dict(dipolelayer='z')             # 1D wire
poissonsolver = PoissonSolver(nn=3, relax='GS')  # FD solver
```

## Parallelization

```python
parallel = dict(kpt=4)           # 4 k-point groups
parallel = dict(band=2)          # 2 band groups
parallel = dict(domain=(2,2,2))  # 2×2×2 domain decomposition
parallel = dict(sl_auto=True)    # automatic ScaLAPACK
```

## PAW Setups

```python
setups = 'paw'                   # default PAW
setups = dict(Cr='14')            # 14-electron Cr (high accuracy)
setups = dict(Ni=':d,6.0')        # DFT+U on Ni d-orbitals, U_eff=6 eV
setups = dict(Mn=':d,4.0', O=':p,2.0')   # multiple corrections
setups = {5: 'ghost'}            # ghost atom at index 5 (int key)
setups = 'sg15'                  # SG15 ONCV pseudopotentials
```

DFT+U syntax: `':orbital,U_eff'` or `':orbital,U_eff,normalize'`
- `normalize=0` → VASP convention (important for p-orbital corrections)
- Semicolons separate multiple orbital corrections: `':d,4.0,0;p,2.0,0'`
