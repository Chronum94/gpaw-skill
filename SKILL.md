---
name: gpaw
description: >
  Expert assistant for GPAW (Grid-based Projector Augmented Wave) DFT calculations.
  Use for DFT calculations including ground-state SCF, band structure, DOS/PDOS, magnetic systems (spin-polarized,
  DFT+U, SOC), NEB/transition states, TDDFT optical spectra, and vibrational modes.
  Generates complete, ready-to-run Python scripts with converged parameters. Diagnoses
  SCF convergence failures, mixer issues, PAW setup pitfalls, and log file anomalies inline.
allowed-tools: Read Write Edit Bash WebFetch
license: MIT
metadata:
  skill-author: Anubhab Haldar
---

# GPAW Calculation Skill

## Overview

Produce complete, executable GPAW Python scripts with correctly converged parameters.
Always diagnose convergence problems inline: identify the root cause, explain it in one
sentence, and provide a corrected snippet. Default functionals are **PBE** or ***R2SCAN***, used as `MGGA_X_R2SCAN+MGGA_C_R2SCAN`; default
plane-wave cutoff starts at **400–500 eV** and tightens if precision demands it.

---

## Trigger Conditions

Use this skill when the user:
- Asks to set up any GPAW calculation (SCF, band structure, DOS, magnetic, NEB, TDDFT, phonons)
- Pastes a `.txt` log file and asks why the calculation did not converge
- Asks about mixer, smearing, k-points, PAW setups, eigensolver, or stress tensor settings
- Asks for a convergence test workflow
- Asks to interpret GPAW output or extract a result from a `.gpw` file

---

## Core Script Template (All Calculations)

Every script follows this skeleton:

```python
from ase import Atoms
from ase.build import bulk, molecule, surface   # as needed
from ase.io import read
# Use from gpaw import ... for old GPAW
# Use gpaw.dft import ... for new GPAW.
# Prefer new GPAW.
from gpaw.new import GPAW

# 1. Build structure (ASE uses Å throughout)
atoms = ...

# 2. Calculator
calc = GPAW(
    mode=dict(name="pw", ecut= 500), # or 'fd' / 'lcao'
    xc='PBE', # MGGA_X_R2SCAN+MGGA_C_R2SCAN for R2SCAN
    kpts=dict(size=(k, k, k), gamma=True),
    occupations=dict(name='fermi-dirac', width=0.05), # see Smearing section
    convergence=dict(
        energy=5e-4,                    # eV/electron
        density=1e-4,                   # electrons/electron
        eigenstates=4e-8,               # eV²/electron
    ),
    mixer=dict(
        backend='msr1',
        beta=0.05,
        nmaxold=8,
        weight=100,
    ),                                  # tune for metals; see Mixer section
    txt='calc.txt',
)
atoms.calc = calc

# 3. Run
atoms.get_potential_energy()
calc.write('calc.gpw') # mode='all' to save wave functions.
```

---

## Calculation Workflows

### 1. Ground State SCF — Bulk (PW mode)

```python
from ase.build import bulk
from gpaw.dft import GPAW

atoms = bulk('Fe', 'bcc') # Not giving lattice constants uses ASE default lattice constants
atoms.set_initial_magnetic_moments([2.2] * len(atoms))  # omit if non-magnetic

calc = GPAW(
    mode=dict(name="pw", ecut=500, dtype="double"), # Double for CPU, single-precision for GPU.
    xc="PBE",
    kpts=dict(size=(8, 8, 8), gamma=True),
    occupations=dict(name='fermi-dirac', width=0.1),            # 0.1 eV for metals; 0 for insulators/molecules
    spinpol=True,                            # omit if non-magnetic
    convergence=dict(energy=5e-4, density=1e-4, eigenstates=4e-8),
    txt='bulk.txt',
)
atoms.calc = calc
e = atoms.get_potential_energy()
calc.write('bulk.gpw')
print(f'Energy: {e:.4f} eV')
```

**Key flags:**
- `spinpol=True` required for magnetic systems; GPAW may auto-detect but being explicit avoids surprises.
- Always set `magmoms` before SCF for open-shell metals — without them the calculation may collapse to the wrong magnetic state.
- For transition metals with incomplete d-shells, include ≥5 extra empty bands: `nbands=-5` (relative to occupied) or an explicit count, or prefer PPCG eigensolver.

---

### 2. Ground State SCF — Slab / Surface

```python
from ase.build import fcc111
from gpaw.dft import GPAW

slab = fcc111('Al', size=(2, 2, 4), vacuum=10.0)

calc = GPAW(
    mode=dict(mode='pw', ecut=500),
    xc='PBE',
    kpts=dict(size=(4, 4, 1), gamma=True),  # no k-points along vacuum axis
    occupations=dict(name='fermi-dirac', width=0.1),
    poissonsolver=dict(dipolelayer='xy'),        # dipole correction for asymmetric slabs
    convergence=dict(energy=5e-4, density=1e-4),
    txt='slab.txt',
)
slab.calc = calc
slab.get_potential_energy()
calc.write('slab.gpw')
```

**Dipole correction:** Required when the two surfaces of a slab are chemically different (adsorbates on one side, or asymmetric terminations). Omit for symmetric slabs.

---

### 3. Ground State SCF — Molecule (FD mode)

```python
from ase.build import molecule
from gpaw import GPAW

mol = molecule('H2O')
mol.center(vacuum=6.0)

calc = GPAW(
    mode='fd',
    h=0.2,                                   # grid spacing Å; 0.15–0.18 for tighter
    xc='PBE',
    occupations=dict(name='fermi-dirac', width=0.0),   # integer occupations for molecule
    convergence=dict(energy=5e-4, density=1e-6),
    txt='mol.txt',
)
mol.calc = calc
mol.get_potential_energy()
calc.write('mol.gpw')
```

**FD grid tips:**
- For open-shell atoms, set `hund=True` to apply Hund's rule initial occupations.

---

### 4. Ground State SCF — Large System (LCAO mode)

```python
from gpaw.dft import GPAW

calc = GPAW(
    mode='lcao',
    basis='dzp',                             # double-zeta polarized; 'szp(dzp)' for cheaper
    xc='PBE',
    kpts=dict(size=(4, 4, 4), gamma=True),
    occupations=dict(name='fermi-dirac', width=0.1),
    convergence=dict(energy=5e-4, density=1e-4),
    parallel=dict(sl_auto=True, use_elpa=True), # use_elpa only if GPAW installed with Elpa.
    txt='lcao.txt',
)
atoms.calc = calc
atoms.get_potential_energy()
calc.write('lcao.gpw')
```

**LCAO basis sets:** `sz` < `dz` < `szp` < `dzp` < `tzdp`. DZP is the standard balance point. Always set `symmetry='off'` for vibrational/Raman post-processing.

---

### 5. Band Structure

Two-step workflow: (1) converged SCF on k-mesh, (2) non-self-consistent along high-symmetry path.

```python
# Step 1 — ground state
from ase.build import bulk
from gpaw.dft import GPAW

atoms = bulk('Si', 'diamond', a=5.43)
calc = GPAW(
    mode=dict(mode='pw', ecut=500),
    xc='PBE',
    kpts=dict(size=(8, 8, 8), gamma=True),
    occupations=dict(name='fermi-dirac', width=0.01),
    convergence=dict(energy=5e-4, density=1e-4, eigenstates=4e-8),
    txt='Si_gs.txt',
)
atoms.calc = calc
atoms.get_potential_energy()
ef = calc.get_fermi_level()
calc.write('Si_gs.gpw')
print(f'Fermi level: {ef:.4f} eV')

# Step 2 — band structure at fixed density
bs_calc = GPAW('Si_gs.gpw').fixed_density(
    nbands=16,
    symmetry='off',                          # must be off for complete k-path
    kpts=dict(path=bs_calc.atoms.cell.bandpath().path, npoints=100),
    convergence=dict(bands=8),                # converge only occupied + a few empty
    txt='Si_bs.txt',
)
bs = bs_calc.band_structure()
bs.plot(filename='bandstructure.png', show=False, emax=10.0)
bs_calc.write('Si_bs.gpw')
```

**Post-processing:**

```python
from gpaw import GPAW
import matplotlib.pyplot as plt

calc = GPAW('Si_bs.gpw', txt=None)
bs = calc.band_structure()
ax = bs.plot(emax=6.0, emin=-6.0)
plt.savefig('band_structure.png', dpi=150)
```

---

### 6. Density of States (DOS / PDOS)

```python
from gpaw import GPAW
from gpaw.dos import DOSCalculator
import numpy as np

calc = GPAW('bulk.gpw', txt=None)

# Total DOS
dos = DOSCalculator.from_calculator(calc)
energies = np.linspace(-10, 5, 500)
total = dos.raw_dos(energies, spin=0, width=0.1)   # width=0 → tetrahedron

# Projected DOS (atom index, angular momentum)
pdos_d = dos.raw_pdos(energies, a=0, l=2, spin=0, width=0.1)   # d-orbital on atom 0

import matplotlib.pyplot as plt
ef = calc.get_fermi_level()
plt.figure()
plt.plot(energies - ef, total, label='Total DOS')
plt.plot(energies - ef, pdos_d, label='d PDOS atom 0')
plt.axvline(0, color='k', ls='--', lw=0.8)
plt.xlabel('Energy − E_F (eV)')
plt.ylabel('DOS (states/eV)')
plt.legend()
plt.savefig('dos.png', dpi=150)
```

**Broadening:** `width=0.0` enables linear tetrahedron interpolation (best for metals); Gaussian broadening (`width=0.1`) is fine for insulators/molecules. For PDOS use `l=0` (s), `l=1` (p), `l=2` (d), `l=3` (f).

---

### 7. Structure Relaxation + Cell Optimization

```python
from ase.build import bulk
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from gpaw import GPAW, PW

atoms = bulk('Si', 'fcc', a=5.5)   # slightly off equilibrium
calc = GPAW(
    mode=PW(500, dedecut='estimate'),     # dedecut corrects Pulay stress
    xc='PBE',
    kpts=dict(size=(6, 6, 6), gamma=True),
    occupations=dict(name='fermi-dirac', width=0.1),
    convergence=dict(eigenstates=1e-10),   # tight eigenstates for stress
    txt='relax.txt',
)
atoms.calc = calc

ucf = FrechetCellFilter(atoms)
opt = BFGS(ucf, logfile='relax.log')
opt.run(fmax=0.01)                         # 0.01 eV/Å for cell; use 0.05 for atoms-only

print('Relaxed a:', atoms.cell.lengths())
calc.write('relaxed.gpw')
```
---

### 8. Equation of State / Lattice Constant Convergence

```python
from ase.build import bulk
from ase.eos import EquationOfState
from gpaw import GPAW, PW
import numpy as np

a0 = 4.05  # starting guess (Å)
volumes, energies = [], []

for scale in np.linspace(0.94, 1.06, 7):
    atoms = bulk('Al', 'fcc', a=a0 * scale)
    calc = GPAW(mode=dict(name='pw', ecut=500), xc='PBE',
                kpts=dict(size=(8, 8, 8), gamma=True),
                convergence=dict(energy=5e-4),
                txt=f'eos_{scale:.3f}.txt')
    atoms.calc = calc
    energies.append(atoms.get_potential_energy())
    volumes.append(atoms.get_volume())

eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()
a_eq = (v0 * 4) ** (1/3)           # FCC: 4 atoms per conventional cell
print(f'a₀ = {a_eq:.4f} Å,  B = {B/units.GPa:.1f} GPa')
eos.plot('eos.png')
```

---

### 9. Magnetic Calculations

#### 9a. Spin-Polarized SCF (ferromagnetic)

```python
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

atoms = bulk('Fe', 'bcc', a=2.87)
atoms.set_initial_magnetic_moments([2.2] * len(atoms))

calc = GPAW(
    mode=dict(name='pw', ecut=500),
    xc='PBE',
    kpts=dict(size=(12, 12, 12), gamma=True),   # metals need dense k-mesh
    occupations=FermiDirac(0.1),
    spinpol=True,
    nbands=-10,                                    # 10 extra empty bands per spin
    mixer=dict(backend='pulay', beta=0.05,
               nmaxold=8, weight=100),
    convergence=dict(energy=5e-4, density=1e-5, eigenstates=4e-8),
    txt='fe_fm.txt',
)
atoms.calc = calc
atoms.get_potential_energy()
mm = atoms.get_magnetic_moments()
print('Magnetic moments:', mm)
calc.write('fe_fm.gpw')
```

#### 9b. Antiferromagnetic

```python
atoms = bulk('Fe', 'bcc', a=2.87, cubic=True)   # supercell to allow AFM
magmoms = [2.2, -2.2] * (len(atoms) // 2)
atoms.set_initial_magnetic_moments(magmoms)
# rest of calculator identical to FM case
```

#### 9c. DFT+U

```python
from gpaw import GPAW, PW, FermiDirac
from ase import Atoms

# NiO antiferromagnet with Hubbard U on Ni d-orbitals
a, b = 4.19, 4.19 / 2**0.5
m = 2.0
nio = Atoms('Ni2O2', pbc=True, cell=(b, b, a),
            positions=[(0,0,0),(b/2,b/2,a/2),(0,0,a/2),(b/2,b/2,0)],
            magmoms=[m, -m, 0, 0])

calc = GPAW(
    mode=PW(600),
    xc='PBE',
    kpts=dict(size=(4, 4, 4), gamma=True),
    occupations=FermiDirac(0.05),
    spinpol=True,
    setups=dict(Ni=':d,6.0'),                  # U_eff = U − J = 6.0 eV on d-orbitals
    # Multi-orbital: setups=dict(Ni=':d,4.0,0;p,2.0,0')
    convergence=dict(energy=5e-4, density=1e-5),
    txt='nio_u.txt',
)
nio.calc = calc
nio.get_potential_energy()
calc.write('nio_u.gpw')
```

**DFT+U syntax:** `':orbital,U_eff'` or `':orbital,U_eff,normalize'` where `normalize=0` matches VASP convention (minimal effect for d/f, significant for p). Separate multiple corrections with semicolons.

#### 9d. Spin-Orbit Coupling (post-SCF)

```python
from gpaw import GPAW
from gpaw.spinorbit import soc_eigenstates
import numpy as np

# Requires a converged spinpol or non-spinpol GPW file
calc = GPAW('band_structure.gpw', txt=None)
soc = soc_eigenstates(
    calc,
    n1=0,          # first band index
    n2=20,         # last band index
    theta=0.0,     # spin quantization axis polar angle (degrees)
    phi=0.0,       # azimuthal angle
)
e_km = soc.eigenvalues()   # shape: (k-points, bands)
# Spin character: ⟨σ_z⟩ per state
sz = soc.spin_projection(v=2)  # v=0(x), 1(y), 2(z)
```

**When to use SOC:** Heavy elements (Pt, W, Bi, WS₂), topological insulators, magnetic anisotropy, Rashba/Dresselhaus splitting. SOC is applied non-self-consistently by default (scalar-relativistic SCF + SOC diagonalization) — accurate for most band-structure purposes.

---

### 10. NEB (Nudged Elastic Band)

```python
# neb.py — run with: mpiexec -np 12 gpaw python neb.py
from ase.io import read
from ase.neb import NEB
from ase.optimize import BFGS
from gpaw import GPAW, PW

# Load initial and final images (pre-relaxed)
initial = read('initial.traj')
final   = read('final.traj')

n_images = 3
images = [initial.copy() for _ in range(n_images + 2)]
images[-1] = final.copy()

# Linear interpolation of intermediate images
neb = NEB(images, climb=True, parallel=True)
neb.interpolate()

# Attach a GPAW calculator to each image (each gets its own MPI sub-communicator)
for image in images[1:-1]:
    image.calc = GPAW(
        mode=dict(name='pw', ecut=500),
        xc='PBE',
        kpts=dict(size=(4, 4, 1), gamma=True),
        convergence=dict(energy=5e-3, eigenstates=1e-7),  # slightly looser for NEB
        txt='-',
    )

opt = BFGS(neb, logfile='neb.log')
opt.run(fmax=0.05)    # eV/Å
```

**Key points:**
- Run with `mpiexec -np N*n_images` where N CPUs per image. GPAW distributes automatically via `parallel=True`.
- Initial and final structures must be independently relaxed to `fmax ≤ 0.05` before NEB.
- Climbing image (`climb=True`) locates the true saddle point; activate only after preliminary convergence.
- Use `fmax=0.05` for activation energies; `0.01` for precise barrier heights.

---

### 11. TDDFT Optical Spectra

#### 11a. Linear-Response TDDFT (finite systems)

```python
# Step 1: ground state with unoccupied bands
from ase.build import molecule
from gpaw import GPAW
from gpaw.lrtddft import LrTDDFT, photoabsorption_spectrum

mol = molecule('Na2')
mol.center(vacuum=5.0)

calc = GPAW(
    mode='fd',
    h=0.25,
    xc='PBE',
    nbands=20,                              # need unoccupied states
    convergence=dict(density=1e-6),
    txt='na2_gs.txt',
)
mol.calc = calc
mol.get_potential_energy()

# Step 2: LR-TDDFT
lr = LrTDDFT(calc, xc='LDA',
             istart=0,                     # first electron-hole pair
             jend=10)                      # last unoccupied index
lr.write('lr_na2.gz')

# Step 3: spectrum
lr2 = LrTDDFT.read('lr_na2.gz')
photoabsorption_spectrum(lr2, 'na2_spectrum.dat',
                         e_min=0.0, e_max=10.0,
                         width=0.1)       # Lorentzian broadening eV
```

#### 11b. Time-Propagation TDDFT (real-time, extended or finite)

```python
# gs_calc.py — ground state
from gpaw import GPAW, PW, FermiDirac
from ase.build import bulk

atoms = bulk('Na', 'bcc', a=4.23)
calc = GPAW(mode=PW(400),
            xc='PBE',
            kpts=dict(size=(4, 4, 4), gamma=True),
            occupations=FermiDirac(0.1),
            txt='na_gs.txt')
atoms.calc = calc
atoms.get_potential_energy()
calc.write('na_gs.gpw')

# tddft_run.py — time propagation
from gpaw.tddft import TDDFT
from gpaw.tddft.laser import GaussianPulse

td_calc = TDDFT('na_gs.gpw',
                td_potential=GaussianPulse(strength=1e-5,
                                           time0=0.0,
                                           frequency=1.0,   # eV
                                           sigma=0.3,       # eV (bandwidth)
                                           polarization=[1, 0, 0]),
                txt='tddft.txt')

td_calc.absorption_kick(strength=1e-5, direction=0)   # or use pulse above
td_calc.propagate(time_step=10.0,        # attoseconds
                  iterations=1500,
                  dipole_moment_file='dm.dat',
                  restart_file='td.gpw',
                  dump_interval=100)
```

**Post-processing (Fourier transform → spectrum):**

```python
from gpaw.tddft.spectrum import photoabsorption_spectrum
photoabsorption_spectrum('dm.dat', 'spectrum.dat',
                         folding='Gauss', width=0.2,
                         e_min=0.0, e_max=20.0,
                         delta_e=0.01)
```

---

### 12. Vibrational Modes

```python
from math import cos, sin, pi
from ase import Atoms
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from gpaw import GPAW

# H2O molecule
d, t = 0.9575, pi/180 * 104.51
h2o = Atoms('H2O',
            positions=[(0,0,0),(d,0,0),(d*cos(t),d*sin(t),0)])
h2o.center(vacuum=4.0)

calc = GPAW(
    mode='lcao',                            # LCAO faster for vibrational displacements
    basis='dzp',
    xc='PBE',
    convergence=dict(density=1e-6),
    symmetry='off',                         # mandatory for vibrational analysis
    txt='h2o.txt',
)
h2o.calc = calc

# First relax
QuasiNewton(h2o, logfile='opt.log').run(fmax=0.01)  # tight relax before vibrations

# Vibrational analysis
vib = Vibrations(h2o, delta=0.02)    # displacement Å; default 0.01
vib.run()
vib.summary(method='frederiksen')   # Frederiksen correction for imaginary modes
vib.write_jmol()                     # for visualization

# Access frequencies
freqs = vib.get_frequencies()
print('Frequencies (cm⁻¹):', freqs.real)
```

**Raman (extended systems via electron-phonon):**

```python
# Requires LCAO supercell + elphraman post-processing (separate workflow)
# See: gpaw.readthedocs.io → Vibrational → Raman spectroscopy for extended systems
```

---

## Convergence Reference

### SCF Convergence Criteria

| Criterion | Default | Tight | Use case |
|-----------|---------|-------|----------|
| `energy` (eV/e⁻) | 5×10⁻⁴ | 5×10⁻⁵ | Force/stress calculations |
| `density` (e⁻/e⁻) | 1×10⁻⁴ | 1×10⁻⁶ | Vibrational, properties |
| `eigenstates` (eV²/e⁻) | 4×10⁻⁸ | 1×10⁻¹⁰ | Stress tensor, hybrids |
| `bands` | 'occupied' | 'all' or int | Unoccupied-state properties |

### K-point Density Guidelines

| System | Starting grid | Converged |
|--------|--------------|-----------|
| Bulk metal (Fe, Al) | (6,6,6) | (12,12,12) or `density≥4.0` |
| Bulk insulator/semiconductor | (4,4,4) | (8,8,8) |
| Surface slab (4×4) | (4,4,1) | (8,8,1) |
| 1D wire | (1,1,8) | (1,1,20) |
| Molecule / cluster | Γ only | Γ only |

Use `kpts=dict(density=3.5)` (Å⁻¹) for automatic density-based sampling.

### Plane-Wave Cutoff Guidelines

| System / XC | Starting ecut | Converged |
|-------------|--------------|-----------|
| PBE bulk metals | 400 eV | 500–600 eV |
| PBE insulators | 400 eV | 500 eV |
| Cohesive energy (high accuracy) | 500 eV | 600–800 eV |
| GW / hybrid inputs | 300 eV | 500 eV |
| Stress tensor | 500 eV + `dedecut='estimate'` | 600 eV |

Convergence test: compute total energy at ecut = 300, 400, 500, 600 eV; plot; choose where ΔE < 5 meV/atom.

### Smearing (Occupations)

| System type | Setting |
|-------------|---------|
| Metal | `dict(name='fermi-dirac', width=0.1)` or `dict(name='methfessel-paxton', width=0.1)` |
| Semiconductor / insulator | `dict(name='fermi-dirac', width=0.01)` |
| Molecule / cluster | `dict(name='fermi-dirac', width=0.0)` |
| Open-shell atom | `hund=True, occupations=dict(name='fixed-uniform')` |

Always verify total energy is converged w.r.t. smearing width. For final energies on metals, extrapolate to T→0 using `MethfesselPaxton`.

---

## Mixer Reference

The mixer controls density update stability. Defaults work for many cases; tune when SCF oscillates.

```python
# Standard (non-magnetic bulk)
mixer=dict(backend='msr1', beta=0.05, nmaxold=8, weight=100)

# Magnetic metals (separate up/down mixing)
mixer=dict(backend='pulay', beta=0.02, nmaxold=8,
           method='separate', weight=100)

# Strongly correlated / hard to converge
mixer=dict(backend='msr1', beta=0.02, nmaxold=10, weight=100)

# Molecules / non-periodic (looser mixing fine)
mixer=dict(backend='pulay', beta=0.25, nmaxold=4, weight=1)
```

**Tuning rules:**
- Lower `beta` (0.02–0.05) for magnetic metals or strongly correlated systems.
- Higher `weight` (50–200) improves stability for metals with sharp Fermi surfaces.
- `nmaxold=3–5` for molecules; `8–16` for periodic systems.
- Switch to `'msr1'` backend if Pulay oscillates even with small `beta`.

---

## PAW Setups Reference

Default setups are correct for most elements. Flag these special cases:

| Element / situation | Recommendation |
|--------------------|----------------|
| Chromium (Cr) | `setups=dict(Cr='14')` — 14-electron setup for high accuracy (d+s+p) |
| Lanthanides | Require version ≥ 24.11.0 setups; check `gpaw install-data` |
| DFT+U elements (Ni, Co, Fe, Mn, U) | Add `setups=dict(El=':d,U_eff')` — do NOT use a different PAW dataset |
| Ghost atoms | `setups={idx: 'ghost'}` where idx is atom index |
| All-electron reconstruction | `setups=dict(H='ae')` — use only for post-processing, not SCF |
| SG15 / HGH pseudopotentials | `setups='sg15'` or `setups='hgh'` — non-PAW; rarely needed |

Environment variable: `GPAW_SETUP_PATH` must point to the unpacked dataset directory (colon-separated for multiple). Verify with `gpaw info`.

---

## Log File Diagnostics

The GPAW `.txt` log has these key sections:

```
Input parameters:           <- echo of all calculator settings
...
Positions:                  <- atomic positions (Å)
...
              iter  time  total energy  log10-change:  |dn|  eigst  bands
iter:   1  09:01  -100.123456  +99      3.2e-02  2.3e-01  ...
iter:   2  09:02  -101.234567  -0.48    1.1e-02  8.4e-02  ...
...
iter:  23  09:15  -101.456789  -4.52    8.3e-05  3.1e-08  ...  <- converged
```

**Column meanings:**

| Column | Meaning | Healthy value |
|--------|---------|---------------|
| `total energy` | DFT total energy (eV) | Decreasing monotonically |
| `log10-change` | log₁₀(ΔE) | Should decrease each iteration |
| `\|dn\|` | Density change | Decreasing; < 1×10⁻⁴ at convergence |
| `eigst` | Eigenstate residuals | < 4×10⁻⁸ at convergence |
| `bands` | Bands not converged | 0 at convergence |

**Diagnose from log:**
- Energy oscillates back and forth → reduce mixer `beta`; switch to `'msr1'`
- `|dn|` stuck at a plateau → wrong initial magnetic moments; increase `nbands`
- Poisson solver error → check cell symmetry; set gpts divisible by 8
- `bands` column never reaches 0 → increase `nbands` or `eigensolver=dict(niter=5)`
- Calculation reached `maxiter=333` → system is intrinsically hard; first try mixer tuning, then increase `maxiter` only as last resort

---

## Common Failure Patterns and Fixes

### Metal without smearing
**Symptom:** Energy oscillates wildly; `|dn|` never decreases.
**Fix:** Add `occupations=FermiDirac(0.1)` and `nbands=-10`.

### Magnetic calculation converges to wrong state
**Symptom:** Final magnetic moment near zero despite expecting ferromagnetism.
**Fix:** Set explicit `magmoms` before calling `get_potential_energy()`. Use `hund=True` for isolated atoms.

### Poisson solver did not converge
**Symptom:** `"Poisson solver did not converge"` in log, usually on molecules.
**Fix:** (1) Break cubic symmetry with `cell[1,1] += 0.0001`; (2) ensure gpts divisible by 8; (3) set low Fermi temperature for isolated systems.

### Stress tensor is wrong / Pulay stress
**Symptom:** Relaxed volume is clearly wrong; energy changes significantly with cell perturbation.
**Fix:** Use `PW(ecut, dedecut='estimate')` and tighten `convergence=dict(eigenstates=1e-10)`.

### G0W0 / GW band gap not converged
**Symptom:** Large changes in quasi-particle gap with ecut or nbands.
**Fix:** Use `ecut_extrapolation=True`; test ecut at 200, 300, 400 eV; use at least 3× more nbands than occupied bands.

### DFT+U insulator stays metallic
**Symptom:** No gap opens despite large U value.
**Fix:** Check initial magmoms are set correctly for AFM ordering; ensure correct U orbital (`:d` vs `:f`); try a different starting density with `random=True`.

### LCAO vibrations have imaginary modes at equilibrium
**Symptom:** Negative frequencies for non-degenerate modes.
**Fix:** Tighten geometry optimization to `fmax=0.005` (not 0.05); tighten `convergence=dict(density=1e-6)`; set `symmetry='off'`.

---

## Quick-Reference Checklists

### Before running any periodic calculation
- [ ] Cell vectors correct (ASE uses Å)
- [ ] `pbc=True` set on Atoms object
- [ ] k-points non-trivial (not Γ-only for metals)
- [ ] Smearing applied for metals
- [ ] `magmoms` set if magnetic system expected
- [ ] Sufficient `nbands` for system (metals: +5 to +10 empty)
- [ ] `txt` log file specified

### Before running a band structure
- [ ] Ground state converged and saved with `calc.write('gs.gpw')`
- [ ] `symmetry='off'` in band calculation
- [ ] `fixed_density()` used (not a fresh SCF)
- [ ] `convergence=dict(bands=N)` set for desired number of bands

### Before running vibrational analysis
- [ ] Structure relaxed to `fmax ≤ 0.01` (not 0.05)
- [ ] `convergence=dict(density=1e-6)` in calculator
- [ ] `symmetry='off'` in calculator
- [ ] LCAO or FD mode (PW can be used but LCAO is faster for displacements)

### Before running NEB
- [ ] Initial and final endpoints independently relaxed
- [ ] Same calculator settings on all images
- [ ] `climb=True` enabled after first convergence pass
- [ ] MPI process count = N_CPUs × N_internal_images

---

## Useful Snippets

### Restart from .gpw file

```python
from gpaw import restart
atoms, calc = restart('calc.gpw', txt='restart.txt')
atoms.get_potential_energy()   # continues SCF
```

### Extract data from converged calculation

```python
from gpaw import GPAW
calc = GPAW('calc.gpw', txt=None)
atoms = calc.get_atoms()
e = calc.get_potential_energy()
ef = calc.get_fermi_level()
forces = calc.get_forces()
density = calc.get_pseudo_density()           # 3D array
charges = calc.get_hirshfeld_charges()        # Hirshfeld charge analysis
```

### Full Hamiltonian diagonalization (for GW input)

```python
calc.diagonalize_full_hamiltonian()
calc.write('gs_allbands.gpw', 'all')          # saves wavefunctions
```

### Parallelization hints

```python
# Band parallelization for many-band calculations
calc = GPAW(..., parallel=dict(band=4))

# k-point parallelization (default for most jobs)
calc = GPAW(..., parallel=dict(kpt=8))

# Domain decomposition
calc = GPAW(..., parallel=dict(domain=(2, 2, 2)))
```
