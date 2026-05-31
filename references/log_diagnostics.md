# GPAW Log File Diagnostics

## Log File Structure

A GPAW `.txt` log has these major sections:

```
  ___ ___ ___ _ _ _
 |   |   |_  | | | |          GPAW xx.x.x
 | | | | | . | | | |          ...
 |_ _|___|___|_____|          ...

Input parameters:
  ...                          ← echo of all settings

System: ...                    ← unit cell, atoms
...

k-point calculation: ...       ← k-point grid info

              iter  time  ...  ← SCF iteration table header
iter:   1  ...
iter:   2  ...
...
Converged after N iterations.

Energy contributions relative to reference atoms: (SHE)
  ...
Free energy:   -XX.XXXX eV
Extrapolated:  -XX.XXXX eV     ← T→0 extrapolation
```

## SCF Iteration Table Columns

```
iter:  N  HH:MM  E_total  log10(ΔE)  |dn|    eigst   bands
```

| Column | Meaning | Healthy |
|--------|---------|---------|
| `iter` | Iteration number | Increasing |
| `time` | Wall clock time | — |
| `E_total` | DFT total energy (eV) | Decreasing then stable |
| `log10(ΔE)` | log₁₀ of energy change | Decreasing monotonically |
| `\|dn\|` | Density change norm | Decreasing; < 1×10⁻⁴ |
| `eigst` | Eigenstate residuals | Decreasing; < 4×10⁻⁸ |
| `bands` | Bands not converged | Reaches 0 |

## Reading Energy Output

```
Energy contributions relative to reference atoms: (SHE)
  Kinetic:       +XX.XX eV
  Potential:     -XX.XX eV
  External:      +XX.XX eV
  XC:            -XX.XX eV
  Entropy (-ST): -XX.XX eV    ← smearing entropy
  Local:         -XX.XX eV
  --------------------------
  Free energy:   -XX.XXXX eV  ← E_DFT (includes -TS)
  Extrapolated:  -XX.XXXX eV  ← E(T→0) via Methfessel-Paxton or FD
```

Use `Extrapolated` energy for total energy differences (reaction energies, adsorption, etc.).
Use `Free energy` only when comparing to thermodynamic quantities at finite T.

## Failure Pattern Diagnostics

### Pattern 1: Energy oscillates without decreasing
```
iter:  5  -101.12  -0.3   3.2e-02  8.1e-02
iter:  6  -100.98  -0.1   4.1e-02  9.3e-02
iter:  7  -101.15  -0.2   3.8e-02  8.7e-02
```
**Cause:** Mixer beta too large, or wrong magnetic state.
**Fix:**
```python
mixer = dict(backend='msr1', beta=0.02, nmaxold=10, weight=100)
```
Also check: are initial `magmoms` physically reasonable?

### Pattern 2: `|dn|` stuck at plateau, eigenstates OK
```
iter: 12  -101.45  -3.2   8.3e-03  1.1e-09
iter: 13  -101.45  -3.2   8.2e-03  9.8e-10
iter: 14  -101.45  -3.2   8.3e-03  1.2e-09
```
**Cause:** Magnetic system converged to metastable non-magnetic state; or insufficient empty bands.
**Fix:**
```python
atoms.set_initial_magnetic_moments([2.2, 2.2, ...])  # reset and restart
nbands = -10  # more empty bands
```

### Pattern 3: `bands` column never reaches 0
```
iter: 20  -101.45  -5.1   3.1e-05  4.0e-09  12
iter: 21  -101.45  -5.2   2.9e-05  3.8e-09  12
```
**Cause:** Not enough bands or eigensolver not given enough iterations.
**Fix:**
```python
nbands = -15                                    # more empty bands
eigensolver = dict(name='ppcg', niter=5)     # more eigensolver iterations
```

### Pattern 4: Poisson solver error
```
Poisson solver did not converge!
```
**Cause A:** Cubic cell with molecular symmetry confuses FFT Poisson.
**Fix A:** Break cubic symmetry:
```python
mol.center(vacuum=6.0)
mol.cell[1,1] += 0.0001
mol.cell[2,2] += 0.0002
```
**Cause B:** gpts not divisible by 8.
**Fix B:** Set `gpts` explicitly: `gpts=(32, 32, 32)`.

### Pattern 5: Very slow convergence (> 100 iterations) for metals
**Cause:** Fermi surface instability; smearing width too small; or too few k-points.
**Fix:**
```python
kpts = dict(size=(12, 12, 12), gamma=True)
occupations = FermiDirac(0.2)              # start with wider smearing
mixer = dict(backend='msr1', beta=0.02, nmaxold=12, weight=150)
nbands = -15
```

### Pattern 6: `maxiter` reached
```
Did not converge!
```
**Cause:** Calculation is genuinely hard OR has a fixable physical issue.
**Approach:** Fix physical issue first (magmoms, smearing, k-points). Only increase `maxiter=500` as last resort after mixer tuning.

## Useful grep Commands

```bash
# Track energy per iteration
grep 'iter:' calc.txt | awk '{print $2, $4}'

# Check convergence
grep -E 'Converged|Did not converge' calc.txt

# Find magnetic moments
grep 'Magnetic' calc.txt

# Find Fermi level
grep 'Fermi' calc.txt

# Check number of k-points used
grep 'k-points' calc.txt

# Check total computation time
grep 'Total' calc.txt | tail -5
```

## Extracting Results from a Converged Log

```python
from gpaw import GPAW

calc = GPAW('calc.gpw', txt=None)

# Energies
e_free = calc.get_potential_energy()           # free energy (includes -TS)
e_extrap = calc.get_potential_energy(force_consistent=False)  # extrapolated

# Electronic structure
ef = calc.get_fermi_level()                    # Fermi level (eV)
eigs = calc.get_eigenvalues(kpt=0, spin=0)     # eigenvalues at k-point 0
occ = calc.get_occupation_numbers(kpt=0, spin=0)

# Magnetic
mm = calc.get_magnetic_moments()               # per-atom moments
total_mm = calc.get_magnetic_moment()          # total moment

# Forces and stress
forces = calc.get_forces()                     # eV/Å, shape (N, 3)
stress = calc.get_stress()                     # Voigt notation, eV/Å³

# Density
rho = calc.get_pseudo_density()                # 3D array on grid
rho_ae = calc.get_all_electron_density()       # all-electron density
```
