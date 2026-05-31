# Beyond-DFT Methods in GPAW

## G0W0 Quasiparticle Corrections

### Workflow

```python
# Step 1: DFT ground state — save ALL wavefunctions
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

atoms = bulk('Si', 'diamond', a=5.43)
calc = GPAW(
    mode=PW(300),
    xc='LDA',                            # LDA recommended as GW starting point
    kpts=dict(size=(4, 4, 4), gamma=True),
    occupations=FermiDirac(0.001),
    convergence=dict(eigenstates=1e-8),
    txt='gs.txt',
)
atoms.calc = calc
atoms.get_potential_energy()
calc.diagonalize_full_hamiltonian()      # compute ALL bands
calc.write('gs_allbands.gpw', 'all')    # 'all' saves wavefunctions

# Step 2: G0W0
from gpaw.response.g0w0 import G0W0

gw = G0W0(
    calc='gs_allbands.gpw',
    bands=(3, 5),             # (HOMO index, LUMO index) → 0-indexed from valence bottom
    nbands=30,                # bands for self-energy sum (converge this)
    ecut=100,                 # local field cutoff eV (converge this)
    ecut_extrapolation=True,  # 1/E^(3/2) extrapolation; provide list for multi-point
    integrate_gamma='WS',     # Wigner-Seitz k→0 correction (recommended)
    domega0=0.02,             # frequency grid spacing (eV)
    omega2=15,                # max frequency (eV)
    ppa=False,                # full frequency integration (ppa=True is ~10× faster)
    filename='Si_gw',
)
result = gw.calculate()
```

### Key Convergence Parameters

| Parameter | Effect | Converged range |
|-----------|--------|-----------------|
| `ecut` | Basis for screened interaction | 100–400 eV; use `ecut_extrapolation` |
| `nbands` | Self-energy sum | 3–5× number of occupied bands |
| k-points | BZ sampling of screened W | 4×4×4 for testing; 8×8×8 production |
| `domega0` | Frequency grid | 0.02 eV |

### Results Structure

```python
# result.pckl contains:
result['qp']    # quasiparticle energies, shape: (spins, k-IBZ, bands)
result['eps']   # KS eigenvalues
result['sigma'] # self-energy
result['Z']     # renormalization factor

# Band gap from result
vbm = result['qp'][0, :, 3].max()   # highest VB across all k
cbm = result['qp'][0, :, 4].min()   # lowest CB across all k
gap = cbm - vbm
```

### Variants

```python
# Plasmon-Pole Approximation (PPA) — 10× faster, ~0.1 eV less accurate
gw = G0W0(..., ppa=True)

# GW with vertex corrections
gw = G0W0(..., xc='rALDA', fxc_mode='GWG')

# 2D systems (MoS₂, graphene)
gw = G0W0(..., truncation='2D', q0_correction=True)
```

## Bethe-Salpeter Equation (BSE)

```python
from gpaw.response.bse import BSE

# Requires: converged GW calculation or DFT groundstate
bse = BSE(
    calc='gs_allbands.gpw',
    valence_bands=range(4, 8),      # occupied band indices
    conduction_bands=range(8, 12),  # unoccupied band indices
    nbands=30,
    ecut=50,
    gw_skn=None,                    # pass GW QP energies to improve starting point
    filename='bse_Si',
)
bse.calculate()
```

**Convergence:** k-points and `ecut` dominate. For 2D: use `truncation='2D'` and high k-density.

## Hybrid Functionals (PBE0 / HSE06)

### Non-Self-Consistent (Recommended for Band Gaps)

```python
from gpaw.hybrids.eigenvalues import non_self_consistent_eigenvalues
from gpaw import GPAW

# PBE ground state
calc = GPAW('pbe_gs.gpw', txt=None)
# PBE0 eigenvalues (non-SCF)
e_pbe0 = non_self_consistent_eigenvalues(calc, 'PBE0')
```

### Self-Consistent Hybrid

```python
calc = GPAW(
    mode=PW(400),
    xc='PBE0',                        # or 'HSE06'
    kpts=dict(size=(4, 4, 4), gamma=True),
    parallel=dict(sl_auto=True),        # ScaLAPACK required for hybrids
    txt='pbe0.txt',
)
```

**Cost:** SCF hybrid is ~100× more expensive than PBE. Always try non-SCF first.

## Linear-Response TDDFT (LR-TDDFT)

For finite systems (molecules, clusters). See main SKILL.md for full workflow.

```python
from gpaw.lrtddft import LrTDDFT, photoabsorption_spectrum

lr = LrTDDFT(
    calc,
    xc='LDA',          # XC kernel for response (can differ from ground state)
    istart=0,          # first electron-hole pair
    jend=20,           # last unoccupied state
    nspins=2,          # 2 = include singlet+triplet; 1 = singlet only
)
```

**Key convergence parameter:** `jend` — include enough unoccupied states to cover the energy range of interest.

## RPA Correlation Energies

```python
from gpaw.response.rpa import RPACorrelation

rpa = RPACorrelation(
    calc='gs_allbands.gpw',
    ecut=[50, 75, 100, 150, 200],    # list → extrapolation
    filename='rpa',
)
e_rpa = rpa.calculate()
```

Used for: van der Waals binding, accurate reaction energies, benchmark-quality results.
Requires: tight DFT groundstate + full wavefunction diagonalization.

## GLLB-SC Band Gap Functional

Orbital-dependent functional giving accurate band gaps at DFT cost:

```python
calc = GPAW(
    mode=PW(500),
    xc='GLLBSC',
    kpts=dict(size=(8, 8, 8), gamma=True),
    txt='gllb.txt',
)
atoms.calc = calc
atoms.get_potential_energy()

# Extract discontinuity correction
response = calc.hamiltonian.xc.response
dxc = response.get_band_gap()
ks_gap = ...  # from eigenvalues
qp_gap = ks_gap + dxc
```

**Note:** GLLB-SC requires the GLLBSC setup files (installed with the default setup package).
