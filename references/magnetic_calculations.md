# GPAW Magnetic Calculations Reference

## Initial Magnetic Moments

Magnetic calculations are sensitive to initialization. Always set physically meaningful initial moments.

```python
# Method 1: on Atoms object (before attaching calculator)
atoms.set_initial_magnetic_moments([2.2, 2.2])          # ferromagnetic Fe BCC
atoms.set_initial_magnetic_moments([2.2, -2.2])         # antiferromagnetic

# Method 2: in Atoms constructor
from ase import Atoms
atoms = Atoms('Fe2', ..., magmoms=[2.2, 2.2])

# Method 3: get reasonable initial values
# Element    spin-only moment (μB)
# Fe (bcc)   2.2
# Co (hcp)   1.6
# Ni (fcc)   0.6
# Mn         3.0–5.0 (depends on phase)
# Cr         antiferromagnetic, ±0.4 in Cr bulk
```

**Rule:** If the expected moment is non-zero, always set it. GPAW may converge to a non-magnetic solution if initialized at zero.

## Spin-Polarized Calculator Settings

```python
from gpaw import GPAW, PW, FermiDirac

calc = GPAW(
    mode=PW(500),
    xc='PBE',
    kpts=dict(size=(12, 12, 12), gamma=True),  # dense k-mesh for metals
    occupations=FermiDirac(0.1),
    spinpol=True,                                 # explicit; don't rely on auto-detect
    nbands=-15,                                   # extra empty bands per spin channel
    mixer=dict(
        backend='msr1',
        beta=0.02,                             # lower beta for magnetic systems
        nmaxold=10,
        weight=100,
        method='separate',                     # separate up/down mixing
    ),
    convergence=dict(energy=5e-4, density=1e-5, eigenstates=4e-8),
    txt='magnetic.txt',
)
```

## Post-Calculation: Reading Magnetic Moments

```python
from gpaw import GPAW

calc = GPAW('magnetic.gpw', txt=None)
atoms = calc.get_atoms()

mm_per_atom = calc.get_magnetic_moments()    # shape (N,), in μB
mm_total    = calc.get_magnetic_moment()     # scalar total moment

print('Per-atom moments:', mm_per_atom)
print('Total moment:    ', mm_total, 'μB')

# Spin-up/down eigenvalues
e_up   = calc.get_eigenvalues(kpt=0, spin=0)
e_down = calc.get_eigenvalues(kpt=0, spin=1)
```

## DFT+U

### When to Use

Apply DFT+U when PBE incorrectly gives a metallic ground state for a known insulator (NiO, MnO, CoO, Fe₂O₃), or when d/f band positions are systematically wrong.

### Setup Syntax

```python
# Single orbital
setups = dict(Ni=':d,6.0')     # U_eff = 6.0 eV on Ni-d
setups = dict(Fe=':d,4.0')     # U_eff = 4.0 eV on Fe-d
setups = dict(Ce=':f,5.0')     # U_eff = 5.0 eV on Ce-f

# Multiple elements
setups = dict(Mn=':d,4.0', O=':p,3.0')

# Multiple orbitals on one element (semicolon separator)
setups = dict(Co=':d,5.0,0;p,1.0,0')   # d and p, VASP normalization

# Normalization: default=GPAW convention; ',0'=VASP convention
# Difference is negligible for d/f; significant for p-orbital U
```

### Choosing U Values

Published U values for common systems (PBE+U, U_eff = U−J):

| System | Element | Orbital | U_eff (eV) | Reference |
|--------|---------|---------|------------|-----------|
| NiO | Ni | d | 5.0–7.0 | Dudarev et al. |
| CoO | Co | d | 5.0–6.0 | |
| MnO | Mn | d | 4.0–6.0 | |
| Fe₂O₃ | Fe | d | 4.0–5.0 | |
| LiFePO₄ | Fe | d | 4.3 | Ong et al. |
| VO₂ | V | d | 2.0–4.0 | |
| CeO₂ | Ce | f | 5.0–6.0 | |
| UO₂ | U | f | 4.0 | |

Always verify: the insulating gap should open; moment should match experiment.

## Spin-Orbit Coupling (SOC)

SOC is applied non-self-consistently (on top of scalar-relativistic SCF):

```python
from gpaw import GPAW
from gpaw.spinorbit import soc_eigenstates

# Requires a converged band structure calculation (symmetry='off')
calc = GPAW('band_structure.gpw', txt=None)

soc = soc_eigenstates(
    calc,
    n1=0,      # band range start
    n2=30,     # band range end (exclusive)
    theta=0.0, # spin quantization polar angle (deg)
    phi=0.0,   # azimuthal angle (deg)
    scale=1.0, # SOC strength scaling (1.0 = physical)
)

e_soc = soc.eigenvalues()         # shape: (k, n2-n1)
sz = soc.spin_projection(v=2)     # ⟨σ_z⟩, shape: (k, n2-n1)
```

**When SOC matters:** Heavy elements (Pt, Au, W, Bi, Pb, Tl, Hg), non-centrosymmetric systems (Rashba/Dresselhaus), topological materials, magnetic anisotropy energy.

**Self-consistent SOC:** Available in new GPAW via `spinors=True` in calculator settings (experimental). Not recommended for production without benchmarking.

## Spin Spiral Calculations

For non-collinear spin configurations:

```python
calc = GPAW(
    mode=PW(500),
    xc='PBE',
    kpts=dict(size=(8, 8, 8), gamma=True),
    experimental=dict(magmoms=[[0, 0, 2.2], [0, 0, -2.2]]),  # non-collinear
    txt='spinspiral.txt',
)
```

See `tutorialsexercises/magnetic/spinspiral/` for full workflow including constraint forces.

## Magnon Dispersion (Magnetic Force Theorem)

```python
from gpaw.response.mft import MFT

mft = MFT(
    calc='fe_gs.gpw',
    kpts=dict(path='GH', npoints=50),
)
magnons = mft.calculate()
```

Requires: converged ferromagnetic ground state with tight `eigenstates` convergence.

## Fixed Magnetic Moment

To constrain the total magnetic moment during SCF:

```python
occupations = dict(name='fermi-dirac', width=0.1, fixmagmom=True)
```

This fixes the total moment to the initial value from `magmoms`. Useful when comparing energies of different magnetic configurations at the same moment.
