# GPAW PAW Setups Reference

## What Are PAW Setups?

PAW (Projector Augmented Wave) setups serve the same role as pseudopotentials but within the all-electron PAW framework. They encode the frozen-core approximation, partial waves, and projector functions for each element.

## Installation and Path

```bash
# Install latest setups (v24.11.0 as of 2025)
gpaw install-data ~/gpaw-setups

# Manual install
tar -xf gpaw-setups-24.11.0.tar.gz -C ~/
export GPAW_SETUP_PATH=~/gpaw-setups-24.11.0

# Verify
gpaw info
```

`GPAW_SETUP_PATH` accepts colon-separated paths; first match wins.

## Default Setup Selection

The default `setups='paw'` is correct for most elements. Only specify explicitly when:
1. A non-default variant is needed for accuracy
2. Applying DFT+U (via the `:orbital,U` syntax)
3. Using ghost atoms or all-electron reconstruction

## Element-Specific Recommendations

| Element | Issue | Recommendation |
|---------|-------|----------------|
| **Cr** | Default 6-electron setup misses [Ar]3p semi-core | `setups=dict(Cr='14')` for accurate d-band |
| **Lanthanides** (La–Lu) | Require v24.11.0+ setups | `gpaw install-data` with latest package |
| **Early 3d TM** (Ti, V, Cr, Mn) | Semi-core 3s3p important for oxides | Consider harder setups or explicit semi-core |
| **Alkali metals** (Li, Na, K) | Semi-core may matter for intercalation | Check convergence; typically default is fine |
| **H in DFT+U** | Avoid U on H | Only apply U to TM d/f orbitals |

## DFT+U Setup Syntax

The DFT+U correction is applied through the `setups` keyword, not a separate parameter:

```python
# Single correction
setups = dict(Ni=':d,6.0')          # U_eff=6 eV on Ni d-orbitals
setups = dict(Fe=':d,4.0')          # U_eff=4 eV on Fe d-orbitals
setups = dict(U=':f,4.0')           # U_eff=4 eV on U f-orbitals

# Multiple orbitals on same element
setups = dict(Ni=':d,4.0,0;p,2.0,0')  # d and p corrections, VASP normalization

# Multiple elements
setups = dict(Mn=':d,4.0', O=':p,3.0')

# With non-default dataset AND U correction
# Not directly supported — use default dataset + DFT+U syntax
```

**Normalization flag:**
- Default (no flag): GPAW normalization — projections normalized to augmentation sphere
- `,0` suffix: VASP convention — no normalization
- Effect is minimal for d/f orbitals (>90% within augmentation sphere)
- Effect is significant for p orbitals — match to published U values accordingly

## Available Non-PAW Options

```python
setups = 'sg15'    # SG15 ONCV norm-conserving pseudopotentials
setups = 'hgh'     # Hartwigsen-Goedecker-Hutter GTH pseudopotentials
setups = 'hgh.sc'  # HGH with semi-core electrons
```

These are rarely needed; use PAW unless a specific method (e.g., SG15 for GW with external codes) requires it.

## Ghost Atoms

Ghost atoms carry basis functions but no electrons/nuclear charge. Used for basis set superposition error (BSSE) correction:

```python
from ase import Atoms
from gpaw import GPAW

# Set ghost at atom index 1 (by integer key)
setups = {1: 'ghost'}

# Set ghost by element (affects ALL atoms of that element)
setups = dict(H='ghost')
```

## All-Electron Reconstruction

For post-processing only (never use in SCF):

```python
# Reconstruct all-electron density for analysis
setups = dict(H='ae', O='ae')
```

## Available Functionals in Setup Files

Each element has PAW datasets generated for:
- LDA
- PBE
- revPBE
- RPBE
- GLLBSC (for band gap calculations via GLLB-SC functional)

The functional in `setups` must match `xc`. GPAW handles this automatically when you set `xc='PBE'`.

## ACWF Benchmark Reference

The GPAW PAW datasets are benchmarked against the ACWF (Absolute Convergence With FLEUR) reference set. At `ecut=1000 eV`, maximum absolute errors across the periodic table are below 1%. At `ecut=500 eV` (production), errors are < 5 meV/atom for most elements.
