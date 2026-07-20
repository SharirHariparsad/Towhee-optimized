# Towhee 8.2.3 Simulation Setup Guide

Towhee is a Monte Carlo molecular simulation engine supporting NVT, NpT, GCMC, and Gibbs ensembles with numerous force fields. This guide walks through the complete workflow.

---

## 1. Build the Software

```bash
cd src/towhee-8.2.3
./configure
make
```

This produces `Source/towhee` (the main executable) and builds all utilities in `Utils/`. For MPI parallel support: `./configure --enable-mpi`.

---

## 2. Understand the Input/Output File System

| File | Direction | Purpose |
|------|-----------|---------|
| `towhee_input` | Input | Main simulation control file (required) |
| `towhee_initial` | Input | Restart file (coordinates, box, etc.) |
| `towhee_coords` | Input | Raw coordinates for `initstyle='coords'` |
| `towhee_par*` | Input | Partial coordinates for CBMC building |
| `towhee_backup` | Output | Periodic restart backup |
| `towhee_final` | Output | Final configuration |
| `box_*_step_*.pdb` | Output | PDB trajectory files |

---

## 3. Structure of `towhee_input`

The input file is keyword-driven, free-format. Keywords are case-insensitive. A full description is in `Manual/projects/towhee/input/towhee_input_v8_2_x.html`.

### 3a. Simulation Control Block

These are the first keywords in the file:

```
ensemble          'nvt' | 'npt' | 'uvt'          ! NVT, NpT, or grand canonical
temperature       <temp>                          ! Kelvin
pressure          <press>                         ! [bar] (NpT only)
chempot           <mu>                            ! Chemical potential (uVT only)
nmolty            <N>                             ! Number of molecule types
nmolectyp         <n1 n2 ...>                     ! Count of each molecule type
numboxes          1 or 2                          ! 2 for Gibbs ensemble
stepstyle         'cycles' or 'moves'             ! How nstep is interpreted
nstep             <N>                             ! Total steps/cycles
printfreq         <N>                             ! Output frequency
blocksize         <N>                             ! Block averaging size
moviefreq         <N>                             ! Trajectory dump freq
backupfreq        <N>                             ! Restart file frequency
```

### 3b. Initialization Block

```
linit             .true. | .false.                ! Generate new config (T) or restart
initboxtype       'dimensions'                    ! Box specification via hmatrix
initstyle         'full cbmc' | 'coords' | etc.  ! How to place molecules
initlattice       'simple cubic' | 'fcc' | etc.   ! Lattice type for initial placement
initmol           <n1 n2 ...>                     ! Molecules per box
inix iniy iniz    <nx ny nz>                      ! Lattice dimensions
hmatrix           <3x3 matrix for 1 box, 6x3 for 2>
```

### 3c. Monte Carlo Move Probabilities

Each MC move type has a keyword and subordinate parameters:

```
pmvol             0.002d0                         ! Volume move probability
  pmvlpr          1.0d0                           ! Per-molecule probability
  rmvol           0.1d0                           ! Max volume displacement
  tavol           0.5d0                           ! Target acceptance

pm2boxcbswap      0.01d0                          ! Gibbs swap (type + box pair)
  pm2cbswmt       1.0d0                           ! Per-molecule-type prob
  pm2cbswpr       1.0d0                           ! Per-box-pair prob

pmuvtcbswap       0.25d0                          ! Grand canonical swap
  pmuvtcbmt       1.0d0

pmcb              0.33d0                          ! Configurational-bias regrowth
  pmcbmt          1.0d0
  pmall           0.0d0                           ! All-atom regrowth fraction

pmback            0.0d0                           ! Crankshaft
pmpivot           0.0d0                           ! Pivot
pmconrot          0.0d0                           ! Concerted rotation
pmtracm           0.75d0                          ! Center-of-mass translation
  pmtcmt          1.0d0
  rmtrac          0.5d0                           ! Max translation (Angstrom)
  tatrac          0.5d0                           ! Target acceptance
pmrotate          1.0d0                           ! Rotation
  pmromt          1.0d0
  rmrot           0.05d0                          ! Max rotation (radians)
  tarot           0.5d0
```

**Rules of thumb:**
- Move probabilities are cumulative -- the difference between consecutive `pm*` values is the actual probability of that move type
- The last move type must have `pm* = 1.0d0`
- For `nmolty > 1`, `pm*mt` arrays contain probabilities per molecule type (cumulative)

### 3d. Force Field Block

```
ffnumber          <N>                             ! Number of force field files
ff_filename       /path/to/towhee_ff_<name>       ! Path(s) to force field file(s)
classical_potential  'Lennard-Jones' | 'Mie' | 'Buckingham' | etc.
classical_mixrule    'Lorentz-Berthelot' | 'Geometric'
lshift            .false.                         ! Shifted potential
ltailc            .true.                          ! Tail corrections
rmin              1.0d0                           ! Inner cutoff
rcut              14.0d0                          ! Cutoff radius (Angstrom)
rcutin            10.0d0                          ! Inner cutoff for CBMC
electrostatic_form  'none' | 'coulomb'
coulombstyle      'ewald_fixed_kmax' | 'ewald_auto' | etc.
  kalp            5.6d0                           ! Ewald alpha
  kmax            5                               ! Ewald k-vectors
  dielect         1.0d0                           ! Dielectric constant
```

**Available force fields** (in `ForceFields/`):

| File | Type |
|------|------|
| `towhee_ff_TraPPE-UA` | TraPPE united-atom (alkanes, alcohols, etc.) |
| `towhee_ff_TraPPE-EH` | TraPPE explicit-hydrogen |
| `towhee_ff_OPLS-aa` | OPLS all-atom |
| `towhee_ff_OPLS-ua` | OPLS united-atom |
| `towhee_ff_Amber96` | AMBER96 |
| `towhee_ff_Charmm27` | CHARMM27 |
| `towhee_ff_UFF` | Universal Force Field |
| `towhee_ff_DREIDING` | DREIDING |
| `towhee_ff_COMPASSv1` | COMPASS |
| `towhee_ff_MMFF94` | MMFF94 |
| `towhee_ff_SPC-E` | SPC/E water |
| `towhee_ff_TIP3P/4P/5P` | TIP3P/4P/5P water |
| `towhee_ff_LJium` | Lennard-Jonesium |
| `towhee_ff_HardSphere` | Hard sphere |
| `towhee_ff_SquareWell` | Square well |

### 3e. Molecule Definition Block

One block per molecule type. Multiple input styles are available:

**Style 0 — Basic Connectivity Map (most common):**

```
input_style       'basic connectivity map'
nunit             <N>                             ! Number of atoms/beads
nmaxcbmc          <N>                             ! Max CBMC chain length
lpdbnames         F                               ! PDB atom naming
forcefield        '<FF name>'                     ! Must match ForceFields file
charge_assignment 'bond increment' | 'manual'

! For each unit (1 to nunit):
unit ntype        <i>  '<atom type name>'         ! From force field file
vibration         <nbonds>
                  <atom1> <atom2> ...
improper torsion  <0 | type>
! If 'manual' charge assignment:
unit ntype qqatom <i>  '<type>'  <charge>
```

**Style 2 — Manual input (full control):**
Used for specifying bonds, angles, torsions explicitly when the force field file doesn't contain sufficient pattern matching.

---

## 4. Step-by-Step Simulation Workflow

### Example 1: NVT Simulation of n-Pentane (TraPPE-UA)

Create `towhee_input`:

```
ensemble
'nvt'
temperature
336.0d0
nmolty
1
nmolectyp
300
numboxes
1
stepstyle
'cycles'
nstep
100000
printfreq
1000
blocksize
5000
moviefreq
10000
backupfreq
5000
runoutput
'full'

linit
.true.
initboxtype
'dimensions'
initstyle
'full cbmc'
initlattice
'simple cubic'
initmol
300
inix iniy iniz
8    8    7
hmatrix
30.0d0 0.0d0 0.0d0
0.0d0 30.0d0 0.0d0
0.0d0 0.0d0 30.0d0

pmvol
0.002d0
  pmvlpr    1.0d0
  rmvol     0.1d0
  tavol     0.5d0
pmcb
0.33d0
  pmcbmt    1.0d0
  pmall     0.0d0
pmtracm
0.67d0
  pmtcmt    1.0d0
  rmtrac    0.5d0
  tatrac    0.5d0
pmrotate
1.0d0
  pmromt    1.0d0
  rmrot     0.05d0
  tarot     0.5d0

ffnumber
1
ff_filename
/your/path/ForceFields/towhee_ff_TraPPE-UA
classical_potential
'Lennard-Jones'
classical_mixrule
'Lorentz-Berthelot'
lshift
.false.
ltailc
.true.
rmin
1.0d0
rcut
14.0d0
rcutin
10.0d0
electrostatic_form
'none'

input_style
'basic connectivity map'
nunit
5
nmaxcbmc
5
lpdbnames
F
forcefield
'TraPPE-UA'
charge_assignment
'bond increment'
unit ntype
1    'CH3*(sp3)'
vibration
1
2
improper torsion
0
unit ntype
2    'CH2**(sp3)'
vibration
2
1 3
improper torsion
0
unit ntype
3    'CH2**(sp3)'
vibration
2
2 4
improper torsion
0
unit ntype
4    'CH2**(sp3)'
vibration
2
3 5
improper torsion
0
unit ntype
5    'CH3*(sp3)'
vibration
1
4
improper torsion
0
```

Run: `./Source/towhee`

### Example 2: Gibbs Ensemble (VLE of n-Pentane)

Use the input shown above (Gibbs_Ensemble/TraPPE_Pentane example). Key differences from NVT:
- `numboxes 2`
- `initmol 260 \n 40` (molecules split between boxes)
- `hmatrix` has two 3x3 blocks
- Swap moves enabled: `pm2boxcbswap 0.01d0`, `pmavb1 0.02d0`

### Example 3: Grand Canonical (Ethane in Amber96)

Key differences:
- `ensemble 'uvt'`
- `chempot 14359.426d0`
- `pmuvtcbswap 0.25d0`
- Molecule definition uses `charge_assignment 'manual'` with explicit `qqatom` values

---

## 5. Output Analysis

Key output files:
- **Standard output** (stdout or `towhee_output`): Energies, pressures, densities, acceptance ratios
- **`towhee_chempotdata`**: Widom insertion chemical potential data
- **`towhee_backup`**: Restart file (rename to `towhee_initial` to restart)
- **`towhee_final`**: Final configuration
- **`box_*_*.pdb`**: PDB trajectories (view in VMD, PyMOL, etc.)
- **`towhee_citations`**: Paper citations for the algorithms used

---

## 6. Utility Programs

| Utility | Purpose |
|---------|---------|
| `pdb2towhee` | Convert PDB to Towhee input |
| `xyz2towhee` | Convert XYZ to Towhee format |
| `analyse_movie` | Analyze trajectory files |
| `analyse_histogram` | TMMC/histogram reweighting analysis |
| `fitcoex` | Fit vapor-liquid coexistence curve |
| `car2towhee` | Convert MSI/Accelrys CAR format |

Build individually: `cd Utils && make pdb2towhee`

---

## 7. Key Tips

1. **Adjust max displacements** (rmtrac, rmrot, rmvol) to achieve target acceptance (~0.3-0.5). The simulation adjusts them automatically when `tatrac`/`tavol` is set, but the adjustment frequency is controlled by `trmaxdispfreq`/`volmaxdispfreq`.

2. **Force field atom type names** must match exactly what is in the `ForceFields/towhee_ff_*` file. Check the Nonbonded Types section of the force field file for valid names.

3. **For charged molecules**, ensure `electrostatic_form 'coulomb'` is set and choose an appropriate `coulombstyle`. For most condensed-phase simulations, `ewald_fixed_kmax` is recommended.

4. **CBMC parameters**: `cbmc_formulation` and `cbmc_setting_style` control how the configurational-bias algorithm works. For TraPPE-UA use `'Martin and Frischknecht 2006'` / `'Martin and Frischknecht'`. For Amber/OPLS use `'Martin and Siepmann 1999 + Martin and Thompson 2004'` / `'default ideal'`.

5. **Restarting**: Copy `towhee_backup` to `towhee_initial` and set `linit .false.` in the input file.

6. **Parallel runs**: Build with `./configure --enable-mpi`, then use `-p jobfarm` flag with `towhee_parallel` input file.
