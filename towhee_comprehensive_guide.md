# Towhee 8.2.3 — Comprehensive Simulation Guide for Beginners

This guide is written for someone new to computational chemistry. It explains not just *how* to set up Towhee simulations, but *why* you make each choice. The concepts here apply to molecular simulation in general, not just Towhee.

---

## Table of Contents

1. [Preface: What Are We Actually Doing?](#1-preface-what-are-we-actually-doing)
2. [Part 1: Monte Carlo vs Molecular Dynamics](#2-part-1-monte-carlo-vs-molecular-dynamics)
3. [Part 2: Ensembles — Which One Do I Need?](#3-part-2-ensembles--which-one-do-i-need)
4. [Part 3: Force Fields — Choosing the Right Physics](#4-part-3-force-fields--choosing-the-right-physics)
5. [Part 4: Monte Carlo Moves — How the Simulation Explores](#5-part-4-monte-carlo-moves--how-the-simulation-explores)
6. [Part 5: Parameter Selection — From Beginner to Confident](#6-part-5-parameter-selection--from-beginner-to-confident)
7. [Part 6: Reading the Output — A Complete Walkthrough](#7-part-6-reading-the-output--a-complete-walkthrough)
8. [Part 7: Interpreting Results — What Does It All Mean?](#8-part-7-interpreting-results--what-does-it-all-mean)
9. [Part 8: Worked Examples from the Codebase](#9-part-8-worked-examples-from-the-codebase)
10. [Part 9: End-to-End Workflow](#10-part-9-end-to-end-workflow)

---

## 1. Preface: What Are We Actually Doing?

**The core problem:** We want to know macroscopic properties of a substance (density, vapor pressure, heat of vaporization, etc.) starting from only our knowledge of the molecules and how they interact.

**The core idea:** If we generate many random configurations of our molecules and measure the energy of each one, the *average* behavior of those configurations will match the *real* behavior of the substance — so long as the configurations are generated correctly. The trick is generating them correctly.

**What Towhee does:** It starts molecules in a box, then repeatedly makes small random changes (moving a molecule, adding one, changing the box size, etc.). It uses a mathematical rule (the Metropolis algorithm) to decide whether to keep each change. Over time, the system naturally explores the configurations that matter most, and we can average properties over those configurations.

**The key insight you need:** The *ensemble* you choose determines what stays constant during the simulation. The *force field* determines the energy when molecules are near each other. The *Monte Carlo moves* determine how the system explores new configurations. Understanding these three choices is 90% of setting up a correct simulation.

---

## 2. Part 1: Monte Carlo vs Molecular Dynamics

There are two main families of molecular simulation:

| Property | Molecular Dynamics (MD) | Monte Carlo (MC) |
|----------|------------------------|-------------------|
| How it moves | Follows physical forces (Newton's laws) | Random moves accepted/rejected by energy |
| Time | Yes — you get dynamics (diffusion, etc.) | No — no time information |
| Time scale limit | ~microsecond (expensive) | No physical limit (good for rare events) |
| Phase equilibria | Difficult | Easy (Gibbs ensemble) |
| New molecules | Cannot add/remove | Can add/remove (grand canonical) |
| Free energies | Complex | Straightforward |

**When to choose MC over MD (Towhee's strengths):**
- You want vapor-liquid coexistence data (Gibbs ensemble)
- You want adsorption isotherms (grand canonical)
- You have long, branched, or cyclic molecules (CBMC makes this efficient)
- You need chemical potentials
- The system has slow dynamics that MD cannot sample in reasonable time

---

## 3. Part 2: Ensembles — Which One Do I Need?

An *ensemble* is the set of conditions that stay constant during the simulation. Each ensemble corresponds to a different real-world experimental situation.

### The Four Main Ensembles in Towhee

#### Canonical Ensemble (NVT) — `ensemble 'nvt'`

**What stays constant:** Number of molecules (N), Volume (V), Temperature (T)

**Real-world analogue:** A sealed, rigid container in a temperature-controlled bath.

**What you can calculate:** Internal energy, pressure (via virial), radial distribution functions, molecular conformations, radius of gyration.

**When to use it:**
- You are testing a force field against known data
- You are equilibrating a system before switching to another ensemble
- You want to study molecular structure (e.g., protein folding, polymer conformations)
- You have a system at a known density and just need energies

**Examples in the codebase:**
- `Canonical_Ensemble/TraPPE_Molecules` — 6 different TraPPE molecules at 400 K
- `Canonical_Ensemble/Charmm27_ubiquitin` — a protein at 300 K
- `Canonical_Ensemble/NaCl_1x1x1` — ionic crystal at 900 K (Fumi-Tosi potential)
- `Canonical_Ensemble/Fris_Walls` — LJ fluid between walls

#### Grand Canonical Ensemble (muVT) — `ensemble 'uvt'`

**What stays constant:** Chemical potential (mu), Volume (V), Temperature (T). The number of molecules *fluctuates*.

**Real-world analogue:** An open container connected to a reservoir of the chemical. Molecules freely enter and leave to maintain equilibrium.

**What you can calculate:** Adsorption isotherms (how many molecules enter the box at a given chemical potential/pressure), pore filling, capillary condensation.

**When to use it:**
- You are studying adsorption in porous materials (zeolites, nanotubes, MOFs)
- You want to compute an adsorption isotherm
- You are studying nucleation or cavitation
- You need the chemical potential-density relationship

**Key parameter:** `chempot` — this is the chemical potential you impose. Finding the right value is the hardest part. If you have an experimental vapor pressure, you can estimate it, or you can run multiple simulations at different mu values.

**Examples in the codebase:**
- `Grand_Canonical_Ensemble/Amber96_Ethane` — ethane adsorption
- `Grand_Canonical_Ensemble/Walt2001_Nanotube` — methane in a carbon nanotube
- `Grand_Canonical_Ensemble/Steele_Wall` — LJ fluid near a wall

#### Gibbs Ensemble (NVT-Gibbs) — `ensemble 'nvt'` with `numboxes 2`

**What stays constant:** Total molecules (N), Total volume (V), Temperature (T). Molecules swap between two boxes.

**This is Towhee's specialty.** The two boxes represent coexisting phases (liquid and vapor). Molecules move between boxes until the chemical potential and pressure are equal in both — which is exactly the condition for phase coexistence.

**What you can calculate:** Vapor-liquid coexistence densities, vapor pressure, heat of vaporization, critical temperature.

**When to use it:**
- You want the vapor-liquid phase diagram of a pure fluid or mixture
- You need boiling points, vapor pressures
- You are studying how mixture composition affects phase behavior

**Key parameters:**
- `nmolectyp` — total molecules across both boxes (usually 200-500)
- `initmol` — how many molecules in each box initially (e.g., 260 in liquid box, 40 in vapor box)
- `hmatrix` — two 3x3 matrices, one per box
- `pm2boxcbswap` — probability of swapping molecules between boxes (critical for this ensemble)
- `pmavb1` — aggregation-volume bias (increases swap acceptance)

**Examples in the codebase:**
- `Gibbs_Ensemble/TraPPE_Pentane` — n-pentane at 336 K (~0.7 Tr)
- `Gibbs_Ensemble/EPM_VLCC` — CO2 at 265 K (with VLCC output for fitting)
- `Gibbs_Ensemble/Cu_VLE` — copper at 4000 K (EAM potential for metals)
- `Gibbs_Ensemble/Gordon` — water/NaCl at 703 K (high T electrolyte)

#### Isobaric-Isothermal Ensemble (NpT) — `ensemble 'npt'`

**What stays constant:** Number of molecules (N), Pressure (P), Temperature (T). The box volume fluctuates.

**Real-world analogue:** A flexible container (like a balloon) in a temperature-controlled bath, with a constant external pressure.

**What you can calculate:** Density at a given T and P, compressibility, liquid structure.

**When to use it:**
- You want to know the density of a liquid at ambient conditions
- You are studying how pressure affects molecular structure
- You need initial configurations for other ensembles

**Key parameters:**
- `pressure` — in units of kPa (1 atm = 101.325 kPa)
- `pmvol` — probability of attempting a volume change
- `rmvol` — maximum fractional volume change

**Examples in the codebase:**
- `Isobaric_Isothermal_Ensemble/Charmm22_Ethanethiol` — ethanethiol at 298 K, 1 atm
- `Isobaric_Isothermal_Ensemble/Solid_LJium` — LJ solid at T=0.6, P=100 (uses cell moves, not just volume)
- `Isobaric_Isothermal_Ensemble/Au_Cu_Switch` — metal alloy with composition switching

### Advanced Ensembles

#### Transition Matrix Monte Carlo (TMMC) — `tmmc_flag .true.`

TMMC is like grand canonical but *much* more efficient. Instead of imposing a chemical potential, it builds a "transition matrix" that records how likely the system is to move between different numbers of molecules. From this matrix, you can extract the entire density of states and properties at *any* chemical potential in a single simulation.

**When to use it:** You need the full phase behavior over a range of conditions. TMMC is especially good for systems with free energy barriers (e.g., nucleation).

**Examples in the codebase:**
- `TMMC/SPC_E_WATER` — water at 500 K with CBMC, autofit Gaussian generation
- `TMMC/LJ` — Lennard-Jonesium at T=1.0

### Decision Tree

```
What do you want to know?
│
├─ "What is the density of this liquid at room conditions?"
│   → Isobaric-Isothermal (NpT)
│
├─ "What is the vapor-liquid coexistence curve?"
│   → Gibbs Ensemble (NVT with 2 boxes)
│
├─ "How much gas adsorbs in this porous material?"
│   → Grand Canonical (muVT)
│
├─ "What is the structure of this molecule in solution?"
│   → Canonical (NVT) — place molecule + solvent at known density
│
├─ "What is the free energy / phase diagram over a wide range?"
│   → TMMC
│
├─ "I'm testing/validating a force field"
│   → Canonical (NVT) — compare to experimental data at known conditions
│
├─ "I need to study a mixture"
│   → Any ensemble works. Use nmolty > 1 and set nmolectyp for each.
│     Gibbs for VLE of mixtures, NpT for liquid mixture densities.
│
└─ "I have a solid crystal"
    → NpT with cell fluctuations (pmcell), or NVT at experimental density
```

---

## 4. Part 3: Force Fields — Choosing the Right Physics

### What is a Force Field?

A force field is a mathematical formula that gives the energy of a system based on the positions of the atoms. It is a *model* — an approximation of reality. The quality of your simulation depends directly on the quality of this model.

In Towhee, force fields are stored in `ForceFields/towhee_ff_*` files. Each file contains:

1. **Nonbonded parameters** — sigma and epsilon for Lennard-Jones, or equivalent for other potentials
2. **Bond parameters** — equilibrium distances and force constants (or fixed bond lengths)
3. **Angle parameters** — equilibrium angles and force constants
4. **Torsion (dihedral) parameters** — energy barriers for rotation around bonds
5. **Charge information** — either as explicit charges or as "bond increments" that automatically assign charges based on bonding patterns

### Energy Formula (Simplified)

The total energy is:

```
E_total = E_bonds + E_angles + E_torsions + E_nonbonded + E_electrostatic

E_nonbonded = sum over atom pairs: 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
E_electrostatic = sum over atom pairs: q1*q2 / (4*pi*epsilon0*r)
```

### How to Choose a Force Field

This is the single most important decision in your simulation. Here is a guide based on what the Towhee examples show:

| System Type | Recommended Force Field | Example(s) |
|---|---|---|
| Alkanes (hydrocarbons) | **TraPPE-UA** (united-atom) | Gibbs_Ensemble/TraPPE_Pentane |
| Alcohols, ethers, carbonyls | **TraPPE-UA** or **OPLS-aa** | Canonical/TraPPE_Molecules |
| CO2 | **EPM2** | Gibbs_Ensemble/EPM_VLCC |
| Water | **SPC-E**, **TIP3P**, **TIP4P**, **TIP5P** | TMMC/SPC_E_WATER |
| Proteins | **Amber96**, **Charmm19/22/27** | Canonical/Amber_Villin |
| Nucleic acids | **Charmm27** | Canonical/Charmm27_Nucleic_Acid |
| Ionic liquids / salts | **FumiTosi** (Buckingham) | Canonical/NaCl_1x1x1 |
| Metals | **EAM** (Embedded Atom Method) | Gibbs_Ensemble/Cu_VLE |
| Semiconductors | **Stillinger-Weber** | Gibbs_Ensemble/Vink2001_Silicon |
| Zeolites | **Catlow** or **Dubb_Catlow** | Gibbs_Ensemble/Catlow_Zeolite_4a |
| Generic organics | **UFF** or **DREIDING** | Canonical/UFF |
| Pharmaceuticals | **MMFF94** | Isobaric_Isothermal/MMFF94_Ethyleneoxide |
| Simple models | **LJium**, **HardSphere**, **SquareWell** | Canonical/Hard_Sphere |
| Explosives | **Dick1994** + H-bond | Isobaric_Isothermal/Dick1994_PETN |
| Clays | **ClayFF** | — |

### United-Atom vs All-Atom

**United-atom (UA):** Hydrogen atoms on carbon are not explicit — they are "merged" into the carbon atom (e.g., CH3, CH2, CH are single interaction sites). This reduces the number of atoms by ~3x, making simulations much faster.

**All-atom (AA):** Every hydrogen is explicit. More accurate for some properties but 3x more expensive.

**Rule of thumb:** Use TraPPE-UA for hydrocarbons and simple organics. Use Amber, CHARMM, or OPLS-aa for systems where hydrogen bonding or detailed geometry matters (proteins, nucleic acids).

### What If I Need a Force Field That Isn't Pre-Installed?

Look in `ForceFields/` for existing files. There are 87 files. If your molecule is not covered:

1. Pick the closest force field (e.g., TraPPE-UA covers most organic functional groups)
2. Write a custom force field file by modifying an existing one
3. Use the `'basic connectivity map'` input style with atom type names that match your force field file

---

## 5. Part 4: Monte Carlo Moves — How the Simulation Explores

In MD, molecules move because forces push them. In MC, we propose random changes and accept/reject them based on energy. The *types* of changes we propose are called "moves."

Each move type has a probability (the `pm*` keywords), and these probabilities are **cumulative** — the actual probability of a move is the difference between its `pm*` and the previous `pm*`. The last move must have `pm* = 1.0`.

### Essential Moves and When to Use Them

#### Translation (`pmtracm`, `pmtraat`)

Picks a molecule and moves its center of mass (COM) or individual atoms by a random displacement.

**When you need it:** Always. This is the basic sampling move. Even in the simplest monatomic fluid, you need translations.

**Sub-parameters:**
- `rmtrac` — maximum COM displacement (Angstrom). Target ~30-50% acceptance.
- `tatrac` — target acceptance ratio. The code auto-adjusts `rmtrac` to hit this.
- `pmtcmt` — cumulative probabilities for each molecule type

**Examples from the codebase:**
- `Canonical_Ensemble/TraPPE_Molecules`: `pmtracm 0.75`, `rmtrac 0.5`, `tatrac 0.5`
- `Gibbs_Ensemble/TraPPE_Pentane`: `pmtracm 0.67`, `rmtrac 0.5`, `tatrac 0.5`

#### Rotation (`pmrotate`)

Rotates a molecule around its center of mass.

**When you need it:** For any non-spherical molecule. Monatomic molecules (like LJ spheres) don't need it.

**Sub-parameters:**
- `rmrot` — maximum rotation angle (radians). Start with ~0.05-0.1 radian.
- `tarot` — target acceptance ratio.
- `pmromt` — cumulative probabilities per molecule type

#### Volume Change (`pmvol`)

Changes the box volume isotropically (all dimensions scaled equally).

**When you need it:** In NpT ensembles to find the equilibrium density.

**Sub-parameters:**
- `rmvol` — maximum fractional volume change. Start with ~0.1-0.2.
- `tavol` — target acceptance (~0.3-0.5).
- `pmvlpr` — cumulative probability per molecule (not per type)

#### Cell Change (`pmcell`)

Changes the box shape non-isotropically (different dimensions independently).

**When you need it:** For solid-state simulations where the unit cell can change shape. Rarely needed for liquids.

**Example:** `Isobaric_Isothermal_Ensemble/Solid_LJium` uses `pmcell 0.005` (0.5% probability) with `pmvol 0.0`.

#### CBMC Swap Moves

##### Two-Box CB Swap (`pm2boxcbswap`)

Grows a molecule in one box while deleting it from the other (Gibbs ensemble).

**When you need it:** This is the core of the Gibbs ensemble. Without it, you cannot exchange molecules between the liquid and vapor boxes.

**Sub-parameters:**
- `pm2cbswmt` — per-molecule-type probabilities
- `pm2cbswpr` — per-box-pair probabilities

**Example:** `Gibbs_Ensemble/TraPPE_Pentane`: `pm2boxcbswap 0.01`

##### Grand Canonical CB Swap (`pmuvtcbswap`)

Creates or deletes a molecule in a single box (grand canonical ensemble).

**When you need it:** For adsorption simulations, cavitation, or any muVT simulation.

**Example:** `Grand_Canonical_Ensemble/Amber96_Ethane`: `pmuvtcbswap 0.25`
**Example:** `TMMC/SPC_E_WATER`: `pmuvtcbswap 0.65` (TMMC needs many swaps)

##### Single-Box CB Swap (`pm1boxcbswap`)

Creates/deletes within one box (used in mixtures where you need to change composition).

#### CBMC Regrowth (`pmcb`)

Regrows a portion of a molecule using Configurational-Bias Monte Carlo.

**When you need it:** For molecules with internal degrees of freedom (rotatable bonds). This is what makes Towhee efficient for complex molecules — a translation move alone would rarely find the right conformation.

**Sub-parameters:**
- `pmall` — fraction of CB moves that regrow the entire molecule (vs partial regrowth)
- `pmcbmt` — per-molecule-type probabilities

**Example:** `Gibbs_Ensemble/TraPPE_Pentane`: `pmcb 0.33` — 33% of moves are CB regrowth

#### Aggregation-Volume Bias (`pmavb1`, `pmavb2`, `pmavb3`)

AVB identifies "holes" in the liquid where there is free volume, and targets molecule insertions to those holes. This dramatically improves acceptance rates for swap moves in dense phases.

**When you need it:** In Gibbs or grand canonical simulations of dense liquids. Without AVB, insertions almost always fail because a randomly placed molecule would overlap with existing ones.

**Sub-parameters:**
- `avb1rad` — radius of the "aggregation volume" sphere (Angstrom). Typically ~5-8 A.
- `pmavb1in` — probability of inserting *in* vs *out* of a cluster

**Example:** `Gibbs_Ensemble/TraPPE_Pentane`: `pmavb1 0.02`, `avb1rad 8.0`

#### Conformational Moves (`pmback`, `pmpivot`, `pmconrot`, `pmcrback`)

These are specialized moves for changing the internal conformation of chain molecules:
- **Crankshaft** (`pmback`): Rotates a segment around the line connecting its two endpoints
- **Pivot** (`pmpivot`): Rotates a terminal segment around a bond
- **Concerted rotation** (`pmconrot`): Changes internal dihedrals in a coordinated way

**When you need them:** For long chain molecules, especially in dense phases where CB regrowth alone may be inefficient. Often not needed for small molecules.

#### Composite Moves (`pmcomposite`, `pm1boxcomswitch`)

Multiple molecules are moved simultaneously as a composite operation.

**When you need it:** In mixtures where swapping identities is useful (e.g., converting molecule type A to type B).

### Which Moves Do You Actually Need?

Here is a minimal set for each common scenario:

| System | Minimum Moves |
|--------|---------------|
| Monatomic fluid, NVT | Translation |
| Monatomic fluid, NpT | Translation + Volume |
| Small molecule (butane), NVT | Translation + Rotation + CB regrowth |
| Small molecule, Gibbs | Translation + Rotation + CB regrowth + Two-box CB swap + AVB1 |
| Small molecule, GCMC | Translation + Rotation + CB regrowth + Grand canonical CB swap + AVB1 |
| Long chain (decane), NVT | Translation + Rotation + CB regrowth (+ Pivot for dense phases) |
| Rigid molecule (water), NVT | Translation + Rotation |
| Solid, NpT | Translation + Cell change (no volume — cell handles it) |
| Protein, NVT | CB regrowth (for backbone + sidechains) + Translation + Rotation |

### Setting Move Probabilities

The probabilities are cumulative. If you have three moves:

```
pmtraat      0.10    ! 10% single-atom translation
pmtracm      0.40    ! 30% COM translation (0.40 - 0.10)
pmrotate     1.00    ! 60% rotation (1.00 - 0.40)
```

**Guidelines for first attempts:**
- Small molecules (1-5 atoms): Translation 40%, Rotation 30%, CB regrowth 30%
- Medium molecules (5-20 atoms): Translation 20%, Rotation 10%, CB regrowth 70%
- Long chains: Translation 10%, Rotation 5%, CB regrowth 85%
- Gibbs ensemble: Same as above but allocate 1-2% to CB swap and 1-2% to AVB1

---

## 6. Part 5: Parameter Selection — From Beginner to Confident

### How Many Molecules?

| System | Typical N | Reason |
|--------|-----------|--------|
| Simple fluid, NVT | 100-500 | Large enough for good statistics |
| Simple fluid, NpT | 200-500 | Need good density statistics |
| Gibbs (VLE) | 200-500 total | Split between boxes, need ~50 minimum in each |
| GCMC | 0-500 | Starts with some, then fluctuates |
| Protein + solvent | 1 protein + 1000-10000 waters | Depends on system size |
| Small molecule in box | 1-100 | Depends on target density |

**Rule of thumb:** Use at least ~200 molecules for bulk liquid properties. Make the box large enough that the minimum dimension is > 2x the cutoff radius.

### Cutoff Selection

```
rmin   1.0d0    ! Inner cutoff — atoms closer than this are not allowed (hard core)
rcut   14.0d0   ! Cutoff radius — interactions beyond this are ignored
rcutin 10.0d0   ! Inner cutoff for CBMC — atoms closer than this use full calculation
```

| System | rcut (A) | Why |
|--------|----------|-----|
| LJ fluid | 2.5-3.0 sigma | Beyond this, LJ forces are negligible |
| TraPPE-UA | 14.0 | Standard for this force field |
| Water (SPC/E) | 9.0-10.0 | With Ewald for electrostatics |
| Charged systems | 10.0-14.0 | Ewald handles long-range, but rcut still matters |

**Critical:** For charged systems, you MUST use Ewald summation (`electrostatic_form 'coulomb'`, `coulombstyle 'ewald_fixed_kmax'`). The cutoff alone is not sufficient for charges — they are too long-ranged.

### CBMC Parameters

Configurational-Bias Monte Carlo grows molecules atom-by-atom, choosing favorable conformations. The key parameters control how many trial positions are attempted:

```
cbmc_dihedral_generation    'ideal + autofit gaussian'
nch_tor     100             ! Number of trial dihedral angles
cbmc_bend_generation        'ideal + autofit gaussian'
nch_bend_a  100             ! Number of trial bending angles
cbmc_bond_generation        'autofit gaussian'
nch_vib     1               ! Number of trial bond lengths (1 = fixed bond)
```

**What these mean:** When growing a molecule, for each step Towhee randomly generates `nch_tor` trial dihedral angles, computes the energy for each, and selects one based on Boltzmann weighting. More trials = better quality but slower.

**Starting values:** `nch_tor 50-100`, `nch_bend_a 1-100`. For simple molecules with few conformations, 50 is fine. For long flexible chains, use 100+.

**The `cbmc_formulation` and `cbmc_setting_style` keywords:**

These control the detailed algorithm. The two most common combinations:

```
! TraPPE-UA (united-atom):
cbmc_formulation      'Martin and Frischknecht 2006'
cbmc_setting_style    'Martin and Frischknecht'

! OPLS, Amber, CHARMM (all-atom):
cbmc_formulation      'Martin and Siepmann 1999 + Martin and Thompson 2004'
cbmc_setting_style    'default ideal'
```

The difference is in how pre-nonbonded and 1-4 interactions are handled during the CBMC growth process. Using the wrong combination can give wrong results — match it to the force field.

### Chemical Potential (Grand Canonical)

The `chempot` parameter is the chemical potential you impose. Finding the right value is iterative:

1. **From ideal gas:** mu_ideal = kT * ln(rho * Lambda^3), where rho is number density, Lambda is thermal de Broglie wavelength. This gives a rough starting point.

2. **From a Gibbs simulation:** Run a Gibbs simulation and note the chemical potential from the output (`u (Gibbs Total)`). Use this value in a GCMC simulation at the same temperature.

3. **Trial and error:** Start with a guess. Run a short simulation. If the box fills up completely (N reaches nmolectyp), mu is too high. If N drops to zero, mu is too low. Adjust by ~500-1000 K units (in units of K, which is mu/kB).

### Equilibration vs Production

A proper simulation has two phases:

**Equilibration (discard these configurations):**
- The system starts from an artificial configuration (e.g., molecules on a lattice)
- Properties drift until they reach steady values
- This period must be discarded

**Production (keep these configurations):**
- Properties fluctuate around stable averages
- These are used for computing final results

**In the input file:**
```
controlstyle  'equilibration'    ! Use this for the first portion
! OR
controlstyle  'production'       ! For the actual data collection
```

Alternatively, just run equilibration manually and discard the first portion of the output (e.g., 10-20% of the total steps).

**How long to run?** There is no universal answer. Check that:
1. Block averages stabilize (see Part 6)
2. Running averages of key properties (density, energy) are flat
3. The standard deviation from block averaging decreases as 1/sqrt(N_blocks)

A minimum production run is typically 10,000-50,000 MC cycles for small molecules, and 50,000-1,000,000 for complex molecules.

### Initial Configuration

**From scratch (`linit .true.`):**
- Towhee places molecules on a lattice (simple cubic, FCC, etc.)
- Uses CBMC to grow the molecules with correct conformations
- The box size and lattice determine the initial density

**From a restart (`linit .false.`):**
- Reads `towhee_initial` (previously saved `towhee_backup` or `towhee_final`)
- Continues from where the previous simulation left off
- Max displacements and box dimensions are read from the restart file

**Always start from scratch for a new system.** Once equilibrated, use the restart for production runs.

---

## 7. Part 6: Reading the Output — A Complete Walkthrough

When Towhee runs, it prints information to standard output. This section walks through a representative output (from the `Canonical_Ensemble/TraPPE_Molecules` example). The output file is the file you would redirect to `towhee_output`.

### Section 1: Header and Input Echo (Lines 1-120)

```
MCCCS Towhee - Version 8.2.3   (July 22 2021)
Reading from towhee_input file: towhee_input
ensemble: nvt
temperature:  400.000
nmolty:            6
nmolectyp:            1            1            1            1            1            1
```

**What to check:** Verify all input parameters were read correctly. This is your first sanity check. If something was misformatted, it may show a default value here.

### Section 2: Force Field and Molecule Setup (Lines 120-410)

```
opening forcefield file:    1
Lorentz-Berthelot Mixing rules
    using the TraPPE-UA  force field
unit:    1 name:CH3*(sp3)
unit:    2 name:Ccccc(sp3)
...
Charges assigned for Molecule Type:    2
Unit:    1 nbname: CH3*(sp3)  Charge:    0.26500
Unit:    2 nbname: Och(sp3)   Charge:   -0.70000
Unit:    3 nbname: Ho         Charge:    0.43500
Total charge for Molecule Type:    2 is:    0.00000
```

**What to check:**
- Are the atom types being recognized from the force field file?
- Are charges assigned correctly? The total charge of each molecule should be zero (or the expected net charge).
- Are all bond, angle, and torsion types being found?

```
Bond Types
Type:      1 Style: Fixed  Length: 1.5400
Angle Types
Type:   1 Style: Standard Harmonic  Angle:    114.000 Constant:    31250.0
Torsion Types
Type:   1 Style: Old UA OPLS Cosine Series
          k0:       0.0 k1:     355.0 k2:     -68.2 k3:     791.3
```

**Check for "no types found" or zero parameters** — that would mean the force field file doesn't have parameters for your molecule's connectivity.

### Section 3: Nonbonded Parameter Table (Lines 333-404)

```
Nonbonded Force Field
Lennard-Jones 12-6 potential
Num. Atom(i)    Num. Atom(j)         sigma    epsilon
   3 CH3*(sp3)     3 CH3*(sp3)      3.7500    98.0000
   3 CH3*(sp3)     4 CH2**(sp3)     3.8500    67.1416
   4 CH2**(sp3)    4 CH2**(sp3)     3.9500    46.0000
```

**What to check:** This cross-referencing table tells you every pair's sigma and epsilon. Verify these against published force field values. The Lorentz-Berthelot mixing rules mean sigma_ij = (sigma_i + sigma_j)/2 and epsilon_ij = sqrt(epsilon_i * epsilon_j).

### Section 4: Initial Energies (Lines 449-494)

```
Initial Energies for Box     1
Total vibration (atoms ij bonded together)                     0.000 [K]
Total angle (atoms ijk bonded in sequence)                 11026.767 [K]
Total torsion (ijkl bonded in sequence)                    40778.674 [K]
Total nonbond (more than 3 bonds apart) vdw                28662.693 [K]
Total nonbond (more than 3 bonds apart) Electric          -33109.591 [K]
Total classical                                             47358.543 [K]
Aggregate van der waals (bond and nonbond)                 28662.693 [K]
Aggregate Electric (bond and nonbond)                       -197.863 [K]
initial virial pressure in box  1 =         512.17
```

**Energy units:** Energies are reported in **Kelvin** (K) and **kcal/mol**. The conversion is: E in K = E in kcal/mol * (kcal/J) / kB. 1 kcal/mol = 503.22 K at T=0 K, but Towhee uses E_K = E_kcal * 1000 / R where R is the gas constant, so 1 kcal/mol ≈ 503.22 K.

**Energy breakdown:**
- **Vibration** (bonds): Energy from bond stretching (typically 0 for rigid bonds — Towhee uses fixed bond lengths by default)
- **Angle**: Energy from angle bending
- **Torsion**: Energy from dihedral rotation + 1-4 nonbonded interactions
- **Nonbond vdw**: van der Waals between atoms separated by >3 bonds + intermolecular
- **Nonbond Electric**: Coulomb interactions (includes real space, self, correction, reciprocal for Ewald)
- **Total classical**: The sum of everything

**What the values tell you:**
- If any energy component is absurdly large (> 10^8), atoms are overlapping — the initial configuration is bad. This would happen if the box is too small or `rmin` is violated.
- For condensed phases at equilibrium, intermolecular vdW is large and negative, while intramolecular terms (angle, torsion) are positive.
- The virial pressure gives an initial estimate. It may be absurd initially (millions of kPa) — that's fine, it will relax during equilibration.

### Section 5: Per-Cycle Output (Lines 496-740, abbreviated)

```
+++++ start of markov chain +++++

Cycle      Box   Energy [K]  Volume [A^3] Press. [kPa] Molecules
      1 B: 1  0.4736E+05  0.6400E+05        512.2     1    1    1    1    1    1
      2 B: 1  0.3993E+05  0.6400E+05        512.2     1    1    1    1    1    1
```

**Each cycle line shows:**
- Cycle number
- Box number (for multi-box simulations)
- Total energy (K)
- Volume (A^3) — constant in NVT
- Pressure (kPa) — from virial
- Count of each molecule type

**What to watch for:**
- Energy should fluctuate around a stable average. If it drifts continuously, the system is not equilibrated.
- Large jumps followed by slow return may indicate rare events (good — you want to sample them).
- In NpT, the volume should fluctuate around an average. If it keeps growing or shrinking, you haven't equilibrated.

### Section 6: Block Averages (Lines 501-738, repeated every `blocksize` cycles)

```
Block Averages (BA) for block     1
BA Box: 1 Volume [A^3] 0.64000000E+05
BA Box: 1 V^2 [A^6] 0.40960000E+10
BA Box: 1 Specific density [g/ml] 0.14448461E-01
BA Box: 1 Virial Pressure      [kPa] 0.00000000E+00
BA Box: 1 Total Classical 0.46126843E+05
BA Box: 1 Inter vdw       -.26917970E+02
BA Box: 1 Angle           0.11051225E+05
BA Box: 1 Torsion         0.40439030E+05
BA Box: 1 Intra vdw       0.27784299E+05
BA Box: 1 Coulomb         -.33120792E+05
BA Box: 1 Tail vdw        -.11028461E+02
BA Box: 1 Number density [nm-3] Type   1 0.15625000E-01
BA Box: 1 Mol Fraction Type   1 0.16666667E+00
BA Box: 1 Stress Tensor Virial S_xx          [kPa] 0.00000000E+00
BA Box: 1 Radius of Gyration Type:  1        3.26317
```

**Block averages are your main diagnostic tool:**

- **Volume**: Should be stable. Standard deviation across blocks tells you the uncertainty.
- **Specific density** (g/ml): Use this to compare with experimental densities.
- **Total Classical**: The total potential energy.
- **Inter vdw**: Intermolecular van der Waals energy. Large and negative for condensed phases.
- **Intra vdw**: Intramolecular van der Waals (atoms in same molecule but >3 bonds apart).
- **Coulomb**: Electrostatic energy. For nonpolar molecules (like alkanes), this should be ~0.
- **Number density** (nm^-3): Molecules per cubic nanometer.
- **Mol Fraction**: Mole fraction of each component.
- **Radius of Gyration**: Average size of each molecule type (important for polymers).

**How to check convergence:** Compare block averages. After enough cycles, the block averages should vary only by their statistical uncertainty (not drift). The standard deviation of the block averages divided by sqrt(N_blocks) gives the standard error of your final average.

### Section 7: Move Statistics (Lines 748-791)

```
* Configurational-Bias REGROWTH Moves *
 Molecule type:    1 Box:    1
     Length  Attempts  Regrown  Accepted  %Regrown   %Accep.
       11         1         1         1    100.00    100.00
       15         1         1         0    100.00      0.00

* COM Translation Moves *
Molecule:  2 Box: 1 Attempts:              4. Accepted:              4. Accepted: 100.000 %

* Rotation Moves *
Molecule:  1 Box: 1 Attempts:              2. Accepted:              2. Accepted: 100.000 %
```

**These are your acceptance rates.** General guidelines:

| Move | Target Acceptance | If too high | If too low |
|------|-------------------|-------------|------------|
| Translation | 30-50% | Increase `rmtrac` | Decrease `rmtrac` |
| Rotation | 30-50% | Increase `rmrot` | Decrease `rmrot` |
| Volume | 20-40% | Increase `rmvol` | Decrease `rmvol` |
| CB regrowth | 20-50% | Complex — adjust CBMC trials | Increase trial counts |
| CB swap | 1-10% | Use AVB | Increase AVB radius |
| AVB insertion | 5-30% | (Good) | Increase `avb1rad` |

**Very low acceptance** for swap moves is the most common problem in Gibbs and GCMC simulations. This is why AVB exists — it can increase swap acceptance 10-100x.

**For CB regrowth, note the "Length" column:**
- "Length 5" means a 5-atom segment was regrown
- 0% acceptance for long regrowths is normal — regrowing part of the molecule is easier than all of it
- `pmall` controls the fraction of regrowth attempts that regrow the entire molecule vs partial regrowth

### Section 8: Final Energies (Lines 792-836)

Same format as initial energies. Compare them:
- Final energy should be lower (more negative) than initial — the system has relaxed
- The energy components should be physically reasonable (no huge numbers)
- The total should be stable if the system is equilibrated

### Section 9: Averages Table (Lines 837-886)

```
Averages               Units Type       Box  1
 Volume                 nm^3       0.64000E+02
 Specific Density       g/ml        0.01444846
 Number Density         nm-3    1      0.01562
 Mole Fraction                  1    0.1666667
 Radius of Gyration        A    1    3.2785054
 Ideal Pressure          kPa       0.51790E+03
 Total Classical           K        0.4124E+05
 Inter vdw                 K       -0.1708E+03
 Angle                     K        0.1120E+05
 ...
 U                    kJ/mol        0.3429E+03
 C_v: T^2<H^2>-<H>^2 J/K mol        0.3898E+03
```

This is the final summary of the entire simulation. Key quantities:
- **Specific Density**: The average mass density. This is what you would compare with experiment.
- **Number Density**: Molecules per nm^3. Useful for calculating molecular volumes.
- **Ideal Pressure**: From ideal gas law p = NkT/V. For dense liquids, the virial pressure will be very different (and correct).
- **U (kJ/mol)**: The average molar internal energy.
- **C_v**: Heat capacity at constant volume, from energy fluctuations.

### Section 10: Block Averages Summary (Lines 887-908+)

```
Block Averages (  5 blocks) Units Type Box Average      Standard Deviation
 Specific Density            g/ml        1  0.14448E-01  0.00000E+00
 Total Classical                K        1  0.41245E+05  0.24443E+04
 Inter vdw                      K        1 -0.17077E+03  0.81999E+02
 ...
```

**This is the most important section for error analysis.** The standard deviation across blocks, divided by sqrt(number of blocks), gives the standard error of the mean. For example:

```
Total Classical Average = 41245 K
Standard Deviation = 2444 K
Number of Blocks = 5
Standard Error = 2444 / sqrt(5) = 1093 K
```

**Rules of thumb:**
- If the standard deviation is > 50% of the average for properties like density, you need more sampling
- For energies, 5-20% relative uncertainty is typical
- The standard deviation should decrease as you run more blocks (it scales as 1/sqrt(N))

### Special Outputs for Gibbs Ensemble

The Gibbs output adds:

```
H_vap (pV/n=RT)      kJ/mol        0.2189E+02
H_vap (thermo p)     kJ/mol        0.1909E+02
```

**Heat of vaporization** computed two ways:
- `pV/n=RT`: Assumes ideal gas for vapor. H_vap = U_vap - U_liq + RT
- `thermo p`: Uses thermodynamic (dU/dV) pressure. More accurate if vapor is not ideal.

These should be close to each other. Large differences indicate the vapor is not ideal (too dense).

**Chemical potentials in Gibbs output:**

```
u (Density)               K    1        0.000   -5333.340
u (NVT Insertion)         K    1        0.000    3536.847
u (Gibbs Total)           K    1        0.000   -1793.727
```

At equilibrium, the total chemical potential should be equal in both boxes. If `u (Gibbs Total)` differs significantly between boxes, the simulation has not equilibrated.

### Special Outputs for Grand Canonical

The GCMC output mainly shows molecule counts fluctuating. The key quantity is the **average number of molecules at the imposed chemical potential**. You would run multiple simulations at different mu values to build an adsorption isotherm.

---

## 8. Part 7: Interpreting Results — What Does It All Mean?

### Checking for Equilibration

Before any analysis, you must confirm the system is equilibrated. Here is a systematic approach:

1. **Plot energy vs cycle number.** The first portion shows a trend (drift) as the system relaxes from the initial configuration. After equilibration, the energy should fluctuate around a stable mean with no systematic drift.

2. **For NpT simulations, plot volume vs cycle.** The volume should reach a stable average. An expanding or contracting volume means the pressure is not yet at the target.

3. **For Gibbs simulations, plot the number of molecules in each box.** The molecules should partition between boxes and fluctuate around stable averages. Equal chemical potential means equilibrium — check the `u (Gibbs Total)` values.

4. **Check block averages.** If block 1 and block 5 give noticeably different averages (beyond statistical fluctuations), the system is not equilibrated.

**How to handle non-equilibrated data:** Discard it. Only use data from after the equilibration point. You can split the simulation into an equilibration phase (manually set `nstep` to cover just equilibration, or run separate equilibrate/production runs).

### What Is a "Good" Result?

Compare with experimental data or literature simulation values:

| Quantity | Typical Accuracy (with good FF) |
|----------|-------------------------------|
| Liquid density | 1-3% |
| Vapor pressure | 5-20% |
| Heat of vaporization | 5-15% |
| Chemical potential | Within kT of reference |
| Radial distribution | Peak positions within 0.1 A |
| Radius of gyration | 5-10% (for polymers) |

**If your results are far off:**
- Check that the force field is appropriate for your system
- Verify that the molecule definition (atom types, connectivity) matches the force field
- Ensure you ran long enough for proper sampling
- Check that electrostatics are handled correctly (Ewald parameters)
- Compare with published benchmark values for the same force field

### Common Red Flags

| What You See | Likely Cause | Fix |
|---|---|---|
| Energy = 10^8+ K | Atoms overlapping | Increase rmin, check initial box size, use smaller initial molecules |
| 0% acceptance for all moves | Box too dense or rmin too large | Reduce system size, adjust rmin |
| N drops to 0 in GCMC | Chemical potential too low | Increase chempot |
| N hits nmolectyp in GCMC | Chemical potential too high | Decrease chempot |
| Box 2 (vapor) in Gibbs has 0 molecules | mu too high in liquid | Redistribute molecules, adjust chempot |
| Gibbs chemical potentials differ | Not equilibrated | Run longer |
| Volume keeps shrinking in NpT | Pressure too low or initial box too large | Check pressure units (kPa!), adjust initial box |
| Volume keeps expanding | Pressure too high | Same fix |
| Acceptance for swaps = 0% | No AVB, or avb1rad too small | Enable AVB, increase radius |
| "No torsion types found" | Connectivity not matching FF file | Check atom type names and vibration connectivity |
| Block averages not converging | Not enough sampling | Run longer, increase blocks |
| Density wildly wrong | Wrong force field or molecule definition | Double-check atom types and FF choice |

---

## 9. Part 8: Worked Examples from the Codebase

### Example 1: NVT of TraPPE Pentane Isomers (`Canonical_Ensemble/TraPPE_Isomers`)

**Goal:** Simulate 5 different pentane isomers at 300 K in the NVT ensemble.

**Key decisions:**
- *Why NVT?* We are testing structure, not phase behavior. We want conformations at a known temperature.
- *Why TraPPE-UA?* United-atom is fast and accurate for alkanes.
- *Why linit = .true.?* Starting from scratch.
- *Why simple cubic lattice?* Simple and works for low-density systems.

**How to interpret the output:**
- Compare the radius of gyration of each isomer — linear pentane (n-pentane) should have the largest Rg, while highly branched isomers (neopentane) should have the smallest.
- The intramolecular energy distribution tells you about conformational preferences (gauche vs trans).

### Example 2: Gibbs Ensemble of TraPPE Pentane (`Gibbs_Ensemble/TraPPE_Pentane`)

**Goal:** Compute vapor-liquid coexistence of n-pentane at 336 K (close to the normal boiling point of 309 K).

**Key decisions:**
- *Why Gibbs?* Need vapor-liquid coexistence directly.
- *Why 300 molecules?* Standard for Gibbs — enough to get good statistics in both phases.
- *Why 260 and 40 initial split?* Box 1 (liquid) needs many molecules to form a dense phase; Box 2 (vapor) needs fewer — it will be mostly empty space.
- *Why init box sizes 30 A and 60 A?* The liquid box is smaller (dense), the vapor box is larger (gas phase needs space).
- *Why AVB with radius 8 A?* Typical insertion radius. Targets the free volume in the liquid.
- *Why linit = .false.?* This example is set up to restart from a pre-equilibrated file. For a fresh start, set it to .true.

**How to interpret the output:**
- Box 1 (liquid) density should be ~0.58 g/ml (close to experimental ~0.61 g/ml at 336 K).
- Box 2 (vapor) density should be much lower (~0.01 g/ml).
- The heat of vaporization should be ~22 kJ/mol (experimental ~24 kJ/mol at boiling point).
- The chemical potentials in both boxes should be equal at equilibrium.

### Example 3: Grand Canonical of Ethane in Amber96 (`Grand_Canonical_Ensemble/Amber96_Ethane`)

**Goal:** Determine how many ethane molecules adsorb at a given chemical potential and temperature.

**Key decisions:**
- *Why Amber96?* Explicit hydrogen — all-atom model for ethane.
- *Why charge_assignment 'manual'?* Amber uses explicit partial charges, not bond increments.
- *Why nstep = 100 moves (not cycles)?* Moved-based stepping gives finer control.
- *Why pmuvtcbswap = 0.25?* 25% of moves attempt insertions/deletions.
- *Why louthist = T?* Output histogram data for post-processing.

**Charge assignment in this example — the manual method:**
```
unit ntype qqatom
1    'CT'   -0.1926d0
vibration
4
2 3 4 5
...
```
Each carbon (CT) gets -0.1926 e, each hydrogen (HC) gets +0.0642 e. The total sums to zero (2 * -0.1926 + 6 * 0.0642 = 0). This is essential — the simulation will check and warn if total charge is non-zero.

### Example 4: NpT of LJ Solid (`Isobaric_Isothermal_Ensemble/Solid_LJium`)

**Goal:** Simulate solid LJium at very low temperature (T=0.6) and high pressure (P=100 in LJ units).

**Key decisions:**
- *Why lshift = .true., ltailc = .false.?* For solids with short cutoffs, shifting the potential avoids discontinuities. Tail corrections are inappropriate here because the system is not isotropic.
- *Why rcut = 2.5?* Standard for LJ systems — the potential is negligible beyond this.
- *Why pmcell = 0.005 (not pmvol)?* Crystals can change shape non-isotropically. Isotropic volume changes miss the correct crystal structure.
- *Why pmplane/pmrow?* Slab/cylinder moves for solid rearrangements — more efficient than COM translations in a dense crystal.
- *Why only 216 atoms (6x6x6)?* Small enough to be fast, large enough for bulk properties.

**What to look for:** The cell dimensions should converge to the equilibrium crystal structure at T=0.6, P=100. The energy should be very stable (negligible fluctuations at low T).

### Example 5: TMMC of SPC/E Water (`TMMC/SPC_E_WATER`)

**Goal:** Use Transition Matrix Monte Carlo to compute the phase behavior of water at 500 K.

**Key decisions:**
- *Why TMMC instead of Gibbs?* TMMC gives the entire density of states, from which you can compute properties at any chemical potential. Also handles the free energy barriers between liquid and vapor.
- *Why random_allow_restart = .false.?* TMMC needs reproducible random numbers for the transition matrix.
- *Why RANLUX random number generator?* High-quality RNG for precise work. The default DX-1597 is also fine.
- *Why initmol = 0?* Start with zero molecules and let them be inserted by GCMC swaps.
- *Why pmuvtcbswap = 0.65?* TMMC needs many insertion/deletion attempts to build the transition matrix.
- *Why detailed CBMC parameters?* Water is a small rigid molecule, but the CBMC parameters for dihedral/bend/bond generation are explicitly set because the ideal + autofit Gaussian method needs them.

**What to look for:** The output includes `towhee_histogram` files with the transition matrix and histogram data. These are processed by `analyse_histogram` (in `Utils/`) to compute the phase diagram.

---

## 10. Part 9: End-to-End Workflow

### Complete Workflow for a New Simulation

#### Step 1: Define Your Goal

Ask yourself:
- What property do I want? (density, phase diagram, adsorption, structure)
- What molecule(s) am I studying?
- What are the experimental conditions? (T, P, composition)

#### Step 2: Choose the Ensemble and Force Field

Use the decision tree from Part 2 and the force field guide from Part 3.

#### Step 3: Determine System Size

- Minimum 200 molecules for bulk properties
- Box size > 2x rcut in each dimension
- For Gibbs: enough molecules so both boxes have > 50 after equilibration
- For GCMC: start with 0 molecules or a small number

#### Step 4: Write the Input File

Start from a similar example in `Examples/`. Modify:
1. Temperature, pressure, chemical potential
2. Molecule types and counts
3. Force field path
4. Move probabilities
5. Output frequencies
6. Initial box size (set density close to expected value)

#### Step 5: Run a Short Test

Start with `nstep = 100-1000` (or 10 cycles). Check:
- No errors during setup
- Energies are reasonable
- Move statistics show some acceptance
- Initial configuration was built successfully
- Charges are balanced

#### Step 6: Equilibrate

Set `linit .true.` and run 10,000-50,000 cycles (depending on system complexity). The goal is to reach a stable state where:
- Energy fluctuates around a constant mean
- Volume (for NpT) or molecule counts (for GCMC/Gibbs) are stable
- Chemical potential (for Gibbs) is equal in both boxes

#### Step 7: Production Run

From the equilibrated system:
- Copy `towhee_backup` to `towhee_initial`
- Set `linit .false.`
- Set `nstep` to your production length (at least 5-10x the equilibration length)
- Choose `blocksize` so you get 5-20 blocks for error analysis

#### Step 8: Analyze Output

- Read the averages and block averages sections
- Compute mean and standard error from block averages
- For Gibbs: extract liquid and vapor densities, heat of vaporization
- For GCMC: average number of molecules at this chemical potential
- For NpT: density at the given T and P

#### Step 9: Validate

- Compare to experimental data
- Compare to literature simulation values (same force field)
- Check that results are converged (longer run gives same answer within error)
- If results are poor, revisit your choices: force field, system size, equilibration time

### Quick Reference: What to Do When Results Are Wrong

| Problem | Check First |
|---------|-------------|
| Density too high | Box too small? rmin too small? Wrong force field? |
| Density too low | Box too large? ltailc off for large rcut? |
| Gibbs VLE curve shifted | Check for finite-size effects (too few molecules). Use tail corrections. |
| No swap acceptance | Enable AVB. Increase avb1rad. Check pm2boxcbswap > 0. |
| Energies huge (10^6+) | Atoms are overlapping. Check initial box size and molecule placement. |
| Non-zero total charge | Check manual charge assignments. Bond increments should sum to zero. |
| Different results on different CPUs | This is normal — you're sampling a distribution. Run longer to converge. |

---

## Appendix: File Format Quick Reference

### `towhee_input` Organization

```
! 1. Simulation control
ensemble 'nvt'
temperature 300.0
nmolty 1
nmolectyp 100
numboxes 1
stepstyle 'cycles'
nstep 10000

! 2. Output control
printfreq 1000
blocksize 5000
moviefreq 10000
backupfreq 5000

! 3. Initialization
linit .true.
initboxtype 'dimensions'
initstyle 'full cbmc'
initlattice 'simple cubic'
initmol 100
inix iniy iniz 5 5 4
hmatrix
30.0 0.0 0.0
0.0 30.0 0.0
0.0 0.0 30.0

! 4. MC moves
pmvol     0.0
pmcb      0.40
  pmcbmt  1.0
  pmall   0.0
pmtracm   0.70
  pmtcmt  1.0
  rmtrac  0.5
  tatrac  0.5
pmrotate  1.0
  pmromt  1.0
  rmrot   0.05
  tarot   0.5

! 5. Force field
ffnumber 1
ff_filename /path/to/towhee_ff_TraPPE-UA
classical_potential 'Lennard-Jones'
classical_mixrule 'Lorentz-Berthelot'
lshift .false.
ltailc .true.
rmin 1.0
rcut 14.0
rcutin 10.0
electrostatic_form 'none'

! 6. Molecule definition
! One block per molecule type
input_style 'basic connectivity map'
nunit 5
nmaxcbmc 5
lpdbnames F
forcefield 'TraPPE-UA'
charge_assignment 'bond increment'
unit ntype   1  'CH3*(sp3)'
vibration 1
2
improper torsion 0
unit ntype   2  'CH2**(sp3)'
vibration 2
1 3
improper torsion 0
! ... continue for atoms 3, 4, 5
```

### Common Mistakes

- **Forgetting cumulative probabilities:** The `pm*` values must increase monotonically, ending at 1.0.
- **Using the wrong units:** Temperature in K, pressure in kPa (not atm or bar), distances in A.
- **Wrong force field path:** The `ff_filename` path must point to the actual file on your system.
- **Electrostatics without Ewald:** For charged molecules, you MUST use Ewald summation or a proper long-range method.
- **Charges not balanced:** The total charge of each molecule type and the entire system must be zero (for neutral systems).
- **Not discarding equilibration:** The first portion of the simulation depends on the artificial initial configuration and must be discarded.
- **Too few molecules:** For bulk properties, < 100 molecules introduces significant finite-size effects.
- **Copying blind from examples:** The examples use `/towheebase/` paths — change these to match your installation.

### Further Reading

- The manual in `Manual/projects/towhee/` has detailed documentation for every input keyword
- Allen & Tildesley, *Computer Simulation of Liquids* — the classic textbook
- Frenkel & Smit, *Understanding Molecular Simulation* — excellent for MC methods
- Each force field has its own published paper with parameter tables; look in `towhee_citations` output for suggested references

---

## Colophon

This guide was produced by analyzing the Towhee 8.2.3 codebase at `src/towhee-8.2.3/`, including all 60+ example directories across 6 ensemble types, the full algorithm documentation, force field files, and standard output examples. The examples used for reference span:

- **Canonical (NVT):** 22 simulations (TraPPE, Amber, CHARMM, COMPASS, DREIDING, UFF, SPC/E, etc.)
- **Gibbs:** 16 simulations (alkanes, CO2, water, metals, zeolites, etc.)
- **Grand Canonical (muVT):** 4 simulations (adsorption, energy biasing, walls, nanotubes)
- **Isobaric-Isothermal (NpT):** 15 simulations (liquids, solids, alloys, explosives, Henry's law)
- **TMMC:** 2 simulations (LJ, SPC/E water)
- **Histogram/VLCC/Convert/DFT:** Specialized utilities
