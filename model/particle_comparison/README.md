# Conditional particle assignments and actual-data comparison

This packet tests a supplied mass attachment to the current interacting
two-site full-Fock field. It also reruns the fixed spectral-alpha/ACT
comparison against actual observations. All targets are already known:
the work is **retrospective**, including the two named candidate assignments.
It supplies no eligible global physical dictionary and adds zero CP11 physical
comparisons. No result confirms that our universe is an OPH universe.

## Exact finite spectrum

For the open `(2,1,1)` geometry, coordinate speed `1/8` and interaction `U=1`,
the hopping magnitude is `t=1/16`. `compare.certify_source()` constructs the
current full-Fock Hamiltonian, requires every entry to be exactly a Gaussian
integer divided by 16, and checks all five number-sector characteristic
polynomials with exact arithmetic. All 16 states are retained. This binds the
following spectrum to the executable action, independently of numerical pole
clustering or visibility thresholds:

| Number | Eigenvalues, with multiplicity |
| --- | --- |
| 0 | 0 |
| 1 | −1/16 (2), +1/16 (2) |
| 2 | 0 (3), 1, 1/2−√17/8, 1/2+√17/8 |
| 3 | 15/16 (2), 17/16 (2) |
| 4 | 2 |

Charge-one transition frequencies lie in the eight-element signed superset
`{7/16−√17/8, −1/16, 9/16−√17/8, 1/16, 15/16,
7/16+√17/8, 17/16, 9/16+√17/8}`. Some entries may be dark for a
particular probe. Including every entry makes an exclusion conservative.
A Hamiltonian scalar offset cancels from transitions; a common chemical
potential shifts each addition frequency by `−μ`. A removal frequency is
its negative, so absolute-energy attachment also covers removal probes.
The parameter conjugate
to these frequencies is Hamiltonian model coordinate, with `ℏ_model=1`.
No physical proper clock or MeV scale is inferred from it.

## Two declared real-particle candidate assignments

Choose the supplied half-filled reference: the unique `N=2` ground of
`K=H−N/2`. It is distinct from the finite-temperature `β=1, μ=0` Gibbs
reference. The lowest hole and particle excitations are doublets:

| Conditional assignment | Relative number charge | Exact gap |
| --- | --- | --- |
| Electron-like hole, `N=1` | −1 | `Δ=√17/8−1/16≈0.4528882032` |
| Positron-like particle, `N=3` | +1 | the same `Δ` |

Assign one common scale `κ=(m_e c²)/Δ` using the
[CODATA 2022 electron mass energy](https://physics.nist.gov/cuu/pdf/wall_2022.pdf),
`m_e c²=0.51099895069(16) MeV`. Both conditional mass energies then equal
`0.51099895069 MeV`; their calibration errors are fully correlated.
**Neither absolute mass is independently predicted.** The one remaining
dimensionless prediction is zero normalized particle–antiparticle mass
difference. It satisfies the retrospective
[PDG 2025 bound](https://pdgweb.lbl.gov/2025/listings/rpp2025-list-electron.pdf)
of `8×10⁻⁹` at 90% confidence, obtained from
[Fee et al.'s positronium spectroscopy](https://doi.org/10.1103/PhysRevA.48.192)
assuming the positronium Rydberg constant is exactly half the hydrogen one.
This is a conditional equality comparison, not a g-factor comparison,
independent positron absolute-mass measurement, or unique OPH signature.
Admission checks enforce the mass observable, MeV calibration units,
dimensionless difference, confidence level, source citations, Rydberg
assumption and calibration/comparison roles. Duplicate JSON keys,
nonfinite values and undeclared conversions reject before comparison.

The physical vacuum, electric-charge dictionary, Lorentz spin, continuum
dispersion, stable physical pole and same-excitation defect/transport/clock
map remain unestablished. Doublet multiplicity alone does not give physical
spin. The internal particle–hole symmetry supplies a useful candidate relation;
it does not identify the excitations with observed electrons and positrons.


## Executable unit and operator dictionary

[`dictionary.py`](dictionary.py) supplies the complete sixteen-state benchmark
map, retaining the same reference and candidate projectors:

| Model quantity | Assigned physical quantity |
| --- | --- |
| Excess number | `Q=e(N−2I)` |
| Internal spin | `S_phys=hbar*S_internal` |
| Ground-relative generator | `H_phys=kappa*(K−E0 I)` |
| Hamiltonian model coordinate `s` | `t_phys=(hbar/kappa)*s` |
| Model lattice coordinate `x` | `x_phys=(c_phys/c_model)*(hbar/kappa)*x` |
| Oriented number current `J` | `I_phys=(e/Tunit)*J` |

The exact SI definitions of `h`, `e` and `c` come from
[BIPM](https://www.bipm.org/en/measurement-units/si-defining-constants).
With the single electron-energy calibration, one model time unit is
`5.83360161713e−22 s` and one lattice unit is `1.39909581423e−12 m`.
The assigned lattice spacing is `3.6231056256` reduced Compton lengths.
This finite benchmark therefore does not resolve that length scale.
The code checks the complete charge/spin/symmetry maps, Heisenberg evolution,
local continuity and Peierls current units. Exact arithmetic enclosures and
fully correlated calibration uncertainty are retained separately from
unquantified physical-dictionary error. The assigned Hamiltonian clock has
no derived intertwiner with native incident or modular time.

## Electromagnetic applicability checks

[`electromagnetism.py`](electromagnetism.py) lifts the existing Peierls law
to all sixteen Fock states and tests the two actual charged projectors.
Both have nonzero electric dipole and current response, with complete
transition weights `1/4` and `1/256`, and first dipole moment `1/32` in model
units. The full operator continuity identity, charge-conjugation signs and
spectral double-commutator sum rule hold. These are meaningful external-field
charge responses.

The open bond has no loop flux. Its phase is removable by a gauge
transformation, so its spectrum cannot measure a magnetic moment or `g`.
The actual bare onsite interaction has constant Fourier kernel; it does not
have the momentum dependence of static Coulomb exchange. This distinguishes
the supplied contact action from a directly identified Maxwell interaction;
it does not exclude an emergent photon in another limit. A dynamical photon,
annihilation photons and physical QED identification are unestablished.


## Dispersion and finite-volume controls

The open dimer has no translation momentum. The diagnostic therefore declares
a separate periodic four-site ring with the same local Weyl/contact action,
retaining all 256 Fock states. It checks the signed fermionic translation,
number blocks, unique half-filled ground, complete Fourier-probe residues,
particle–hole relation and the free symbol with both massless nodes.

The lowest charged gap at `k=0,pi` is `0.4175457854`; at `k=±pi/2` the two
spin-resolved probes have lowest visible gaps `0.4723570583` and
`0.5944851107`, with substantial residues. These disagree with
`E(k)^2=m^2+c_model^2*sin(k)^2` under the declared dimer mass calibration.
The ring/dimer rest-gap ratio `0.9219621585` is a finite-volume/reference
change; even using the ring's own rest gap does not repair the bare-speed
benchmark. Independent entrywise Jordan–Wigner reconstruction reproduces the
numbers, periodic signs and hopping normalization.

Fitting only the lowest envelope gives `c_eff=0.2208544942`, versus supplied
`c_model=0.125`. There are only two distinct `sin(k)^2` values, so fitting an
intercept and slope cannot test the functional form or recover a physical
velocity. Renormalized Lorentz behavior in another limit is not excluded.
The result rejects the stated finite dictionary benchmark; it establishes
neither electrons/positrons nor a universal OPH no-go. These diagnostics do
not add source execution on the ring or a physical photon clock.

## Exhaustive three-lepton exclusion, including energy-origin freedom

Test `m_i=κ|ω_i−μ|` with arbitrary **one common** `κ>0` and real `μ`,
for all assignments to the already declared broad ratio windows
`m_μ/m_e∈[206,207]`, `m_τ/m_e∈[3470,3485]` in
[`chiral_candidate`](../common_source/chiral_candidate/README.md).
These are broad empirical diagnostics, not a joint confidence region.
There are `8×7×6=336` ordered distinct frequency assignments and two
electron signs, hence **672 exact cases, none feasible**. Equal-frequency
assignments cannot fit distinct target masses. No probe or assignment is
selected after testing and retained as a success.

The reduction is elementary. For frequencies `a,b,c` assigned respectively
to electron, muon and tau, set `x=m_μ/m_e`, `y=m_τ/m_e`,
`q=(c−a)/(b−a)` and let `s_e=±1`. Because `x,y>1` and `κ>0`, the
other two signs must be `sign(b−a)` and `sign(c−a)`. Eliminating the
scale and chemical potential yields

`y=|q|x+sign(c−a)s_e(1−q)`.

Its slope is positive, so each attainable interval follows exactly from
the two rational endpoints of `x`. Every interval is either below the tau
window (largest upper endpoint approximately `2163.0075`) or above it
(smallest lower endpoint approximately `11655.6566`). Computation and
ordering take place in `Q(√17)`; outward rational enclosures only format
results. Independent linear programs retain scale and offset as unknowns
and check every case. A synthetic passing case verifies that chemical
potential freedom has not accidentally been suppressed.

This rejects the bounded three-generation attachment even with a common
chemical potential. It does not reject OPH, a different volume/action,
other observables, sector-dependent energy shifts or a nonlinear physical
dictionary. At fixed `μ=0`, the seven nonzero absolute frequencies already
have dynamic range only `22.8794`, far below the electron–muon ratio.

## Actual-data replay and eligibility

`cross_scale.build()` recomputes the fixed spectral expression and ACT
posterior weights, rather than trusting its retained receipt. The inverse
coupling is `137.035999189566…`, `0.598388…` reported CODATA standard
uncertainties above its observed value. The declared primordial map lies
at lower-tail posterior weight `3114/508985≈0.006118` in the retained
P-ACT-LB projection. This marginal posterior weight is **not a p-value**
or a joint model likelihood. The standard cosmology and calibration models,
missing physical photon-residue map and missing common error contract remain
explicit in [`empirical_diagnostics`](../empirical_diagnostics/README.md).
Earlier oscillator and chiral histories are not reexecuted by this replay.

The source prediction is hashed before observations are reopened. Because
the targets were already exposed, this is a reproducibility boundary, not
a prospective validation claim. CP11 remains open until CP07/CP10 supply
one eligible physical dictionary and a complete, independent comparison.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest model.particle_comparison.test_compare
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_comparison.compare --check model/particle_comparison/results.json
```

`results.json` binds the exact certificate, all-case digest, calibration,
observations, actual-data replay and current implementation/documentation
hashes. To regenerate it after a reviewed change, remove only that generated
receipt and run the same command with `--write`.
