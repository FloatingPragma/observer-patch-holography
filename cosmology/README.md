# OPH Cosmology And Dark Gravity

This directory contains publicly released research companions, not journal-published, on the cosmological continuation of Observer Patch Holography. They form layer 7 of the one OPH tower, cosmology and tests: the same executable process that must realize observers, clocks, fields, matter and gravity in layers 0 to 6 must also fix the vacuum, the dark sector, the primordial response and black-hole records, and pass independent prospective tests.

Each paper opens with the shared statement of the [executable-model hypothesis](../docs/OPH_EXECUTABLE_HYPOTHESIS.md) and the tower, followed by its own hypothesis, its finite theorems, the joins the executable model must realize before its conclusions apply, and the result that would refute the branch. Supplied spatial domains, FLRW backgrounds, species content, record-production profiles, radial lifts and clock calibrations are premises of the corresponding branch. Their use establishes no preferred clock of the underlying observer network. Modular parameters, record precedence and calibrated local duration have distinct roles; see the [working terms](../docs/INTERPRETATION_AND_TERMS.md).

The executable model supplies native repair, records, local event clocks and a two-dimensional screen continuum. It supplies no covariant cosmological source, physical scale, radial lift, expanding background, photon sector or quantum vacuum, so every cosmological conclusion here is a conditional implication until those joins are made on the same process.

## Paper Map

| Paper | Hypothesis and place in the tower |
| --- | --- |
| [Observer Patch Holography and the Dark Sector](oph_dark_matter_paper.pdf) ([source](oph_dark_matter_paper.tex)) | Dark matter is non-luminous modular charge carried on unrecovered overlap collars. Machine-checked classification, saturation and profile theorems; under deep-regime scale covariance and quadrature composition, the deep limit of Milgrom's modified dynamics with a fitted `a_0`, whose Tully–Fisher exponent 4 sits 2.66 standard deviations above the seen sample's fit. Refuted by native repair that recovers every collar or by a calibrated galaxy likelihood excluding the deep law |
| [Observer-Screen Synchronization as a Candidate Inflation-Free Branch](oph_inflation_without_inflaton_observer_screen_synchronization.pdf) ([source](oph_inflation_without_inflaton_observer_screen_synchronization.tex)) | Observer-screen synchronization does the jobs of inflation. Record-density expansion is exact on a supplied population; the edge-center tilt `n_s = 1 - P_star/48` is conditional on a declared generator density; flatness is identified, not selected; horizon coherence is unsupplied; under a bounded form-congruent identification the native density field forces zero tilt |
| [Finite-Source Cosmology](oph_cosmology_finite_source_cmb_program.pdf) ([source](oph_cosmology_finite_source_cmb_program.tex)) | The primordial sky is the image of one finite observer-screen source map. Conditional source-screen spectrum, radial uniqueness routes, exact golden-population Fourier control, and the source, transfer and likelihood conditions the executable model must meet |
| [Cosmological Vacuum And Structure Formation](oph_cosmological_vacuum_and_structure_formation.pdf) ([source](oph_cosmological_vacuum_and_structure_formation.tex)) | The vacuum is the refinement limit of a source-derived finite quotient ensemble, and the vacuum term is a screen capacity. The native density law is reflection positive but its generator is `A`, not `sqrt(A)`, so no quantum vacuum is identified; the fixed-capacity equation of state is in tension with DESI DR2 |
| [Cosmology Data And Likelihood Contracts](oph_cosmology_data_likelihood_contracts.pdf) ([source](oph_cosmology_data_likelihood_contracts.tex)) | The test interface of layer 7: comparison classes, model and data contracts, the no-data-use rule, and the retrospective ledger of every receipt-conditional value on seen data |
| [Standard Photon–Baryon Transport as the Transfer Interface for OPH Sources](oph_boltzmann_transport_derivation.pdf) ([source](oph_boltzmann_transport_derivation.tex)) | A conditional reduction: on a branch that supplies geometry, photon convergence, vertex and collision premises, transport is the Ma–Bertschinger system plus two OPH interfaces (repair-exchange rate and collar stress) that vanish on the standard branch |
| [A Finite Black-Hole Information Ledger](oph_black_hole_information_ledger.pdf) ([source](oph_black_hole_information_ledger.tex)) | Black-hole information as a finite record ledger: machine-checked unitary dynamics and sign-unitary record channels, known finite entropy theorems, stated conjectures for thermality and relocation, and the declared crossing-read area law |

The [formal radial-lift theorem fragment](../paper/tex_fragments/RADIAL_LIFT_THEOREMS_330.tex)
contains the one-shell non-identifiability proof, physical source-dilation
theorem, radial cross-covariance tomography, exact amplitude conversion, and
finite-window bound. The
[simulator contract](../code/cosmology/SIMULATOR_RADIAL_CONTRACT_330.md)
specifies the fail-closed evidence split used by finite runs.

The [physical CMB theorem program](physical_cmb_theorem_program.md) specifies the supplied and unsupplied physical-event, order-faithful-placement, count--volume, continuum, source, lift, stress, abundance, transfer, and likelihood inputs in one place.

## Expansion, Area Law, and Public Data

On the golden source-record family of the spacetime construction, the expansion history is a record density: the comoving read law generates the same finite causal order for every scale-factor profile, the scale factor enters only as the count weight `a^4` per comoving cell and conformal tick, the count clock converges to the proper-time ratio for shrinking diamonds, and the redshift is the fourth root of a record-count ratio for congruent comoving windows with constant scale factor within each (or their shrinking-window limit) ([fragment](../paper/tex_fragments/SOURCE_NET_FLRW_RECORD_DENSITY.tex), [Lean](../Lean/Geometry/SourceNetConformalRecordDensity.lean)). The number of unordered read pairs crossing a surface per layer converges to an area law, `(pi/4) n^2 a^4` per unit area, with the Planck-area dictionary `n^2 a^4 = 1/(pi l_P^2)`; this count is a declared one-nat-per-read quantity on a supplied three-dimensional population, distinct from native regional entropy, which carries no boundary-seam coefficient on the two-dimensional source screens ([entropy audit](../code/geometry/ENTROPY_FIRST_LAW_AUDIT.md)) ([fragment](../paper/tex_fragments/CROSSING_READ_AREA_LAW.tex), [Lean](../Lean/Geometry/CrossingReadAreaLaw.lean), [receipt](../code/causal_refinement/crossing_read_area_law_receipt.json)).

The conditional primordial, background, capacity and dark-sector values are compared with public measurements in the [cosmology postdiction ledger](../code/cosmology/postdiction_ledger/COSMOLOGY_POSTDICTION_LEDGER.md), produced deterministically and checked by an independent verifier, with its brackets machine-checked in [Lean](../Lean/ObserverPatchHolography/EinsteinBranch/CosmologyLedgerBrackets.lean); every row is on seen data and carries its class. The tilt FZ-15 and fixed-capacity FZ-13 protocols were frozen on 2026-10-04, with status `frozen_stamped_upgrade_pending` in the [register](../claims/frozen_prediction_register.json); their complete protocols and anchoring requirements govern eligible comparisons. The [primordial registration proposal](../code/cosmology/primordial_freeze/PRIMORDIAL_REGISTRATION_PROPOSAL.md) holds the unadopted tensor FZ-16 and running/isocurvature FZ-17 proposals. The [operative FZ-15 target](../evidence/custody/falsification/frozen_targets/fz15_2026-10-04/target.md) specifies the tilt comparison.

## Dark-Gravity Structure

Under the collar-stress attachment, recovery defects bound a candidate dark source. Conserved comoving charge gives the familiar inverse-volume dilution. Two explicit scaling and composition premises select a spherical deep-acceleration law with baryonic Tully–Fisher scaling: the deep limit of Milgrom's modified dynamics, with `a_0` fitted.

A physical interpretation must link this source to relativistic stress, abundance, lensing, cosmic expansion and calibrated readouts on the same process. The dark-sector paper states the hypotheses and compares the resulting formulas with galaxy data.

## Reproducibility

Companion radial-lift and clock-certificate code lives in
[`../code/cosmology/`](../code/cosmology/). The larger simulator and
visualization surfaces are:

- [OPH physics simulator](https://github.com/muellerberndt/oph-physics-sim)
- [Interactive simulation](https://simulation.floatingpragma.io)

The cosmology papers are focused research companions and are not part of the core release bundle.
