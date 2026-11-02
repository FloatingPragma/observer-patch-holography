# Spectrum, refinement and geometric matter inputs — 2 November 2026

This campaign advances five existing completion issues. The admitted inputs
include the simulator, observer, scheduler, repair and memory specifications.
Every result below identifies its own state, algebra, operation, clock and
refinement. A changed clock or coupling law is derived explicitly, with its
own reference; finite source energies are not renamed physical masses.

| Issue | Why this work matters | Result and scope |
| --- | --- | --- |
| [#1138 — reference and access](https://github.com/FloatingPragma/observer-patch-holography/issues/1138) | A continuing detector needs a prepared reference and controlled time estimates. | FA84 gives explicit fixed-selector trace-norm convergence to a faithful number-sector reference. FA85 resolves the distinct centered ground's local spectrum. Neither supplies a general adaptive joint stationary state or a centered native detector. |
| [#1139 — protected particles](https://github.com/FloatingPragma/observer-patch-holography/issues/1139) | Protection must survive the actual complete repair records. | FA84 derives the sharp native fixed-number coherence factor `5/6`, including retained classical histories and untouched quantum references. Coherent memories and active recovery require their specified additional laws. |
| [#1142 — interacting continuum](https://github.com/FloatingPragma/observer-patch-holography/issues/1142) | The previous bounded-coupling continuum is free. | FA86 constructs controlled interacting spin dynamics at fixed finite graph, and proves the charged continuum remains free throughout `U=o((log c)^(3/2))`. It explicitly exposes the missing volume-uniform estimate. |
| [#1143 — measured species](https://github.com/FloatingPragma/observer-patch-holography/issues/1143) | A full spectral measure is needed before any species dictionary can be tested. | FA85 certifies all 228 local positive charged lines, weights and spatial residue matrices of the centered cube, retaining the unobserved many-body complement. This is a source-clock energy readout, not a physical species assignment. |
| [#1145 — Einstein inputs](https://github.com/FloatingPragma/observer-patch-holography/issues/1145) | The existing geometry theorem needs matter and modular data on the same state. | FA87 supplies the same Weyl vacuum's conserved stress, geometric ball modular generator, finite relative first law and small-ball remainder. Area normalization, entropy stationarity and responsive geometry remain open. |

## What is now calculated

The [repair proof](../repair_coherence_lifetime/README.md) checks 983,040 native
source columns and all 531,441 ternary occupation differences. For a connected
weighted graph with normalized positive seam weights and fair endpoint coins,
the minimum cut `lambda` gives fixed-number coherence decay `(1-lambda)^n`.
For the twelve-port uniform selector, `lambda=1/6`. The diagonal population
gap is `(5-sqrt(5))/60`, using the explicitly imported exclusion-gap theorem.
Population mixing traces out the record tape; the coherence theorem retains
every unnormalized classical history. These are different statements.

The [centered spectrum](../centered_particle_spectrum/README.md) concerns
`K=H-N/2`, `c=1/8`, `U=1`, and its unique half-filled ground. It does not
replace the original uncentered ground's 78-line spectrum or its four-bin
native experiment. Exact symmetry projectors yield a 32-dimensional ground
block and a 228-dimensional local charged block. Integer eigensystem
witnesses are checked with rational residual, norm and spectral bounds;
numerical diagonalization only proposes these witnesses. Both charge
branches and all certified tiny weights are retained.

| Centered quantity | Certified source-clock interval |
| --- | --- |
| Ground energy | `[-4.12405873337830, -4.12405873337828]` |
| Smallest visible excitation | `[0.36168271180597, 0.36168271180599]` |
| Largest visible excitation | `[3.88643475495059, 3.88643475495061]` |

The visible energy ratio is less than eleven. No choice of distinct visible
lines can therefore reproduce the fixed charged-lepton mass ratios with one
positive multiplicative scale. This rejects that direct finite-gap dictionary;
it does not rule out particles arising from other states or a continuum limit.
The [complete table](../centered_particle_spectrum/lines.csv) retains every
energy and weight interval, and the receipt includes the full spatial residues.

The [strong-coupling proof](../strong_coupling_refinement/README.md) constructs
the low-band embedding into the actual Fock space. With `tau=c|E|` and
`U>=4tau`, its complete-output error is at most
`min(2,4tau/U+3|t|tau^3/U^2)`. It includes intermediate doublon–hole states
and output leakage. At fixed graph, `U=c^2/J` yields interacting Heisenberg
evolution. Refining a physical three-dimensional cube makes the bound
extensive; this calculation does not prove a fixed-`J` relativistic limit.
The separate free-limit extension supplies a new compactness argument for
unbounded `U=o((log c)^(3/2))`. Its eight Weyl branches still have
`E=v|p|` and zero rest mass.

The [modular/stress proof](../weyl_geometric_modular/README.md) applies the
free-fermion modular theorem to that same complete Weyl vacuum. It checks
the spin frames, the necessary half-gradient spinor term, the filled-sea
projector identity, stress conservation and the normalization `2pi/v`.
An explicit class of normal state perturbations has a finite Araki-relative
first law. A stress variation with spatial Hessian norm at most `M` has
ball remainder at most `4pi^2 M R^6/(35v)`. The result supplies sector-relative
matter parts of G2–G5 in the existing gravity stack; it supplies no area law
or generalized-entropy stationarity assumption by definition.

## Reproduction and custody

From the repository root, run the strict tests for the four packages:

```sh
python -m pytest -q -W error model/repair_coherence_lifetime model/centered_particle_spectrum model/strong_coupling_refinement model/weyl_geometric_modular
python -m model.repair_coherence_lifetime --check
python -m model.centered_particle_spectrum --check
python -m model.strong_coupling_refinement --check
python -m model.weyl_geometric_modular --check
```

The work is isolated from an independent receipt rebuild in the main
checkout. The starting commit `fc2aab8c0` has pre-existing stale source
bindings: both the completion-program and finite-map validators reject
them before this campaign's ledger edits. New package receipts distinguish
their direct checks from inherited custody. In particular a deterministic
`--check` can pass while the receipt explicitly reports
`dependency_custody_pass=false`; it verifies that the pending status and the
new calculations reproduce, not that old experiments were rerun.
The modular input's historical FA82 checkpoint
`6a66874bb72384e3dfb3fddaeb97b5df089fac18` was separately checked against all
512 bound Git source blobs with no mismatch.

The [scoped validation record](spectrum_refinement_geometry_review.json)
records **140 passing strict tests**, four matching cold replays, and
successful auxiliary-only TeX compilation of all eight primary papers
and all four new proof fragments. Claim and premise views validate.
The full completion-program and attachment-map checks remain blocked by
older source bindings; their failures are retained in that record.
No released PDF is replaced in this campaign.

No old experiment is silently recertified, no accepted gate is broadened,
and no measured-particle comparison is promoted. Full interacting
relativistic matter with the observed massive and massless species, its
native detector and common units, and a responsive geometry satisfying the
remaining Einstein premises still require construction.
