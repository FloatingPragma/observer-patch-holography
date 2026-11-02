# Ordered spin sector, native observers and spatial particle channels

2 November 2026. This pass develops five existing construction issues:
[#1138](https://github.com/FloatingPragma/observer-patch-holography/issues/1138),
[#1139](https://github.com/FloatingPragma/observer-patch-holography/issues/1139),
[#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142),
[#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143), and
[#1145](https://github.com/FloatingPragma/observer-patch-holography/issues/1145).
They connect the observer/reference, protected sector, interacting limit,
particle dictionary and existing Einstein theorem stack. None is closed by
the scoped results below. The simulator, scheduler, repair, memory and
observer controls are explicit admitted specifications; no physical
identification is inferred merely from being able to inspect them.

## Constructed results

| Attachment | Result | Domain |
| --- | --- | --- |
| FA92 | Ordered ground sector and positive excitation spectrum accumulating at zero | Infinite cubic spin-1/2 antiferromagnet of the protected source limit; fixed lattice spacing |
| FA93 | Complete native spin gates, field repair and probes, with a finite local-history resource construction | Supplied input marginals or explicitly prepared Neel product; declared controller and repair schedule |
| FA94 | Nine invertible finite-lattice stress responses converging in a common lattice/radius limit | U=0 free filled-band ground; specified shell family and lab-slice stress probes |
| FA95 | Actual spatial cubic representation of every certified local particle line | Unchanged finite centered cube; full CAR lift and exact coherent central filters |

## Genuine gapless spectrum in the interacting spin branch

[The ordered-sector proof](../spin_goldstone_sector/README.md) matches the
existing source-derived interaction to the precise reflection-positive
antiferromagnet hypotheses. It imports the stated infrared and ground-selection
theorems explicitly. A rational 64-cubed integral enclosure removes the rounded
numerical step from the order bound and supplies `m_s>1/1000`.

For an L2-normalized smooth momentum packet supported on
`p/2<=|k|<=p<=1`, with any known `0<m<=m_s`, a smooth energy-filtered transverse
spin observable has squared ground-state norm at least `m^2 p/96`.
Its positive energy support lies within
`(m J p/16,64 J p/m^2)`, and its mean energy is at most `48 J p/m^2`.
It is orthogonal to the entire zero-energy subspace. Taking `m=1/1000` and
letting p decrease proves actual gaplessness, beyond an eigenvalue spacing
that might only count ground degeneracy.

The selected spin states embed in the same large-coupling source family and
its complete recorded repair limit. This is an interacting fixed-spacing
theory; no isolated particle pole, Lorentz symmetry or spatial continuum is
deduced. Native preparation of the selected ordered ground is not claimed.

## Native observations with explicit precision and cost

[The native observer construction](../native_spin_observer/README.md) uses
five pulses to encode one spin and 62 to implement an exact ideal Heisenberg
bond word. Complete adjacent C and E field-repair dilations cost 69,635 and
94,207 pulses, retaining their off-code action, raw records and first failures.
Nonadjacent endpoints use fermionic routing with its CAR signs.

Six finite executions compare complete output states, including entangled
reference columns. The mathematical local-history cutoff has trace error
`(2m+1) B_R(T)` for m interventions, uniformly in an untouched reference.
The finite-box product formula adds `2 A J^2 T delta`, where A counts
overlapping bond pairs. A rational command algorithm derives dyadic integrated
port coefficients from the exact source law; its full-channel error is at
most `58 N 2^-p` for N pulses and p coefficient bits. Executed binary64 words
retain their separate fixed coefficient errors.

The finite resource algorithm concerns its specified preparation and
repair/probe/rotation schedule. Arbitrary adaptive intervention synthesis
requires its own counted words. The initial ground marginals are not supplied
by the locality theorem, and repairs are scheduled between completed gate
words. These are visible interfaces, not hidden assumptions.

## Lattice stress responses for the existing gravity stack

[The stress refinement proof](../lattice_weyl_response/README.md) controls
the free negative-energy lattice projector without assuming a spectral gap
or a momentum cutoff. Exact cell averaging and the H4 shell bounds give both
field and first-derivative convergence. Symmetric, finite-site quadratic
stress observables then converge as smeared forms; a pointwise estimate is
not substituted for a derivative estimate.

Writing `rho=R/L` and `q=a/L`, every normalized matrix entry differs from
the explicit continuum response by at most
`C_lattice q/rho^4+C_ball rho^2`, with rational constants retained in the
receipt. The sufficient trajectory `rho<=2^-28`,
`q<=rho^7/(1+C_lattice)` preserves invertibility with inverse infinity norm
at most `272/15`. All eight source valleys remain in the algebra. The
specified perturbation excites one limiting valley and leaves seven spectators.

These are lab-slice stress measurements. At rest their continuum target is
the modular-ball energy; finite lattice detectors are not asserted to be
modular generators. No lattice modular-log convergence is claimed.
The result supplies a matter-response refinement input to the existing G5
target. Area normalization, stationarity, geometric null balance, responsive
metric and universal positive Newton coupling remain separate obligations.

## Spatial labels and retained energies

[The cubic-channel calculation](../cubic_particle_channels/README.md) applies
the actual Weyl-spinor spatial rotation, including the full CAR lift. Its
48-element character table, signed ground invariance, commutants and central
filters are verified exactly. It classifies all 228 charged and all 270
positive local neutral lines from the existing spectral certificates.

| Lowest visible channel | Source-clock energy enclosure | Actual cubic representation |
| --- | --- | --- |
| Full-Fock neutral gap | `[0.01502125417612,0.01502125417614]` | A1 + E |
| Local charged line, either charge | `[0.36168271180597,0.36168271180599]` | G1 + G2 + H |
| Local density line | `[0.81180294910720,0.81180294910723]` | T1 |

The lowest neutral triplet therefore cannot be called a spatial vector solely
from its hidden-spin-one label. Every visible charged line contains G1, G2
and H, with complete-probe trace fractions 1/4, 1/4 and 1/2. Cubic labels have
nonunique continuum spin extensions; the table does not assign electron,
photon or other measured species, or turn finite gaps into calibrated rest masses.

## Ledger and verification boundary

The [campaign contract](gapless_native_spatial_response_campaign.json) and
[validation review](gapless_native_spatial_response_review.json) retain exact
evidence and checks. The finite attachment map grows to 95 reviewed domains,
90 realized within their stated domains. The 346 canonical claim identities,
84 canonical premise rows, existing Lean library and historical physical
acceptance statuses are preserved. This adds no physical comparison or
canonical physical discharge.

The flagship is updated by replacing its current construction summary;
detailed proofs live in the focused packages and supplementary fragments.
No CI job or workflow is added. New direct receipts distinguish their passing
checks from inherited source-hash failures. Historical receipts and accepted
evidence are preserved for the separate CI cleanup.

The remaining completion work is a common interacting relativistic spatial
continuum with the required massive and massless sectors, a tested physical
species and unit dictionary, and the Einstein stack's geometric response
premises on that same family. The present pass does not declare these achieved.
