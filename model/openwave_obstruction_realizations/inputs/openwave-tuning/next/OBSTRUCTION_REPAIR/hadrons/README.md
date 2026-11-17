# Hadron obstruction repair

The missing-colour claim is construction-specific: the frozen finite native
colour-chain source already carries SU(3) Gauss links. Its physical scale and
its transfer to the common electromagnetic source remain open. Merely copying
the electron unit onto it is not a transfer proof.

This package supplies a concrete positive-current amendment that opens the
previously forbidden neutron-like beta threshold. It does **not** reproduce
physical nucleon masses, the physical neutron-proton splitting, a beta decay
record, the nuclear spectrum or a physical confinement scale.

| Problem | Finding and possible repair | Executed result |
|---|---|---|
| Row 14 ordering and row 25 threshold | The inherited vectorlike number prescription has 0<a<b<1, so its scalar flavour gap is always b-a<1 electron unit. A new local interaction must change the relative binding. | The positive native-current square described below gives gap [1.1131736031415502,1.1131736031706540]. Its ground beta threshold is open; old and reflected-current controls fail. |
| Row 15 nucleon masses | The original local singlet trial has energy at most 57/16, independent of the colour string coefficient. Thus its electron-normalized proton candidate stays below 1.821 MeV; increasing string tension cannot repair the >515-fold discrepancy. | 48 exact endpoint certificates on L=2,3,4,5 and kappa=1/2,1,2, both charge members. A new source mass-generation or rest-energy identification is required. The current-square amendment leaves this scale problem unresolved. |
| Row 17 helion and alpha | An isospin-scalar interaction cannot repair the mirror gap. The existing four-site oscillator vanishes on helion and lowers alpha by at most 63.875 keV, whereas the remaining alpha excess exceeds 11.69 MeV. | The frozen rational full-source enclosures exclude both held-out targets. An independently fixed charge-dependent response and substantially different attractive dynamics are required; a larger fitted coefficient is not a source-derived solution. Halo outgoing channels remain absent. |
| Row 22 colour and tension | Independent electric and colour registers remove the one-Abelian-field conflict. The existing tree SU(3) bridge supplies the colour object, but its kappa and physical clock/length mapping are independent of the electron anchor. | No new 1 GeV/fm prediction. The common scale needs an explicit source-level transfer, and the colour coefficient must retain its ancestry. The existing integer-flux 10,000-fold failure is retained. |

## New local interaction

On each cell let B_c=(N_quark-N_antiquark)/3 and let Q_c be the actual
electric current restricted to the colourful internal module. Declare

    H_new = H_colour + H_vectorlike_number + U sum_x (B_c,x-Q_c,x)^2,
    U=1,

using the already declared colour coefficient U=1. On the three-quark sector
the added current is (-N_up+2N_down)/3. No hadron mass or measured mass
difference chooses U. The new operator was selected retrospectively to test
the known flavour-scalar obstruction; it is an architectural proposal, not a
claim of a unique OPH law or an independent physical prediction.

The current is a colour and spin scalar; its square is positive, commutes
with electric charge, and has finite cell support. It is invariant under
charge conjugation because both B_c and Q_c change sign. It vanishes on the
colourless electron sector. Thus it preserves the old positive full-Fock
lower bound instead of creating an unstable attractive vacuum. It breaks the
old global weak flavour symmetry, which is precisely why the scalar-shift
no-go no longer applies. All charge, spin and occupation complements remain
in the full source definition; the executed spectra concern the complete
14-dimensional fixed projections of the two-cell three-quark singlet.

This is a new join of the native colour-chain Hamiltonian and the original
vectorlike free number-energy intervals. The original baryon-ordering package
also had a weak Casimir term; its full source is not asserted identical to
this amended chain.

## Exact tests and physical comparison

`audit.py` transforms the normalized bosonic occupation matrix by the diagonal
Gram factors prod(n_i!). Every transformed entry is rational, including the
new interaction. Exact congruence inertia bounds the lowest eigenvalue at
both charge assignments. The inherited up/down number intervals are then
propagated outward as scalar shifts. No numerical eigensolver supplies a
certificate.

The old gap is [0.11663300395906927,0.11663300398817310]. A separate positive
grounded-chain electric-field term at the existing coefficient 1/16 lowers
it to [0.09583275835342729,0.09583275838253112]; this tested alternative fails.
The current-square source raises it above one; swapping the two flavour
components of the current reverses it to approximately -0.8799075952.
The amended gap is approximately 0.568831 MeV, with about 0.057832 MeV
available above a zero-mass-neutral beta threshold. The comparison gap is
1.293333 MeV. Opening a channel is therefore a repaired structural
obstruction, not a successful physical mass prediction or executed decay.

`native_extension.py` executes the new controlled diagonal on all 32 code
and control columns for each charge member. Its 652 actual M6 pulses retain
every first-export complement. The maximum diagnostic operator error is
2.91e-14 and complete Gram error is 5.82e-14; the exact rational native error
budgets accompany the output. The false-control and padded columns are
identities. The full amended filter, preparation, long-time evolution and
decay were not executed. In particular, the old neutron filter cutoff must
be rebuilt: its old low band cannot be transferred after the energy changes.

Run with:

    python openwave-tuning/next/OBSTRUCTION_REPAIR/hadrons/audit.py
    python openwave-tuning/next/OBSTRUCTION_REPAIR/hadrons/native_extension.py

`INPUTS.json` binds local snapshots of the reviewed source and inherited
receipts. `NATIVE_RESULT.json` also binds the compiler helpers used live.
These are producer results requiring an independent review.
