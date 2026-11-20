# Empirical status

The [executable model](../model/README.md) reproduces particle
stability, an observable of nature. Bound objects persist in free
evolution, read on the observers' agreed records. A pair of loads written on
one seam by a local emitter is read as a pair on some seam at every reading,
with probability at least one minus a proved bound
([binding theorem](../model/contrast_memory_spherical_v1/PROOF.md)); it is read on seams
other than the first; nothing holds it but the contrast-memory law; and the
same law without its agreement energy loses it. At coupling `kappa = 8` the
largest probability of reading no pair is `0.0096` on one cell and on
federations of 20, 80 and 320 cells, against `0.80` for the exchange law
alone on 20 cells. Three loads on a face are read as an agreed face object in
the same way. The coupling is supplied.

Under a second candidate law on the same carrier, the
[scalar channel](../model/scalar_channel/README.md), an object on a ring of
four seams is kept for every interval, and a local detector reads its
excitation from counts at five waits. After a pulse its seam record swings
with no drive at one line, which an observer reads at `4.8990` of its own
reference and finds equal to the level difference of the object read
separately. A changed law with the same alphabet, routing and symmetry loses
both ([cage](../model/scalar_channel/cage/README.md)). The
[observation ledger](registers/OBSERVATION_LEDGER_V3.md) records these as
particle persistence and an internal clock attained on this construction.

The target is agreement between observers' measurement statistics and their
effective interpretations. Every test of the model compares distributions of
admissible observer records with what observers in our world record. The
state of the simulation is no physical space, and a difference between
microscopic state encodings is no failure of the model. The theory proves that
observers who inhabit the simulation can describe their experiences as
relativistic 3+1D continuum spacetimes (the flagship's proposition on smooth
effective field descriptions). A comparison with a laboratory record requires
identifying the same preparation, intervention and readout in the model and
the experiment before comparing their probabilities.

## What the model has to supply next

The other observables of the theory on this model are work in progress: a
selection of the couplings, a record of currents that later events read, the
class of an object as a record, the twelve port currents from the law's own
services on the bound object, one limiting speed for all objects, one law
that carries the objects of both candidate laws, a clock that starts with no
supplied phase, a refinement map that preserves seams, an observer clock
normalized to proper time, and the regional response that the Einstein
theorem consumes.
No laboratory dataset is compared with the model's records.

## Comparisons of the theory with public data

These comparisons use the theory's conditional maps, without the executable
model's records. The linked ledgers keep the individual rows, inputs, less
successful outcomes and reference editions; their numerical agreement is not
automatically a likelihood or an independent prediction. Reusing records,
overlapping cosmological likelihoods or reproducing a fitted reference does
not multiply independent evidence.

| Family and authoritative evidence | Established comparison and boundary |
| --- | --- |
| [Fine structure](../code/P_derivation/README.md) | Structured inverse-alpha output `136.994835177413`, a `300.388218 ppm` deficit against CODATA 2022; the gauge-width variant gives `137.035660136947`, a `2.474095 ppm` deficit. The declared bare map's hadronic part is the quark-screening ansatz `S = 1 - N_c alpha_3/pi` (about `0.887`). Replacing the ansatz with the measured KNT19 payload gives `136.3827548` on `[136.3670, 136.3985]`, which excludes the measured value; closing the gap with source anchor and transport fixed would need the payload `0.032244`, outside the measured `0.027609 +- 0.000112`, or an anchor shift of `0.635`. Closer agreement of a variant does not select the physical continuation. |
| [Cosmology postdictions](../code/cosmology/postdiction_ledger/COSMOLOGY_POSTDICTION_LEDGER.md) | Primordial tilt/running, tensor, isocurvature, non-Gaussianity, curvature, capacity and dark-energy comparisons retain their premise and reference ancestry. Edge tilt near `0.966021` is compatible with several subsets. For the ACT DR6 P-ACT-LB combination (Planck, ACT DR6, CMB lensing and DESI 2024 (DR1) BAO) the [posterior-chain audit](../code/cosmology/tilt_posterior/README.md) places it at the lower `0.612` percent tail position (`3114/508985`) under baseline Lambda-CDM. That retrospective tension is a posterior position, not a p-value or a frozen verdict. Seen-data diagnostics are not frozen scores. |
| [WMAP angular-template exposure](../evidence/custody/falsification/frozen_targets/fz11_2026-07-31/frozen_target_spin_six_primitive_port_prediction_2026-07-31.md#prior-related-exposure) | The retained custody document reports a preregistered A5 search with null family-wide `p=0.64`. This is a documented exposure, not an independently replayed raw-data campaign in the current comparison packages. Its data and CMB likelihood class are excluded from an FZ-11 verdict; the frozen dispersion-coefficient relations have no qualifying physical comparison. |
| [Fixed-capacity dark-energy comparison](../code/cosmology/fixed_capacity_wlaw/README.md) | Retained DESI DR2 Gaussian diagnostics for `(-1,0)` range from `2.60` to `3.93` standardized units across overlapping combinations. These are not frequentist exclusions; chain-based interpretation and eligibility are governed by the declared protocol. |
| [SPARC and dark-sector diagnostics](../code/cosmology/rar_deep_regime/README.md) | Deep-regime acceleration/BTFR comparisons fit an acceleration constant and expose tension under declared weighting and mass-to-light choices. Full-RAR calibration reuses the same data and a published interpolation; `a0` is fitted. Lensing/quadrupole calculations with uncertified integration or supplied field laws are not native physical confirmation. |
| [Frozen physical targets](FROZEN_PREDICTION_LADDER.md) | Cosmology, angular-invariant, dispersion, record-count and Kerr targets have their own frozen, pending, deferred or void status. A registered prediction, custody proof or synthetic sensitivity test is not an acquired empirical verdict. Their immutable contracts are authoritative. |

The [public-data audit](research/PUBLIC_DATA_PREDICTION_AUDIT.md) gives the
primordial tilt, coupling and photon-dispersion comparisons in detail. The
[falsification programme](OPH_FALSIFICATION_PROGRAM.md) and
[observation ledger](registers/OBSERVATION_LEDGER_V3.md) distinguish these
bounded exclusions and conditional comparisons from programme-wide claims.
This page is a navigation and interpretation summary, not a replacement
prediction register or a statistical combination of dependent results.

## Reproduction

From the repository root:

```sh
python -m model
python -B -O -m model.scalar_channel.cage
python -m pytest -q -W error model/contrast_memory_spherical_v1 model/scalar_channel
python -m pytest -q code/cosmology/tilt_posterior
```

The [reproduction guide](../REPRODUCE.md) gives the pinned environment and
the theorem, source and publication replays. Passing those checks establishes
reproducibility of the stated analysis; it does not turn synthetic records
into acquired data or authenticate a physical apparatus.
