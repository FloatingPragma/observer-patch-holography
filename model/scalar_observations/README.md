# Shared scalar observations across spatial refinements

This packet gives a fixed spatially smeared field detector the same meaning
across an **interacting** spatial family. The physical box, field bandwidth,
preparation, action parameters, detector and Hamiltonian model times stay fixed.
The source executes both measurements and the feedback on the disturbed field.
Every classical record remains in the comparison.

It extends the [canonical scalar regulator](../relativistic_scalar/README.md)
and the [smooth instrument](../common_source/SMOOTH_REFINEMENT.md). The change
is a common physical observation contract across grids, with interacting
spatial error bounds. The earlier eight-site packet remains intact.

## The observation contract

The supplied scalar action is restricted to seven real Fourier modes on the
fixed torus `[0,2*pi)^3`: the constant and the first sine/cosine modes along
each axis. The code compares grids with 5, 9, 17, 33 and 65 points per axis.
Those grids sample the **same seven modes**; they do not supply that many
independent field degrees of freedom. The bare mass and speed are one,
`g=lambda/V=1`, and the oscillator basis frequency is two. No physical mass,
coupling, vacuum, geometry or clock is inferred from OPH in this packet.

The measured quantity is the fixed smearing

    Q = integral sqrt(2/V) cos(x) phi(x) d^3x.

The principal protocol is:

1. Execute the first native repair on the uniform binary 66-state carrier.
   Its current selects `exp(+i Q/4)` when positive and `exp(-i Q/4)` otherwise,
   acting on the product oscillator blank. This is a supplied preparation,
   not the interacting vacuum.
2. Evolve for Hamiltonian model time `1/8`. Read both Kraus outcomes
   `K0=cos((Q+pi/2)/2)` and `K1=sin((Q+pi/2)/2)`.
3. First bit one applies `exp(i Q/4)` to that disturbed field. Evolve for
   another `1/8` and repeat the same measurement.
4. Last bit one selects native seam F; bit zero selects E. Retain its coin,
   current and final port load along with the first repair and both field bits.

For the finite source execution, every `Q` in the preparation, detector and
kick is `Q_d=P_d Q P_d`. These functions of a finite matrix are complete
instruments, but differ from projecting the unrestricted operations. The
Hamiltonian instead compresses **after** forming its quadratic and quartic
polynomials. This distinction is tested explicitly.

There are 1,056 complete native/quantum classical histories at each code.
The field is neither reset nor postselected in the principal history.
Separately labelled controls delete the interaction, delete the first-bit
field kick, or replace the disturbed state by a fresh prepared state.
Disconnecting the last-bit native actuator is a separate feedback control.

## What is bounded

The [derivation](DERIVATION.md) proves two different spatial comparisons.

| Comparison | Target and scope |
| --- | --- |
| Source execution at fixed occupation | The same finite-code observations with exact spatial gradients. Rational bounds cover temporal splitting and spatial discretization; numerical and physical service errors are separate. |
| Seven unrestricted interacting oscillators | The same fixed-band observations with exact gradients. A finite-energy bound includes measurement backreaction and conditional feedback and decreases as `O(n^-2)`. |
| Finite occupation to unrestricted observations | The [occupation certificate](OCCUPATION_DERIVATION.md) bounds the entire principal history, including preparation, both instruments and conditional feedback. Exact residuals and counted rounding give an ideal source-to-unrestricted record TV bound below 0.054 at the finest retained grid. |

The quartic spatial integral is exact already at grid size five. Only the
gradient coefficient `kappa_n=[sin(pi/n)/(pi/n)]^2` changes. A four-point
counterexample in the tests shows why smaller grids can alias the quartic.

At four levels per mode (16,384 field amplitudes), the finest grid with
64 symmetric steps per interval has combined ideal classical-history TV
bound **below 0.001427**, relative to the exact-gradient **four-level** target.
The temporal contribution is below 0.00004588. This is an all-input finite-code
bound; the actual compiler is also checked against independent sparse evolution.

For the unrestricted fixed-band theorem, positivity gives
`W(Q) >= (sum Q_j^2)^2`, so the spatial comparison only needs a finite initial
energy. The full first measurement raises unconditional energy by `eta^2/8`;
the feedback's energy contribution is also bounded. Neither this theorem nor
the finite-code certificate removes the seven-mode bandwidth restriction or
establishes relativistic locality.

## Occupation adequacy and observations

The zero-duration detector checks compare two source read settings with the
analytic characteristic function of the supplied Gaussian blank. No parameter
is fitted, and these checks are not experimental calibration.

A stronger exploratory protocol, with detector strength two and kicks one,
exposed a failure of the four-level approximation: raising the measured mode
from four to eight levels changes the positive-kick history by about **0.0129**
in TV, larger than its roughly 0.0047 interaction-deletion contrast. Its small
blank-detector error did not predict the error of the continued history.
This stress comparison is retained explicitly in the receipt.

The principal protocol uses detector strength one and kicks one quarter.
The same cutoff comparison changes its positive-kick history by about
**0.000151**, while deleting the interaction changes it by about **0.00240**.
The receipt includes both preparation signs and a separate reference with
six levels on every mode. These are numerical diagnostics of the entire
finite approximation, including preparation, dynamics, detector and feedback.
They are not certified distances to unrestricted oscillators. The settings
were chosen after internal numerical exploration, without experimental data;
this is not a preregistered empirical prediction.

The [occupation receipt](occupation_results.json) certifies this principal
history. Finite Hermite polynomial
curves are checked against the full unbounded Hamiltonian and coordinate
operators using exact rational residuals. This bounds all occupation boundary
terms, temporal polynomial error and actual coefficient rounding, without
assuming convergence of a quartic Taylor series. The complete-output norm
bound is conservative: it also counts changes in the surviving field that
the two detector bits need not resolve.

A larger auxiliary calculation resolves the specific event **the two field
bits differ**, conditional on a positive first native current. It uses six
levels on the spectator modes and twelve on the measured mode for the
interacting action. For the free comparison, the six spectator evolutions
factor out of every record probability exactly. Unrestricted terminal effects
are evaluated with rational enclosures for their Gaussian matrix elements.
After both error budgets, the event probability is greater under `g=1` than
under `g=0` by **more than 0.00008**. The first native current is positive with
exact probability `5/66`, so a positive separation also holds for the full
native record law. Every branch remains in this comparison; conditioning here
specifies a retained event and does not alter the principal protocol.

These auxiliary codes are proof calculations, not additional source executions.
The four-level source-to-unrestricted bound is larger than the interaction
signal: it does not by itself certify that signal. The separate event
calculation proves it directly for the unrestricted target. The torus, seven
modes and action are still supplied, and no measured physical data is used.

The [published oscillator comparison](experimental_readout/README.md) provides a physical readout test: 21 complex axis observations
calibrate the Gaussian response relation, which predicts 100 withheld complex
observations with per-spin-mean RMSE 0.07639. The original data and alternate
coordinate conventions are retained. This retrospective, postselected
preparation/readout test does not identify the seven-mode dynamics or clock.

The separate [executable comparison](experimental_readout/EXECUTABLE_READOUT.md)
prepares a calibrated displaced thermal state and emits those reserved
readouts through native source operations. A two-parameter axis fit supplies
the preparation; the full 64-level detector supplies the predicted binary
probabilities. Its occupation certificate covers finite preparation and both
measurement outcomes, while statistical calibration and apparatus errors
are separate. This direct execution has a separate factorization
analysis and its distinct uncertainty calculation.

The next dynamical observation obligation is parameter identification on
declared calibration observations and prediction of withheld physical
responses, with a justified action and clock mapping to those data.
The [interacting band extension](../interacting_band_extension/README.md)
performs one controlled band enlargement. It adds the `cos(3x)`, `sin(3x)`
pair at momentum `(3,0,0)`, derives the induced low-mode mass-squared shift
`g/4` for the declared blank, and bounds the complete low-mode observations of
the continued, disturbed history on its finite occupation code. The remaining
omitted modes, repeated enlargement with an unrestricted occupation limit, and
relativistic locality are open.
The [native source-selection obligation](../SOURCE_SELECTION.md) also remains:
compiling a supplied scalar action does not select its physical law or values.

## Reproduce

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.scalar_observations --check model/scalar_observations/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.scalar_observations.occupation --check model/scalar_observations/occupation_results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/scalar_observations
python3 -m model.verify --hashes-only
```

[results.json](results.json) retains the complete classical laws, local
full-phase source residuals, rational error budgets and labelled controls.
The [original receipt bytes](data/occupation_source_receipt.json.gz) are retained
in a compressed snapshot; their decompressed SHA-256 and the original protocol
are frozen and hash-bound by the separate [occupation receipt](occupation_results.json).
The live receipt can therefore refresh source custody without changing this
certificate's input. Their numerical cutoff diagnostics retain their original
scope. The occupation receipt retains exact
complete-history budgets and the unrestricted interaction event intervals.
The larger occupation diagnostics use independent sparse reference evolution;
they are explicitly not source executions. The primary two- and four-level
codes both execute through the twelve-port compiler. At two levels the
compressed squared-coordinate and onsite matrices are diagonal, so temporal
splitting is exact; that degeneracy does not establish occupation adequacy.

The resource ledger counts quantum compiler words, including their code
load/export events, along the worst recorded branch. Two surrounding native
repairs are listed separately; classical routing/control overhead is uncosted.
Matrix caching accelerates host simulation without reducing the declared
source pulse count. Floating diagnostics are not outward interval certificates
or physical gate-error budgets. `all_pass` reports the stated software and
finite-model checks; no experimental dataset is used and no empirical
agreement or completion of OPH is claimed.
