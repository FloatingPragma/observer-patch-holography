# Passive energy uptake from the same phase-wall daughters

Two supplied local receivers can gain energy from the daughters of the existing
phase-wall source. At source exposure 32 and receiver delay 32, both receivers
are excited on **36,025 of 67,108,864 original wall launches**, against **7** on
the ordered sham. Each receiver's excitation energy is 32.5153809758 exchange
units. These are simulated terminal energy records from the native material
flight followed by a declared detector interaction. They identify energy
delivery by the daughter material; OW-08's radiation clause remains open.
[Complete records and work accounts](results.json) · [Executable](measure.py)

## Source, calibration and intervention

The source is unchanged: a four-dimer ring with two opposite phase walls and
two three-dimer wings in the existing boundary-conditioned ground preparation.
The dimer-number guards, native contrast 8, both central ordered phases and
source/sham launch definitions come from [daughter_spectrum](../daughter_spectrum/README.md).
Every source is freely evolved until exposure 8 or 32. A sharp complete ring
read then holds its observed orientation during the daughter assay. Other
ring words remain in the source-miss bin. Nothing is enrolled by a receiver
click. The source-order probability at exposure 32 is 0.0307646345 for walls
and 0.9994893094 for the sham.

Each receiver is an added two-level apparatus, initially in its ground state.
Its fixed Hamiltonian is

\[
H_{\rm receiver}=H_{\rm wing}-E_{\rm wing,0}
 +\omega |1\rangle\langle1| +g X_{\rm far}X_{\rm detector},
\qquad \omega=32.5153809758,\quad g=1/16.
\]

Here the identity tensor factors are implicit. The only interaction with the
wing is its native far-dimer swap. The apparatus coupling, two-level energy
gap, local addressing, source hold and terminal receiver-energy measurement
are supplied. This is a static reciprocal interaction during the assay.
There is no inverse wing evolution or frequency-selective time-bin
postselection. It does not establish that the detector coupling is itself a
native Interval Echo process.

The frequency is calibrated from the middle of the three known first-band
wing gaps before evaluating source responses. The retained receiver delays
are 0, 16, 32 and 64. On either native boundary, at delay 32 the middle
first-band eigenstate excites the receiver with probability 0.9704883, while
the ground gives 0.000001174. All eight input levels, every delay, both
boundaries and each control Hamiltonian are retained in the calibration.
The finite source and spectrum already motivated these settings; this is
apparatus development with no prospective confirmation claim.

For each setting the record has nine outcomes: neither, left only, right
only and both receiver excitations for each ordered ring phase, followed by
source miss. There are 2^26 original launches per cohort, seed 81168. The
measured native wall-minus-sham double-click difference is 0.0005367100,
above the simultaneous 99% Hoeffding difference half-width 0.0003553980
across all 24 arm/exposure/delay comparisons. The inferred exact difference
is 0.0005314767. The sample and the derived probability are distinct fields.

## Where the receiver energy comes from

The unmeasured native source flight conserves the original wall-minus-sham
mean surplus of 65 units. The receipt separately records the work of the
complete source read and hold, including missed ring words. At exposure 32
it changes the full wall ensemble's mean energy by +0.0263774383 and the
sham's by +0.0367737094. A conditional source-energy balance is not obtained
by subtracting these unconditioned work values from a heralded branch.
The earlier [bond measurements](../energy_budget.json) remain the
unconditioned account of energy leaving the freely evolving source.

The detector coupling is then switched on with zero mean work because each
detector starts in its energy ground and has zero X expectation. For every
heralded phase the static coupled Hamiltonian conserves energy. At exposure
32 and delay 32, the wall ensemble's expectations per original launch are:

| Quantity after the source read | Energy |
| --- | ---: |
| Initial wing excitation in ordered-source branches | 1.99607822475 |
| Wing energy decrease during receiver interaction | 0.986568420875 |
| Detector energy increase | 0.986565373596 |
| Final wing-detector interaction energy | 0.000003047279 |

Wing decrease equals detector gain plus the final interaction energy.
Switching the coupling off costs the negative of that final interaction
energy, which is retained explicitly. The inferred detector gain per
ordered wall source is 32.068165 units across both receivers. A both-click
record reads 2ω = 65.0307619516 units of receiver energy. This selected
energy is not an estimate of the full ensemble's original formation energy.
Ground-state virtual excitations and coupling work are included in the sham
and calibration, rather than being called emitted energy. The ledger's
coupled-energy drift is below 2×10^-9 per original launch for every record.

## Controls and remaining criterion

With both attachments cut at release, both receivers give zero double
counts at the quoted setting. The original wing input is kept, so the cut
is a quench with its own small background. Normalizing by source-order
probability gives equal detector distributions for wall and sham to
10^-10, as source-guide factorization requires. The held wing Hamiltonian
removes the attachment term and retains both internal contacts. The test
independently matches it to each corresponding block of the full source
generator.

With contrast zero, the same source, routing, guards, receiver frequency
and detector interaction also give zero double counts. Its exact
wall-source probability is 2.95×10^-11; the changed law moves the guide
spectrum away from this calibrated hard line. With receiver coupling zero,
the calibrated excited input transfers no energy. The controls would falsify
an interpretation in which the reported native coincidence is merely an
initial receiver population or a source-independent guide background.

The result supplies a passive, energy-accounted read of the same daughters
whose material-contact displacement and downward transitions are measured
in the parent experiments. It derives no propagation speed or massless
dispersion, transverse polarization, freely propagating vacuum radiation,
photon number, electron/positron identity or positronium channels. The
material guides and the protected source remain essential scope. The
[literal criterion and accepted example](../daughter_spectrum/PROTOCOL.md)
therefore leave the full OW-08 claim partial and on HOLD.

## Replay

From the RER root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_08_annihilation.passive_receiver.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_08_annihilation/passive_receiver
```

The default replay only compares the receipt. `--output PATH` creates a new
record exclusively and serializes before opening it. Four tests independently
construct the guarded wing from port loads, match the complete ring blocks,
compare the coupled propagator with a direct matrix exponential, test
disconnection, retain every cohort and miss, and reconstruct detector energy
from the terminal click probabilities. Raw normalization residuals and the
sampling-complement correction remain in each record. Input and executable
hashes bind this finite apparatus calculation.
