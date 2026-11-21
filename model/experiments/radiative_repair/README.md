# Opposite charges and a neutral field

Two opposite charges prepared at their lowest frequency can convert into a pair of neutral quanta. On the
20-cell carrier the observer reads 2.553 units of neutral-field energy after 32 clock units. The initial
contact costs 0.090 units. After allowing for that cost and sampling error, at least **2.436 units** of the
initial charged energy reach the neutral field. The source also has preparations that never convert. On
one cell their space has dimension 97 out of 144 charged-pair states.

## Local law and source map

This experimental construction takes the labeled ports, internal seams, external matching and emitter
addresses from the [A-contrast-memory-spherical-v1 carrier](../../executable_model/carrier.py). It places two fermion modes
a and b and one neutral oscillator c at each port. The oscillator retains occupations 0, 1 and 2. The
charges of a and b are opposite by design. The source has one coefficient J > 0 and no species-specific
rest-energy term. With ℏ = 1, its Hamiltonian is the sum of positive local squares

    H/J = Σ_edges [(aᵢ+aⱼ)†(aᵢ+aⱼ) + (bᵢ+bⱼ)†(bᵢ+bⱼ) + (cᵢ−cⱼ)†(cᵢ−cⱼ)]
          + Σ_ports (bᵢaᵢ−cᵢ²/√2)†(bᵢaᵢ−cᵢ²/√2).

Each edge contributes all three displayed terms. The factor 1/√2 gives the normalized two-boson state unit
matrix element. Equal strengths, the relative edge signs, the local alphabets and the conversion channel
are supplied. The half-edge convention of the contrast-memory law of A-contrast-memory-spherical-v1 and the seam/wedge
fermion of the [scalar channel](../../candidate_laws/scalar_channel/README.md) have different generators. Only their declared port carrier
is used here.

The vacuum has zero energy. The conserved signed charge is Nₐ−Nᵦ and the conserved total occupation is
K = Nₐ+Nᵦ+N_c. A local conversion consumes one of each charge and releases two neutral quanta. The number two
is built into the term. Transfer probabilities, the released energy and the response to a later intervention
are calculated. At K ≤ 2 the contact is 2J times the projector onto the difference of the local charged pair
and neutral pair, normalized to one. A local rotation to that basis, a bit write, phase 2JΔt and the inverse
write and rotation execute its individual propagator. Its coherent sum with the edge terms defines the
continuous-time source used here. Digital Trotter error and autonomous laboratory switches are outside this implementation.

For total charge zero and K = 2, Xᵢⱼ is the amplitude for a positive charge at i and a negative charge at j.
The symmetric matrix Y represents two identical neutral quanta, with total norm Σ|Yᵢⱼ|². The probability
of a distinct unordered outcome {i,j} is 2|Yᵢⱼ|². The diagonal outcome has probability |Yᵢᵢ|². For adjacency
A and degree matrix D, put Q = D+A and L = D−A. In units of J the generator sends

    X ↦ QX + XQ + diag(diag X − diag Y)
    Y ↦ LY + YL + diag(diag Y − diag X).

This includes every state of the sector. Even fermion operations are primitives on the graded local registers.
Implementing them with the contrast-memory law on ordinary port bits requires a separate compiler.
The Fock-space tests reconstruct all local generators with fermion
signs and truncated oscillator ladder factors. The blocked-conversion control retains both diagonal contact
costs and removes the cross term. Parking the neutral field retains its onsite costs and removes neutral hopping.
Both controls use the same position measurement. The energy settings follow each Hamiltonian. The parked
field uses onsite number detectors and the blocked contact uses a coincidence projector. One-cell records
resolve individual ports. Federation records resolve each cell, every ordered pair of signed detector
regions and every unordered neutral coincidence.

## A rest frequency from the port adjacency

With one excitation of either charge sign, the contact vanishes and the one-body generator is JQ. Its
one-cell characteristic polynomial in J units is (x−10)(x−4)⁵((x−5)²−5)³. The observer addresses one labeled port and compares a
freely evolving arm with a parked copy. Two analyzer phases read the real and imaginary parts of its
return amplitude. An isolated neutral seam with the same coefficient J supplies the clock. Its two local
occupation settings give cos(2Jt) and sin(2Jt). Unwrapping uses a supplied bound on the interval between
samples. The reader selects its line count from the data, fits the first 384 clock readings and predicts
the remaining 128.

The positive-charge lowest line is **2.763965** in clock units and the negative-charge line is **2.763970**,
against the exact value 5−√5 = 2.763932. Both signs share one one-body generator, so the two readings differ
by sampling alone; the first clock tick is read 4.07 × 10⁻⁴ high, which the fit carries into the line, and
with a free phase the readings fall within 10⁻⁵ of the exact value. All four charged lines appear, with
weights 3/12, 5/12, 3/12 and 1/12 in order of frequency. The neutral reference gives lines 0, 5−√5, 6 and
5+√5. The shared number 5−√5 follows here from the port adjacency.

On a federation, every external seam adds a positive square to the sum of one-cell charged operators.
Equality at the one-cell lower bound requires a lowest-band vector in each cell and amplitudes summing to
zero across each external seam. Two adjacent cells supply six coefficients and at most five incident
constraints. An exact rest packet can therefore have support confined to those two cells, including when all
other cells are present. The prepared profiles satisfy the full generator equation on the 20-, 80- and
320-cell carriers. Their counted lines agree with the point-source one-cell line.

The rest-band multiplicities on these three carriers are exactly **30, 120 and 480**. To count the constraints,
the verifier multiplies the other three factors of the cell minimal polynomial, leaving a rank-three
multiple of the lowest-band projector over ℚ(√5). The three external-port rows of a cell are either
independent, or two antipodal rows cancel and the third is independent. Reduction at √5 = 45 modulo 101
certifies the nonzero minors. The cancellations are verified in both rational coefficients. In a dependence
of external seam equations, an independent tripod forces all three seam weights to zero. An antipodal
tripod equates two weights and forces the third to zero. Propagating these equalities on each actual carrier
gives zero free dependencies. Subtracting the number of external seams from three times the number of cells
gives the stated multiplicities. The extensive flat band permits stationary packets.

With the one-cell operator fixed, any static positive rank-one term per external seam leaves an extensive
lowest charged band at the same energy. For every choice of strengths and phases, a preparation in that
band keeps constant detector probabilities during free evolution. The exact multiplicities persist when
both endpoint factors of every seam are nonzero. [Proof and integer verifier](FLAT_BAND.md)

## Conversion, transport and the part that cannot radiate

A localized positive and negative charge begin on different ends of an internal seam. On 20 cells, up to
6.270% of all trials return two neutral quanta, and up to 4.921% return neutral quanta in two distinct cells
over seven specified waits. The blocked control has zero neutral counts. The parked control permits conversion
and has zero split coincidences. Starting with a neutral pair yields signed pairs in up to 27.80% of all trials.
The experiments are finite reversible evolutions.

Exchanging the two charged addresses transposes X. An antisymmetric X has zero diagonal and remains
antisymmetric under QX+XQ. It can never convert. The even and odd superpositions of the two initial
assignments have identical initial signed-position probabilities. Only the even preparation radiates.
An ordered preparation produces exactly half its neutral probability at every wait.

There is a stricter obstruction on one cell. Let H₀ = Q⊗I + I⊗Q and let B insert the twelve coincident
charged-port states. A charged vector stays nonradiating at every time precisely when it is orthogonal
to every H₀ᵏB. Vanishing of the neutral amplitude forces the contact source B†X to vanish at every time,
which leaves free charged evolution. Conversely, all powers vanishing imply this condition by the
exponential series. The numerical decomposition proposes a projector P onto this span. Integer arithmetic
checks symmetry, P²=P, PH₀=H₀P and PB=B with common denominator 480. A nonzero rank-47 minor of the
integer powers through order nine modulo 1009 matches tr P. This certifies the entire span.

Of the 97 dark directions, 66 are antisymmetric and 31 are symmetric. The exact upper bounds on **total**
neutral probability for an ordered point source at port-graph distances 0, 1, 2 and 3 are 1, 29/120, 29/120
and 1/2, respectively. They are projections onto the converting span and need not be attained at a finite
wait. Graph distance labels the preparation and gives no physical length measurement.

## A receiver that converts the field back

The relay starts from charges localized at the two ends of the emitter seam. The energy-release experiment
uses a separate preparation in the lowest band. During relay emission the only active converters are in
cell 0. The local meter in the neighboring cell 6 reads neutral occupation 0, 1 or 2 near clock reading
8.00014. Its three projectors preserve coherence within each outcome. The receiver converter is then
enabled or kept blocked. For every meter outcome, every
terminal signed and neutral region pair is counted. Observing two arrivals guarantees that emission took place
at the sender and both neutral excitations reached the receiver before the second intervention.

One clock unit later, **528 of 16,777,216** trials with the receiver enabled have two neutral arrivals followed
by both opposite charges at the receiver; the record law gives 3.02 × 10⁻⁵, an expected 506. The blocked arm
gives **0 of 16,777,216**, a sampled zero: its record law gives 4.65 × 10⁻⁹, an expected 0.08. The
denominators include every meter outcome. A simultaneous one-sided lower bound for the excess across the four waits is
2.59 × 10⁻⁵ at this wait at nominal 99.9% coverage. At the longest delay the blocked arm gives 79 such records.
Returning neutral quanta can reconvert at the sender and charged excitations can travel back. The complete
late-time control is part of the result.

Enabling the receiver immediately after the meter supplies exactly zero mean switch work, even conditional
on any meter outcome of positive probability. The change in the Hamiltonian connects only the local charged
coincidence (meter value 0) and the local neutral coincidence (meter value 2) at each receiver port. For
every meter projector Pₕ, the diagonal block Pₕ(H_on − H_off)Pₕ therefore vanishes. The replay checks all
matrix elements changed by the switch. Before the meter, a coherent superposition of the two coincidences
can have nonzero switch work. The meter can also exchange energy with the system. The identity concerns
the instantaneous switch of the modeled Hamiltonian. The cost of operating the external apparatus is outside this model.

## Where the outgoing energy comes from

Balanced local analyzers read occupations of each charged symmetric bond mode and each neutral
antisymmetric bond mode. Twice their mean occupations in J units give the bond energies. A local rotation
and detector read the contact projector. Each setting uses a separately initialized ensemble. Summing
their expectations gives the Hamiltonian energy, including the contact. The relay also reads the energy
before and after the local meter and the converter switch. Their measured differences enter the work
budget. Source-on free intervals conserve the expectation. The apparatus that prepares, measures and switches
the source is supplied.

Two opposite charges can occupy the same lowest-band orbital v. Their free initial energy is exactly
2(5−√5)J. Bringing their profiles together contributes J Σ|vᵢ|⁴ in contact energy. Disjoint compact
lowest-band packets have zero contact and are exact stationary states of the full generator. For the
overlap preparation, subtracting the entire initial contact from the outgoing neutral energy bounds the
drop in charged energy from below. By conservation, that drop equals the final neutral energy plus the
final contact energy minus the initial contact energy. The final contact energy is nonnegative. There is
no switch or intermediate read during this interval.

| At clock reading 32 | Initial charged energy | Initial contact | Final neutral energy | Neutral probability | Lower bound from charged energy |
| --- | ---: | ---: | ---: | ---: | ---: |
| One cell | 5.52838 | 0.15002 | 3.70987 | 0.66330 | 3.53365 |
| 20 cells | 5.52736 | 0.09014 | 2.55299 | 0.45403 | 2.43662 |

Energies are in units of J. On 20 cells the neutral energy on internal bonds outside both initially occupied
cells is 1.18043. Disjoint packets and the blocked-conversion controls give zero neutral counts. The
experiment compares eight waits with a simultaneous error allowance. The charged line is a frequency of
this construction. A physical identification of its neutral excitations as light requires transverse
polarizations, physical propagation and electromagnetic coupling. The fermionic signs, the number of
neutral daughters and the common coupling are architectural inputs. The construction supplies one charged
gap.

## Replay, costs and error conventions

Install the [model dependencies](../../requirements.txt) and set `OPENBLAS_NUM_THREADS=1` and
`OMP_NUM_THREADS=1`. From the repository root run

```sh
python -B -O -m model.experiments.radiative_repair --cell-only
python -B -O -m model.experiments.radiative_repair
python -B -O -m model.experiments.radiative_repair.relay
python -B -O -m model.experiments.radiative_repair.rest
python -B -m pytest -q model/experiments/radiative_repair
```

The short replay covers one-cell conversion and all four compact-profile spectra. Each full replay
rebuilds the source, regenerates all raw counts, decodes them and checks the archived result. The three
protocols fix the preparations, interventions, domain and gates. The source manifests bind the executable
files and carriers. `--record` requires a fresh directory. A false scientific gate is retained.

There are 12 local states per port. The exact neutral, K=2 sector has n²+n(n+1)/2 basis states,
222 on one cell and 86,520 on 20 cells. Two ordered n×n arrays store it, with Bose symmetry checked.
A hop involves two neighboring ports. A product-form digital implementation also needs a chosen schedule.
The simulator directly propagates the sparse sector Hamiltonian. A phase read uses a parked arm and a
controlled spatial recombination, the clock uses a separate dimer, and a coherent rest preparation needs
at most s−1 adjacent rotations per charge on a connected support of s ports, with ideal relative phases.
The relay needs a three-valued local meter record per trial. Port and cell labels describe the computational
hardware. The pointers supply the observables.

The three experiments store 32, 10 and 11 integer arrays. Respectively they simulate 415,923,961,856,
254,611,030,016 and 894,326,276,096 trials aggregated by setting. The rest-pair experiment uses 514 clock
settings, 40 position settings and 52,752 energy settings, each with 16,777,216 preparations. Every
setting retains complementary outcomes. The counter digitizes cumulative probabilities on a 2⁻³⁰ grid and
samples conditional binomials with 52-bit uniform midpoints. Energy error bounds use independence of
separately prepared ensembles. In the rest-pair experiment the largest detector reading is b=4 in J
units. Each one-body operator has norm at most twice the maximum
degree. Two excitations on degree at most six contribute at most 24 and the contact contributes at most two.
With total expected energy at most 26, the sum of variances of the local
mean estimates is bounded by b(26+ε)/N for N shots per setting and total digitization bias ε. The code
uses Bernstein's inequality, a union bound for at most 256 energy estimates and the worst rounding bias.
It subtracts two error allowances in the table. For the relay it inverts the exact binomial tails in both
arms and reserves error probability for digitization of the conditional draws. These sampling statements
use the ideal random-draw interpretation of the seeded simulation. Hamiltonian evolution, analyzer
identities and the two spectral dimension certificates have separate numerical and algebraic checks.
A second solver applies the Hamiltonian directly to the two amplitude matrices and expands the propagator
in Chebyshev polynomials. On the 20-cell carrier it independently checks the receiver intervention and the
released energy at clock reading 32. Its series truncation bound is below 10⁻¹³ in state norm.

The test the two experiments address is the annihilation of an opposite-charge pair into a neutral field,
read on the observers' records, and the reverse conversion of that field into a charged pair at a receiver.
