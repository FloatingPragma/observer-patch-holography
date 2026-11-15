# Choose observations that can identify the action

This experiment tests the [calibration ambiguity](../identifiability_checks/README.md#action-and-geometry-ambiguities)
with an executed intervention. The four-level calibration cannot
see `B(N)=binom(N,4)` at any precision. A coherent preparation involving level
four, followed by a phase-sensitive read, can see it. The retained native
witness uses **41 pulses per trial** and separates its two supplied candidate
laws by total variation **0.7071067812**. No empirical outcomes enter this
packet, and no physical Hamiltonian is selected from the source specification here.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.identification_design --check model/identification_design/results.json
python -m pytest -q -W error model/identification_design
```

## The response map and its exact blind directions

Consider a specified finite background `H0` and the enlarged action family

\[
 H(c)=H_0+\sum_{k=0}^{5}c_k\binom Nk.
\]

This is a declared family of diagonal corrections, including a cubic number
term `binom(N,3)` and the positive local `B(N)`. It is not a
parameterization of every interaction or a source derivation of `H0`.
The number basis and allowed coherent preparation/read services are trusted.

Prepare `psi_n=(|0>+|n>)/sqrt(2)` and read the projector onto
`(|0>-i|n>)/sqrt(2)`, retaining its complement. For any Hermitian `H`,

\[
 2p'_n(0)=H_{nn}-H_{00}.
\]

Indeed `rho_0n=1/2` gives
`dot(rho)_0n=i(H_nn-H_00)/2`; the other occupied diagonal contributions to
the read cancel. No commutation of `H0` with the corrections is needed.
Subtract the specified background's diagonal gap. The response matrix for
the correction coefficients is exactly

\[
 A_{nk}=\binom nk-\binom 0k.
\]

The implemented exact rational rank calculation gives:

| Allowed probe occupations | Response rank | Invisible coefficients |
| --- | ---: | --- |
| 1, 2, 3 | 3 | global energy `c0`, `c4`, `c5` |
| 1, 2, 3, 4 | 4 | global energy `c0`, `c5` |
| 1, 2, 3, 4, 5 | 5 | global energy `c0` |

The greedy designer chooses the smallest occupation that increases rank:
first four, then five. The nonconstant square response is triangular with
diagonal one. Its exact inverse is the finite-difference transform

\[
 c_k=\sum_{n=1}^{k}(-1)^{k-n}\binom kn g_n,
 \qquad g_n=2p'_n(0)-(H_{0,nn}-H_{0,00}).
\]

It identifies this finite family modulo its energy gauge if the exact slopes
are available. It does not identify an unrestricted action. Adding
`binom(N,6)` leaves another exact blind direction until an appropriate higher
occupation is accessed. The global-energy gauge is relative to this
unconditionally evolved preparation/read interface; coherent control over
an entire evolution versus a bypass would be a different interface.

The parent packet's actual probes are not replaced or claimed to have this
rank. Rather, even the stronger hypothetical complete access to its entire
four-level code cannot see the compressed `c4,c5` terms: both operators are
exactly zero there. The probes here identify a concrete way out of that
obstruction. A local noncommuting `N+Q^4/24` background also passes the
derivative identity in an independent commutator check. That check executes
no complete coupled field and supplies no long-time error guarantee.

## What the native observer actually executes

The locked finite-time witness compares the following **component** pair on
the eight-level code:

\[
 H_b=N+\tfrac14\binom N2+\tfrac18\binom N3+b\binom N4,
 \qquad b\in\{0,1\}.
\]

Its time is the supplied model coordinate `t=pi/2`. It is not seconds, host
runtime, or a source-derived observer clock. The complete computation is:

1. Initialize the three-qubit code in `|0>` and execute one native `Y` pulse
   on bit two, preparing `(|0>+|4>)/sqrt(2)`.
2. Compile and execute the **complete eight-level** diagonal
   `exp(-i t H_b)`, including its global phase, in 38 native pulses.
   No target matrix or eigensystem is inserted into evolution.
3. Execute a fixed midpoint phase compensation and an `X` read pulse on bit
   two. Both candidates receive exactly the same preparation and detector.
4. Read all three bits. Retain every one of the eight output probabilities;
   no outcome is discarded or used for postselection.

The compiler builds these pulses from the existing twelve-port source
controls and coherent proper-code load/export services. Its blank register,
analogue control angles and final computational read are supplied
apparatus. The receipt binds the actual word hashes, source inputs, resource
counts, and full-phase evolution residuals. The latter are floating numerical
checks, not outward arithmetic certificates.

The source produces

| Candidate | `P(read code 0)` | `P(read code 4)` |
| --- | ---: | ---: |
| `b=0` | 0.1464466094 | 0.8535533906 |
| `b=1` | 0.8535533906 | 0.1464466094 |

The remaining six codes have zero ideal probability and remain in the
numerical record. An independent dense Hamiltonian exponential and explicit
measurement matrix reproduce the full law. The producer requires every
full-law discrepancy to be finite and below `2e-12`. The exact ideal probability difference
is `1/sqrt(2)`. This measurement attains the trace distance between the two
ideal output pure states in this component experiment.

In the negative control, occupation is committed before the same final read.
Both candidates then give `P(0)=P(4)=1/2`. The receipt retains the full joint
law of the earlier occupation and later read code. This shows precisely
which coherence the proposed measurement needs; a population read alone does
not expose the added diagonal action.

## Sampling and what rank does not guarantee

For fresh independent trials, decide between the locked candidates by whether
the fraction of code-zero outcomes exceeds `1/2`. Under either ideal law,
the [Hoeffding bound](https://www.tandfonline.com/doi/abs/10.1080/01621459.1963.10500830)
gives

\[
 P(\text{wrong decision}\mid b)\le e^{-m/4}.
\]

Thus **19 ideal independent trials** suffice for this conservative bound to
fall below 1% for either candidate. These are a proposed sampling budget;
the packet computes exact-model probabilities and does not generate or claim
laboratory shots. If an independently justified absolute probability error
`epsilon` is available for both candidate laws, the corresponding conditional
bound is

\[
 \exp[-2m\max(0,1/(2\sqrt2)-\epsilon)^2].
\]

For example, an **assumed** error of 0.01 requires 20 trials for the same bound.
The floating matrix residual is not promoted to this physical error bound.
Preparation error, coherent calibration error, memory between repetitions and
detector bias need their own control before applying it to an experiment.

Exact response rank is also not a small-shot recovery theorem. The inverse's
row absolute sums are `1,3,7,15,31`, amplifying gap errors. A direct slope
estimator `2*(p_hat-1/2)/t` has sampling radius
`2 sqrt(log(2/alpha)/(2m))/t` with confidence at least `1-alpha`, before any
finite-time derivative bias. At `t=1/128`, 200 trials and `alpha=0.01`, this
radius exceeds 29. Longer-time fits can be much more informative, but a
single finite time can alias phase differences: the witness loses all
separation at `t=2*pi`. Independent times or bounded parameter intervals are
required for broader identification.

## Identification scope and remaining requirements

Use exact response nullspaces to decide which interventions to add before
running large coupled reconstructions. Explore candidate laws with a checked
semantic evaluator, then execute source words for locked discriminators.
The present witness needs neither a large eigenbasis compilation nor another
high-precision replay of probes whose response is identically zero.

Coupled-field identification requires applying the coherent local preparations
and phase reads while the full supplied field evolves, and compare candidate
laws at times with controlled bias and sampling cost. The 19-trial guarantee
above belongs only to the displayed commuting component pair. It cannot be
transferred to the noncommuting field without its own predicted laws and
error budget.

This makes the action-selection obligation testable; it does not resolve it
from A1--A3. A source-derived result needs a source/interface condition
that excludes an alternative, with a counterexample when that condition is
removed. The [native source-selection tests](../SOURCE_SELECTION.md) are the
appropriate lane for that attachment. The planar cube realization
is also a valid metric negative control: the energy probes do not
derive three spatial dimensions or identify the physical meaning of graph
distance. The counterexample and detector laws have distinct scopes.
