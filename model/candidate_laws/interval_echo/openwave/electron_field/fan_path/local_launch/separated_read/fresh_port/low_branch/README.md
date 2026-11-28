# Low energy on the four-fan charged guide

The narrowest resolved pair below the low-energy cutoff carries very little
material current. Its 1010 charged-face read is separable for every
wait, even though the broader 1011 detector can recover a small
nonseparable comparison after a coherent energy admission. The neighboring
pair has a nonseparable 1010 record, with a full-rank local channel that
rules out exact recovery of an arbitrary qubit at that port. These
results use the unchanged 3,520-row written-field comparison law on the
320-cell Interval Echo board. [Protocol](PROTOCOL.md) · [Records](records.json)

## Energy count and material motion

The dyadic cut $81/16$ has exactly 14 eigenvalues below it. Disjoint
integer residual enclosures give one simple line per interval. This
extends the archived two-line count to the dense low part of the spectrum.
The first six line centers, with the lowest numbered zero, are

| Line | Energy in native generator units |
| ---: | ---: |
| 0 | 3.5471620199 |
| 1 | 3.9054660831 |
| 2 | 4.5749255577 |
| 3 | 4.5787679077 |
| 4 | 4.6003260107 |
| 5 | 4.6086276679 |

The pair of lines 2 and 3 is split by approximately $0.00384235$.
A projection onto these lines is an ideal global diagnostic. In any
normalized input admitted to this band, the probability of an actual
charge in the source cell is at most **11.5454%**. For a cut after fan
$j$, let $P_j$ be the signed-load location flag on the later part of
the path. The model continuity equation uses $i[H,P_j]$. The directly
reconstructed bridge head gives the same signed change and flips when
material crosses. Inside a two-line band the current has zero diagonal
and off-diagonal magnitude
$|E_1-E_0|\,|\langle e_0,P_je_1\rangle|$ in the band labeling.
The three maximal current magnitudes are approximately
$2.4603\times10^{-6}$, $2.0404\times10^{-6}$ and
$1.3264\times10^{-6}$ in material units per native clock unit.
At any fixed admitted input, the load probability on the far side of
each cut can swing by at most 0.1281, 0.1063 and 0.0691 percentage points
respectively over all waits. Each individual energy line is stationary and carries zero directed
bridge current. Conjugating the charge complements the entire head
record and reverses the signed current convention.

With the original localized 98 source there is **zero** raw remote arrival
at zero wait. Inserting the ideal 2–3 energy pass before the 1010 read makes
a remote click possible even at that zero wait, at about **0.17714%** of
all starts. Applying the same global energy pass after a 1010-first
charge read gives exactly zero support at zero wait. These are distinct
measurement orders. No low-energy probability produced by the global
projection is counted as autonomous transport from a sharp sender.

## Local receiver and reference

Each local detector accepts the named signed material address and retains
all compatible local words. All other cells are traced. The two words
at 1010 are the same cell-owned face compared in the charged timing test.
The 1011 read is its four-word two-face detector with both local analyzer
bits. The sender is the same local two-column preparation at 98 with an
untouched reference. The band Bell columns used for channel comparison
span the two resolved energies and are not assumed to be a local
preparation. [Charged receiver instruments](../README.md)

For lines 2 and 3, the band Bell bank has 1010 acceptance about **13.86%**
and a strictly positive two-qubit partial-transpose minimum. The source
to band matrix is invertible. Free evolution is diagonal in the band,
so the source and any free wait act as an invertible congruence on the
reference in the local comparison. Every wait gives the same separability
verdict at 1010. For the 1011 two-face detector the corresponding
partial-transpose minimum is around $-1.1852\times10^{-7}$ per band Bell
input. Its sign also survives the reference congruence at every wait
after a successful ideal global band admission. Its per-original-source
negative magnitude exceeds $3.3802\times10^{-12}$. The receipt gives
local Pauli moments and integer trial witnesses for the band-Bell
comparisons, with the source reference analyzer last.

The source loading for this pair has average pass probability **1.26024%**
of original starts. The two singular-value-squared source weights are
about 0.02518 and 0.00002859. The population that misses the global
diagnostic and the population that passes it but misses the named receiver
remain in separate bins. At 1010, independently inputting the two fine
energy lines gives conditional local head distributions with trace
distance less than $8.20\times10^{-5}$ after registering arrival.
The corresponding 1011 bound is below $9.32\times10^{-4}$.
These bounds limit the best binary guess from any receiver-only head
measurement with equal priors in the accepted cohort. Resolving the
energies themselves would erase the two-line reference coherence.

For the adjacent pair 4 and 5 the 1010 two-word detector is NPT after the
ideal coherent band admission, with a band-Bell minimum near
$-5.5075\times10^{-4}$. Average admission on the original source is
**0.243610%**. Even this favorable local result cannot give an exact
reversible transmitted qubit. The accepted local Choi matrix is
positive definite. For any nonzero receiver recovery map, its positive
lower operator bound implies a post-recovery Choi rank of at least two.
The identity qubit channel would require rank one. The argument permits
local heralded filters. It uses no remote read. The same full-rank obstruction applies to the closer pair and to the
original lowest pair at 1010. The 4–5 candidate has at most 18.765%
source-cell loading within the normalized band, and the largest of its
three maximal cut currents is about $2.3295\times10^{-4}$.

Admitting all fourteen lines would keep **9.69562%** of the
original launches under the ideal global diagnostic. At the held
late clock 31.0625 it gives **0.33529%** joint original-start
acceptance at 1010 with a positive partial transpose, and **0.25353%**
at 1011 with a negative partial transpose. The bypassed original
1011 read accepts **3.17714%** at that clock. The subtraction of
records before and after global spectral admission is not a
decomposition into physical arrival histories.

## Energies and turns

A turn preserving every resolved energy line can only change its phase,
since all fourteen slots have rank one. The two orders of any pair of
such turns agree. In cell 84 the only proper pose preserving both
active physical aperture positions is the identity. The boundary
conclusion does not exclude driven operations or changes of apparatus.

At 1010 the cell-owned face interchange and the head-phase operation
are valid detector actions. Applying them to the whole guide on the
admitted lines changes its energy reports. On lines 2 and 3 the
face operation sends about $8.14\times10^{-5}$ and
$9.16\times10^{-5}$ outside the pair, and the head phase sends about
26.69% and 26.01% outside. Using both to read a local face therefore
supplies no degenerate, energy-preserving spin rotation. The band
calculation places a joint restriction on energy, local information
and actual material current. It has no external rest-mass or angular
calibration.

The executable checks the inertia on the exact half-integer generator
using a factorization and inverse residual. It encloses the normalized
integer trial vectors by their gap to the rest of the spectrum and
propagates those errors through the local read and phase evolution.
A separate sparse evolution reproduces the raw source's late record.
All charged-sector statements have an exact opposite-sign check.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.fan_path.local_launch.separated_read.fresh_port.low_branch.measure
~~~
