# A native measurement that leaves all six outcome spaces intact

One excess occupies one of twelve data ports on the 20-cell Interval Echo board.
A vacancy at port 239 balances it. Four neutral pointer dimers start in their
zero orientation. The first reader publishes which of six adjacent address pairs
contains the excess. Within each pair, either coordinate preparation and their
uniform superposition return to their respective preparation with probability
one. A reader that also writes the low address bit gives the same six effects
and the same coordinate returns. Its uniform-state return probability is one
half in every pair. [Measurement and source bindings](measure.py) · [Recorded counts](results.json)

## The two instruments

The native material copy preserves the control load, restores its host and
toggles a neutral pointer if the control load differs from one. Each of the six
control/pointer input words acquires the same phase
\(h=\exp(-133\pi i/260)\). Routing conjugates this operation by ternary swaps.
Every intervening port returns to its input value and the vacancy stays at port
239. The exact local pulse identity and the numerical checks on all 81 output
coordinates are part of the address reader. [Material copy and routing](../full_address/README.md) · [Symbolic local identity](../full_address/test_address.py)

The first instrument executes the 14 copies for the three high bits. For an
excess at address \(a\) it leaves the data word intact and writes
\(2\lfloor a/2\rfloor\) on the four physical pointers. Let
\(Q_a=|a\rangle\langle a|\) and \(P_j=Q_{2j}+Q_{2j+1}\).
The complete outcome map is \(\rho\mapsto P_j\rho P_j\), since its only
nonzero branch is \(h^{14}P_j\). All six projectors have rank two and sum
to identity on the source code. The 17 raw labels retain unused words and the
malformed-pointer outcome, which have zero probability on the enrolled domain.
These zero-effect outcomes have no normalized outcome-certain input state.

The control executes all 20 native copies. The public reader again publishes
only the pair. Summing over the fine record gives
\(\rho\mapsto Q_{2j}\rho Q_{2j}+Q_{2j+1}\rho Q_{2j+1}\).
The phase from the 20 copies is common to every path. Each outcome effect is
still \(P_j\). This controller uses the same board, pulse library, initial
pointers and reported outcome alphabet.

## Preparation, return measurement and finite-probe theorem

The source starts with the excess at port zero. A full exchange prepares the
other coordinate when requested. For the coherent input, a half exchange on
(0,1) and a degree-reference pulse of duration \(\pi/16\) at port zero
prepare \((|0\rangle+|1\rangle)/\sqrt2\) up to a common phase. Routed
transpositions carry those two addresses into the selected pair. After the
instrument, inverse routing returns them for analysis. The coordinate reader
measures position. The coherent reader uses the independent X-analyzer sequence
of a reference, half exchange and reference with durations
\(\pi/16,\pi/2,3\pi/16\). Its port-zero effect is the uniform projector.
The exact two-coordinate Hamiltonian has symmetric and antisymmetric
eigenvalues 8 and 9. The reference costs on its two coordinates are 32 and 8.
These spectra give both pulse identities by direct multiplication. [Exact return tests](test_return.py) · [Native analyzer](../phase_table/acquire.py)

The terminal full-address read acts on the same four pointers. If the first
read left word \(r\), the second one reads \(r\mathbin{\mathrm{XOR}}a\)
for a terminal address \(a\). The detector recovers \(a\) by XOR of the
two local readings. This recovers the same address effect for every first
pointer word in the experiment, including the control's fine records.
All 294,912 simulated starts are included in the 36 joint histograms. Each
joint law is rounded to integer weights on the \(2^{-40}\) grid, with the most
probable cell taking the complement, and every start compares one PCG64
integer below \(2^{40}\) with those weights. The rounding stays inside the
\(3\times10^{-12}\) norm budget of the propagation, and the counts replay
identically on every platform. The tested instrument returned all 8,192 starts
at every setting. Uniform returns under the control ranged from 4,066 to 4,162
per setting.

On each two-dimensional outcome space the admitted map has sure effect.
The two coordinate returns and the uniform return therefore discharge the
three exact probability hypotheses of
`finite_pure_return_probabilities_select_lueders_one` with \(n=2\). The
theorem selects the identity channel on the block for every matrix, including
every mixed state. The source pulse identities supply the exact premises.
The finite counts check the apparatus implementation. This realizes PR-131
for every populated outcome of one native projective instrument on the
specified protected code. The fine-read control violates the same premise
under an equal-effect implementation. The calibrated pulses, isolation,
sharp pointer readings and Born sampling rule belong to the apparatus.
[Lean finite-probe theorem](../../../../../Lean/EventAlgebra/FiniteProbeInstrumentSelection.lean)

From the repository root, replay with one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
.venv/bin/python -B -O -W error -m model.experiments.interval_echo_consensus.signed_effects.pair_return.measure
.venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/signed_effects/pair_return
```
