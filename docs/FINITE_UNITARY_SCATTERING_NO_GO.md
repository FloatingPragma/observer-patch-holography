# Fixed-cutoff unitary asymptotics: a no-go and its boundary

## Result

Let $G$ be a Hausdorff topological group. If the natural powers
$g^n$ converge in $G$, then $g=1$. Consequently, for a nonidentity
finite-dimensional unitary matrix $U$, the sequence

\[
U,\;U^2,\;U^3,\ldots
\]

has no ordinary limit in the inherited matrix topology. In finite dimension,
entrywise and matrix-norm topologies agree, so changing between the usual
finite-dimensional norm choices does not evade the obstruction. Nor does
asking for full weak-operator convergence at that same fixed finite dimension:
the standard finite-dimensional linear-operator topologies coincide.

The kernel-checked ambient-matrix statement is
`finite_unitary_ambient_powers_have_no_limit` in
`Lean/QFT/FiniteUnitaryScatteringNoGo.lean`. Its abstract carrier theorem,
`tendsto_powers_forces_identity`, applies to every Hausdorff topological
group. It uses no OPH axiom, physical calibration, target value, or numerical
fit.

## Proof mechanism

Assume $g^n \to a$. Removing the first term from a convergent sequence does
not change its limit, hence $g^{n+1} \to a$. Continuity of multiplication
also gives

\[
g^{n+1}=g^n g\longrightarrow ag.
\]

Hausdorff uniqueness yields $a=ag$. Group cancellation then gives $g=1$.
The proof depends only on the topology and group law.

## Consequence for a scattering architecture

An exact nontrivial unitary update at a fixed finite cutoff cannot produce a
scattering carrier by letting that update itself settle to a limiting
operator. An architecture using such an update needs some additional
operation, for example comparison with a reference or free evolution, a
selected projected scalar or observable readout, an infinite-dimensional weak
limit reached through a continuum or infinite-volume limiting procedure, an
open-system channel, or a finite-time recurrence-aware operational readout.
The theorem identifies this requirement without choosing among those routes.

This is a no-go for direct convergence of one sequence $U^n$. It is not a
no-go for scattering. In particular, wave operators use relative evolutions
of the form $U(-t)U_0(t)$. The Lean theorem
`identical_relative_evolution_tendsto` records an exact scope control:
$(U^n)^{-1}U^n=1$ converges for every $U$, including nonidentity $U$
whose own powers cannot converge. Comparison limits for distinct dynamics
require additional hypotheses and are outside this result.

The theorem also makes no claim about asymptotic states, convergence of a
selected scalar or projected observable, infinite-dimensional weak limits,
subsequential or Cesàro limits, resonance poles, cross sections, the optical
theorem, renormalization-group flow, or the existence of a continuum quantum
field theory. It does not discharge the interacting QFT, renormalization, or
scattering lane by itself.

## Independent exact replay

`code/qft/finite_unitary_scattering_no_go.py` performs a dependency-free,
integer-arithmetic replay on the identity, a swap matrix, and a quarter-turn
matrix. For each orthogonal step it checks

\[
U^{n+1}-U^n=U^n(U-I)
\]

and verifies that the squared Frobenius displacement remains exactly
constant. The two nonidentity controls therefore fail even the Cauchy
condition. The replay also checks that identical relative evolutions are the
constant identity and rejects a nonorthogonal shear. These finite controls
exercise the proof mechanism and its scope; the general proof is the Lean
theorem.

Run the checks from the repository root:

```bash
python3 code/qft/finite_unitary_scattering_no_go.py
python3 -m pytest -q code/qft/test_finite_unitary_scattering_no_go.py
cd Lean && lake env lean QFT/FiniteUnitaryScatteringNoGo.lean
```
