# Sparse native-coarse path window with an error certificate

This prototype runs the local-coarse modular one-load instrument on the actual
canonical ring geometry without enumerating its complete modular Gauss sector.
The source and complete output are `oph_openwave/sparse_native.py`; a complete certificate can be emitted with
`python -m oph_openwave.sparse_native --full-ledger --output receipt.json`. Run:

```sh
python -m oph_openwave.sparse_native
python -m pytest tests/test_sparse_native.py
```

Only an explicitly requested output file is written. No Fable dependency is imported.
No reserved acceptance or integer-cutoff configuration is executed.

## Declared finite instrument

Enumerate the radius-R ball from port 0 with zero records in the half-hop state
graph, with registers modulo N. Project H_m onto this ball. Keep its hopping
matrix sparse and evolve the initial vector with `expm_multiply`.

For every seam, C/E law and fair coin, group sources with the same native local
outcome into a single partial-permutation Kraus operator. C null outcomes retain
coherence between all states with that outcome, including transit states. This
is the named local-coarse quantum lift, a supplied extension of the classical
record rule. It does not measure global position or all records.

A native jump whose target is absent from the ball receives a separate retained
`boundary_rejected` outcome and leaves the source state unchanged. Its native
record is null and its applied current is zero. This rejection branch is a
regulator, not a native outcome. Splitting it from accepted jumps avoids a
noninjective sum of identity and translation in one Kraus operator. Every source
appears exactly 4M times with squared amplitude 1/(4M), and every branch map is
injective. Thus sum K†K=I. Removing the rejection branches loses trace, which a
persisted test detects.

The readout reports native rejection probability, outer-shell population and
the norm of the omitted coherent derivative QHP ψ at the end of the interval.
The latter two are diagnostics, not finite-time coherent loss probabilities.
Coherent projection error is bounded separately below.

## One-interval bound

The reference is the COMPLETE modular Gauss sector at the same N and the same
initial basis state. No finite ring serves as an infinite reference, and no
quadrature or unknown reference state enters this bound. Let d_max be maximum
port degree and

    h = sqrt(m² + 2 d_max),       x = |t| h,
    T_j(x) = exp(x) x^j / j!.

Both complete and projected generators have norm at most h. The bipartite
incidence block has absolute row sum at most d_max and absolute column sum at
most 2, and its mass term anticommutes with the off-diagonal block. Powers of
complete and projected H applied to the initial basis state agree through R.
Two Taylor tails therefore bound their coherent state-vector difference by
2 T_(R+1). Their pure-state trace distance is at most

    delta = min(1, 2 T_(R+1)).

The projected coherent vector has norm at most T_(R-1) outside radius R−2:
its Taylor polynomial through degree R−2 has no support there. A native jump
spans at most two half-hops. Consequently the complete native channel and the
window channel agree on every vector supported inside radius R−2. If psi is
the projected coherent state and chi its normalized projection to that inner
ball, their pure-state trace distance is the outside norm. Two applications of
channel contractivity and a triangle inequality give

    D(C_full(psi), C_window(psi)) <= 2 min(1, T_(R-1)).

For commit probability p, the complete coherent-then-event instrument and the
window instrument thus satisfy

    D(output_full, output_window)
       <= min(1, delta + 2 p min(1, T_(R-1))).

For R<2 use outside norm 1. The same argument applies to the classical-outcome
extension of the instrument, with zero reference weight assigned to the
additional rejection label. Here D is half the trace norm. Every final event
probability differs by at most the bound.

The implementation evaluates an upper bound on h using an exact rational dyadic
square-root majorant. It upper-bounds exp(x) with a rational Taylor sum and a
geometric tail: after term k, all subsequent term ratios are at most x/(k+2)<1.
All arithmetic in the bound is `Fraction` arithmetic for the exact binary input
parameters. The final float is rounded upward. This certifies the mathematical
projection/instrument error majorant. It does not certify floating-point error
in SciPy's propagator; numerical comparisons and trace controls are reported
separately.

This is an initial-basis-state, one-interval certificate. Repeated commits change
the preparation, so repeatedly adding this particular initial-state bound would
be invalid. A multi-interval state-dependent or support-expansion estimate is
required. The Bernoulli event probability is supplied; no native clock derivation
or equality to the continuous per-seam GKSL evolution is inferred.

## Executed checks

The complete triangle modulus-6 model is the exact finite control, with the
previously explored m=.23, t=.37, p=.29. Radius 2 uses five states and differs
by 0.01078392192 in trace distance, below bound 0.8701541577. Radius 4 uses
nine states and differs by 0.00007431167405, below bound 0.09219419153.

The actual `flat_ring(1,1)` has 24 ports and 63 seams. Its radius-4 modular
window uses 55 states. One interval at m=.23, t=.1, p=.29 gives a mathematical
trace-distance bound 0.005842811201 against its complete modular sector.
The retained native-jump rejection probability is 3.80651e-12. The small
rejection probability alone would not justify a small coherent error.

The tests cover: Kraus completeness with rejecting lost-branch control,
coarse null-outcome coherence, exact-full comparison and actual-ring execution,
and invalid-bound-parameter/zero-interval controls. No physical row is completed.

## Certified numerical output without normalization

`certified_interval()` evaluates the sparse Taylor action with exact complex
rational coefficients. The mass, time and probability inputs are treated as
their exact binary dyadic values. The discarded state-vector tail epsilon is
bounded by the same rational exponential majorant. On the actual 55-state ring
window, order 11 gives epsilon ≤9.05707e-15.

The polynomial vector phi is not renormalized. Its exact squared norm differs
from one by −3.58519e-16; this defect is retained. Applying the finite Kraus maps
to phi gives an exactly trace-accounted subnormalized positive output. Every
position and ledger probability is a rational sum of squared amplitudes, with
its float conversion error reported separately.

Let A=|phi><phi|−|psi><psi|, where psi is the exact projected unitary vector.
Then ||A||_1 ≤ epsilon(2+epsilon) and Tr A is the exact norm defect d. For a
probability effect, including after a trace-preserving quantum instrument,

    |Tr(E A)| <= (||A||_1 + |Tr A|)/2
              <= [epsilon(2+epsilon)+|d|]/2.

The |d| term is necessary: half the trace norm alone does not bound effects for
a trace-changing approximation. A persisted scalar counterexample rejects
that incorrect rule. The one-interval regulator bound is added to this numerical
effect bound. Exact rational endpoint calculations and outward float conversion
produce intervals for every full-modular position and outcome probability.
For the ring control, the numerical effect error is ≤9.23633e-15 and the
combined error is ≤0.005842811200167. These are complete mathematical-map and
arithmetic certificates; no normalized-CPTP property is asserted for the
subnormalized polynomial itself.

The Taylor helper is private: its norm certificate requires the generated
projected incidence Hamiltonian and matching mass/degree metadata. Public
`certified_interval()` supplies those from its own declared construction; the
certificate presumes its Hamiltonian/state maps are unmodified. Source hashes
bind the implementation. The default runner returns compact controls and
certificate summaries; `--full-ledger` retains every position/outcome interval.
