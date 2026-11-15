# Pure persistence in the finite interacting source

The source is the C10 six-orbital CAR triangle with six modular record
fibres per occupation mask. Its Hamiltonian and native instruments are
specified in [README.md](README.md). The classification under review is
Codex C12 `pure_persistence.py`; the independent executable
[persistence_audit.py](persistence_audit.py) imports only the auditor's
direct Gauss/CAR construction. It compares the reconstructed matrix with
the actual saved Hamiltonian before making a rank claim.

For a unit vector ψ, positivity implies that
`sum_r K_r |ψ><ψ| K_r* = |ψ><ψ|` only if each `K_r ψ` is a scalar multiple
of ψ. A moving C or E branch is an injective signed partial permutation
with square zero: after its endpoint transfer, its original directional
occupation condition fails. Therefore its eigenvalue on ψ must be zero.
There is no cancellation between different input columns because this
partial map is injective. Across both directions and all three seams,
the moving domains cover exactly the configurations with unequal port
occupations. A common native pure fixed state is consequently supported
on the subspace S where the three ports are all empty or all filled.
Conversely every native branch is scalar on S. The independent integer
branch reconstruction checks these statements on all 384 coordinates,
including zero outcomes. Coin, seam and MIX weights only supply the
corresponding positive scalar normalization.

Let P project onto S, `H_S = PHP|_S`, and `C = (I-P)HP|_S`. Any H-invariant
subspace within S is contained in both `ker C` and `ker(C H_S)`.
The independent audit represents every frozen binary64 entry as its
exact rational value and reduces it modulo the prime 65521. Denominators
are powers of two and hence invertible. A nonzero maximal minor modulo
this prime is also a nonzero rational minor. Full column rank of the
stacked matrix `[C; C H_S]` therefore proves that its rational kernel is
zero, without a floating singular-value threshold or a repeated-nullspace
algorithm.

| Source number | Dimension of S | Rank of C modulo 65521 | Rank of stacked constraints |
| --- | ---: | ---: | ---: |
| 0 | 6 | 0 | 0 |
| 1 | 18 | 17 | 18 |
| 2 | 18 | 18 | 18 |
| 3 | 12 | 12 | 12 |
| 4 | 18 | 18 | 18 |
| 5 | 18 | 17 | 18 |
| 6 | 6 | 0 | 0 |

The endpoint sectors are their entire number blocks; number conservation
proves H-invariance there. Thus the maximal H-invariant subspace of S has
dimension twelve, entirely in N=0 and N=6. Exact rational LDL congruence
counts one eigenvalue in each of twelve mutually disjoint rational
intervals. This independently proves that all endpoint energies are
simple and that the two number sectors share no energy. Pure states
stationary under H and the complete native instrument are precisely
these twelve eigenrays. Pure nonstationary unitary orbits can also use
superpositions within the twelve-dimensional invariant subspace.

The former uniform one-fermion transit circulation is a bare-Hm
eigenvector. Under the interacting H its energy variance is
`6.449511111111113`; its port-occupation leakage at time `.37` is
`0.004122491766362765`. These are independent complete-H propagation
comparisons against C12's stored result. They do not measure an absorbing
survival time.

## Algebraic-time extension

Suppose a finite Hermitian H has algebraic entries and t is a nonzero
real algebraic number. Its eigenvalues are algebraic. If two distinct
eigenvalues had the same phase under `U = exp(-itH)`, then
`t(λ-μ) = 2πk` for a nonzero integer k. The left side is algebraic and
the right side transcendental, a contradiction. Thus U and H have
exactly the same spectral projectors. Polynomial interpolation on the
distinct phases expresses H as a polynomial in U. Every U-invariant
subspace is therefore H-invariant.

This applies to the frozen rational C10 H at `.37`, whether this time
denotes the rational 37/100 or the exact binary64 value. It does not
assert equality of those two numerical times, and it does not transfer
to arbitrary transcendental time or an arbitrary externally phased H.

For the selected step with a retained no-commit branch
`sqrt(1-p) U` and native branches `sqrt(p) K_r U`, where `0<p<1`, a pure
stationary output requires Uψ to be collinear with ψ and every K_rψ
to be collinear with ψ. The native argument puts ψ in S; the algebraic
time argument makes ψ an H eigenvector. Hence the same twelve-ray
classification applies to this specified discrete step. A pure orbit
that stays pure at every step must have every U iterate in S. Its finite
U-cyclic span is U-invariant and therefore H-invariant, so it also lies
in the twelve-dimensional endpoint subspace. This is a repeatability
statement, not a claim about a single pure output.

These conclusions concern one exact finite source and its complete
instrument. Mixed stationary states, metastable behaviour, approximate
lifetimes, propagating packets, number conservation and physical
particle identity are separate questions. In particular, the N=3
maximally mixed state remains stationary, as the producer states. No
physical row status follows from this classification.
