# Cell-observer entropy converges to the full regional entropy

This theorem concerns the actual finite-entropy continuum projector P in
`../gravity/continuum_source/`. It also applies to the flattened time
generator because that generator has the same projector. It supplies a
previously missing qualitative observer-refinement limit. It does not
give a numerical error for the present cell size or commute this limit
with a shape derivative or a large-ball limit.

## Assumptions and their source

Let chi be multiplication by a fixed bounded measurable spatial region B.
Write P_beta for the source's local spin reference and D=P-P_beta. The
continuum-source proof, section 4, establishes the one-sided localization

    chi D belongs to S_q for every q>1.

It uses D:L2 -> H3, a compactly supported smooth cutoff equal to one on B,
and the localized Sobolev embedding. This is stronger than only assuming
chi D chi is in S_q. Its adjoint D chi is in the same ideal. On a compact
source the global estimate already suffices.

Let Q_n be finite-rank orthogonal observer projections such that

    Q_n <= chi,   Q_n -> chi strongly,   [Q_n,P_beta]=0.

Nested normalized scalar cell-indicator modes, with all spin components,
satisfy these conditions when their union is dense in L2(B). For a fixed
union of cubes one can repeatedly subdivide every cube. For a fixed ideal
ball one can use nested dyadic cubes wholly contained in that ball; their
union covers it up to its measure-zero boundary. Merely selecting cell
centres near a sphere can include cells outside B and does not establish
these hypotheses for the ideal ball. The digital-shape error is separate.

## Off-diagonal convergence in a Schatten ideal

The actual finite observer covariance is C_n=Q_n P Q_n on Q_n H. Put

    A_n = Q_n P (1-Q_n) = Q_n D (1-Q_n),
    A   = chi P (1-chi) = chi D (1-chi).

The reference terms vanish because the observer projections commute with
P_beta. For any fixed 1<q<2,

    A_n-A = (Q_n-chi) D (1-Q_n) + chi D (chi-Q_n).

In the first term, (Q_n-chi)D=(Q_n-chi)chi D. A uniformly bounded operator
sequence converging strongly to zero multiplies a fixed S_q operator to
zero in S_q norm. This follows first for finite-rank operators and then
by their S_q density. The right-multiplication version follows by taking
adjoints, since chi-Q_n is self-adjoint. Both terms therefore tend to zero:

    ||A_n-A||_(S_q) -> 0.

Also ||A_n||_(S_q)<=||chi D||_(S_q). This retains the contribution from
unobserved modes within B as well as modes outside B; no cell complement
has been discarded from the covariance.

## Entropy is continuous along these particular approximants

Because P is a projection,

    C_n(1-C_n)=A_n A_n^dagger,
    C_B(1-C_B)=A A^dagger.

Every singular value s of A_n or A lies between zero and one half. Define

    g(s)=h((1-sqrt(1-4s^2))/2),
    h(t)=-t log(t)-(1-t)log(1-t),   g(0)=0.

The symmetry h(t)=h(1-t) gives exactly

    S(Q_n)=sum_j g(s_j(A_n)),   S(B)=sum_j g(s_j(A)).

This identity includes multiplicities. It is not an approximate
replacement of entropy by a quadratic functional. As s approaches zero,
g(s)=O(s^2 |log s|)=o(s^q) for q<2.

For completeness, choose M bounding all the S_q norms. Then
s_j(A_n)<=M j^(-1/q), uniformly in n, and the same holds for A. Given
epsilon>0, choose delta so that g(s)<=epsilon s^q for 0<=s<=delta.
Choose a finite J with M J^(-1/q)<=delta. The tail after J is at most
epsilon M^q for every n and for the limit. Each of the finitely many
leading singular values converges because S_q convergence implies
operator-norm convergence. Continuity of g gives convergence of their
finite sum. Letting epsilon tend to zero proves

    lim_n S(Q_n) = S(B) < infinity.

No monotonicity of S(Q_n) is asserted or needed. Entropy can decrease when
additional observed modes recover correlations with previously traced-out
ones. The finite-observer bulk contribution is included and converges to
the full-region value at fixed B.

## Exact scope of the result

This establishes the observer-cell limit for a fixed source and fixed
region under the stated projections. The existing source proof supplies
the needed one-sided localized ideal estimate. The gravity worker gave a
second mathematical reading and confirmed both that input and the need
to distinguish a fixed digital cube union from an ideal sphere.

The theorem supplies no computable S_q tail constant for a particular
cell size, no uniform large-region bound, no interchange of derivative
and observer limit, no area coefficient, and no matter/geometry selection.
`check.py` tests the exact projection/entropy identities on a complete
finite source; it is not a numerical proof of the continuum convergence
rate. These distinctions are retained when using this result in the
Einstein construction.
