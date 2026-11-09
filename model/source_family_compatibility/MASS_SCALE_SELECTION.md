# Fixed-source scaling and the massive-spectrum boundary

## Actual source and reference

For the centered full-CAR source, keep c-normalized coupling U/c=8 and
one chosen normalized thermodynamic ground prescription fixed. Use its
actual signed original CAR spectral measure, including all addition and
removal weights. Let nu_k be the canonical 16x16 positive matrix measure
of total mass I constructed in SPECTRAL_RESIDUE/NODAL_CHANNELS. Take
the positive scalar trace measure tr(nu_k)/16, denoted nu_k below.
Any nonzero positive original matrix response is visible in this trace.

Choose one node k_* and one consistent local Bloch chart, with the specified
unitary boundary gluing. Set c=v/a, U=8v/a, with v>0 fixed; physical
momentum relative to the node is p=(k-k_*)/a and physical signed energy
is E=v epsilon/a. This exact homothety uses a fixed reference; it does
not cover an independently varying selected state, coupling or clock.

Define positive measures, on every compact momentum region for all
sufficiently small a,

 M_a(f)=integral dp integral nu_(k_*+ap)(d epsilon)
                          f(p,v epsilon/a),

for f in C_c(R^3 x R). The energy measure is signed in its energy
coordinate, not in its weights. Every weight is nonnegative.
The trace makes the construction independent of an internal Bloch gauge.
Matrix-valued convergence gives the same scalar conclusion; no presumed
quasifree or one-particle dynamics is used.

The all-volume source identity is literal:

    K(c,8c)=c K(1,8), c>0.

It includes the centered chemical term. Each finite full-ground projector
and its normalized density matrix are unchanged by positive multiplication.
Choose the same thermodynamic ground subsequence, exterior convention and
symmetry average for every c. The resulting GNS dynamics has the same
reference and generator multiplied by c; all original CAR weights stay
unchanged, while excitation energies are multiplied by c. This consumes
`charge_scaling.py` and the canonical momentum measure in
`SPECTRAL_RESIDUE.md`; it does not select a different ground at each scale.

## Exact dilation and joint vague limits

Let D_lambda(p,E)=(lambda p,lambda E), lambda>0. Changing p to lambda p
in the integral gives

 M_a(f o D_lambda) = lambda^(-3) M_(a/lambda)(f).     (1)

Equivalently for Borel sets contained in the local chart domains,

 mu_a(D_lambda A) = lambda^3 mu_(lambda a)(A).        (2)

The a/lambda in (1) and lambda a in (2) are intentional. A convention
using the inverse dilation interchanges them. The pointwise fibre identity
is F_a(lambda p,lambda B)=F_(lambda a)(p,B).

For a compact set C contained in a bounded momentum box P,
mu_a(C)<=volume(P), uniformly in a, because every nu_k has mass one.
Thus the measures are locally uniformly bounded, even when part of the
original spectral probability escapes to infinite physical energy.

Assume the FULL joint vague limit M_a -> M as a decreases to zero exists.
It is a positive Radon measure, possibly zero and possibly deficient in
total energy weight. Passing (1) for compact continuous tests gives

 M(f o D_lambda)=lambda^(-3) M(f),
 mu(D_lambda A)=lambda^3 mu(A).                      (3)

The latter follows as equality of Radon measures, so it applies to Borel
sets, including a singular shell. No continuity of individual energy
fibres or of a pole residue is required. No exchange of a small-frequency
limit with a momentum limit occurs in this argument.

## Positive-mass shell weight

Fix any proposed limiting relativistic speed w>0. For any M>0 define

 H_M={ (p,E): E^2-w^2|p|^2=M^2 }.

Here M is the invariant rest-energy parameter, expressed in energy units.
If w and these coordinates were physically identified, the corresponding
rest mass would be M/w^2. No physical unit conversion is supplied here.

Suppose mu(H_M)>0. Radon inner regularity gives a compact subset
A of H_M with 0<mu(A)<infinity. Choose countably infinitely many
distinct lambda_j in [1,2]. The sets D_lambda_j A lie on distinct
hyperboloids H_(lambda_j M), hence are pairwise disjoint. They are all
inside the compact continuous image [1,2]xA. But (3) gives

 mu(D_lambda_j A)=lambda_j^3 mu(A)>=mu(A)>0,

contradicting local finiteness on that common compact set. Consequently

             mu(H_M)=0 for every fixed M>0.          (4)

This excludes any positive original-weight discrete relativistic massive
particle shell in the full limit, not merely an isolated one. Continuous
timelike spectral weight is compatible with (3); so are massless shells,
zero-energy responses, or the zero limit. The statement does not infer
which of these the actual source realizes. In particular it does not
prove Z positive or zero and it does not identify a physical massless
particle, Lorentz symmetry or an interacting QFT.

For a complementary threshold statement, suppose a nonzero limiting
channel measure is itself dilation covariant and its support lies in the
future causal cone E>=w|p|. The support is dilation invariant. For any
support point, its invariant mass sqrt(E^2-w^2|p|^2) is finite and is
multiplied by lambda. Scaling toward zero proves that the infimum of
invariant mass over this FULL channel support is zero. A positive lower
threshold for that channel is impossible under these additional premises.
If the point is massless the conclusion is immediate. This argument does
not identify an actual causal-cone source channel. A fixed physical-energy
filter may select a channel that is not dilation invariant; the conclusion
then does not apply to that selected submeasure.

Continuous timelike spectral weight can remain, with invariant masses
approaching zero. Neither (4) nor the full-channel threshold result is a
theorem excluding every unstable resonance, an analytically continued
resonance pole, or every dressed charged excitation. Those objects require
additional definitions and source estimates. No delta-function pole is
assumed as a necessary description of every physical particle.

## Discrete refinement and subsequences

Suppose a_n decreases to zero and a_(n+1)/a_n tends to one, and suppose
M_(a_n) converges vaguely along the WHOLE sequence. For arbitrary b->0,
choose a_n>=b>a_(n+1) and put lambda=b/a_n. Then lambda->1, and (1)
gives

 M_b(f)=lambda^(-3) M_(a_n)(f o D_(1/lambda)).        (5)

The transformed tests converge uniformly to f, with supports in one
compact set. The preceding volume(P) bound controls their error. Thus
M_b -> M along arbitrary b->0. Full convergence on such a discrete
family implies (3)-(4). This includes a_n=1/n and a_n=1/n^2.
No additional source current or energy convergence rate is needed for
this interpolation: it is only the exact simultaneous rescaling (1).

For a geometric sequence a_n=r^n with 0<r<1, or a sequence with
r_n=a_(n+1)/a_n -> r in (0,1), whole-sequence convergence gives a
precise intermediate conclusion. Exact rescaling gives

    M_(a_(n+1))(f)=r_n^(-3) M_(a_n)(f o D_(1/r_n)).

Uniform compact mass bounds permit r_n->r in the tests. Consequently
mu(D_r A)=r^3 mu(A), and iteration gives this covariance for every
integer power r^j. This is discrete covariance, not (3) for arbitrary
positive lambda. A positive-weight mass-M shell then requires positive
weight on every shell r^j M, j in Z. Such a countable orbit can have
summable weights on a compact region, so the continuous-dilation
contradiction does not apply. Any nonzero invariant future-cone channel
has invariant-mass infimum zero, by taking j->+infinity. An
isolated single massive species without this orbit is inconsistent with
the discrete covariance. This conclusion assumes convergence along the
WHOLE geometric sequence, not merely a sparse subsequence of it.

For an arbitrary sparse subsequence a_n, its limit need NOT satisfy (3).
The sequence lambda a_n can have a different limit; exact scale covariance
then relates different tangent measures rather than making one tangent
measure dilation invariant. An assumed existence of one subsequential
effective continuum therefore does not by itself yield (4). A restriction
only to isolated finite regulators also has no scaling-limit conclusion.
The source bounds do not prove the full joint
vague limit, so this assumption must remain explicit.

## Source-family selection

For this exact fixed-ratio, fixed-state, unrenormalized original CAR
readout, deciding the nodal atom Z is sufficient to decide bounded-window
survival for unit packets supported O(a) from the node via NODAL_LOCALITY, but it cannot produce a positive discrete
mass spectrum in a full ballistic limit. If that full limit exists and
Z>0, any surviving original massive shell is excluded by (4).

A proposed path to finite positive masses must therefore justify at least
one departure from these hypotheses: source-selected running U/c or other
dimensionless couplings, a changing normalized reference with controlled
limits, a different source-selected clock/scale dictionary, a justified
finite physical cutoff, or a separately reconstructed sector with its own
operational map. Merely multiplying masses or inserting another reference
does not provide such a derivation. A sparse tangent limit would require
its own stability/universality argument and is not a unique full continuum;
a whole geometric limit must retain its discrete mass-dilation orbit.

An arbitrary a-dependent field or probability renormalization is not
included: the present measures retain original complete probabilities.
Even a conventional homogeneous field normalization changes the degree
of (3), rather than creating a fixed positive discrete mass scale.
General nonhomogeneous redefinitions require a separate source derivation.

This theorem concerns one defined family and observable class. It does
not exclude other justified running families or collective even sectors.
No physical species or mass is assigned.


## Executable controls and remaining premises

`mass_scale_selection.py` reuses the actual full 16-state default dimer
at c=1/8,U=1 and reconstructs its positively rescaled source. It verifies
the full matrix identity, the common exact ground projector, every source
spectral projector and all charge-energy probabilities. The two creation
branches a* and a form a complete instrument by the CAR, and all energy
records, including zero weights, are retained. This finite control does
not turn a dimer into a three-dimensional momentum mode and does not
claim execution of a large native detector tape. Additional rational
full-ground controls reuse `charge_scaling.homogeneous_dimer`.

Exact rational controls also bracket n^(-power) refinements without
floating roots, verify Lorentz quadratic-form scaling, and lower the
invariant mass of a hypothetical future-cone support point. The separate
regressions integrate a sample joint measure to check the Jacobian and
index a/lambda in (1), rather than only repeating the formula. Such
implication controls are not replacements for the actual ground.

The source action homogeneity and the conditional measure theorem are
proved. Existence of its full joint continuum measure, its nonzero
weight, its causal support and its physical field dictionary remain
unproved. Sparse tangent limits, finite cutoffs and source-selected
running parameters must retain their stated scope. The result explains
why deciding the nodal atom alone cannot finish a positive
massive spectrum on this fixed normalized source path.
