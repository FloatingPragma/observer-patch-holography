# Fixed-source quantum geometry refinement

## Source and observation domain

The carrier is the full tensor product of two native eight-level cells,
two geometry oscillators with occupations0 through4, and the field vacuum
plus all49 modes of the seven-by-seven Dirichlet scalar grid. Its dimension
is80000. Native energies, the rational harmonic interpolation F, induced
matrices S=FᵀLF and M=FᵀF, stiffness100 and oscillator frequencies are the
literal inputs of `openwave_dynamic_wave/source.py`. No measured gravity
parameter enters this construction.

Write X=a+a† on the complete five-dimensional oscillator and
Q₀=c₊X₊+c₋X₋, Q₁=c₊X₊−c₋X₋. The source is

H = Hgeometry + exp(Q₀)K₀ + exp(Q₁)K₁
    + [0 direct-sum sqrt(exp(2FQ)Lexp(2FQ))].

The harmonic oscillator boundary completion is Hgeometry=ω₊N₊+ω₋N₋.
All64 native sectors and both field number sectors reduce H. Each term is
positive. The exact unique global ground is native00, geometry vacuum and
field vacuum. A wave-sector ground certificate concerns its stated native
sector, rather than the global vacuum. The scalar extension is number
conserving; canonical time-dependent covariant field quantization is not
asserted.

The coordinate roots are0 and ±sqrt(5±sqrt(10)). Their eigenvectors follow
the Hermite recurrence, normalized with positive first coordinate.
Interval arithmetic gives |Q₀|,|Q₁|<0.19. Positivity of F and row sums at
most one imply the same bound at each grid site. The resulting bounds are

Hgeometry <9, Hnative <8, Hwave <9/2, hence0≤H<32.

A complete source read can therefore change mean source energy by at most32.
The preparation, setting registers and record-copying apparatus are
specified finite services. This bound does not supply their physical energy
or time calibration.

For c=0,1,2,3,4, let P_c retain both geometry occupations at most c and
leave native and field coordinates unchanged. The comparison generator is
P_c H P_c. These are compressions of one fixed construction. They are not
obtained by replacing Q with P_c Q P_c inside the exponential or square root.
The distinction is strict: the native10 geometry-vacuum diagonal element
of the master compression exceeds that of the independently exponentiated
cap-one construction by at least

E₁(c₊⁴+c₋⁴)/12 >0.

The two additional closed four-step oscillator walks through occupation2
give the displayed coefficient. The exponential series has nonnegative
closed-walk contributions here, so higher terms cannot cancel it. This is
a source-identity countercontrol, without a claim that either declared
architecture is inadmissible.

## Positive roots and finite arithmetic

Each commuting geometry coordinate block has K=exp(2q)Lexp(2q).
The grid bound λmin(L)≥8/64 and |q|<0.19 give λmin(K)>1/25.
For the emitted dyadic symmetric matrix A, exact integer multiplication
and interval K entries enclose r=||A²−K||. A positive Gershgorin lower bound
α for A gives

||A−sqrt(K)|| ≤ r/(α+1/5).

This follows by solving A(A−sqrt(K))+(A−sqrt(K))sqrt(K)=A²−K,
using the positive Sylvester inverse. All25 root blocks are retained.

The one-mode eigenframe is enclosed entrywise. Its tensor-square frame error
is bounded by2δ+δ² plus the explicitly bounded Kronecker rounding.
For a frame error ε and root norm below9/2, reconstruction error is at most
10ε. The analogous native contribution has norm below8 and error at most20ε.
The executable checks the small-error hypotheses. Each reconstructed entry
sums25 products of three factors. The conservative gamma_200 bound controls
that sum. A matrix of dimension d with entry errors at most e has operator
error at most de. The source assembly and its final2⁻⁴² dyadic rounding are
therefore bounded in operator norm. The primary wave-block bound is below
3.15×10⁻¹⁰.

The arithmetic contract is binary64 round-to-nearest without overflow.
For n scalar operations, gamma_n=nu/(1−nu), u=2⁻⁵³. Matrix products,
Frobenius sums and complex dot products have explicit operation-count bounds
in the producer. Numerical proposals do not certify themselves: root
residuals are checked with exact integer products, and spectral proposals
carry either exact integer residuals or a gamma_n enclosure of those products.

## Complete event propagation

One evolution event is exp(−iH/2). A normalized dyadic matrix B=Hhat/32
has absolute row sums below one, so ||B||≤1. The executable uses the
Chebyshev expansion of exp(−16iB) through degree64. Each Bessel coefficient
J_n(16) is enclosed by its exact rational power series through k=128 and
its first omitted alternating term. The omitted Chebyshev tail is bounded by

6·8⁶⁵/65! /(1−8/66) <10⁻³⁰.

For the computed recurrence, let
ζ=2 gamma_(2d)||abs(B)||+5u and
G=64·65·ζ/2. The exact second-kind Chebyshev bound |U_k|≤k+1 gives recurrence
error at most G/(1−G). The producer checks G<1/1000 and includes coefficient,
summation and source-Hamiltonian errors. The resulting per-event state error
is independent of the supplied normalized vector. Sixteen events have error
at most16ε(1+ε)¹⁵, plus the certified initial normalization error.

The packet is the exact complex integer preparation defined by
`packet_integers`, normalized by its exact squared integer norm. The
geometry starts in its vacuum. The complete experiment retains a half-weight
field-vacuum branch and half-weight packet branch. At each of sixteen events,
receiver x≥5 is read. Every receiver-site/geometry outcome and every final
no-hit outcome is retained. Field-vacuum and unoccupied native-sector outcomes
have their explicit probabilities. No postselection is removed; the reported
mean arrival is labelled as conditional on a hit.

Orthogonal event histories form a Stinespring isometry. Telescoping the
sixteen approximate gates bounds the whole history vector, without adding
probability errors independently over every pixel. A state-vector error δ
gives complete probability L1 error at most2δ+δ². Scalar summation and norm
rounding have explicit additional guards. This controls transverse means,
hit probabilities and conditional arrival means.

## A reference-stable finite compression theorem

Embed each compressed geometry output into the master carrier. The receipt
bounds the norm difference of the complete arrival-history vectors and of
the final field-vacuum vectors. The admitted initial subspace is

span{native00,native10} tensor geometry vacuum
  tensor span{field vacuum, the exact packet}.

The four corresponding columns have orthogonal native/field labels before
any subsequent readout. Their difference columns therefore remain orthogonal
in the common history dilation. The maximum column-norm error bounds the
isometry operator norm on this entire subspace. Tensoring an arbitrary
reference and applying any common complete output instrument preserves the
trace-distance bound. In particular it preserves joint native clock reads;
the comparison is stronger than a match of position marginals.

The complete arrival-instrument total-variation upper bounds are
approximately0.166912,0.022479,0.003021 and0.000358 for nested caps0,1,2 and3.
They are finite comparisons to cap4. They do not prove convergence as the
master cutoff tends to infinity. The source has conditional delay in
[0.0105311,0.0105598] in its event-parameter convention and displacement in
[−0.001776299,−0.001775262] sites at event8.

## Finite Schur and native spectral bounds

Write the exact master wave block as [[A,B*],[B,D]] for P_c and its
complement. Positivity and the oscillator energy imply

D≥bI, b=(c+1)ωmin+(Eleft+Eright)exp(−0.19)+1/5.

The code uses a downward rational enclosure of b and an upward bound β on
||B||, obtained from the exact dyadic row/column sums plus source error.
If a₂ is a certified lower bound for the second eigenvalue of A, the
quadratic-form inequality with a positive auxiliary scalar gives

λ₂(H) ≥ [a₂+b−sqrt((b−a₂)²+4β²)]/2.

This is the finite Schur/min-max shortcut. It verifies the isolation premise
without diagonalizing the1225-dimensional master matrix. The embedded
compressed ground proposal is an exact integer vector with exact rational
Rayleigh quotient and squared residual for Hhat. Source assembly error
widens its Rayleigh and residual bounds. When the upper Rayleigh bound is
below the displayed second-level lower bound, Temple's inequality encloses
the exact ground. For native10 and cap-three trial the interval is
[1.35873578062,1.35873578145]. The native00 one-wave interval is
[0.55176591113,0.55176591179].

Every one of the64 native vacuum blocks has dimension25. Dyadic eigenframe
Gram and residual products give its full spectrum. If the orthogonality
error is η<1/2 and the source/eigenvalue norm bounds are32, polar correction
gives the conservative spectral enclosure2R+128η. These actual spectra,
not a transferred static-lapse ratio, are available to clock comparisons.

The construction supplies no physical G, seconds, photon identity, Einstein
regional premises or infinite geometry/spatial limit. The event parameter
is a convenient label for a finite executable word. Native clock records
supply operational comparisons in `openwave_relational_clock`.
