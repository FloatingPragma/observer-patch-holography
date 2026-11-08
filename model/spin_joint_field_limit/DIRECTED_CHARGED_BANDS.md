# Directed charged bands on the joint source

The positive-field construction supplies source-even collective charged measurements with
arbitrarily narrow *spatial momentum* support and a uniform original
probability. Their energy intervals remain broad. The finite control below
resolves actual order-three translations on one periodic source and measures
a strictly positive energy width in each charge and direction. Neither
statement identifies a particle mass or a sharp continuum dispersion.
Its internal U(1) charge is not identified with an original odd CAR charge
or a measured Standard Model quantum number.

## 1. Site inversion and the same source reference

Use the state and clock of `PROOF.md`: `p=n^-2`, `B=J n^-5`, followed by
the prescribed positive-field infinite ground and a joint history
subsequence. Inversion of site coordinates `x -> -x` preserves the cubic
exchange bonds and `epsilon_x=(-1)^(x1+x2+x3)`. It therefore preserves
`H_B`, its unique positive-field ground, and the finite local ground
approximations on inversion-symmetric boxes. This is a literal site
permutation; no relativistic parity or transformation law of physical
angular momentum is assumed.

The limiting GNS inversion I fixes the vacuum and commutes with the source
clock and internal U(1). It takes each configuration field with profile F
to the same field with profile `F(-X)`. It reverses the common spatial
translations, obtained from the even microscopic translation subgroup.
Consequently it maps the joint spectral event `(E,k)` to `(E,-k)`.
The canonical field called P in the parent proofs is not spatial momentum.

Let `rho=|k0|>0`, `0<delta<=1/4`, and `r=3rho/2`. Choose a real, even,
smooth Fourier profile supported inside the two disjoint caps

    C_+={k:|k-k0|<=delta rho}, C_-=-C_+,

with L2 norm one. Its inverse Fourier transform F is real, even and
Schwartz. Both caps lie in the parent annulus `r/2<=|k|<=r`. Such profiles
exist for every direction k0 and every positive delta. Shrinking delta
changes F; this is not a single state approaching a plane-wave atom.

## 2. Original weight in each directed record

Use precisely the bounded arctan creation instrument of
`ARCTAN_CHARGED_BANDS.md`, formed with this F. Its two creators B_+,B_-
and positive complete complement have norm bounds and charge covariance
on the same state and clock. Since F is even, inversion fixes each creator.
Let

    a=3/387200, b=387200,
    Pi_sigma=1_[ar,br](L) 1_C_sigma(P_space), sigma=+,-.       (1)

The original proof of the band bound works with `Pi_++Pi_-` in place of
its full annular projection. In fact, the *linear canonical test vector*
P(F)Omega has spatial spectral support in `C_+ union C_-`. Thus the
excluded spatial part contributes zero to its norm; the same first and
inverse energy moment bounds control the excluded energy part. The weak
commutator and charge-pairing argument are otherwise unchanged. It gives

    ||(Pi_++Pi_-) B_q Omega||^2 >= 1/387200, q=+1,-1.        (2)

This does not assert that the nonlinear B_q itself has only those two
momentum supports. Its other momenta and energies are retained among the
complementary outcomes.

Inversion exchanges Pi_+ and Pi_- while fixing B_q Omega. The disjoint
events therefore have equal original probabilities, and (2) implies

    ||Pi_sigma B_q Omega||^2 >= 1/774400                    (3)

for all four choices of q and sigma. The instrument has seven records:
each charge followed by either cap in (1) or the rest of its complete
energy-momentum measure, and the original creation complement. No success
distribution is normalized to one.

The momentum norm in each successful record lies in
`[(1-delta)rho,(1+delta)rho]`. Its energy interval is

    [9rho/774400, 580800rho].                              (4)

Equation (3) is uniform in delta and rho. The spatial profile, cutoff
radius and preparation cost need not be uniform when delta shrinks. The
existing finite-source compiler and cutoff bounds apply to each fixed
profile. They do not provide an effective source-to-continuum convergence
rate or a fixed spatially local apparatus for all delta.

Smooth complete detectors require the parent boundary-weight and
resolution allowances. A complete-output half-trace error eta changes
each record probability by at most eta. A sharp cap is a spectral statement
until that resolution step is supplied; it is not an assertion of an exact
finite-time momentum PVM in the continuum.

The cap has a direct translation readout. For the unnormalized successful
vector `psi=Pi_sigma B_q Omega`, put `w=||psi||^2`. Its actual spatial
translation U(y) satisfies

    ||(U(y)-exp(i sigma k0.y))psi||^2
       <=min(4,delta^2 rho^2 |y|^2) w,
    w-Re(exp(-i sigma k0.y)<psi,U(y)psi>)
       <=min(2,delta^2 rho^2 |y|^2/2) w.                   (5)

These follow directly from `|exp(it)-1|<=min(2,|t|)` on the cap. The
second line is an original, unnormalized controlled-translation
interference score, with the unsuccessful records retained. It supplies
a quantitative spatial phase while making no assertion about a time
phase `exp(-i omega(k0)t)`. The energy width is exactly the remaining
distinction.

## 3. An actual periodic source with directed translation outcomes

The finite control uses six spin-one-half sites on the periodic cycle
`0-1-2-3-4-5-0`, with `J=1`, `n=2`, `p=1/4` and `B=1/32`. Its full
64-dimensional Hamiltonian is produced by
`native_ordered_ground.source.matrix`. The physical Hamiltonian has the
same exchange and staggered-field terms as `PROOF.md`. All matrices in
the calculation are expressed together in the native Marshall frame.
This finite cycle is a source instance, not an exhaustion of the
three-dimensional joint ground.

Write `F_x=cos(2pi x/3)`, and define the contraction

    A=(1/4) sum_x epsilon_x F_x S_x^3,
    B_q=(1/8) sum_x epsilon_x F_x(S_x^3-q i S_x^2).          (6)

Here `||A||<=1/2`; B_q is its charge-q Fourier component, so
`||B_q||<=1/2`. In the Marshall frame (6) becomes

    B_q=(1/32) sum_x (2F_x)(epsilon_x Z_x-q i Y_x).

It obeys `[Q,B_q]=q B_q` for the actual conserved charge
`Q=(1/2)sum epsilon_x X_x` in that frame. The creation complement
`I-B_+^*B_+-B_-^*B_-` is positive; its square root supplies the third
Kraus operator. The integer producer independently checks its positivity
by a row-norm bound.

The even translation T sends each site x to x+2 and has order three. It
commutes with H_B and Q. Its three actual spectral projections are

    P_l=(I+exp(-2pi i l/3)T+exp(-4pi i l/3)T^2)/3,
    l=0,+1,-1.                                           (7)

These are translation phases, not the two-element reflection labels of
the eight-site cube. An odd translation in the Marshall frame exchanges
the source charges and cannot replace T in the joint charge-momentum
readout.

Perron positivity makes the unique ground invariant under each site
permutation symmetry. Translation by three sites exchanges B_+ with
`-B_-` and reverses Q, while commuting with T and H_B. It forces the
ground charge to be zero and pairs the two charged spectral measures.
Site inversion fixes each B_q and exchanges P_+ with P_-. Also
`B_q+T B_q T^*+T^*B_q T=0`. Hence each zero-momentum preparation record
has probability zero, and the four nonzero charge/direction records have
identical complete energy measures.

The same-clock energy interval `[8/3,18]` contains the entire positive
spectrum of `L=(H_B-E_B)/p`; its complement includes the ground. The
instrument retains all two charge choices, all three translation phases,
both energy events and the original creation failure: thirteen records.
Eight records have zero probability on this ground, four are the directed
successes, and one is the preparation complement. Every effect and zero
label remains present. The inherited finite native compiler can implement
these finite contractions and controls; no complete pulse tape is executed
by this certificate.

## 4. Certified source values and the remaining energy width

The floating eigensolver proposes a complete integer-quantized basis at
48 bits. Exact integer Gram and Hamiltonian residual calculations admit
it, identify the isolated ground and bound its vector error. No solver
residual or eigenvalue is accepted on floating evidence. The producer
requires ground-vector error below `10^-8`, ground energy within
`[-4.309116,-4.309115]` and gap at least `2/3`. It checks all 64 levels
before deriving the displayed full positive energy interval.

For a real image v=B_+Omega and any power j of the real Hamiltonian,
commutation with T gives the exact nontrivial-character moment

    <v,H_B^j P_+ v>=<v,H_B^j P_- v>
      =(<v,H_B^j v>-<v,H_B^j T v>)/3.                     (8)

The first three moments j=0,1,2 use only integer products and rational
normalization of the proposed ground. Ground-vector error epsilon gives
an error at most `epsilon ||H_B||^j/2` in each original moment, because
`||B_+||<=1/2`. Interval division uses the strictly positive original
probability. The variance of H_B equals the variance of H_B-E_B, so the
energy-reference interval enters the mean but not the variance.

Every fresh calculation must lie within these fixed admitted intervals:

| Quantity for each of the four records | Certified interval |
| --- | --- |
| Original probability | [0.00481981, 0.00481982] |
| Mean energy on L | [6.097426, 6.097427] |
| Energy variance on L | [0.048956, 0.048957] |
| Original preparation complement | [0.98072072, 0.98072076] |

Each original directed probability exceeds `1/208`. Each relative energy
variance is greater than `1/800`. Thus this literal momentum-resolved,
charge-resolved source response is not a single energy line. The finite
variance does not imply a continuum obstruction: the ring and its source
parameters are fixed, rather than the prescribed growing-volume joint
limit. Conversely, the continuum probability (3) is proved from its
source theorem and is not inferred from this finite example.

A shrinking continuum energy shell requires the actual susceptibilities
and force slacks of `WARD_DYNAMICS.md`, or an independently reducing
positive-weight spectral branch. Narrow spatial momentum, exact charge
conservation and the broad energy bounds (4) do not supply that missing
estimate. The construction makes directed measurements available while
keeping their measurable energy width explicit.
