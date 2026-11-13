# Native proton–pion channel support

The declared CAR source already has a nonzero weak output in its complete
four-quark/one-antiquark sector. Particle identification needs more than
that total norm: the output must actually couple to the native composite
fields. This package computes that projection without a mass-fit vector.

The native colour-epsilon contraction supplies five proton fields with
spin and isospin one half. Include the translated onsite field and all
four positions of a spin-zero down/anti-up pion. Keeping total occupied
span at most one and fixing its leftmost site gives 23 composite columns.
Literal fermionic creation words antisymmetrize each product. Their Gram
matrix is positive definite: all 23 columns are independent. The union
with the actual weak output has 2736 occupation states. Thus

    Pi_C = Y (Y*Y)^(-1) Y*

is an exact orthogonal projector, including the exchange/Pauli relations;
its complement is retained. With P the seven-coordinate Lambda inclusion,
G=P*P and c=sin(theta)cos(theta), direct integer CAR contraction gives

    M = P* W* Pi_C W P / c².

The generalized eigenvalues of (M,G) are 6/5, 241/150 (twice), and
(277±sqrt(54409))/200 (twice each). In particular

    (1/5)G < M < (13/5)G.

The displayed rational eigenvalue bounds use an outward square-root
bound. Both M and the complementary compact Gram are positive; the
projection does not hide another zero or replace the current by a fitted
transition matrix.

The local vertex leaves the occupied spatial span unchanged: it replaces
the strange quark and creates a quark–antiquark pair at that same site.
Consequently noncompact Lambda inputs have no component in this compact
23-field output. The complete rest-pole certificate bounds its noncompact
weight by q. Every normalized full Lambda rest pole therefore obeys

    c² a(1-q) <= ||Pi_C W Psi||² <= c² b,

where a,b are the exact rational bounds on M relative to G. This remains a
k=0 fibre statement. No normalizable Lambda packet is asserted. The
normalizable nucleon packet in the separate package does not change that
domain automatically.

The component probe W/36 and sqrt(I-W*W/1296) form a complete instrument on
the entire uds source and Fock complement. Split its success into the
23-field output, other compact output, and every noncompact output. The
receipt gives bounds for all four terminal events. The actual prepared
seven-coordinate trial has success 1/7, failure 6/7 and its certified
rest-pole overlap; the original-input event lower bound includes both
factors. Nonselective trial preparation has zero mean energy change on
P7/7. Neither the weak component probe nor the composite-channel readout
is asserted to conserve energy; their complete mean work is retained.

## Exact support control and the open scattering premise

Every literal weak word creates an anti-up at the site where it converts
strange to up. In a proton plus negative-pion product both up quarks must
belong to the proton and the anti-up to the pion. Therefore direct
projection of W onto bare proton/pion supports that are spatially
disjoint is identically zero, on the complete initial uds space. This
selection statement is checked for all 72 words; it does not depend on a
compact initial approximation.

That zero is not a decay amplitude. Local creation followed by strong
propagation can produce separated outgoing particles. Conversely a
nonzero compact response and an open sum-of-cluster energy shell do not
prove that the interacting spectral density is nonzero at the Lambda
energy. Pi_C generally fails to commute with H5 and is not its spectral
projector. The missing calculation is the outgoing wave operator or full
H5 resolvent boundary value applied to this actual weak output, with the
same finite-range strong action. Full weak vacuum and other channels must
then be included. No autonomous decay, physical lifetime or new empirical
validation is inferred from this source result.

## Actual kinetic propagation into disjoint composite supports

The complete native nearest-neighbour hopping component supplies a positive
constructive step beyond the bare support zero. Contract all 96 directed
quark/antiquark hopping words against the actual weak output and the two
onsite proton/pion products at distinct neighbouring cells. Their Gram is
108 times the identity. The exact 2-by-7 contraction, before the common
kinetic and Cabibbo factors and the source epsilon normalization, is

    [0, -54, 18,  0, 36,  0,  0]
    [0,   0,  0, 36,  0, 18, 54].

For the normalized second Lambda coordinate, the leading separated
amplitude is (3/2)tc while ||W psi||=sqrt(6)c. This input is an actual
spatially localized finite-support state; it is not the Lambda rest pole.
On the complete infinite N5 carrier the off-diagonal kinetic generator K
has norm at most 10t. The diagonal kinetic term is the constant 10t and
adds only a phase. The spectral-theorem Taylor remainder gives

    ||(exp(-i tau K)-I+i tau K) W psi||
       <= tau² (10t)² sqrt(6)c/2.

At tau=1/(1000t), direct separated overlap starts at zero and the evolved
amplitude is bounded below by

    c [3/2000 - sqrt(6)/20000] > 0.

The receipt uses an outward rational square-root upper bound and includes
the W/36 probe probability. Probe failure and the entire orthogonal N5
success outcome are retained. This establishes finite-time native kinetic
propagation into disjoint one-cell-separated composite supports on the
full infinite kinetic carrier. It neither omits those supports from the
strong interacting space nor claims that the full autonomous strong
Hamiltonian has this derivative: correlated propagation and every other
strong term must still be included in the actual scattering calculation.
The separated-field projection generally performs work; no zero-work
claim is made for that readout or for switching the kinetic substep.
