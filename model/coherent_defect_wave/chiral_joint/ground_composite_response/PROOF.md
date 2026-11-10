# Native composite response on the full physical local ground

## Source, reference and probe

The source is the thirty-mode native fermion Fock space tensor the
five-state scalar. Its Hamiltonian contains all three quantum Yukawa
channels and the scalar excitation projector, each with coefficient
one. The clock is `exp(-it H_local)`, with hbar one. The physical domain
has all twelve local gauge generators zero and no compensating boundary
charge. No intersite Hamiltonian is added.

The [complete local reference theorem](../local_reference/PROOF.md)
proves a global physical ground doublet, with raw L=1, B3=0 and spin
one-half. Its reference is `rho=P_g/2`. This normalized projector is
stationary and rotation invariant. Its selection is the declared source
preparation law; this package supplies no probability of preparing it
from an arbitrary input. Its odd fermion parity is retained when the
probe apparatus is attached.

The actual composite is

`J_alpha=phi_0 c_(alpha,11)-phi_1 c_(alpha,10)`.

The four recorded settings are removal by J_alpha and addition by its
actual adjoint, for alpha up and down, each with probability one-quarter.
The successful Kraus operator is `K=J_kind,alpha/2`. No success branch
is renormalized. This is the response of the interacting ground; the
blank-state preparation `J^dagger Omega` is a different input.

## Complete source and target spaces

The particle-hole and representation construction in the local
reference theorem places the ground in the unique quark C=54 module
of dimension 140. J acts only on leptons and the scalar, so it preserves
that complete quark module. It cannot reach any of the excluded quark
representations. The quark module has transformed occupation twelve;
placing all quark modes before the six lepton modes in the fermion
tensor convention therefore contributes parity plus one to the odd
lepton operator. The native particle-hole transform is applied to the
operator polynomial before its six-mode action is evaluated.

The source contains the entire sixteen-dimensional highest-spin-one-half
block at raw L=1 and its exact spin-lowered partner. Their direct sum has
dimension 32. This is sufficient for the complete ground projector,
whose rank two is independently checked from the projector polynomial.
The source metric is the pullback of the actual positive Fock metric.

Removal reaches raw L=0 and addition raw L=2. In either case the entire
physical target in the maximal quark module has dimension twenty-six:
eleven spin singlets and five spin triplets. The implementation builds
the complete weak/hypercharge-singlet kernels at magnetic weights
`2S_z=-2,0,2`, of dimensions five, sixteen and five. Exact kernel
dimensions agree with an independent weight-character count. Every
native composite image is checked against the full target inclusion,
and the reduced Hamiltonians satisfy the full invariant intertwiner
and metric-Hermiticity equations.

If G_s and G_t are the two positive metrics, the adjoint in reduced
coordinates is

`J^sharp=G_s^(-1) J^dagger G_t`.

It is also checked by applying the actual adjoint CAR/scalar polynomial
before taking its source matrix element. That full adjoint image can
contain spin-three-halves states, so an unproved reverse-image closure
is not assumed. The Euclidean transpose of the coordinate matrix is
not used as its physical adjoint.

Set `y=E(E-1)` and define

`p_o(y)=y^3-313y^2+31708y-1041468`,
`p_v(y)=y^3-313y^2+31738y-1041088`.

The exact squarefree minimal polynomials of the complete source and
target matrices are

`p_s(E)=(E-1)(y-105)p_o(y)`,
`p_t(E)=(E-1)(y-103)p_v(y)`.

Both have degree nine. They are checked by evaluating them on every
matrix coordinate. The target has six simple energy roots from p_v,
two roots of y=103 with rank three each, and E=1 with rank fourteen.
Thus all twenty-six target states, including every energy-one state,
are retained. The corresponding source ground root e_g is the smallest
root of `p_o(E(E-1))`; its interval is contained in the independent
global-ground interval from the local reference certificate.

## Exact spectral weights

For a finite self-adjoint matrix H annihilated by a squarefree polynomial
p, the projector at a root e is

`P_e = [p(H)/(H-e)]/p'(e)`.

Here the quotient denotes the polynomial quotient in the formal
variable H, evaluated at the root e. The implementation uses exact
coefficients

`c_i(e)=sum_(j=i+1)^degree(p) p_j e^(j-i-1)`

and reduces them modulo the relevant irreducible energy factor. It
checks the source projector trace two and each complete target
projector rank algebraically. The distinct real roots and the full
matrix annihilator give a complete orthogonal projector resolution.

All transition weights follow from the rational cross moments

`M_(j,i)=Tr(H_s^i J^sharp H_t^j J)`, for `0<=i,j<=8`.

For target root e_f the joint recorded probability is

`w(e_f)=Tr(P_g J^sharp P_(e_f) J)/32`.

The factor 1/32 consists of setting probability 1/4, the successful
Kraus factor 1/4 and ground-density normalization 1/2. The receipt
contains every cross moment and each reduced bivariate numerator,
with denominator `p_s'(e_g) p_t'(e_f)`. These formulas specify exact
algebraic weights independently of the printed intervals.

Rational root isolation to width below `10^(-110)`, rational interval
evaluation and outward dyadic rounding provide the recorded bounds.
All nine weights are strictly positive in each of the four settings.
Every excitation `e_f-e_g` is positive. The lowest line is the global
local physical gap, with a positive joint probability near 0.017786
in each setting. The total joint success per setting is near 0.04113.
The exact matrix moments give equal spectra and weights for both spin
settings and for the addition/removal pair. These are source
conjugation and rotation relations, with no measured particle or
antiparticle species assigned.

Summing the exact target projector resolution gives the success
probability directly as `Tr(P_g J^sharp J)/32`. Its independent
one-variable interval overlaps the sum of all nine weight intervals.
The completeness equality follows from the exact projector identity;
interval containment alone is not used to prove normalization.

## Failure, fermionic reference and complete records

The [full native neutral apparatus](../neutral_composite/PROOF.md)
proves the effect spectrum of `A=J^dagger J` is contained in `{0,1,2}`
on the entire fermion/scalar carrier. Its exact polynomial projectors
are P0, P1 and P2. The positive failure operator is

`F=P0+(sqrt(3)/2)P1+(1/sqrt(2))P2`,

so `K^dagger K+F^dagger F=I` as a full normal-ordered operator
identity. This equality is recomputed for every setting. All twelve
gauge currents commute with both J and F. The code also verifies the
full native spinor commutator identities relating both probe spins.

The full forty-dimensional unitary dilation acts on the twenty-state
active two-mode/scalar space and a two-state fermionic reference. Its
block form is

`U=[[sqrt(I-K^dagger K),-K^dagger], [K,sqrt(I-K K^dagger)]]`.

Both unitary identities and source-plus-reference parity and number
conservation are replayed. Spin rotation gives the other probe-spin
control. The reference begins in its even branch. On success it has
odd parity and the complementary raw lepton charge: plus one for
removal, minus one for addition. Its S_z similarly compensates the
selected source spin component. The declared controller and its
retained reference therefore implement an even total operation;
the odd system Kraus map is not treated as an autonomous operation.

Failure is the full operator F on the actual source. It can enter
spin-three-halves states outside the 32-dimensional initial source
space. It is never replaced by a square root of a compressed source
effect. Its probability is evaluated by the full completeness
identity, and no energy selection is imposed on that branch.

For each setting the complete record alphabet contains nine successful
energy outcomes, one failure, and the orthogonal complement of the
complete target space. The last has exactly zero probability for this
input because every composite image was checked. The setting register,
failure and reference are retained. There are forty-four joint records,
forty with positive probability; their exact sum is one. Orthogonal
record pointers admit the usual classical record copying without
copying the underlying quantum state.

## Domain and verification

This is a finite, gauge-neutral composite response on the full physical
local ground with the declared preparation law. It supplies no
interacting many-cell vacuum, momentum-dependent propagation,
relativistic pole or observed particle identification. Charged local
sectors compensated by gauge links belong to a different reference
problem. The probe controller is a finite apparatus intervention,
rather than part of the idle local Hamiltonian.

The certificate recomputes all rational source and target matrices,
native operator images, full metric adjoints, exact moments, projector
ranks and weights. It uses no numerical eigenbasis. A separate test
compares generalized-eigenvector amplitudes with every exact weight
under a declared absolute tolerance `2*10^(-10)`; that diagnostic does
not widen any certified bound. Hostile tests reject an omitted dark
energy factor, an incomplete composite image, the wrong adjoint charge
sector, a quark-changing probe and missing source custody. All declared
inputs, the local-ground receipt and the repository-local import
closure are pinned, with missing files rejected.
