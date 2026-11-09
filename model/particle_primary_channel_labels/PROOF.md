# Exact channels of the physical plaquette response

## Source and five conserved occupation counts

Use the unchanged [primary member](../particle_common_contract/PROOF.md),
including its open Wilson Hamiltonian H, full product link cutoff,
source Gauss projector, polynomial density rho, idle canonical reference
and source clock. The [finite response](../particle_primary_survival/PROOF.md)
uses the Hermitian contraction `A_p=ReTr_M15(U_p)/15`, with the entire
coordinate product compressed after multiplication.

The actual native central generator splits M15 into five inequivalent
blocks. In the integer convention Q6=6Y their `(q,dimension)` pairs are
`(-4,3),(-3,2),(1,6),(2,3),(6,1)`. Let P_q be their one-body projectors
and define on the complete CAR carrier

```
N_q = sum_x dGamma(I_4 tensor P_q),
N_F = sum_q N_q,                 Q6_total = sum_q q N_q.       (1)
```

The source has no interblock term or pairing term. Every native link
matrix commutes with P_q, because P_q is a central spectral projector
and G is connected. The actual onsite `3 beta_D/a` acts on the Dirac
factor. Each hopping is a fermion bilinear with coefficient
`(-beta_D+beta_D gamma_j)R(U_j)/(2a)`, together with its adjoint.
The CAR identity `[dGamma(P),c*Mc]=c*[P,M]c`, applied to the global
one-body block projector, proves `[H,N_q]=0`. Electric and magnetic
terms act only on links. Each endpoint Gauss action preserves P_q, so
`[P_G,N_q]=0`; the scalar preparation polynomial then gives
`[rho,N_q]=0` for every q.

The two internal basis charts are checked separately in `channels.py`.
The exterior-mask currents determine the central blocks; the identified
packet chart gives the same five charges and dimensions and is tested
on a non-diagonal native link. Indices from the two charts are not equated.

The simultaneous central gauge transformation at every owned and ghost
vertex fixes every link coordinate. Its remaining action is
`exp(i theta Q6_total)` on matter. Full zero-external-charge Gauss
invariance therefore gives `Q6_total P_G=0`. This restricts the allowed
joint occupations to `sum_q q N_q=0`; it does not set every N_q to zero.
The five occupation counts are source CAR labels. No baryon, lepton or
observed-particle interpretation is inferred from (1).

A_p acts as identity on matter, is gauge invariant and is CAR even.
Consequently its exact time translates, real triangular average and
finite positive quadrature/Cayley filter B have

```
[B,N_q]=0,   [B,P_G]=0,   [B,(-1)^N_F]=0.                      (2)
```

In a simultaneous H and occupation basis a nonzero B_ji has identical
five occupation vectors on i and j. Its original positive spectral
measure is a sum over those reducing sectors with their original rho
weights. No sector weight is divided out. The neutral dilation and its
complementary square root, source echoes and all complete Ramsey Kraus
operators inherit (2). Their neutral records carry none of these five
charges. The reference is retained and has Hamiltonian zero.

## The physical cube action and the open boundary

The allocated ambient link set is asymmetric: each cell owns three
positive outgoing links, including dangling positive-face links. A
naive geometric cube action on that whole allocated set is invalid.
Each outside ghost vertex has exactly one incoming link and no matter.
Its Haar Gauss projection retains the right-invariant functions of that
one link. In L2(G), and hence in the invariant finite word subspace,
these are precisely the constants. The ghosts are distinct one-link
vertices. The physical source therefore factors as the Gauss-invariant
ordinary open interior graph times fixed constant dangling-link states.
Their electric energy is zero, and they occur in no retained matter
hopping or plaquette.

Proper cube rotations about `(side-1)/2` act on the interior graph.
When a rotated edge reverses orientation its coordinate becomes U*.
Haar measure and the per-link word space are invariant under inversion:
matrix entries become transposed conjugate entries, which are in the
same word alphabet. The tensor cutoff consequently intertwines the
rotations. This argument keeps the complete parent plaquette products.

The [binary cube enumeration](../cubic_particle_channels/group.py) is
reused only for its exact 48 unit quaternions and group characters.
Its comparison Hamiltonian and particle spectra are not transferred.
For q=(w,x,y,z), the spin lift in the actual primary Clifford chart is

```
S(q)=wI-x gamma_2 gamma_3-y gamma_3 gamma_1-z gamma_1 gamma_2,
S(q) beta_D S(q)*=beta_D,
S(q) gamma_j S(q)*=sum_k R(q)_(k,j) gamma_k.                  (3)
```

Here the displayed gamma indices run from 1 through 3. The code checks
(3) exactly for every q. A positive hopping coefficient is mapped to
the matching positive coefficient or its adjoint when the edge reverses.
Second quantization of the site permutation and S(q), with the link
pullback, gives a unitary physical representation commuting with H.
The electric trace-metric Casimir is invariant under endpoint exchange;
the real plaquette traces are permuted. Local Gauss constraints are
permuted as well. Thus this representation also commutes with rho.
Its central binary element is CAR parity; adjoint action on the even
plaquette algebra factors through the 24 proper spatial cube rotations.

The source marginal is rotation invariant. The canonical purification
is invariant under `R tensor conjugate(R)`, not in general under R on
the source alone. This covariance is a comparison of state descriptors;
the reference is not counter-evolved or rotated during the readout.
Continuous rotations, Lorentz spin, and a charge-conjugation quantum
number require additional attachments. In particular, complex
conjugation of link coordinates alone does not establish a symmetry of
the entire matter Hamiltonian. No such charge-conjugation label is used.

## A distinct cubic scalar observable with original finite weight

For every active selected member n>=1, the authenticated plaquette has
anchor `(side/2,side/2,side/2)` and plane (1,2). Its proper cube orbit
contains exactly 24 distinct geometric plaquettes. Its coefficient
permutation representation is regular:

```
A1 + A2 + 2 E + 3 T1 + 3 T2.                              (4)
```

The original single plaquette is therefore not itself a cubic scalar.
To see that (4) is an operator statement, act on the physical empty-Fock,
constant-link blank b. Distinct plaquettes have orthogonal character
loops: Haar integration over a link present in just one loop vanishes.
The ten native/conjugate loops at each plaquette are also orthogonal.
Thus the coefficient-to-operator map is injective on b. These statements
hold for the product cutoff because all those single-plaquette vectors
are retained. The blank is only a matrix-element witness, never the
selected preparation.

The concrete scalar average is

```
A_A1=(1/24) sum_(p in orbit) A_p,       ||A_A1||<=1.           (5)
```

It is a distinct finite observable on the same source/state/clock. Its
records are not identified with records of the unaveraged plaquette.
The [survival commutator calculation](../particle_primary_survival/PROOF.md)
gives orthogonal vectors `[H,A_p]b` at different plaquettes, each with
squared norm `K/a²`, K=4949/20. Linearity and orthogonality give

```
||[H,A_A1]b||²=K_A1/a²,       K_A1=K/24=4949/480.            (6)
```

All matter terms annihilate b before and after each A_p. The compressed
magnetic commutator annihilates b by the same retained-vector argument;
no global commutation of compressed coordinate operators is assumed.

Reuse the finite response theorem with rho>=c P_G, c=2^-N, source width
W, filter half-width T=1/W and Ramsey time t=1/W. Let r=(23/24)^4.
The original scalar contrast obeys

```
D_A1 >= d_A1=(11/24) r K_A1 c/(a² W²)>0.                    (7)
```

Polynomial passivity and Hermiticity are unchanged. Choosing the
smaller band endpoint `delta=(1/a)2^-ceil(N/2)` in that theorem gives

```
mu_ideal([delta,W]) >= (r K_A1-1)c/(2a² W²),
mu_finite([delta,W]) >= [(37/96)r K_A1-1/2]c/(a² W²)>0.       (8)
```

The second inequality uses the inherited finite quadrature/Cayley plan
with d_A1: its band-weight error is at most d_A1/4. It is the actual
finite Hermitian contraction's original transfer measure. The same
four-outcome [direct compact service](../particle_compact_mass_probe/PROOF.md)
can be applied to that contraction, with the same conditional compiler
and resource boundaries. No spatial matrix, full response probability
or efficient implementation is asserted by the factored sufficient plan.

Since final excitation is at least positive transfer and at most W,
the [mass-readout coupling](../particle_primary_mass_readout/PROOF.md)
transports (8) to the same full-width final-excitation band without
weight loss. This requires no numerical ground-energy subtraction.
It does not produce a rest-frame observable or an isolated spectral line.

The selected rho is a mixture of initial energy and occupation sectors.
A scalar operator connects equal cubic irreducible representations,
including their distinct multiplicity spaces. It does not force the
initial or final source state to be an A1 vector. A scalar vacuum would
make that further selection, but no unique scalar vacuum is assumed.
Even in such a vacuum, continuum spin 0 and spin 4 both contain A1 on
restriction to the cube group. Scalar gauge loops may overlap neutral
composites and multiparticle states without selecting any named species.

Equations (7)-(8) retain shrinking times, a lower endpoint tending to
zero, a growing upper endpoint and a nonuniform lower weight. They
exclude neither gapless spectral weight nor massive weight. The concrete
additional input for a finite-capacity mass comparison is a resolved
branch, same-state momentum and initial-momentum control, and controlled
finite-domain dispersion and readout errors, as stated in the mass-readout
proof. A continuum mass claim additionally requires uniform survival
in a fixed physical band. Declared SI units can be attached to either
comparison; these channel labels alone supply no source-derived GeV mass.
