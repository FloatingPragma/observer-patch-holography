# A sitewise chiral source with actual native currents

The declared program `coherent-sitewise-chiral-family-v1` has one local
fermion space `K=C^2 tensor R`, where
`R=Lambda^2(C^5) direct_sum Lambda^4(C^5)`. Its two-component memory is a
left-Weyl Cauchy factor, ordered `spin*15+internal`. Its native internal
basis is exactly the fifteen subsets in
[the actual-current attachment](../field_identity_audit/PROOF.md).
The code supplies all `2^30` occupations at each spatial cell. Grading,
binding, routing memory and the Hamiltonian are declared construction data.
The native idle rule fixes the selected even cloud and supplies no binding.

## 1. A local object and every occupation

Fix one canonical odd paired-flux orbit. An occupied object carries that
literal inside/exterior pair, the two-component memory and one actual native
dark cloud. For native component j the dark orbital is

`d_j^dagger=(c_(port2,j,repair0,address0)^dagger
            -c_(port2,j,repair0,address1)^dagger)/sqrt(2)`.

The adjacent source proof executes the full native wave and link controls:
each orbital changes sign, every selected even-degree cloud is fixed, and
all twelve actual native currents act by the exterior R matrices. Combining
both star-related words gives the positive interior Gauss code. Denote this
one-object isometry by B. A vacant whole-object slot is the native-empty
flat Gauss-positive orbit, orthogonal to every occupied B state. The
occupied object is declared odd; its native cloud alone has even parity.
The selected gauge algebra is the trace-free subgroup of the native block
algebra. The additional central native degree current is separate.

At each cell x supply thirty independent whole-object slots. For ascending
orthonormal modes `i_1<...<i_n`, define

`E_x |i_1...i_n> = (1/sqrt(n!)) sum_pi sign(pi)
  B|i_pi(1)> tensor ... tensor B|i_pi(n)>
  tensor |vac>^(30-n)`.

The occupied prefix stays inside that spatial cell. The `n!` distinct terms
are orthogonal; different occupation sets and different n are orthogonal.
Thus `E_x^* E_x=I` on the entire local Fock space, including n=0 and n=30.
The ambient selected slot space has dimension `31^30`; no claim of minimal
storage is made. The code admits all local simultaneous occupations, rather
than projecting to zero or one composite per cell. If `P_x=E_x E_x^*`, then
`E_x a_i E_x^*` obeys CAR with identity `P_x` on its code. These odd operators
are only claimed on that admitted invariant code.

Each occupied source object contains two nontrivial spherical face holonomies.
Consequently the actual total face-count operator over that cell's slots
obeys `D_x E_x=E_x 2N_x`. This is a counting identity. It supplies neither a
chiral mass term nor an energy coefficient. A local native current is the
sum over the thirty actual cloud registers. Its action commutes with the
antisymmetrizer and preserves the occupied prefix, so for every n,

`Q_x(T) E_x=E_x dGamma(I_2 tensor T)`.

This proof uses the actual one-cloud intertwiner. It does not replace a
cloud current by a nominal CAR-slot transition. The receipt checks all
zero-, one- and two-object columns on the full thirty-mode site and three
three-object controls, with every one of the twelve native currents. The
all-occupation statement follows from the displayed multilinear identity.
The native source tick is the identity on each even cloud, hence on their
products. Its `I-Re U_F` energy contributes zero on this code.

## 2. Spatial support without a global particle-slot sort

Use the graded tensor product of these cell codes and set
`E=graded_tensor_x E_x`. Each E_x preserves parity. For an even logical
operator O supported on cells S, its source representative on the code is

`O_source=(E_S O E_S^*) tensor I_(outside S)`.

It has exactly that set of spatial cells as support in the graded source
algebra. For a Hermitian local Hamiltonian this is Hermitian; exponentiation
acts as identity on the orthogonal complement of the local code. For a
unitary local gate W its complete extension is

`W_source=E_S W E_S^*+I-P_S`.

The two orthogonal blocks prove unitarity. Such operators commute with
outside even local observables and preserve the tensor product code.
Even operators supported on disjoint cell sets commute. This supplies the
spatial support transfer absent from a single global particle-slot encoding.
Within a cell, arbitrary access to its thirty slots is an allowed local
controller. Gate complexity and physical implementation costs are not bounded
by this identity.

This is locality in a declared graded architecture. A representation by
ordinary qubit matrices uses parity strings across cells. A theorem that
those strings can all be compiled into bounded-range bosonic gates is not
claimed. An even hopping term between two cells is local in the graded
algebra, regardless of their positions in an auxiliary ordering. Changing
the ordered matrix representation does not change that algebraic support.

## 3. Kinetic and gauge operations on that same code

Any regulated one-body left-Weyl symbol on `C^2` may be attached tensor R,
using the same Fock law and site encoding. The declared continuum principal
symbol is `h_L(p)=-sigma.p tensor I_15`. It has thirty Cauchy rows and one
family of fifteen left-Weyl internal components. There is no four-component
Dirac factor with both spacetime chiralities in the same R. Its negative
energy holes carry the conjugate gauge representation; they do not add a
second independent Weyl family. Choosing this principal symbol alone does
not remove regulator mirrors or prove an interacting chiral continuum.
Those properties belong to the specified regulator and state.

For an external unitary link U_xy in the R representation, the oriented
nearest-neighbor coefficient `i sigma_j tensor U_xy/(2a)` and its adjoint
are Hermitian as a sum. Under site frames G_x,G_y it transforms using
`U_xy -> G_x U_xy G_y^*`. The full exterior functor transfers the one-body
covariance to every occupation. The receipt executes an exact non-Abelian
fundamental frame, its actual exterior R representation, and all three
Weyl axes. This supplied-background coefficient is distinct from the
operator-valued link below; in particular the local sine regulator has its
own additional modes.

A finite dynamical gauge attachment reuses the construction in
[the quantum-link proof](../../standard_model_source_map/GAUGE_ATTACHMENT.md)
with the actual dark-module basis. The five R blocks have dimensions
`3,6,1,2,3`. On each spatial edge introduce

`H_link=C|0> direct_sum direct_sum_r(V_r^* tensor V_r)`,

of dimension `1+sum_r dim(V_r)^2=60`, with labels 0 and `(a,b)` inside each
irreducible block. Set `U_ab=|a,b><0|`. For a Hermitian current T the left
link generator acts as `-T^T` on a and the right generator as T on b. The
matter-plus-link Gauss operators are

`G_x(T)=dGamma(I_2 tensor T)_x+L(T)`,
`G_y(T)=dGamma(I_2 tensor T)_y+R(T)`.

Define the same-code edge Hamiltonian

`H_xy=sum_(spin u,v;a,b) [i sigma_j[u,v]/(2a)]
   c_(x,u,a)^dagger c_(y,v,b) tensor U_ab + adjoint`

plus any real coefficient times the excited-link projector. The operator
is even and supported on the two cells and their edge. Direct coefficient
identities give `[G_x(T),H_xy]=[G_y(T),H_xy]=0`. Independence of normally
ordered CAR bilinears proves these on the complete `2^60` matter Fock
space, not just a one-particle truncation. The receipt checks 72 kinetic
Ward identities and 24 electric-projector identities. Reversing either
endpoint orientation gives a nonzero hypercharge residual. Sitewise E
transfers the matter currents into actual native source currents, so the
same exact Ward identities hold on the encoded source.

The quantum-link registers are declared source architecture. Their raising
operators are nilpotent, `U_ab^2=0`; none has a nonzero eigenvalue. They
cannot be replaced by identity links in a vacuum to import a free Weyl
spectrum. Their interacting reference, electric/magnetic scaling and chiral
regulator require their own calculation. Exact finite gauge covariance alone
does not eliminate mirror fields or construct a chiral quantum measure.
The established R anomaly checks, including four weak doublets, are reused.

## 4. Local Higgs interactions and mass scope

A source-compatible local Higgs vertex is even in the Weyl fields. For a
symmetric R-valued kernel C(H) and `epsilon=i sigma_2`, it has the form

`(1/2) sum_(alpha beta i j) epsilon[alpha,beta] C_ij(H)
  c_(alpha,i)c_(beta,j) + adjoint`.

The operator changes local number by two but preserves parity, so it lifts
through the same E_x without enlarging the spatial support or discarding any
occupation. Native cloud number conservation of the idle rule does not
imply this vertex: the vertex is an explicitly chosen interaction. Its
scalar charge, current identities, apparatus, scalar state and coefficients
must be supplied on that declared source. The invariant tensors are owned
by [the chiral mass proof](../physical_attachment/PROOF.md).

The source construction therefore supplies actual chiral currents, a local
full-occupation code and compatible finite gauge/Higgs operator interfaces.
It does not fix a positive interacting vacuum, a mirror-free local regulator,
Yukawa values or a physical mass scale. The absence of an unbroken constant
Weyl mass does not forbid Higgs-generated masses; it requires their actual
interaction and reference. No Standard Model species or GeV mass follows
from the family weights alone.
