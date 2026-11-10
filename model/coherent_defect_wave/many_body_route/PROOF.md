# A graded flux source with a positive free-field reference

The declared program is `coherent-graded-flux-wilson-v1`. Its degree-one
space is the same `cell × five flux faces × four Clifford components` code
as the [paired-flux route](../mobile_flux/ROUTING_COMPLETION.md). Its full
occupation law is fermionic and its clock is a specified local Hamiltonian.
Both are explicit architectural choices. Neither follows from the bosonic
two-lane SWAP law or from choosing a branch of a Floquet logarithm.

## 1. Every occupation has a literal source encoding

At each cell introduce twenty independent carrier slots, indexed by `(f,c)`.
The five f values refer to the five canonical odd link orbits in the mobile
source; c is the four-component Clifford memory. Slot `(f,c)` admits the
four literal words `0,STAR,g_f,g_f xor STAR`. On its Gauss-positive code,

`|0> -> (|0>+|STAR>)/sqrt2`,

`|1> -> (|g_f>+|g_f xor STAR>)/sqrt2`.

These are orthogonal. The first has no nontrivial face holonomy; the second
has the actual inside/exterior pair, with regional charge one and total
spherical face count two. The tensor encoding of all occupation bits is an
isometry. It preserves every slot's interior Gauss constraint and satisfies
`E* Q_j E=n_j` and `E* D_j E=2n_j`. The source has twenty carriers per cell,
rather than a single carrier that somehow holds twenty simultaneous objects.
All `2^20` local occupation states are retained. Multiple objects can occupy
one spatial cell in distinct slots; double occupation of a single mode is
excluded by the declared fermionic algebra.

Declare the slot code to have even vacuum and odd occupied grading, and use
the canonical anticommutation relations
`{a_i,a_j*}=delta_ij`, ` {a_i,a_j}=0`. In a fixed occupation ordering the
creation matrix carries `(-1)^(number of occupied preceding modes)`. The
literal encoded source gate includes this sign. A plain bosonic swap of
link registers does not supply it. In an ordinary qubit implementation the
same prescription uses parity strings; ordinary tensor-gate locality is not
inferred from locality in the graded CAR algebra. The declared operational
field algebra and its even local observables use the grading.

For a finite logical unitary W, its source lift is `E W E* + I-E E*`.
It is unitary and commutes with all admitted interior Gauss constraints.
Odd field operators are defined on the invariant source code by the same
isometry. Their CAR identity is the identity of that code. The code has no
native wave occupations in this construction. The full literal native
face-changing update need not preserve a slot whose admitted odd orbit is
fixed to f; the [native-routing interaction](../native_routing_join/PROOF.md)
requires a separate transfer to this occupation architecture.

The restriction to one occupied slot agrees isometrically with the paired-
flux source code: position is the occupied carrier's cell, f is its literal
flux orbit, and c is the routing component. Coin changes and spatial hopping
move the actual paired-flux code between slots. The exterior lift `Gamma(R)`
of the one-object routing unitary is also well defined on all occupations,
with degree-N matrix elements `det R[I,J]`. It supplies an exact fermionic
extension of that particular tick. The Hamiltonian below specifies a
different clock; its filled-sea state is not asserted stationary for R.

## 2. A finite-range Hamiltonian and its complete branches

Use the owned Hermitian Clifford matrices `alpha_i,beta`, satisfying
`{alpha_i,alpha_j}=2delta_ij`, `{alpha_i,beta}=0`, `beta^2=I`.
Set lattice spacing a, source time in the same units, Wilson coefficient one
and a unit coefficient converting the literal face count `D=2` into the
mass parameter. The one-body symbol at physical momentum p is

`h_a(p)=sum_i alpha_i sin(a p_i)/a + beta M_a(p)`,

`M_a(p)=2 + sum_i(1-cos(a p_i))/a`.

Its onsite block is `(2+3/a)beta`; its forward-neighbor block is
`-(beta+i alpha_i)/(2a)`, with Hermitian conjugate on the reverse bond.
The action is diagonal in the five face channels. This is an executable
nearest-neighbor source operator on every finite periodic cube.
Second quantization defines the complete clock
`H_a=sum_ij a_i* (h_a)_ij a_j`, including every occupation sector.

The Clifford identity gives exactly

`h_a(p)^2=E_a(p)^2 I`,

`E_a(p)^2=sum_i sin(a p_i)^2/a^2 + M_a(p)^2`.

There are two positive and two negative eigenvalues per face/momentum.
Since `M_a>=2`, `E_a>=2` everywhere, with equality at zero momentum. The
spectral projectors are `P_±=(I±h_a/E_a)/2`. At each corner with n components
equal to pi in dimensionless momentum, `E=2+2n/a`. Only the zero corner has
finite low energy as a decreases. More generally an energy bound `E<=B`
forces `sum(1-cos k_i)<=a(B-2)`, so k approaches zero modulo `2pi`;
then the sine term bounds physical momentum. Thus no additional finite-energy
valleys survive. Removing the Wilson term would retain all eight corner
Dirac copies. Their removal here is a regulator choice, not an erasure of
detector alternatives. Finite cutoff states and their probabilities remain
in the full finite Hilbert space.

The ordered Floquet route has a different branch structure. At its eight
corners its symbol is `(-1)^n exp(-2ia beta)`: four even-n corners share its
small eigenphases, and four odd-n corners lie near the quasienergy boundary.
This is a reason not to silently identify a principal Floquet logarithm with
the chosen Wilson generator. The receipt checks these corner identities.

## 3. Filled sea, particle and hole preparations

On V finite cells the one-body dimension is `20V`. Fill every negative
eigenmode and leave every positive eigenmode empty. This defines one pure
Slater determinant Omega with source pair number `N_sea=10V`. Its finite
energy is `E_sea=-10 sum_p E_a(p)`. The many-body spectrum follows directly
by diagonalizing the one-body operator and retaining every occupation:

`H_a-E_sea=sum_(positive modes) E b* b + sum_(negative modes) E d* d`.

Here a positive-band source annihilator is b, whereas a negative-band source
creator is the hole annihilator d. All terms are nonnegative. Omega is the
unique finite ground, the full-Fock excitation gap is exactly two, and the
nonvacuum sector of relative source charge zero has gap four. Both gaps use
periodic boundary conditions admitting zero momentum. They refer to this
massive free source, not a photon-containing theory or a Yang-Mills claim.

The relative conserved charge is `Q=N-10V=N_particles-N_holes`. A particle
adds one actual paired-flux object to the reference; a hole removes one.
Both have positive excitation energy. This gives a precise particle/hole
interpretation of the two signs of the one-body Dirac spectrum. The charge
is the declared source occupation charge, without electromagnetic or
Standard Model identification. In particular, a hole is not a negative
number of literal source defects. The underlying regional flux is Z2;
the signed charge counts deviations from the chosen many-body reference.

Actual local source charge is `Q_x=sum_(f,c)n_(x,f,c)`, whose reference mean
is ten. It is not a zero-or-one particle-position projector in the filled
sea. Excitation readouts use the same spectral projectors and reference.
The actual native `U(3)×U(2)` currents act trivially on these native-empty
link slots. The five degenerate face channels acquire no color or weak
charges from their count.

On one four-mode zero-momentum block, a finite preparation uses a fermionic
reservoir mode with its initial occupation recorded and probabilities
one half for each setting. A fermionic SWAP with source mode0 conserves
the total source-plus-reservoir number. Reading the reservoir gives four
complete Kraus operators:

`a_0*/sqrt2`, `a_0/sqrt2`, `(I-n_0)/sqrt2`, `-n_0/sqrt2`.

Their squared operators sum to I. Each preparation record has probability
one quarter on the block's filled sea. The add and remove records prepare
a particle or hole with energy two and relative charges +1 or -1. The
no-change records are retained; each has ground and particle-hole-pair
energy components. The complete particle/hole PVM in the receipt retains
all 36 preparation/occupation records, including 30 zero probabilities.
No success-conditioned distribution substitutes for the complete instrument.
The two occupation operators each have spectrum `{0,1,2}` on the entire
sixteen-dimensional block, and the nine joint projectors sum to its identity.
On the full lattice this detector is tensored with the identity on every
other momentum/face block. Completeness does not depend on the prepared
state's support. The reservoir SWAP is a driven preparation service here;
its number conservation does not assert conservation of the chosen source
Hamiltonian without an apparatus work register.

Charge conjugation is also exact for this free symbol. With
`B=alpha_2 beta`, using coordinate indices 1,2,3, the owned matrices obey
`B B*=I`, `B alpha_i* B^dagger=alpha_i` and
`B beta* B^dagger=-beta`. Thus
`B h_a(p)* B^dagger=-h_a(-p)`. The induced finite Bogoliubov particle-hole
map exchanges particle and hole records, reverses Q and preserves excitation
energy. The same identity holds in the continuum. This defines antiparticle
partners within this declared free Dirac field; it identifies no measured
lepton, quark or gauge representation.

In these matrix identities a star means entrywise complex conjugation;
the Hermitian adjoint is written dagger.

## 4. Controlled massive free-field limit

Let `h_D(p)=alpha.p+2beta`. On the physical momentum window `|p_i|<=P`,
the elementary sine and cosine remainder estimates give

`||h_a-h_D|| <= b(a,P)=a^2 P^3/2+3a P^2/2`.

Consequently `|E_a-sqrt(p^2+4)|<=b`, and Duhamel gives
`||exp(-it h_a)-exp(-it h_D)||<=|t|b`. Since both positive energies are at
least two, the normalized-symbol formula implies
`||P_-,a-P_-,D||<=b/2`. This estimate also controls the limiting reference,
rather than proving a one-particle equation on an unrelated state.

For a fixed finite physical torus, identify its common band-limited Fourier
modes at each sufficiently fine a. The filled-sea covariance is P_-, and
the CAR norm of a smeared field is its one-body L2 norm. Thus every normalized
two-point function on a fixed time window differs by at most
`b(1/2+|t-s|)`. Wick expansion bounds a normal-ordered 2n-field determinant
by `n n!` times the largest two-point error. All finite ordered words follow
by CAR reordering. These are convergence statements for the complete free
field and vacuum correlations, not just for a few plane-wave frequencies.

Smooth compact momentum smearings exhaust the one-body space. The CAR norm
bound extends convergence by density; the finite physical torus can then
grow, with the smooth Fourier sums converging to their R3 integrals.
The resulting quasi-free field has the massive Dirac equation, the filled
negative-band covariance and positive particle/hole Fock representation.
The mechanisms are those used by the owned
[free continuum construction](../../continuum_free_limit/README.md): vacuum
annihilation conditions, Wick uniqueness, graded locality and the positive
particle/hole spectrum. Its specific eight-valley centered source and its
interaction bounds are not transferred here. In this quadratic source the
vacuum projectors and dynamics converge directly with the displayed bound.

The continuum current obeys `partial_t rho+div j=0`, with
`rho=psi*psi` and `j_i=psi*alpha_i psi`; `|j.n|<=rho` follows from
`(alpha.n)^2=I`. Integrating outside an expanding ball proves propagation
speed at most one for compact initial data. The distributional Dirac
propagator therefore has causal support. Smeared CAR fields anticommute at
spacelike separation, and their even observables commute. The continuum
one-body mass shell is `E^2-p^2=4`; sums of positive particle/hole energies
satisfy the spectrum condition. The owned
[Clifford spin/current identities](../physical_attachment/PROOF.md) give
spin one half and Lorentz transformation of the continuum Dirac current.
The conserved current inner product on Cauchy surfaces transforms with the
surface as well as the spinor. It defines the unitary one-body Lorentz action;
a pointwise boost alone is not an equal-time unitary. The positive and
negative mass shells are separately invariant under the proper orthochronous
Lorentz group, so their vacuum covariance and particle/hole Fock lift are
covariant with this surface transformation.
The finite lattice itself has its cubic regulator symmetry, not exact boosts.

The finite-source claim is a local graded free field with controlled
massive continuum readouts. The choices of CAR grading, Wilson regulator,
filled-sea preparation and face-count coin coefficient are architecture.
No interaction survives this quadratic law; it supplies neither chiral
Standard Model matter nor its measured masses, coupling constants or GeV
scale. Coupling it to a different source and reference requires an explicit
operator and state transfer.
