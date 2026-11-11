# Scalar readouts on the native quark and lepton source

The source is the full thirty-mode native fermion algebra tensor the
five-state scalar memory. Its three literal quantum interaction tensors
and scalar occupation cost all have coefficient one. The reference is
the actual local physical Gauss-ground doublet divided by two, with clock
`exp(-it H)`. The full physical ground theorem and its complete interacting
probe images are the owning results in `chiral_joint/local_reference`,
`native_quark_response` and `chiral_joint/interacting_charged_response`.

This reference is the ground in the source-only Gauss-singlet domain.
Charged insertions require retained compensating references and use the
full Hamiltonian in their charged targets. No unconstrained-ground or
nonnegative charged-excitation hypothesis is substituted. Ground selection
is the declared preparation law; no preparation probability from arbitrary
input is supplied here. Probabilities below are on that specified input.

## Complete scalar readout and normalization

Write `H=N+V`, `U=2N-I`, `lambda=e(e-1)` and `Delta=1-2e`, where `e` is
the actual physical ground root. The exact source projector gives

`n=<N>=(1-1/Delta)/2`,
`rho_scalar=diag(1-n,n/4,n/4,n/4,n/4)`.

The code computes every scalar population from the actual source basis
and positive Gram metric. The four charged scalar states are retained.
For `phi_k=|0><p_k|+|a_k><0|`, the scalar Gram is

`<phi_i phi_j^dagger>=<phi_j^dagger phi_i>=C delta_ij`,
`C=1-3n/4`.

Each highest-weight scalar orientation has a complete target of dimension
26 inside a 130-coordinate weight space. The scalar map, its adjoint in
the actual metrics, the full Hamiltonian image and the complete weak and
spin kernels are checked exactly. The thirteen-root target annihilator
is evaluated on the entire target. Ground-field moments reconstruct the
complete rational resolvent. Exact gcds retain zero residues; rational
root isolation and outward number-field evaluation give every residue
with its original weight. All moments agree between scalar orientations;
the actual weak current relates the two doublet components.

For every complete scalar target energy `epsilon`, retain addition and
removal weights `C+(epsilon)` and `C-(epsilon)`. Their retarded sums give

`sum(C+-C-)=0`,
`Z_H=sum (epsilon-e)(C++C-) = 2-(7/4)n-(3/4)e`.

The source value is approximately `9.58523463744285`. This number is
computed on the full interacting reference, without importing the
lepton-only source's normalization. The two complete charged spectral
measures independently reproduce this formula.

The full CAR/scalar polynomial calculation also proves

`D+D^dagger=8I-10N-3V`,
`D-D^dagger=[H,sum_k[phi_k,phi_k^dagger]]`,
`D=sum_k[phi_k,[H,phi_k^dagger]]`.

Ground stationarity eliminates the second commutator in expectation;
weak invariance divides the first identity equally between the two
scalar components. This independent computation gives the same `Z_H`.
Thus the full retarded scalar response begins with `Z_H/z^2`.
The finite spectral convention is `phi_canonical=phi/sqrt(Z_H)`;
no relativistic kinetic coefficient or selected pole residue is inferred.

The actual scalar probe uses Kraus operators `phi_k/2` and
`phi_k^dagger/2`, with their charged quantum references. The full complement
is `sqrt(3) N/2`. Each component/orientation branch has total probability
`C/4`; the no-operation probability is `3n/4`, so all five branches sum
to one. Every charged energy residue carries the factor one quarter.
The no-operation branch has weights `3n^2/4` at energy e and
`3n(1-n)/4` at energy `1-e`; all other energy weights are zero.
This follows because N preserves the complete scalar pair containing
the ground. No successful energy branch is normalized separately.

## Full contact tensors

For all thirty native modes, calculate the actual anomalous insertion

`K_ij={ [c_i,H],c_j }`.

Exact normal ordering checks all 900 pairs on the entire Fock algebra.
There are 56 nonzero components, including every color and both scalar
components. Each insertion lies in the four-dimensional scalar span
`(phi_0,phi_1,phi_0^dagger,phi_1^dagger)`. The code computes the mixed
three-point traces from the literal insertion and the independently
computed scalar density. Dividing by the full scalar Gram `C I4`
reconstructs K exactly as a full-domain operator. All nonzero primitive
coefficients have modulus one.

The complete fermion addition/removal normalization is `I30` by CAR.
Consequently each nonzero contact component has modulus `sqrt(Z_H)`,
approximately `3.09600300992148`, in the stated scalar convention.
This source gives a universal contact modulus. Channel-dependent
dynamical responses below use the same primitive tensors and same state.
The contact insertion is the exact equation-of-motion vertex; a general
frequency-dependent amputated vertex is outside this calculation.

## An actual charged scalar pulse

Adjoin a blank and a dual scalar doublet, with reference creation maps
`R_k^dagger=|k dual><0|`. The pulse generator is

`X=sum_k(phi_k^dagger R_k^dagger+phi_k R_k)`.

For channel r, write its literal interaction as
`V_r=sum_k(P_rk phi_k^dagger+P_rk^dagger phi_k)`. Depending on the actual
Higgs orientation, P contains a pair annihilator or its adjoint. These
orientations and signs are extracted from the full source polynomial.
For color multiplicity `d_r=3,3,1`, each `||P_rk||<=2d_r` follows from its
literal CAR words. Define

`B_r=sum_k P_rk R_k`, `A_r=i(B_r-B_r^dagger)/(8d_r)`.

The norm bound is one for every channel. The pointer effects
`(I+A_r)/2` and `(I-A_r)/2` give a complete positive instrument with a
retained source and reference. All twelve global gauge-current commutators
vanish for the pulse and each pointer. These are covariance identities
on the finite cell and its reference; no dynamical link constraint is
inferred. Direct full-polynomial multiplication
proves the blank-reference identity

`<0|i[X,A_r]|0>=V_r/(8d_r)`.

This is a charged Higgs creation/annihilation intervention. Neither scalar
occupation cost nor a coefficient in the interaction is differentiated.
The two pointer probabilities are one half at zero pulse, and their
derivatives are the opposite halves of the calculated expectation slope.

Let `M1_r=sum_j w_rj(epsilon_rj-e)` be the full unconditioned elementary
singlet removal moment, with its one-half failure weight retained. Exact
normal ordering gives

`V_r=-sum_(spin,color in r)(c^dagger[H,c]+h.c.)`,
`<V_r>=-4d_r M1_r`.

The complete authenticated native quark and lepton spectra therefore give

`d< A_r >/dtheta = -M1_r/2`.

The source values are `-0.210891457185039` for each quark channel and
`-0.198310442473788` for the charged-lepton channel. The color average is
part of the same pointer normalization. The sum of all three derived
interaction energies equals `<V>=e-n` exactly. Their difference comes
from the native interacting records, without a per-flavor parameter table.

The complete pulse also has an exact finite-angle readout. Let `N_p` and
`N_a` select the two particle and two antiparticle scalar states. Set
`V_rp={N_p,V_r}`, `V_ra={N_a,V_r}`, so `V_rp+V_ra=V_r`. The literal
blank-reference blocks of `exp(-i theta X)` are

`K_0=diag(cos(sqrt(2)theta),1,1,cos(theta),cos(theta))`,
`K_k=-i sin(sqrt(2)theta)|p_k><0|/sqrt(2)-i sin(theta)|0><a_k|`.

The code verifies their differential equation, initial value and complete
Kraus sum. Full operator multiplication gives

`<A_r>_theta=[sin(sqrt(2)theta)<V_rp>/sqrt(2)+cos(sqrt(2)theta)sin(theta)<V_ra>]/(8d_r)`.

Both channel energies are evaluated on the actual ground using the
complete native removal maps. For a representative singlet,
`K=c^dagger H c-c^dagger c H` is formed with the target Hamiltonian and
the actual metric adjoint. Averaging over the two lepton spin removals,
or using the scalar quark spin and color symmetry, gives
`<V_rp>=-2d_r <{N_p,K+K^dagger}>` and the analogous antiparticle expression.
Their sum agrees exactly with the independent full spectral moment.
At the declared finite angle `theta=pi/(2sqrt(2))`, the receipt records
both probabilities `(1 +/- <V_rp>/(8d_r sqrt(2)))/2` with outward rational
bounds. The small-angle slope and complete finite pulse share one source,
reference and pointer.

## Neutral channels and the full32 spectator extension

The minimal R15 source has no gauge-singlet right-handed neutrino.
On a declared real neutral scalar orientation, its actual antisymmetric
contact matrix has rank 28 and two neutral Weyl coordinates in its
kernel. This is an algebraic slice; the invariant reference has zero
charged scalar expectation and does not select that background.

The native gauge-singlet operator
`J_spin=phi_0 c_(spin,11)-phi_1 c_(spin,10)` is also evaluated on this
same interacting ground. Its complete retained receipt supplies all nine
energy outcomes, all four probe settings, the positive failure operator
and the zero outside-target record. The actual neutral local gap is
approximately `0.0895553` in the source clock. The full CAR equation for
J contains nonzero cubic terms; its gauge-singlet label therefore does
not imply a free or massless neutrino. The package authenticates the
owning receipt and derives the complete addition/removal normalization
from its exact cross moments, with every factor of setting probability
and success probability retained.

The empty five-level occupation supplies a gauge-singlet weight. A
separate declared spectator lift takes `R16=R15 direct_sum 1` and
`H32=H30 tensor I_sterile_Fock`, with spin-major sterile modes 15 and 31.
Its reference is `(P_ground30/2) tensor I4/4`. The two sterile annihilators
commute with H, every sterile contact vanishes, and each complete
addition/removal measure has weight one half at zero energy transfer.
The scalar response is unchanged by the tensor factor. This is the
explicit decoupled extension, without an interaction inferred from a
shared charge label. An interacting Dirac or Weinberg extension requires
its own source law; the optional coefficient-bearing Weinberg artifact
is excluded from this calculation.

These are local observer responses and exact source laws. No generation
copies, observed neutrino flavors, physical masses or oscillation
parameters are assigned by the calculation.
