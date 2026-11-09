# Spatial Fourier records and the physical momentum premise

## The source and the surveyed results

Use the [primary construction](../particle_common_contract/PROOF.md), its
physical polynomial density, all sixty CAR modes per cell, quantum links,
open boundaries, source Hamiltonian and idle conjugate reference. The
[regular spatial algebra](../particle_primary_spatial_algebra/PROOF.md)
has a point-norm continuous action beta of R3 commuting with its source
time action. It does not prove that the selected limiting state is
spatially invariant or that beta has a unitary implementation on its GNS
space. This distinction is material for a particle momentum readout.

The [periodic joint readout](../particle_joint_translation_readout/PROOF.md)
constructs conserved Fourier projectors of actual periodic translations,
including the conjugate reference. Its source and state differ from the
primary open family. `Lean/QFT/PeriodicMomentum.lean` propagates a supplied
commutation relation [H,T]=0 to Fourier polynomials; it does not construct
T on the primary source. The variational translation witness in
`Lean/Variational/TranslationInvariantComposedInstance.lean` is a separate
two-state source. Neither theorem's conclusion transfers by relabelling
primary positions. The construction below instead reuses the actual
[Hermitian probe dilation and finite filter](../particle_filtered_probe/PROOF.md).

## A complete finite spatial Fourier instrument

Choose M distinct admitted centered source positions x_j along one axis,
with spacing h, and the same physical plaquette probe at each position.
The selected instance has M=4 and h=1/4. The authenticated extra position
x_M supplies a boundary comparison. Member2 contains all five positions.
Their probes commute with every Gauss constraint. They can be the exact
unsmeared contractions, ideal fixed-time triangular filters, or their
declared finite positive quadrature/Cayley versions. The executable
primary plan uses the last, with one common width and error schedule.
Write these actual Hermitian contractions as B_j and set

```
C_(0,j)=B_j,       C_(1,j)=sqrt(I-B_j²).
```

The exact algebraic square root and the neutral dilation
`[[B_j,C_(1,j)],[C_(1,j),-B_j]]` are the declared apparatus services.
Prepare a neutral M-dimensional address in the uniform vector, apply the
corresponding controlled dilation, and perform the discrete Fourier
transform on the address. Retain the dilation bit b, every Fourier label
k and the full unchanged reference. The source Kraus operators are

```
K_(b,k)=(1/M) sum_(j=0)^(M-1) exp(-2pi i k j/M) C_(b,j).       (1)
```

Fourier orthogonality gives

```
sum_(b,k) K_(b,k)* K_(b,k)
 =(1/M) sum_j [B_j²+(I-B_j²)]=I.                              (2)
```

These are all 2M fine records, with original probabilities
`p_(b,k)=Tr rho K_(b,k)* K_(b,k)`. Every operator is physical; functional
calculus and neutral address operations preserve all source constraints.
Completeness is on the ambient carrier before physical restriction, so
arbitrary source/reference entanglement is retained. No translation of
the source or the reference is applied. An exact direct-service invocation
allocates an address of dimension M and one neutral bit, giving joint
carrier `2M D_source²`. Operator descriptions, square-root computation,
controller precision, compiler workspace and repeated records require
separate resources; no efficient full spatial execution is asserted.

For the blank branch define the actual spatial Gram matrix
`G_jl=Tr rho B_j* B_l`. Its complete original Fourier law is

```
p_(0,k)=M^-2 sum_(j,l) exp[2pi i k(j-l)/M] G_jl.               (3)
```

This measures coherent spatial phase dependence on the specified source.
A classical random choice of address removes off-diagonal interference
and generally gives a different output. The effects K* K need not be
projections. Formula(3) is therefore a finite-aperture wavevector readout,
not a proved spectral measure of conserved primary translation momentum.
It is valid without a translation-invariant preparation.

There is also a strictly positive result for the separate, unsmeared
version of this actual primary instrument. Let b be the physical
empty-Fock/constant-link vector and take B_j=A_j, the normalized real
plaquette trace. The [native character calculation](../particle_primary_survival/PROOF.md)
gives A_j b=(1/30) times the sum of ten orthonormal native and conjugate
character loops, so ||A_j b||²=1/90. Distinct plaquettes have orthogonal
blank images: an edge present in only one loop has a nontrivial matrix
coefficient whose independent Haar integral is zero. This remains true
for plaquettes sharing some links. Equation(1) therefore gives

```
||K_(0,k)b||²=1/(90M),
Tr rho K_(0,k)*K_(0,k) >= 2^-N/(90M)>0 for every k.           (3a)
```

Here rho>=2^-N P_G is the actual physical density floor and b is only
a matrix-element witness. It does not replace rho. This bound proves
original finite wavevector outcomes, without a momentum eigenstate or
a surviving limiting distribution. The fixed-time filtered instrument
has its own probabilities and errors; (3a) is not transferred to it
without a sufficient closeness estimate.

If every finite filter differs from its ideal contraction by norm at
most delta, square-root continuity gives a two-column difference at most
`sqrt(delta²+2delta)`. The controlled address operation acts on a uniform
orthogonal direct sum, so its column error is bounded by the maximum of
the individual errors, with no factor M. Fourier transformation is unitary.
The complete output half-trace difference is consequently at most
`min(1,sqrt(delta²+2delta))` for every retained reference. Additional
compiler and preparation errors add explicitly. This comparison uses the
same direct two-column instrument on both sides; it does not equate its
failure channel with that of a different selector apparatus.

## The actual open aperture and box boundary

Set q_k=2pi k/(Mh). Form the shifted aperture from B_1,...,B_M in the
same source. Direct reindexing, without any source symmetry assumption,
gives for the blank branch

```
K'_k-exp(i q_k h)K_k
 =exp(i q_k h)(B_M-B_0)/M,       norm<=2/M.                    (4)
```

Shifting r steps has error at most min(2,2r/M). At fixed physical
aperture width R=Mh, the one-step error divided by h is at most 2/R.
Thus making h small alone does not make a fixed-aperture observable an
exact momentum-transfer eigenoperator. Increasing the aperture demands
its own survival, record and error control. For probes in its admitted
spatially smeared regular menu, the spatial quotient relates ideal
translations to reindexing through its finite/infinite dynamics comparison.
This does not supply a continuous spatial action on every unsmeared point
probe, or remove (4)'s aperture boundary term.

The actual finite carrier is open. For L>2, a cyclic owner-cell relabelling sends
the native bond from (L-2,0,0) to (L-1,0,0) to the pair
(L-1,0,0),(0,0,0), which is not a native positive-neighbor bond.
The corresponding owned outgoing link reaches a ghost endpoint L,
not the interior endpoint0. Local incidence also differs: the origin
has no incoming links, while its shifted interior neighbor has one.
This permutation does not preserve the source graph or its Gauss
incidence. A periodic completion changes the program. No cyclic
translation momentum is assigned to this open Hamiltonian.

## Which spatial implementation is actually necessary?

Let (pi,H,Omega) be the GNS representation of the same limiting state
omega. Spatial invariance is sufficient for the canonical construction
`V_y pi(A)Omega=pi(beta_y(A))Omega`. Equality of the inner products makes
this map unitary; point-norm continuity gives strong continuity. It fixes
Omega and commutes with the canonical time group. Conversely, any
implementation fixing Omega implies omega beta_y=omega. Thus invariance
is exactly the premise for this vacuum-fixing construction.

It is not necessary for every possible implementation. A finite faithful
state rho on M_d has GNS space identified with Hilbert-Schmidt matrices
by `[A]->A sqrt(rho)`. For any strongly continuous unitary group U_y
on C^d, beta_y=Ad U_y is implemented there by
`V_y(X)=U_y X U_y*`. This is a genuine unitary group and preserves the
left observable representation, whether or not rho is invariant. Its
prepared vector sqrt(rho) is fixed precisely when [rho,U_y]=0.

The exact control uses rho=diag(3/4,1/4), U_y=exp(i y sigma_x) and
H_source=0. Its standard generator on column-vectorized matrices is
`P=I tensor sigma_x-sigma_x^T tensor I`. The mean ground excess is zero,
while the prepared vector has `Var(P)=2-sqrt(3)>0`. Therefore perfect
energy cooling alone does not establish momentum concentration, even
when an implementation exists. This control is not the primary source.

For the actual limiting representation, an alternative implementation
requires an explicit strongly continuous unitary representation V_y
implementing beta on the same pi, with the required commutation with
source time. Equivalences or normal extensions must include the actual
representation multiplicities and coherent group law; quasi-equivalence
of individual translated states alone is not used as a construction.
No such primary implementation is proved by the finite faithfulness
bound: its lower density floor vanishes with regulator, and finite
faithfulness does not imply faithful or separating limiting GNS vectors.
The existing evaluation-state counterexample in the spatial proof shows
that a continuous algebra action can fail to act on the selected GNS
null space altogether.

Even if an implementation is given, `exp(i k0.y)V_y` implements the
same algebra action and shifts its generator by k0. A vacuum or other
operational momentum origin must fix that freedom. For mixed input,
an invariant joint energy/momentum state gives response transfers;
separate initial momentum control is needed to infer final excitation
momentum. No concentration estimate for the actual primary preparation
is supplied by the Fourier apparatus or the cooling certificate.

The constructive result is the complete same-source finite instrument
(1)-(3), its reference-uniform approximation bound and the actual edge
identity(4). Turning these records into a conserved particle momentum
requires a same-representation implementation and an operational
attachment to its spectral measure, followed by a surviving common
energy-momentum branch. No Lorentz shell or physical mass follows here.

## Independent executed controls

For a two-level control take the stationary pair H=sigma_y,
rho=(I+sigma_y/2)/2 and four probes
(sigma_z,sigma_x,-sigma_z,-sigma_x)/2. The eight original probabilities,
ordered by branch and Fourier label, are
`(0,3/16,0,1/16,3/4,0,0,0)`. Each classical-address blank probability
would instead be1/16. Exact Kraus completeness and the entire unchanged
reference marginal are checked. The nonzero blank effects have
eigenvalue1/4 and are not projections. Additional controls check a
nonzero aperture edge, state-dependent interference, zero probes and
failed source/scope authentication. The small matrices are not a
numerical evaluation of the primary spatial density.
