# A Gauss-constrained ground family on the spatial Wilson source

## 1. One action, full carrier and explicit boundary condition

Fix the spatial Wilson source of
[particle_spatial_ground_state](../particle_spatial_ground_state/PROOF.md):
sixty fermionic modes per cell, the native fifteen-state representation
of `G=S(U(3)×U(2))`, gauge-coordinate links, product per-link word cutoff
`ell`, spacing `a>0`, and fixed finite electric/magnetic coefficients.
The same Hamiltonian and source clock are

```
H = P_ell [ e0/(2a) sum_e C_e + dGamma(beta_D D_W)
                     - kappa0/a sum_p Re Tr_M15(U_p)/15 ] P_ell,
U(t)=exp(-itH).                                             (1)
```

All coordinate products are formed before compression. No Hamiltonian
coefficient is changed by imposing Gauss constraints or choosing a state.
The full Fock carrier and every admitted link word remain in the ambient
space. The constrained state has support in its physical subspace.

For a finite open cube, each cell owns its three outgoing links. Retain
Hamiltonian interactions whose entire support lies in the cube, and
retain the electric terms on all outgoing links. Impose local Gauss
constraints at every endpoint of these links. Endpoints outside the cube
have no matter register and have zero declared external charge. This is
a zero-external-charge boundary condition. It is not the assertion that
a subregion of a larger physical state has zero boundary charge.

The native full M15 representation has zero trace for all twelve Lie
algebra generators. Four Dirac copies therefore have no normal-ordering
background charge. Empty Fock matter and constant Haar link functions
form a nonzero physical vector. The fully occupied sixty-mode local
Slater state is also invariant, since its determinant character is
trivial. This reuses the explicit seed in
[FIELD_ROUTE](../particle_observer_outputs/FIELD_ROUTE.md)
and the exact source trace check. The seed is a finite input, without a
claim that it is the interacting ground.

## 2. Local Gauss projectors and their exact construction

Let `U_x(g)` act on matter at x and by the corresponding endpoint action
on every incident coordinate link. Actions at distinct sites commute:
when they share a link they act by left and right multiplication, which
commute. Matter factors have disjoint support, and their charge operators
are even fermion bilinears. The Lie generators satisfy

```
[G_x^a,G_y^b]=0   if x!=y,
[G_x^a,G_x^b]=i f^(ab)_c G_x^c.                         (2)
```

The components at one nonabelian site do not commute. The commuting
constraint objects are the site Haar projectors

```
P_x=integral_G U_x(g) dg,
P_x*=P_x=P_x²,     [P_x,P_y]=0,
P_G=product_x P_x.                                     (3)
```

With any positive definite Lie-generator metric, the kernel of the
positive operator `sum_(a,b) g^(ab) G_x^a G_x^b` equals the common
invariant kernel. In an orthonormal metric this is the positive Gauss
Casimir. Equivalently, for an arbitrary complete anti-Hermitian generator
basis T, the positive Gram operator `sum_a T_a* T_a` has exactly this
kernel; the Gram operator need not be called a central Casimir.
Since G is connected, Lie invariance is full group invariance. This
computes (3) through finite matrix nullspaces and Gram inverses.

The per-link native-word spaces are invariant under both endpoint
actions. Their Haar-orthogonal projectors commute with these actions,
so the constraints restrict exactly to the finite product cutoff.
Finite word inner products and endpoint generators are computed by the
invariant-tensor and Gram construction in
[the preparation proof, section6](../particle_field_preparation/PROOF.md).
`exact_kernel_projector` implements the kernel step for exact finite
matrices; `joint_projector` checks and multiplies commuting projectors.
No sampled Haar integral or acceptance filtering defines the state.

A reproducible ambient basis can be fixed without choosing a physical
basis: order sites lexicographically, then Dirac components and native
internal packet indices; use that order for the Fock occupation basis.
Order links by their owner and direction. For each link, order the
identity and native matrix entries/conjugates, then words up to `ell`
by length followed by lexicographic order.
Retain the first linearly independent words in this order and apply Haar
Gram-Schmidt with positive norms. All projector and preparation matrices
below are specified in this ambient basis. Any change of basis must be
applied to all operators and the conjugate reference together.

The executable nonabelian control uses the actual native weak doublet,
its coordinate-link coefficient space, one left matter mode and one
right hole in the determinant-one filled sixty-mode carrier. Its gauge
subrepresentation has factors
`V_left tensor V*_right tensor V*_link,left tensor V_link,right`,
of dimension16. The two computed site kernels have ranks4 and4; their
commuting product has rank1. The invariant vector is
`(1/2) sum_(a,b=1)^2 |a,b,a,b>`. All twelve native generators are used.
Omitting the right constraint leaves charged states. This checks the
actual endpoint algebra; it does not claim that this small subspace
reduces the spatial Hamiltonian.

## 3. The same Hamiltonian preserves the physical carrier

Electric Casimirs commute with left and right endpoint actions.
Coordinate hopping transforms the left and right matter indices
contragrediently to the link indices, so each endpoint charge cancels.
A closed plaquette trace cancels every internal endpoint action.
These are the local Ward cancellations in
[native_quantum_gauss](../native_quantum_gauss/PROOF.md), applied to the
coordinate representation in (1). That package's raising-only link
Hamiltonian is not substituted for the coordinate action.

Consequently every term of (1) commutes with every `U_x(g)`, and

```
[H,P_x]=[H,P_G]=0,
H_phys = H restricted to ran P_G.                      (4)
```

The physical carrier is nonempty by the explicit blank. Let its dimension
be `d_phys` and its lowest energy `E_phys`. The ambient source bound
`D<=2^(60V)451^(3V ell)<=2^b`, `b=60V+27V ell`, gives
`1<=d_phys<=D`, so `log d_phys<=b`. The same certified interval `[L,U]`
for H bounds its physical restriction. No new spectral estimates or
physical mass data are needed.

This physical space contains a nontrivial magnetic response. For a
plaquette with four distinct links and `ell>=1`, let
`chi_p=Tr_M15(U_p)` and `phi_p=Re chi_p/sqrt(5/2)`. Both the constant
blank and phi_p are Gauss invariant and belong to the admitted product
word cutoff. The five native irreducible characters and their conjugates
are pairwise inequivalent, as verified by the source charge menu and
the character construction in the preparation proof. Product Haar
integration therefore gives `<1,phi_p>=0` and `||phi_p||=1`.
On empty matter, the actual source matrix element is
`<phi_p,H 1>=-(kappa0/(15a))sqrt(5/2)`: electric and matter terms
annihilate the blank, and distinct plaquette terms are orthogonal by
integration over a link occurring in only one loop. Thus Gauss
projection retains nonzero field dynamics even at per-link degree1.
This is a holonomy excitation, with no protected-defect or particle
identification inferred from that matrix element.

## 4. Deterministic physical preparation and gap-free cooling

On the physical carrier define

```
rho_phys,beta = P_G exp(-beta H)/Tr(P_G exp(-beta H)),
Q = P_G [T_m(beta(U I-H)/(2s))]^s,
T_m(X)=sum_(j=0)^m X^j/j!,
rho_Q=Q²/Tr(Q²),          Psi_Q=vec(Q)/||Q||_F.         (5)
```

Here `s` is a power of two with `s>=max(1,beta(U-L))`. This is a direct
state specification. It does not mean prepare a full-carrier Gibbs run,
measure P_G and retain only accepted runs. The full source Hamiltonian
and the clock in (1) remain unchanged.

Both source and reference in `Psi_Q` lie in their separate Gauss kernels,
since `P_G Q=Q=Q P_G`. The reference transforms in the conjugate ambient
representation and has zero idle Hamiltonian. No reference charge is
needed to compensate a nonsinglet source preparation. The purification
has even joint fermion parity; source and reference can each contain odd
physical sectors. Source density stationarity follows exactly from
`[Q,H]=0`; the joint purification need not be stationary under `H tensor I`.

Let b0 be the source physical blank and its conjugate reference. The
coefficient `<b0,Psi_Q>` is real positive because the restriction of Q
to the physical carrier is positive definite. If the two vectors differ,

```
v=(b0-Psi_Q)/||b0-Psi_Q||,       R=I-2|v><v|.           (6)
```

Then `R b0=Psi_Q`. If they coincide use identity. Since b0 and Psi_Q are
separately source/reference Gauss invariant, R commutes with both gauge
actions independently. It preserves total parity and extends by identity
on the entire orthogonal complement. This is a complete deterministic
preparation unitary, with no rejected input branch. It may be an admitted
finite apparatus primitive or compiled with a separately stated output
error under the existing apparatus contract. Preparation work is not
asserted to vanish, and autonomous C/E cooling is not inferred.

The entropy proof in
[`GapFreeCooling.lean`](../../Lean/Thermodynamics/GapFreeCooling.lean)
applies to the eigenvalues of H_phys. The same scaled Taylor estimate as
[the full-carrier cooling proof](../particle_spatial_ground_state/PROOF.md)
therefore gives

```
Tr(rho_phys,beta H)-E_phys <= log(d_phys)/beta,
y=beta(U-L)/(2s)<=1/2,
r_m=y^(m+1)/(m+1)! / (1-y/(m+2)),
eta=min(1,2s r_m),
D_half(Psi_Q,Psi_phys,beta)<=eta,
Tr(rho_Q H)-E_phys <= b/beta+(U-L)eta.                 (7)
```

For a rational tolerance epsilon, choose `beta=2b/epsilon`, the least
admitted dyadic s, and the first Taylor order with
`eta<=min(1/4,epsilon/(2(U-L)))`. The scalar-width case is immediate.
This terminates without a gap oracle and provides the requested physical
mean-energy tolerance. Finite matrices can be exponentially large in
volume; neither the projector computation nor generic preparation is
claimed efficient. No dense full spatial matrix is allocated by the
plan calculator. Exact matrix prescriptions and finite-error compiled
preparations have distinct stationarity claims.

## 5. Physical ground inequality and the charged-operator boundary

For a stationary density supported on P_G with mean excess at most
epsilon, take a bounded A preserving the physical subspace. In particular,
every gauge-invariant observable qualifies. On that domain,

```
omega(A* [H,A])
 = Tr rho A*(H-E_phys)A - Tr rho A*A(H-E_phys)
 >= -||A||² epsilon.                                  (8)
```

The first term is positive because `A ran P_G` stays in the physical
carrier. Stationarity makes `rho(H-E_phys)` positive on its support,
which bounds the second term. The inequality does not apply to a bare
charged local creation or annihilation operator merely because its
input state is physical.

A charged experiment can attach an explicit boundary reference carrying
the opposite charge and matching fermion parity. Its total instrument
then preserves the combined constraint. However, the accessible combined
carrier may have a ground floor E_access below E_phys. The correct
energy excess is

```
Tr rho H-E_access
 = (Tr rho H-E_phys)+(E_phys-E_access).                (9)
```

A zero-energy compensating reference does not justify deleting the last
term. The native control in section7 has E_access=E_phys by an independent
full-sector spectral proof. No such equality is assumed for the general
spatial source. This preserves the distinction between physical neutral
observables and referenced charged-field probes.

## 6. A constrained spatial ground subsequence

Fix a, ell and the coefficients in (1), and take growing open cubes with
the zero-external-charge boundary condition. Use (5) with physical mean
excess epsilon_R tending to zero. Extend each finite state by empty matter
and constant link states outside its cube. Every global local Gauss
constraint annihilates this extension: interior constraints are included
in P_G; exterior factors are invariant blanks; ghost endpoints have the
explicit zero external charge. This extension is a state construction,
not a claim about the marginal of a larger ground state.

The local operator algebras are finite dimensional. A countable dense
set of matrix entries admits a diagonal weak-star subsequence, giving a
compatible limiting state. The fixed-cutoff Wilson interaction bounds
and local dynamics construction of
[particle_spatial_ground_state](../particle_spatial_ground_state/PROOF.md#6-a-fixed-cutoff-spatial-algebraic-ground-subsequence)
apply unchanged. The constrained finite densities are stationary for the
same finite H, so dynamics convergence passes stationarity to the limit.
For each fixed local positive Gauss operator C_x, its expectation is zero
in every sufficiently large constrained state and hence in the limit.
The limiting GNS cyclic vector is locally Gauss physical.

For a fixed local gauge-invariant observable A, all relevant interaction
and constraint supports lie inside sufficiently large boxes. Equation(8)
then converges to `omega(A*[H,A])>=0`. Thus the restriction to the invariant
observable algebra is a ground state of the actual Wilson coordinate
interaction. Spatial averaging and the same bounded boundary/volume
estimate give a translation-invariant choice if desired. No computable
subsequence, unique state or convergence rate is supplied.

Taking a region of this state preserves its boundary representation
labels and their correlations with the complement. It does not project
the region to the kernel of a box-only Gauss operator with its exterior
links erased. The full-region and source-plus-boundary algebras must be
specified before attaching regional modular or charged detector claims.
This is the same boundary-charge distinction used by
[the finite-code refinement proof](../particle_field_preparation/REFINEMENT.md#the-complete-finite-code-comparison).

This theorem fixes spacing and local word cutoff. It supplies no joint
spacing/cutoff limit, interacting Lorentz covariance, charged scattering,
physical particle masses or Einstein dynamics.

## 7. Complete native physical preparation and charged records

The actual directional control has208 ambient states and
`G=6(N-2)`. Its explicit neutralizing background is-12. This differs from
the zero-trace full M15 matter source; an empty native control is not a
Gauss-zero state. The exact source physical projector is

```
P_native = 1_(N=2) tensor I_flux,     rank P_native=6*13=78.
```

It commutes with the unchanged native Hamiltonian and retains all thirteen
flux labels. The orthogonal complement has130 dimensions. The complete
Gauss monitor `{P_native,I-P_native}` retains its outside-Gauss record,
which has zero probability on the exact direct preparation.

`prepared_control(beta=4,order=8)` evaluates (5) on every physical
reducing block and returns zero amplitudes on the ambient complement.
It returns all82 original block labels, so detector code can retain the
full ambient domain. The polynomial degree is1024, dyadic scale128 and
canonical half-trace error at most
`10604499373/47847722383835136`. Exact rational source arithmetic and
independent ground brackets certify

```
Tr rho_Q H-E_phys <= 114578301787/1099511627776 < 0.105.
```

The native complete compensated instrument is
`K_(+,i)=c_i*/2`, `K_(-,i)=c_i/2`, for all four Dirac components.
Its physical-input effects sum to P_native. Every bare charged output
lies outside P_native. Attach an orthogonal boundary-reference label
of charge `-6 sign` and odd parity to each branch; every combined output
then has total Gauss charge zero and even total parity. The source and
reference charge assignments are checked exactly. Both sign classes
have total probability1/2, and Pauli zeros remain records.

The boundary reference's idle Hamiltonian is declared zero. It can be
implemented as a finite reference register with blank charge0 and eight
charged output labels. The resulting charge-preserving even Stinespring
isometry extends within its finite symmetry blocks to a full unitary.
Its primitive status and work budget have the same apparatus contract
as (6); no physical free reference particle is inferred.

The full-sector native spectral brackets place the absolute minimum
in N=2, so E_access=E_phys for this control. The complete sixteen-output
Cayley/Fourier apparatus of
[particle_ground_spectral_records](../particle_ground_spectral_records/PROOF.md)
is executed on the direct physical preparation. It retains all mode,
sign and phase outputs, the preparation reference and the compensating
boundary reference. Its unnormalized probabilities sum exactly to one.
Two phase bins are not individual energy-line or mass measurements.
