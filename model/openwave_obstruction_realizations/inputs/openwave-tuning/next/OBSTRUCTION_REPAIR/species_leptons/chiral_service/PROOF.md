# Chiral protection, electric dynamics and the full wall carrier

## Wilson symmetry audit

The frozen fractional service uses
`alpha_i=X_D tensor sigma_i`, `beta=Z_D tensor I`, and
`gamma5=X_D tensor I`. Its free cell-clock Hamiltonian is

```
H_W(k)=sum alpha_i sin(k_i)+beta[h m+sum(1−cos(k_i))], h=1/8.
```

The physical massless axial charge must commute with the Hamiltonian.
Here `[alpha_i,gamma5]=0`, but `[beta,gamma5]` has norm 2. At a single
cutoff corner and m=0, `[H_W,gamma5]` has norm 4. The anticommuting
grading `Y_D tensor I` reflects energy signs; it permits both beta m and
the Wilson beta term. It does not forbid an additive beta coefficient.
Colour, vector U1, cubic rotations and the energy-reflection grading all
permit that coefficient. No exact source symmetry protects m=0 against
the declared gauge interactions.

This is a statement about allowed operators, not a calculation of their
renormalized coefficients. Expanding a Wilson link produces a quadratic
gauge seagull with a beta component, but exchange contributions must be
included before inferring a net self-energy. The separate frozen
`mass_response` package rules out the tested leading bare-register
transplant by its positive compact-shell lower bound. Neither result
determines the full interacting pole.

The coframe source in `inputs/gravity_source.py` uses
`H_sigma=D_sigma+beta R^3 Q_sigma^2`. With the corresponding physical
axial gamma5, the second term also breaks naive axial conservation.
Its free gaplessness and coframe response supply no interacting lepton
protection. There is no gauge/current/clock intertwiner between that
source and this wall carrier.

## An exact electric counterexample for spatial overlap

The repository spatial-overlap construction has a modified q5 commuting
with its fixed-background fermion Hamiltonian. Restrict its momentum
symbol to one U1 holonomy theta. Put `gamma=−i beta alpha_1` and

```
V=−cos(theta) I+i gamma sin(theta),
H=beta(1+V), q5=gamma5(1−V)/2.
```

Then `[H,q5]=0` and `q5^2=I−H^2/4`. Its eigenvalues have magnitude
`|cos(theta/2)|`. On `L2(U1)`, let `E=−i d/dtheta` and let f be any
normalized constant four-spinor. Direct differentiation gives

```
[E^2,q5] f=−q5'' f,
q5''=−gamma5 cos(theta)/2+i gamma5 gamma sin(theta)/2,
||[E^2,q5] f||=1/2.
```

The physical electric coefficient a/2=1/4 therefore gives norm 1/8,
or 1/64 in the h=1/8 cell clock. Only fluxes −1,0,1 enter this witness.
It is an exact continuum-rotor calculation, with finite Fourier support.
A matrix Berry connection cannot make q5 covariantly constant because
its eigenvalues vary with theta; a commutator changes eigenvectors and
preserves the eigenvalue derivative. This agrees with the explicit
electric-term qualification in [Creutz, Horváth and Neuberger](https://arxiv.org/abs/hep-lat/0110009).
It excludes this attempted transfer of an axial symmetry. It does not
exclude every interacting overlap or domain-wall construction.

## The separate Euclidean action component

The 448-dimensional control has sixteen four-dimensional sites,
four Euclidean spin components and seven internal components. Its U3
representations are exactly

```
R_u(W)=det(W)^−1 W, R_d(W)=W, R_e(W)=det(W).
```

For `W=exp(−i theta/3) I`, their physical charges are +2/3, −1/3 and −1.
No relabeling of the electron unit occurs. The gamma matrices are
`gamma_E,i=−i beta alpha_i`, `gamma_E,4=beta`, and `gamma5=X_D tensor I`.
With the ordinary four-dimensional Wilson action and rho=1, define

```
H_kernel=gamma5(D_W−I), S=sign(H_kernel),
D=I+gamma5 S, hat_gamma5=−S.
```

On a nonzero kernel gap, `S^2=I` proves
`D gamma5+gamma5 D=D gamma5 D`. The independent variations
`delta psi=i epsilon hat_gamma5 psi`,
`delta barpsi=i epsilon barpsi gamma5` give an exact massless action
Ward identity. The mass insertion is
`D_m=(1−m/2)D+m I`; its variation is exactly
`2 i epsilon m barpsi gamma5(1−D/2)psi`.
This is the modified lattice symmetry of
[Lüscher](https://arxiv.org/abs/hep-lat/9802011).

The gauge current is the derivative of this complete D with respect to
the represented link, including the sign function. If H has eigenvalues
lambda_a and eigenvectors V, its exact derivative is

```
delta sign(H)=V [f_ab (V* delta H V)_ab] V*,
f_ab=(sign(lambda_a)−sign(lambda_b))/(lambda_a−lambda_b),
f_ab=0 for equal-sign eigenvalues.
```

This defines every link vertex, without replacing it by a local Wilson
vertex. The control applies a local electromagnetic gauge variation and
checks `delta D=i[Q_local,D]`. Thus the current divergence, U3 covariance
and axial Ward identity refer to the same Euclidean operator.

The sampled represented generators have norm at most 5/4096. The
telescoping plaquette bound is therefore 20/4096<1/40. The repository
admissibility theorem gives `H_kernel^2` in [1/4,64] and a degree-m
inverse-root approximation with range at most 2m+1 links and tail
`16(15/17)^(m+1)`. For the declared degree 192 the range is 385 links.
This reuses the hypotheses and polynomial estimate of the repository
overlap component and the smooth-background locality result of
[Hernández, Jansen and Lüscher](https://arxiv.org/abs/hep-lat/9808010).
Floating-point residuals and that analytic tail are different errors.

Every internal component here is a full Dirac field. Gauge anomalies
cancel between its left and right parts. The singlet axial Jacobian
contains `2 Tr gamma5(1−D/2)`, the index; it must not be discarded.
For three identical charged representations, a traceless family axial
generator has zero Jacobian because the family trace factors to zero.
The executed smooth background has index zero, to numerical precision;
that observation is not a restriction on all topological sectors.

The action has a valid finite chiral Ward component. It supplies no
reflection-positivity, Hamiltonian state, physical-clock or fluctuating
gauge-vacuum transfer to the fractional service. Full compact gauge
dynamics does not confine links to the admissible subset merely because
its magnetic potential penalizes rough configurations. These facts prevent
using the action control as a proof of the physical charged hierarchy.

## A domain-wall Hamiltonian on the declared clock

The source has S complete Dirac slices per physical cell, with
S=2,4,8,16 all retained. The extra coordinate is an internal regulator
coordinate. Spatial links and gauge transformations are identical on
every slice. Put `alpha4=Y_D tensor I`, `T=−(beta+i alpha4)/2` and
`P±=(I±gamma5)/2`. The massless Hamiltonian is the massless spatial
Wilson operator on each slice, plus T from slice s+1 to s and T* back.
The only input mass term is

```
H_boundary(m)=h m |0><S−1| tensor P+ beta P− + Hermitian conjugate.
```

The coefficients h=1/8, r=1 and the wall height implicit in T are
architectural choices. The boundary masses are the actual432 register
coefficients, with their existing ancestry. No observed mass, critical
mass, target residual or selected wall length enters this construction.
The use of a finite extra coordinate follows the mechanism of
[Kaplan](https://arxiv.org/abs/hep-lat/9206013) and
[Shamir](https://arxiv.org/abs/hep-lat/9303005). Its specific real-time
Hamiltonian is defined by `source.py`; equivalence to their full
interacting Euclidean formulations is not asserted.

At zero spatial momentum define
`E psi=|0>P+ psi+|S−1>P− psi`. Since `T*P+=TP−=0`,

```
E*E=I, H_DW(0,m)E=E(h m beta),
E* dH_DW/dk_i(0,0) E=alpha_i.
```

The remaining free zero-momentum states have energies ±1 in cell units.
The actual electron, muon and tau register poles therefore embed on
exactly the h clock, without changing their ratio. A mass discrepancy
in the register is not repaired by this free transfer.

Each inherited vector and axial probe, replicated on all slices,
compresses through E to its original four-component matrix. This includes
`V=(I,alpha_i)` and `A=(gamma5,2 spin_i)` and hence their left combination.
The same statement holds for the local physical colour-spin auxiliary,
since its spin/colour matrices commute with P±. It is an operator
embedding and supplies no interacting baryon state.

Under a local U3 transformation, every slice transforms with the same
R_u and R_d; the fifth hops and boundary mass are gauge scalars. The
electron uses det(W) with the same physical theta and the same h. The
quark charged weak probe retains its determinant character, and the
electron factor compensates it in the neutral weak contact. Spatial
gauge currents are derivatives of all replicated link hoppings,
including their Wilson contribution. The displayed axial probes are
specified observables; they are not substituted for a conserved Noether
current at finite S.

For a free momentum k, the kinetic alpha term anticommutes with all
fifth-coordinate and Wilson beta terms, so a zero requires every
sin(k_i)=0. At the origin there are four zero states. At any other
corner `W=sum(1−cos(k_i))>=2`; the fifth-hop operator has norm 1, hence
the gap is at least W−1>=1. There is one low vectorlike Dirac field per
internal representation. Every bulk mode and both boundary chiralities
are kept. The complete 64-momentum spectra and boundary probabilities
for each S are in `ALL_BANDS.json`.

## The finite-wall Ward identity and residual bound

Define the axial wall charge by +1 on the first S/2 slices and −1 on
the rest, tensored with identity on spin and internal coordinates.
It compresses to gamma5 through E. It is independent of gauge links,
so electric energy, all spatial gauge hoppings and the local colour-spin
auxiliary commute with it exactly. Its only violation is the fifth hop
across the midpoint, together with the explicit boundary mass.
Thus `dQ5/dt=i[H,Q5]` is an explicit wall-flux operator. The electric
commutator that obstructs spatial overlap is absent. The midpoint flux
is retained in the complete Hilbert space.

For a controlled scalar beta shift delta at zero momentum, with
0<=delta<1 and zero boundary input mass, the normalized trial map is

```
E_delta psi=N sum_s |s>[delta^s P++delta^(S−1−s)P−]psi,
N=sqrt((1−delta^2)/(1−delta^(2S))).
```

The recursion cancels at every interior slice, leaving

```
E_delta* E_delta=I,
E_delta* H^2 E_delta=N^2 delta^(2S) I,
E_delta* H E_delta=N^2 delta^S beta.
```

The min-max principle gives four eigenvalues with absolute value at
most `rho=N delta^S`. The other eigenvalues have absolute value at least
1−delta by the perturbation bound from the delta=0 bulk gap. The complete
spectra and bare-boundary weights are retained for delta=1/4 and 1/2
and every declared S. This is a source-defined suppression mechanism
for a scalar wall shift. A fluctuating gauge vacuum includes noncommuting
spatial kernels and rough links, so it requires its own bulk/mobility-gap
and boundary matrix-element estimates. This scalar bound cannot be
applied to an unknown gauge self-energy as though it were that estimate.

## Complete small quantum-rotor control

The control retains four CAR modes: two slices, two Dirac coordinates,
one physical spin and one unit-charge representation. It couples a loop
holonomy to the spatial sine and Wilson terms and uses electric energy
E^2/4. All sixteen occupations and all fluxes in [−K,K] are retained,
with K=4,8,16. The initial state has flux zero and occupation word 3;
it is a specified two-particle preparation, not a vacuum computation.
Each fixed-number block admits a matching static reference charge for
the one-vertex loop Gauss law. This abelian loop is a control, not a
proved reducing sector of the full U3 link vacuum.

The exact matrix Ward identity equals the lifted midpoint flux and the
electric-axial commutator vanishes. At physical time 1/64 the axial mean
changes from 2 to approximately 1.969113280953057. All occupation/flux
probabilities, including zeros, are in `ROTOR_RECORDS.json`. Flux-changing
terms have total norm at most 64, while the other terms preserve flux.
Leaving the cutoff from zero needs K+1 insertions. The full recorded
law therefore has total-variation error at most
`2 sum_(n>=K+1) 1/n! < 6/(K+1)!`, capped at one. The K=16 bound is
approximately 1.69e−14. Tiny reported boundary probabilities are
diagnostics; the factorial estimate is the stated truncation control.

## Resources, anomalies and physical limits

The quark source uses 24 S complex CAR modes per cell. Adding all three
charged-lepton families uses another 12 S. At eight physical cells the
complete Fock dimension is `2^(288 S)`, before the 25-state auxiliary
at each cell and all 24 full `L2(U3)` spatial links. Neutrino fields and
a dynamical chiral electroweak gauge sector are additional services.
The primitive remains local in physical space and nearest-neighbour in
the regulator coordinate, apart from the explicit mass contact between
its two boundaries at the same physical cell. No full native compiler
or vacuum preparation is supplied for this resource contract.

For each finite S the gauge electric operator on the compact link
configuration space has compact resolvent. Fermion hoppings, auxiliary
couplings and magnetic potentials are bounded operators on the finite
matter/auxiliary factors. Their sum is self-adjoint with compact
resolvent. The Gauss subspace is nonempty: filling two Dirac-spin
coordinates in each slice and both quark flavours yields the local
character det(W)^−2S, cancelled by the declared reference det(W)^2S.
Three equally filled electron families add det(W)^6S, so the combined
reference is det(W)^−4S. These are physical charge references, with
Gauss generators and the electric metric inherited from the fractional
service. They do not establish a physical ground-state particle spectrum.

Every finite carrier is vectorlike and has no gauge anomaly from removing
one chirality. The finite wall-flux Ward term is part of the complete
axial accounting. A continuum axial anomaly or a chiral Standard Model
measure requires its own limit and multiplicity argument; neither is
inferred from a boundary compression. The same physical bare electric
metric is used, but fermion vacuum polarization can change its measured
coupling. No transfer of a Maxwell vacuum or an observed fine-structure
constant is asserted.

Chiral protection, even when established, only forbids the relevant
additive terms within its Ward hypotheses. It does not determine the
finite mass-dependent corrections to charged-particle pole ratios.
The muon residual, its common state/reference/clock pole map, the joined
coframe regulator and the interacting wall residual bound are separate
requirements. `RESULT.json` marks the corresponding physical acceptance
fields false.
