# Complete finite charge classifications and a scoped defect catalogue

Three classifications have exhaustive domains. They do not enumerate all
protected defects of all OPH source realizations, nor identify physical
particles. Additional catalogue entries retain their owning construction's
strictly smaller completeness claim.

## Bare and ordered-cover sources

`model.completion.defect_algebra.report()` reconstructs the actual full C/E
branch maps, including zero and abort slots. Separating effects generate
individual source projectors. Connected transition trees then generate every
matrix unit. The branch star algebra is M78 on the full occupation carrier
and M66 on its binary restriction. The complete commutant is scalar in each
case. Thus no nontrivial source-only Hermitian charge commutes with every
admitted branch. This excludes common reducing sectors, not every invariant
subspace or a controller/history-dependent charge. The faithful-reference
argument for H66 additionally classifies source-only Heisenberg conservation;
it is not extended to the nonfaithful H78 stationary support.

`model.defect_protection.protection.algebra_certificate()` independently
reconstructs the ordered-cover maps. Effect separation, the transition tree
and one odd exchange loop generate exactly

```text
M66 tensor span{I,X} = M66 direct_sum M66.
```

Every branch commutes with Q=I66 tensor X. The full commutant is exactly
span{I,Q}, of complex dimension two. Its sharp sectors have Q=+1 and Q=-1;
the classical 132-state cover itself is connected. The spatial qutrit
realization is an intertwining realization of this same charge, not another
species. Its strictly local preparation carries a nonsharp excess, while
its sharp-charge exterior obstruction is retained. A different admitted
memory law can destroy this charge.

Both classifications are replayed from the source by the catalogue producer.
Other parent packages are bound and authenticated; their complete numerical
or pulse producers are not silently claimed to have been rerun.

## All Pauli-product charges of the selected spin source

Let G be a finite simple graph. Each edge has nonzero XX, YY and ZZ exchange
coefficients, with either the physical or Marshall-frame signs. Fields, if
present, are nonzero single-site X terms. Write a Pauli word modulo scalar
phase as P(x,z), with binary site vectors x,z. Conjugation by P changes only
the sign of each Hamiltonian Pauli monomial. Distinct monomials are linearly
independent, so [P,H]=0 is equivalent to commuting with every nonzero term.
For an edge (i,j), XX gives z_i=z_j and ZZ gives x_i=x_j. The YY condition
is their sum and adds no constraint. An X field at i imposes z_i=0.

Consequently each connected component contributes one global X generator
and, exactly when it contains no field, one global Z generator. These
generators exhaust the Pauli-product charges. Their number modulo phases is
2^(c+c0), where c is the number of components and c0 the number without a
field. This does not classify arbitrary sums of Pauli words commuting with
H: spectral functions and non-Pauli symmetry generators can also commute.

The actual source is reconstructed from
`model.native_ordered_ground.source.matrix()`, not a fitted spin surrogate.
In its Marshall frame it has twelve exchange edges and

```text
H_B = sum_edges(-XX-YY+ZZ-I)/4 - sum_sites X/8.
```

Every entry of its integer numerator is compared with that exact Pauli
polynomial. Exhaustion of all 4^8=65536 words agrees with the graph theorem.
At B=1/4 its only Pauli charges modulo phase are identity and global X.
Removing the field gives identity, global X, global Z and global Y. Thus
no nonidentity Pauli charge supported on a proper subset of the connected
cube is conserved. The statement does not exclude general local conserved
operators, metastable excitations or topology of a different order parameter.

On the one-occupation-per-site code, every complete inherited C/E repair
branch is identity times its recorded coin amplitude. Hence repair alone
protects the entire M256 algebra. Exchange evolution reduces the conserved
Pauli charges to the classification above. A local Z anticommutes with
global X and creates the opposite parity sector; it is not itself a
conserved charge. Its initial exterior equality, energy and complete
detector lifetime follow from `model/defect_localization_readout/PROOF.md`.
The source full-Fock lift includes the CAR sign on doubly occupied sites.

## The continuous charge and two measured creation channels

The Pauli-product census is not the full symmetry group. In the same actual
Marshall frame define Q=sum_x epsilon_x X_x/2, with epsilon=+1 on one cube
sublattice and -1 on the other. Each exchange term commutes with Q because
its two endpoint signs are opposite; the X selection terms commute as well.
`charge.py` checks [Q,H_B]=0 on the complete integer source matrix. It also
constructs the tensor-X spectral basis with every one of its 256 records.
The charge eigenvalues are -4,-3,...,4 with multiplicities
1,8,28,56,70,56,28,8,1. Global X is exactly exp(i*pi*Q), so its Z2 parity
is a subgroup of this U1 symmetry, not an extra independent defect type.

The charge is pulled from the signed two-flavor full-Fock source, where it
is sum_x(c_up^dagger c_down+c_down^dagger c_up)/2. Spin-independent signed
hopping and onsite occupation energy commute with this generator. The
complete C/E source branches also commute with it: empty and doubly occupied
sites are spin singlets, and admitted moves transport the internal spin.
The producer verifies all 6,976 full-Fock input columns for both laws, both
coins and every seam on two, three and four sites, retaining every raw
record. The even local generator commutes with spectators, so the two-site
identity extends to larger carriers. An independent full 256-state signed
Hamiltonian check and the exact spin/Marshall inclusion fix the source-to-Q
map. Local X selection terms also commute with this same charge. General
coherent creation or measurement controls need not preserve it.

A cube automorphism exchanging the sublattices preserves H_B and reverses
Q. The connected negative off-diagonal source matrix has a unique positive
ground, also certified by the complete eigenproposal bounds. Since Q
commutes with H_B, that ground has a sharp charge. The automorphism fixes
the positive ground and reverses its charge, forcing Q=0 exactly.

At site0, the local creation Z decomposes as

```text
A_+ = (Z - i epsilon_0 Y)/2,
A_- = (Z + i epsilon_0 Y)/2,
[Q,A_+] = A_+,  [Q,A_-] = -A_-,  Z=A_++A_-.
```

The two created vectors have sharp charges +1 and -1 and are orthogonal.
Their original probabilities are (1 minus/plus epsilon_0 <X_0>)/2. Both
exceed 4/25 on the certified actual ground. The local instrument

```text
A_q = Pi_(q epsilon_0 X_0) Z_0 = Z_0 Pi_(-q epsilon_0 X_0),
sum_(q=+/-1) A_q^dagger A_q = I,
```

prepares these coherent branches. Apply Z_0 and then measure only X_0,
retaining both outcomes and their original weights. This is the complete
1,092-pulse created experiment in `defect_localization_readout/execution.py`,
with its certified preparation and channel error. Normalizing a selected
branch requires its explicit herald probability. Each nonzero ideal branch
is orthogonal to the ground and has energy at least the same finite ground
gap. That bound is not a particle mass.

The numerical ground proposal is checked against the exact source. If its
normalized vector error is e, its probability outside the exact Q=0 sector,
or outside Q=+/-1 after Z, is at most e². The producer verifies these bounds
with exact rational sums and encloses both nonzero creation probabilities.
An independent charge-probability audit measures all eight commuting local
X operators and groups their 256 fine records by total Q. It has the correct
charge effects, including every zero-weight record. Its state update
dephases within each charge sector: grouping fine records is not the
coherent Lueders projection Pi_Q and does not prepare A_q Omega. Composing
the declared complete native local-pointer channels gives error at most
the sum of their eight channel errors, with all failures retained. This
package does not execute that full fine-readout tape. The coherent local
two-outcome preparation is the separately retained parent experiment.

These are internal source spin charges. A physical electromagnetic charge,
Lorentz spin, stable mass shell and Standard Model identity require the
same-source field and detector attachment. Two charge channels are not
asserted to be two physical particles.

## The same parity acts on the reconstructed field

Use the actual ordered source, state prescription and clock of
`model/spin_canonical_field/README.md`, with the positive-energy statements
in `model/spin_canonical_energy/README.md`. The native pi rotation about
internal axis1 preserves the zero-field source dynamics and the selected
state. It sends both existing scaled transverse quadratures Q_p(F),P_p(G)
to their negatives. It acts exactly on every finite bounded history and
preserves its expectation. Thus on the common subsequential history GNS
space, the rule R[A]Omega=[alpha_pi(A)]Omega is a well-defined unitary:
state invariance preserves all inner products. It obeys R²=I, R Omega=Omega,
and commutes with the actual time and spatial translations. The continuum
Weyl operators transform as R W(F,G) R=W(-F,-G).

Strong regularity gives R Q(F)Omega=-Q(F)Omega and the same identity for
P(G)Omega. The inherited source energy theorem proves these Stone vectors
nonzero with positive energy spectrum and explicit positive spectral weight
in admitted annuli. Thus the conserved native parity has an actual odd
field sector on the same state and clock. It is not merely an internal
label assigned to a separately supplied free field. A sharp particle pole,
Lorentz spin and a physical species do not follow from an invariant sector.
This parity acts on the full inherited canonical history net. The continuous
charge has the configuration-subnet construction below; no second canonical
pair is needed for that smaller algebra.

There is a complete bounded preparation of an odd field vector with a
uniformly positive original probability. Let F=G be a real unit-L2 annular
test, 0<m0<=m<=1/2, s=2/m0 and t=m0/4. Put
U=exp(i*2s*Q(F)), V=exp(i*t*P(F)), and u=<U>. The Weyl relation gives
||[U,V]Omega||=2 sin(m/2). Also

```text
[U,V]=(U-uI)(V-I)-(V-I)(U-uI),
||U-uI||<=2,  ||(V-I)Omega||<=t ||P(F)Omega||<=m0/8.
```

Since sin(m/2)>=m/3, the unitary variance satisfies
||(U-u)Omega||>=5m0/24. Moreover
1-|u|²<=2(1-Re u)=4||sin(sQ(F))Omega||². Consequently

```text
n0 = 25m0²/2304 <= ||sin(sQ(F))Omega||².
```

The complete instrument consists of sin(sQ) and cos(sQ), with effects
summing to identity. They create odd and even outputs respectively from
the invariant reference. No unknown vacuum expectation is needed in this
instrument. Each source finite version is a difference or sum of the
existing source Weyl controls; all outcomes and arbitrary references are
retained. The source energy-kick bound and the triangle inequality in the
energy form give limsup <sin(sQ_p) L_p sin(sQ_p)> <= h0=3s²=12/m0²,
where L_p=K_source/(Jp).

For a fixed positive scaled resolvent parameter lambda, the filter
lambda/(L_p+lambda) preserves parity. A measure of mass n and first moment h
obeys the Jensen bound

```text
integral (E+lambda)^(-2) dmu(E) >= n³/(h+lambda*n)².
```

The complete original probability of odd creation and filter success
therefore has liminf at least lambda²*n0³/(h0+lambda*n0)²>0. Replacing
n0 by n0/2 and h0 by2h0 gives a uniform lower bound for sufficiently small
p, without an effective threshold. Retain the even preparation outcome,
odd filter failure and odd filter success as three distinct complete
records. The filter's microscopic parameter is Jp*lambda, so its clock is
the existing source clock. No normalized successful state is substituted
for these original probabilities.

Orthogonality to the chosen even vacuum does not prove that every other
zero-energy vector is absent. All zero-energy detector outputs are retained
at the fixed sine amplitude above. A separate positive-energy existence
statement follows without that assumption: sin(rQ)Omega/r converges strongly
to Q Omega as rational r tends to zero. Applying any positive-energy
spectral projection with nonzero inherited Q weight proves that some
rational-amplitude odd sine preparation has nonzero weight in that same
positive-energy window. This supplies a source field excitation and a
complete readout criterion, not its numerical mass or a single energy shell.

## The continuous charge and conjugate configuration fields

The symmetry-breaking selection fixes internal axis1 and preserves rotations
about that axis. In the same zero-field ordered state, define

```text
Q3_p(F)=p^2 sum_x epsilon_x F(px) Sx^3,
Q2_p(F)=p^2 sum_x epsilon_x F(px) Sx^2,
P2_p(G)=-p sum_x G(px) Sx^2.
```

Q3 and P2 are the inherited canonical pair. Add the bounded exponentials
of Q2 to the common history compactness construction. Rotation invariance
gives Q2 precisely the same source variance and energy-kick bounds as Q3.
At equal time,

```text
||[Q2_p(F),Q3_p(G)]|| <= (p^4/2) sum_x |F(px)G(px)| = O(p),
[Q2_p(F),P2_p(G)]=0.
```

For bounded selfadjoint A,B the unitary Duhamel bound
`||exp(i(A+B))-exp(iA)exp(iB)||<=||[A,B]||/2` therefore controls
the linear combination of transverse configurations. The operator-norm
commutator also controls their group reordering. The inherited energy-form
transport proof propagates amplitude and test continuity through every
finite time history, including the added Q2 factors. After a further common
subsequence, Q2 and Q3 have strongly commuting equal-time unitary groups.
Their joint spectral calculus is defined. This does not require a common
domain for every unbounded field product. The original canonical history
representation embeds without changing any of its correlation functions.

Take the closed subspace generated from Omega by time histories of these
two configuration groups. It is reducing for the actual positive scaled
generator L and continuum spatial translations. The finite source rotation
alpha_theta, induced by `exp(i theta sum_x Sx^1)`, acts as

```text
alpha_theta(Q3)=cos(theta) Q3+sin(theta) Q2,
alpha_theta(Q2)=cos(theta) Q2-sin(theta) Q3.
```

The O(p) product bound passes this action to the limiting configuration
histories. Source-state invariance preserves their Gram inner products,
so `R_theta A Omega=alpha_theta(A) Omega` defines a unitary on this subspace.
Test continuity makes R strongly continuous in theta, and the exact source
group law gives R_(theta+phi)=R_theta R_phi and R_(2pi)=I. It fixes Omega
and commutes with L and continuum translations. This proves U1 covariance
of the configuration subnet. It does not assert U1 covariance of the entire
original Q3,P2 history net: rotating P2 introduces the additional uniform
momentum whose joint commutator requires a separate fluctuation estimate.

The actual configuration Stone vectors in this common representation obey

```text
Phi_+(F)Omega=(Q3(F)-i Q2(F))Omega,   R_theta Phi_+Omega=e^(i theta) Phi_+Omega,
Phi_-(F)Omega=(Q3(F)+i Q2(F))Omega,   R_theta Phi_-Omega=e^(-i theta) Phi_-Omega.
```

Their domains on Omega follow from the source second moments and the
inherited weak-differentiation proof. They carry charges +1 and -1 of the
same total internal spin symmetry as the local A_+ and A_- source channels.
This connects the source charge to field excitations, without identifying
the finite selected-cube state with an infinite-volume pole state.

There is also a source-induced charge conjugation. For a one-cell odd
translation e, the state-preserving automorphism beta_e=tau_e gamma2 combines
translation with the pi rotation about internal axis2. It sends Q2_p(F) to
-Q2_p(F shifted by pe) and Q3_p(F) to Q3_p(F shifted by pe). Test continuity
removes the vanishing envelope shift. On the configuration history subspace
the resulting unitary C satisfies

```text
C Omega=Omega,  C^2=I,  C R_theta C=R_(-theta),
C Q2 C=-Q2,    C Q3 C=Q3,    C Phi_+ Omega=Phi_- Omega.
```

It commutes with L and the continuum translations. C^2 follows from the
finite beta_e^2=tau_(2e) and the same vanishing-shift bound. It is an internal
charge conjugation of this limiting subnet, not a reversal of continuum
spatial momentum or of time.

Let mu_v be the joint positive-energy and momentum spectral measure of v.
The pi/2 rotation gives mu_(Q2 Omega)=mu_(Q3 Omega). Conjugation C makes the
mixed spectral measure between Q2 Omega and Q3 Omega zero. Consequently

```text
mu_(Phi_+ Omega) = mu_(Phi_- Omega) = 2 mu_(Q3 Omega).
```

Both charged vectors are nonzero, have no zero-energy atom, and inherit the
positive annulus weights and gapless accumulation proved for Q3 in
`spin_canonical_energy/README.md`. For a unit-L2 test supported in
`r/2<=|k|<=r`, each charged vector has weight at least `m0^2/r` in
`m0^2*r/16<=E<=12*r/m0^2`. These are equal spectral measures with explicit
nonzero bounds. They do not assert a sharp relativistic dispersion, Lorentz
spin, a scattering pole or a physical particle-antiparticle identification.
The source's broken internal symmetry supplies these collective charged
excitations; no winding number is assigned to them.

## Catalogue semantics and the particle join

The table separates exhaustive exclusions, memory charges, joint noiseless
subsystems, protected sectors, and prescribed topological backgrounds.
The hidden-spin qutrit and eta range are exhibited sectors, without an
asserted maximal joint noiseless algebra. The compact texture has genuine
opposite shell degrees for its declared order-parameter field; the source
does not autonomously select that background. Four band modes or eight
response lines are not a family or species census.

A complete catalogue of physically relevant defects requires a specified
refining source family and admitted dynamics, a definition of equivalence
and locality, and a classification on that domain. For each candidate the
same dynamics must provide preparation, charge readout, propagation,
spectral residue and one energy convention. A relativistic field attachment
must then establish physical spin/statistics, interactions and experimental
units. No catalogue row bypasses those conditions, and the physical particle
map and physical mass prediction count are empty and zero respectively.
