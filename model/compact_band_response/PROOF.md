# Hopping work on the moving compact interacting band

The source is the compact opposite-degree texture with the complete native
internal current and memory construction of
[compact interacting response](../compact_interacting_response/PROOF.md).
Its four spatial modes carry 180 localized annihilators. The hopping drive
below is a specified simulator control; it is not identified with a physical
metric variation. All source masses, memory coefficients and the quartic
coupling are fixed while this control varies.

## A uniformly isolated moving band

Write the actual one-body operator as

```
h_s=s T+M V,             1/2 <= s <= 3/2, M >= 8.
```

Here T is the same open-cube centered nearest-neighbor hopping, V is the
compact Pauli texture with its two zero cores, and `||T||<=sqrt(3)<2`.
The parameter s changes hopping amplitude. It does not move the core
positions or change the shell degrees. The four zero eigenvalues of M V
give exactly four eigenvalues of h_s in `[-3,3]`; the complementary
eigenvalues have absolute value at least `M-3`. Thus the band projector
P_s has rank four and gap at least `g=M-6>0` from its complement, uniformly
in the containing finite cube.

In an orthonormal eigenbasis let p range over the four band states and q
over every complementary state. Differentiating `P^2=P` and `[h,P]=0`
gives

```
P'=P'P+PP',
<q|P'|p>=<q|T|p>/(E_p-E_q),
P P' P=Q P' Q=0.
```

The opposite off-diagonal block is its adjoint. This formula retains both
signs of the complementary spectrum. In particular,

```
||P'||^2 <= ||P'||_HS^2
 <= (2/g^2) sum_p ||Q T p||^2
 <= 8||T||^2/g^2 < 36/g^2.
```

The volume-independent bound `||P'||<=6/g` follows. It concerns variation
of the finite spectral projector, not spatial lattice refinement. The
complete numerical control at M=64 uses all 686 one-body rows, including
all 682 complementary modes. Its projector derivative norm is approximately
0.0206766, below the rational majorant `3/29`.

## Transport, pairing and the work operator

Set `L_s=[P'_s,P_s]`. Then L is anti-Hermitian, `||L||=||P'||` and
`[L,P]=P'`. The finite ordinary differential equation `U'=L U`, `U(1)=I`
has a unitary solution with `P_s=U_s P_1 U_s*`. Let `W_s=U_s W_1`.
This is a comparison of the source's spectral spaces. Evolution under a
time-dependent h_s alone does not automatically implement this transporter.

The source particle-hole matrix C obeys `h_s C=-C conjugate(h_s)` and
`C conjugate(P_s)=P_s C`. Differentiation gives the same compatibility for
P' and L, hence `C conjugate(U)=U C`. Consequently

```
W_s* C conjugate(W_s)=W_1* C conjugate(W_1),
A_s=W_s* h_s W_s,
A'_s=W_s* T W_s.
```

For the last equality, differentiating A produces
`W_s*(T+[h_s,L_s])W_s`. Its commutator term vanishes because h_s is block
diagonal and L_s off diagonal relative to P_s. Thus Kato transport keeps
the pairing form constant, while deriving the change of the compact kinetic
operator from the actual hopping. Every uniform native current commutes
with `U_s tensor I45`.

On the full fermion algebra let

```
N_s=dGamma(P_s tensor I45),
H_(r,s)=dGamma(h_s tensor I45)+pair(C tensor conjugate(K_r))
                         +lambda(N_s-90)^2.
```

K_r is the unchanged mass/feedback kernel on the complete memory record r.
The derivative of this same Hamiltonian on every record is

```
F_s=partial_s H_(r,s)
   =dGamma(T tensor I45)
     +lambda{N_s-90, dGamma(P'_s tensor I45)}.           (1)
```

The anticommutator in (1) is necessary. The interaction's localization
projector changes when the hopping changes. It cannot be frozen when
computing the work of the chosen source law. An exact eight-state CAR
control rotates one occupied band direction into its complement: the
omitted term has rank four and squared Hilbert--Schmidt norm four at unit
coupling. Direct differentiation of the noncommuting square and the
second-quantized rotation commutator agree exactly.

The positive and negative eigenvalues of P' occur in pairs and at most
four of each are nonzero. Therefore

```
||dGamma(P'_s tensor I45)|| <= 180 ||P'_s||,
||lambda{N_s-90,dGamma(P'_s tensor I45)}||
                         <= (194400 lambda)/(M-6).     (2)
```

These bounds are independent of the exterior volume. The complete hopping
force can be extensive; (2) bounds the added moving-interaction force.

On the transported localized algebra, N becomes the fixed 180-mode number
operator, the pairing form is fixed, and the exact intrinsic derivative is

```
J_s=dGamma((W_s* T W_s) tensor I45),       ||J_s||<=360. (3)
```

The localized/complement factorization and parity make all off-diagonal
one-body terms have zero expectation in an even product state. In a
comparison of two localized states tensored with the same even complement
reference, the complementary term cancels and (1) gives exactly their
contrast of (3). Thus the intrinsic work derivative and the full-source
derivative agree on this reference-subtracted comparison. No equality of
unsubtracted full and core energies is claimed.

For a differentiable drive s(t), unitary evolution under H_(r,s(t)) gives
`d<H>/dt=dot(s)<F_s>` on each fixed record Hamiltonian between switches.
Controller switches, instrument backaction and memory resets require their
own work terms. The Kato transporter is not an extra unrecorded evolution.

For an experiment at fixed s, if the complete-memory theorem bounds the
unread localized state error by delta, (3) bounds its hopping-force
expectation error by `360 delta`.
The complete-history bound gives `360 sqrt(2 delta)` for arbitrary final
history-dependent readouts with operator norm at most 360. These estimates
retain every record; they imply no uniform accuracy on a rare conditioned
history. The relevant time is the declared source clock. No averaging
estimate for an arbitrary changing s(t) is inferred from the fixed-s bound.

## Thermal response of the same driven Hamiltonian

For a supplied positive inverse temperature beta, the full localized Gibbs
reference is `rho_s=exp(-beta H_s)/Z_s`, using the same mean Hamiltonian,
compact kinetic A_s, memory pairing and quartic interaction. It is faithful
on the finite 180-mode Fock algebra. In the transported frame, Duhamel
differentiation and the trace give

```
F_thermal'=Tr(rho_s J_s),
S'=beta[E'-Tr(rho_s J_s)],
F_thermal''=Tr(rho_s J'_s)-beta C_(rho_s)(J_s,J_s),       (4)
```

where C is the noncommuting Bogoliubov--Kubo--Mori covariance. The explicit
second-derivative term matters for a nonlinear Hamiltonian path; a thermal
covariance alone is not its full mechanical response. A source Gibbs reset
and the temperature are preparation inputs, and finite memory records need
not individually preserve this reference.

The numerical control uses the exact reducing sixteen-state slice

```
H4=e1(n0-n1)+e2(n2-n3)+m(pair01+pair23)+u(N4-2)^2.
```

Both positive energies e1,e2 and their derivatives are calculated from
the same complete compact matrix, without dropping its spatial splitting.
At the inherited third lepton coefficient, `u=1/128` and `beta=64`, the
mean hopping force is approximately -0.00215762. The entropy derivative
is approximately -0.13837536. Centered differences of freshly rebuilt
source matrices check the free-energy, energy and entropy derivatives
independently. This is a slice of the mean Hamiltonian, not its complete
180-mode thermal state or the physical vacuum.

## Gravity domain

Equations (1)--(4) give an actual source deformation, its work, and its
thermal entropy relation on the compact interacting carrier. The parameter
s is an admitted hopping drive. A spacetime metric, physical clock and
definition of local stress have not been reconstructed from it. In
particular, neither a scalar mechanical force nor a stationary Gibbs
reference supplies a symmetric conserved Lorentz stress tensor, regional
geometric modular flow, an area coefficient or fixed-volume generalized
entropy stationarity. These results therefore supply finite source work
data without instantiating the Einstein theorem's joint physical premises.
