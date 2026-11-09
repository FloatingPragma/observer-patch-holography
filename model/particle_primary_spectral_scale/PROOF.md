# Spectral scale of the primary Wilson family

## Source and the two comparison domains

Use the [primary contract](../particle_common_contract/PROOF.md), including
its open box, full sixty-mode matter carrier, product link-word cutoff,
physical Gauss constraints, polynomial density Q²/Tr(Q²), unchanged H
and idle conjugate reference. Its spacing is a=2^-n, side L=2^(2n+1),
volume V=L³ and link degree ell=n+1. Every source term has a factor 1/a:

```
Hhat_n=a_n H_n,       epsilon_n=2^(-n-4),
Tr rho_n(Hhat_n-Ehat_0,n)<=a_n epsilon_n=2^(-2n-4).     (1)
```

Hhat is notation for scaling the same operator, not a replacement clock.
Source time t corresponds to dimensionless time t/a. The code verifies
the actual member before producing any source-bound result.

A finite particle comparison may use one specified member/capacity, with
its actual response, original weight, controlled finite-domain dispersion,
particle-channel identification and stated cutoff/approximation errors.
A declared SI comparison scale can convert such an admitted rest-energy
interval to GeV. A full continuum construction is not a prerequisite for
that finite, error-certified comparison. A claim that the same observable
survives the continuum family additionally needs uniform original weight
in a fixed physical band and the corresponding refinement estimates.
Neither kind of rest-mass result is supplied here.

## Exact decoupled anchor on the actual carrier

Separate the dimensionless source into

```
Hhat(lambda)=H0+lambda Vpert,
H0=3 sum_x dGamma(beta_D)_x+(105/2) sum_links Casimir_link,
Vpert=dimensionless Wilson hopping minus (1/7) sum_p A_p,
A_p=P_ell ReTr_M15(U_p) P_ell/15.                     (2)
```

Only lambda=1 is the primary action. The interpolation is a diagnostic
on its actual carrier; the selected polynomial density is not replaced
by its lambda=0 vacuum. Fill all fifteen internal modes in both negative
beta_D rows at each cell and take constant link functions. Call the
normalized vector Omega0. The determinant singlet argument in the
[primary survival proof](../particle_primary_survival/PROOF.md) makes
Omega0 physical, including the exterior constraints. It is the unique
ambient ground of H0 and has energy E00=-90V.

The [actual trace metric](../particle_field_preparation/preparation.py)
gives, for an irreducible representation with SU3 highest weights p,q,
SU2 spin j and integer Q6 charge z,

```
C=(p²+q²+pq+3p+3q)/6+j(j+1)/2+z²/120.              (3)
```

A nontrivial SU3 representation contributes at least 2/3, and a
nontrivial SU2 representation contributes at least 3/8. If both are
trivial, a character of S(U(3) x U(2)) is an integer power of det U3.
On diag(exp(-2it)I3,exp(3it)I2) this character has z in 6Z, so its
nontrivial Casimir is at least 3/10. These statements cover every
Peter-Weyl summand, not only the degree-one list. The product word space
is invariant under the group and its Casimir, so compression retains
these eigenvalues. The native one-dimensional z=6 block attains 3/10
and is present at every selected degree. Thus every nonconstant link
costs at least 63/4 in H0.

Each matter occupation defect above the filled negative-beta state costs
3. The ambient gap is consequently 3. A physical vector of energy less
than 6 would have constant links and exactly one matter defect. Every
native block is a nontrivial local gauge representation, so such a vector
is not Gauss invariant. A particle-hole singlet in one native block at
one site, summed over its internal index between a positive and negative
Dirac row, is physical and costs exactly 6. Therefore

```
gap(H0 on the physical carrier)=6.                    (4)
```

This is the decoupled anchor gap. It is not the gap of Hhat(1).

## Exact local hopping and plaquette denominators

For a bond in direction mu the one-body hopping block is
`S_mu=(-beta_D+beta_D gamma_mu)/2`, with its adjoint on the reverse
orientation. Let P_-=(I-beta_D)/2. On the two-site filled negative-beta
Slater vector only the block `(I-P_-) S_mu P_-` and the reverse directed
block create particle-hole states. Terms preserving the occupied
subspace annihilate this vector by exclusion. The actual Clifford
matrices give a two-site particle-hole norm squared of one per internal
column, with onsite excitation 6. This is checked by the literal
one-body Slater contraction and an independent finite CAR action.

The native link is unitary pointwise and P_- acts only on Dirac indices.
The Slater contraction is therefore independent of the link coordinate.
Integrating the Haar variable preserves that norm. In native block r,
the link matrix coefficients have Casimir C_r. They lie inside every
selected word cutoff before the final compression. Thus the contribution
of one bond to Vpert Omega0 has the exact orthogonal packets

| Q6 | d_r | C_r | Squared norm | H0 excitation Delta_r |
| --- | --- | --- | --- | --- |
| -4 | 3 | 4/5 | 3 | 48 |
| -3 | 2 | 9/20 | 2 | 237/8 |
| 1 | 6 | 21/20 | 6 | 489/8 |
| 2 | 3 | 7/10 | 3 | 171/4 |
| 6 | 1 | 3/10 | 1 | 87/4 |

In particular Delta_r=6+(105/2)C_r. Both bond directions and all Dirac
and native columns are included; the total squared norm is 15.
This uses the entire number-preserving Fock action, rather than a
one-particle replacement for the Hamiltonian.

The plaquette character identity in the survival proof applies with
this filled physical Fock vector as well as with its empty matrix-element
witness. The magnetic contribution on each plaquette is
`-(1/210) sum_(r,sign) L_(r,sign)`. Its ten normalized orthogonal loop
vectors have H0 excitation `210 C_r`. A bond packet changes local
occupations, while a loop packet does not, so these two classes are
orthogonal. Different bond packets have distinct particle/hole positions.
Different plaquettes have a link appearing in only one loop, whose
nontrivial Haar matrix coefficient has mean zero. Consequently their
packets are orthogonal as well, including adjacent plaquettes.

## Finite perturbation coefficients and the volume obstruction

Let R0 be the inverse of H0-E00 on Omega0-perpendicular and zero on
Omega0. Its use here requires only the displayed packets. Define

```
b1=sum_r d_r/Delta_r,       b2=sum_r d_r/Delta_r²,
p1=sum_r 2/(210²*210 C_r),  p2=sum_r 2/(210²*(210 C_r)²).
```

The actual open box has B=3L²(L-1) bonds and P=3L(L-1)² plaquettes.
All coefficients in the following identities are exact source numbers:

```
<Omega0,Vpert Omega0>=0,
||Vpert Omega0||²=15 B+P/4410,
<Vpert Omega0,R0 Vpert Omega0>=B b1+P p1,
||R0 Vpert Omega0||²=B b2+P p2.                       (5)
```

For each finite carrier the nondegenerate anchor eigenvalue has a local
analytic perturbation branch. With an eigenvector phase fixing its
first derivative perpendicular to Omega0,

```
Ehat(lambda)=-90V-lambda²(B b1+P p1)+O(lambda³),
Omega'(0)=-R0 Vpert Omega0.                           (6)
```

The O(lambda³) term in (6) has no supplied volume-uniform bound, and
(6) is not evaluated as an energy prediction at lambda=1. In particular
the first-order dressing norm squared in (5) grows proportionally to V.
Its bounded per-cell value is not a bound on the global state distance.
It also does not establish the convergence or failure of a linked-cluster
expansion. The actual primary expectation of H on Omega0 is -90V/a,
and its variance is `(15B+P/4410)/a²`; Omega0 is not an eigenvector of H.

The [same-source locality proof](../particle_primary_locality/PROOF.md)
gives a bond norm at most 30 and a plaquette norm at most 1/7 in
dimensionless units. Hence J=30B+P/7 bounds ||Vpert||. The elementary
min-max perturbation estimate gives

```
gap(Hhat(lambda)) >= 6-2|lambda|J.                    (7)
```

It guarantees a positive gap only for |lambda|<3/J, a range shrinking
as 1/V. At the primary endpoint its lower bound is negative for every
selected member. The local incident norm is at most `1272/7`, whose
ratio to the anchor physical gap is `212/7`. These large norm bounds
do not disprove stability; they do not certify a small perturbation of
the anchor at lambda=1. Electric packet denominators in the table reduce
the vacuum dressing amplitudes, but they do not control every excited
configuration or overlapping sequence of hops.

A useful strong-coupling continuation must supply a uniform estimate
on those interacting configurations and a remainder or gap-stability
bound reaching lambda=1, uniformly in the box and link cutoff. The
[native color-link calculation](../source_strong_coupling/PROOF.md)
has a different link carrier and action. The
[Hubbard spin limit](../strong_coupling_refinement/README.md) has a
different Fock interaction, coupling trajectory and time scaling.
Neither proves that missing estimate for (2).

## A uniform lattice gap would exclude finite-band primary survival

Suppose, as an additional hypothesis, the actual Hhat_n(1) has a full
physical ground gap at least g>0 for all admitted n. Ground degeneracy
is permitted. Then H_n has gap at least g/a_n. For any physical probe
A_n with norm at most C, its original stationary transfer measure is

```
mu_n(S)=sum_(i,j) p_i |(A_n)_ji|² 1_S(E_j-E_i).
```

Fix 0<delta<R<g/a_n. No initial ground state can contribute a transfer
in [delta,R]: a final ground state has zero transfer and any other final
state has energy at least g/a_n above it. All such records therefore
start in the excited subspace. The actual cooling bound (1) gives

```
mu_n([delta,R]) <= C² Tr(rho_n P_exc)
               <= C² a_n epsilon_n/g
                = C² 2^(-2n-4)/g.                   (8)
```

The strict upper-band inequality matters: a ground-to-gap transition
can lie at the endpoint R=g/a_n. The code returns no bound (8) when
this admission fails. It marks g as a supplied, unproved premise for
the primary endpoint. A generic transition between two high states can
have a fixed small difference, but its original weight is subject to (8).
A three-level passive Q² control executes this possibility; conditioning
on its rare successful records removes exactly the probability relevant
to survival. The control is not the primary action or preparation.
It uses `H=diag(0,64/a,64/a+1)` and a Hermitian operator exchanging
the two excited levels, with a complete ground-level rejection outcome.
For `U=64/a+1`, `beta=64/a`, take
`Q(E)=1+beta(U-E)/2`. This shares the positive-polynomial functional
form, without claiming the primary preparation plan. Its excited
probabilities are faithful and passive. With 0<a<=1, bounding the
unnormalized ground weight below by `2048²/a⁴` and the two excited
weights above by `33²/a²` and `1/a²` proves excess at most
`69761 a/4194304<a/16`. Its fixed positive transfer at frequency one
has vanishing original weight despite a conditional positive fraction
greater than 99/100.

Independently of a gap hypothesis, the exact identity

```
mu_n([delta,R])=muhat_n([a_n delta,a_n R])            (9)
```

uses the same matrix elements and original probabilities. Uniform
fixed-band survival therefore requires nonvanishing response weight
in a dimensionless window shrinking towards zero. The
[mass-readout transport theorem](../particle_primary_mass_readout/PROOF.md)
then places corresponding final excitation weight in
`[a_n delta,a_n(R+r)]`, losing at most `C² epsilon_n/r` in original
weight, for any fixed r>0. Low dimensionless levels without their
operator weights do not satisfy this requirement.

The fixed e0=105 and kappa=1/7 do not themselves prove an O(1) lattice
gap or correlation length. Hhat_n changes volume and word cutoff.
The local estimates above provide a candidate gapped-phase expansion
and quantify its missing uniform control. If such a gap is proved,
(8) rules out fixed-band survival for every uniformly bounded physical
probe on the selected cooling family. A different bounded probe,
fixed-time autocorrelation, or postselection cannot create the absent
ground-relative energy levels. A nontrivial finite-energy limit of this
bounded original-record construction would then require a declared
different trajectory approaching a closing
dimensionless spectral scale, or a different explicit unit policy.
Running couplings are allowed architecture choices with their own
action, state, reference and error analysis. They are not outputs of
the fixed primary family. A change of units alone supplies no continuum
construction. With a fixed declared source speed v matched to c, a
source time unit tau_n proportional to 1/a_n keeps an O(1) dimensionless
gap finite in GeV but also makes metres per source coordinate
`c tau_n/v` proportional to 1/a_n. The physical lattice spacing then
stays constant. This is a different physical scaling domain; spatial
cutoff removal would need its own proof. A spatial correlation length or Lorentz mass
identification additionally needs its state, spatial and field premises;
it is not inferred from (9) alone.

## Consequence for the mass calculation

The existing primary plaquette witness proves a tiny finite contrast at
T_n=t_n=1/W_n, with a band whose endpoints approach zero and infinity.
Its shrinking-time evidence supplies no nonvanishing weight in (9).
The source-derived diagnostic (5) selects a concrete calculation:
control the local perturbation remainder and physical excited branches
at the specified finite member, or prove uniform stability of the
fixed-coefficient phase before interpreting its continuum mass scale.
No endpoint gap, fixed-band collapse, surviving continuum field or
physical mass in GeV is asserted by this package. The finite comparison
route is available once its response, branch and finite errors are
established, even when a continuum survival theorem is absent.
