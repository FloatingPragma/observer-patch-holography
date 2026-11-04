# Ground and thermal references of the complete source

Use the same signed two-flavor cubic source and the five exact local
transformations of [the local spin theorem](../bipartite_spin_limit/LOCAL_LIMIT.md).
The finite domains are periodic cubic boxes with each side a multiple of
four and at least four. The signed hopping is periodic under translation
by two sites. The degree, interaction norms and support counting in that
theorem apply to these tori: distance spheres contain no more sites than
their infinite-cubic counterparts. No term is removed at a charge boundary.

Write

```
G=nu N+T,  nu=U/2,  N=sum_x Q_x,  Q_x=(n_x-1)^2,
B=144c,   J=c^2/U,  nu>=2^20 B,
Y G Y*=nu N+D+Voff,
||D||_(1/2)<=2B,    ||Voff||_(1/2)<=r=B(4096B/nu)^5,
P D P=Hspin+delta,  ||delta||_(1/2)<=eta=2^25 B^3/nu^2,
Hspin=J sum_edges(S_x.S_y-1/4).
```

Every local term of D commutes with the sum of Q_x on its support;
`P Voff P=0`. The local dressing satisfies
`||Y* Q_x Y-Q_x||<=8B/nu`. These are exact local series bounds,
independent of the number V of lattice sites. The input states below
are actual states of G, with all `4^V` occupation states present.

## Charge density from source energy

For a state invariant under translations by two sites, put
`q=Tr(rho N)/V` and `e=Tr(rho G)/V`. The eight translation sublattices
have equal size, so `Tr(rho Q_x)<=8q` for every x. If a Hermitian operator
A obeys `P_Z A P_Z=0`, its diagonal Q block and its two off-diagonal
blocks give

```
|Tr(rho A)| <= ||A|| [p+2sqrt(p(1-p))] <=3||A||sqrt(p),
p=Tr(rho(1-P_Z)),  P_Z=product_(x in Z)(1-Q_x).
```

This follows from Cauchy-Schwarz for the positive state, without a
classical-mixture assumption. Each hopping bond has norm at most c,
zero P compression and defect probability at most 16q. There are three
bonds per site. Consequently

```
|Tr(rho T)|/V<=36c sqrt(q),
nu q<=e+36c sqrt(q),
sqrt(q)<=36c/nu+sqrt(max(e,0)/nu).                       (1)
```

For a ground state, average its density matrix over the even translations.
The average is an exact ground state. A vector in P has source energy
zero, so its actual energy density obeys `e<=0`. Equation (1) therefore
applies with the second term zero. This does not place the finite-c ground
state exactly in P or assume a volume-independent spectral gap.

For the full Gibbs state `rho_beta=exp(-beta G)/Z_G`, with fixed `beta>0`,
faithfulness and translation symmetry follow directly. Its free energy
is no greater than that of a pure P trial state, which is zero. Since
`S(rho_beta)<=V log4`, its energy density satisfies

```
e<=log4/beta<2/beta.                                    (2)
```

Let m denote the right side of (1), using zero for the ground state
and `2/beta` as the thermal energy upper bound. In the rotated state
`rho'=Y rho Y*`, the projection identity and the triangle inequality
in the state Hilbert norm give

```
sqrt(Tr(rho' Q_x))
 <=sqrt(Tr(rho Q_x))+||Y*Q_xY-Q_x||
 <=sqrt(8)m+8B/nu=:d,
p=d^2.                                                 (3)
```

At fixed J and `U=c^2/J`, both m and d are `O(c^-1)`.
The executable majorants replace `sqrt(8)` by 3 and enclose the other
square root above by an exact rational upper endpoint. They require
`p<=1/4` for the thermal entropy bound below.

## A unital spin algebra and the charged complement

At each site identify the complete four-state carrier as charge times
spin, in the order

```
|0,0>=|up>, |0,1>=|down>,
|1,0>=|empty>, |1,1>=|up down>.
```

Then `Q=|1><1| tensor I` and the unital spin embedding is
`j(A)=I_charge tensor A`. It preserves products and adjoints and commutes
with every Q_x. Its operators preserve fermion parity. The extension on
the charge block can change particle number by two; no particle-number
superselection restriction is imposed on this comparison algebra. On P
it is precisely the original spin algebra. For every source state,
`omega_c(A)=Tr(rho' j(A))` is a positive normalized spin state, obtained
by a partial trace over the charge factors. There is no conditioning on P.

For every local D_Z, set `h_Z=P_Z D_Z P_Z`, interpreted on its spin
carrier, and `C_Z=D_Z-j(h_Z)`. Commutation with the local penalty makes
P_Z reducing, so

```
C_Z=(1-P_Z) C_Z (1-P_Z),   ||C_Z||<=2||D_Z||,
D=j(Hspin+delta)+C,        ||C||_(1/2)<=4B.              (4)
```

For a local spin A on s sites, j(A) commutes with every P_Z. Terms
disjoint from its support commute with A. The union bound gives
`Tr(rho'(1-P_Z))<=|Z|p`. Since `n exp(-n/2)<=1` for positive integers,
summing the remaining terms in (4) yields

```
|Tr(rho' j(A)*[C,j(A)])|<=8s B p ||A||^2.               (5)
```

The analogous bounds for delta and Voff are
`2s eta ||A||^2` and `2s r ||A||^2`. For the mean energy, the same support
counting gives `|Tr(rho' C)|/V<=4Bp`. The factor p, rather than sqrt(p),
in these estimates uses both reduction by P_Z and the unital spin algebra.

## Actual source grounds select interacting spin grounds

For any exact source ground density matrix, positivity above its ground
energy gives
`Tr(rho' j(A)*[YGY*,j(A)])>=0`. The penalty commutator vanishes exactly.
Equations (4)--(5) imply

```
Re omega_c(A*[Hspin,A])
 >=-s ||A||^2 [8Bp+2eta+2r].                            (6)
```

The imaginary part has the same absolute error bound; the source
quadratic form itself is real. The right side tends to zero as `c->infinity`
at fixed J, independently of volume. Thus, for any sequence of these
tori with boundary receding from every fixed local support and any sequence
`c->infinity`, weak local accumulation states satisfy

```
omega(A*[Hspin,A])>=0
```

for every local A. This is the algebraic ground-state condition for the
same interacting spin dynamics. It follows by taking the local limit of
(6), since `[Hspin,A]` involves only finitely many neighboring bonds.
Taking `A=I+z B` first shows invariance; the ground quadratic form then
gives a nonnegative implementing Hamiltonian in its GNS representation.
The argument identifies the ground-state class, not a unique ordered
phase, sharp excitation shell or a relativistic vacuum.

## Thermal free energy and finite-volume reference convergence

Let `gamma_beta=exp(-beta Hspin)/Z_spin`, and take rho to be the actual
source Gibbs state. The charge marginal of rho' has entropy at most
`sum_x h(Tr(rho'Q_x))<=V h(p)` when `p<=1/4`, where h is binary entropy
with natural logarithms. The first inequality is subadditivity followed
by dephasing each charge qubit. Araki-Lieb and subadditivity then imply

```
|S(rho')-S(omega_c)|<=V h(p).                           (7)
```

These finite-dimensional entropy inequalities require no tensor product
assumption about the state. The tensor product is the complete onsite
charge/spin dictionary. Dropping the positive penalty from the energy
and using (4) gives

```
F_G(rho)>=F_spin(omega_c)-V[eta+r+4Bp+h(p)/beta].         (8)
```

The source trial state `Y*(|0...0><0...0| tensor gamma_beta)Y`
has the same entropy as gamma_beta. Its penalty and Voff expectations
are zero. Its free energy is at most `F_spin(gamma_beta)+V eta`.
The Gibbs variational identity
`F_H(sigma)-F_H(gamma_H)=D(sigma||gamma_H)/beta>=0`
therefore proves, with `epsilon=2eta+r+4Bp+h(p)/beta`,

```
0<=F_spin(omega_c)-F_spin(gamma_beta)<=V epsilon,
D(omega_c||gamma_beta)<=beta V epsilon.                 (9)
```

For a rational executable bound use `h(p)<=3sqrt(p)<=3d`. For example,
`-p log p<=sqrt(p)` follows by maximizing `-sqrt(p)log p`, and
`-(1-p)log(1-p)<=p<=sqrt(p)`. Thus epsilon is `O(c^-1)` uniformly in
volume. The source and spin free-energy densities differ by at most
`eta+r+4Bp+h(p)/beta` below and eta above. This gives pressure convergence
without choosing a thermodynamic spin phase.

At each fixed finite volume, the quantum Pinsker bound gives

```
||omega_c-gamma_beta||_1<=sqrt(2 beta V epsilon).         (10)
```

The coefficient and logarithm convention are those of
[Carlen and Lieb, equation (2.6)](https://arxiv.org/html/1402.3840v2).
Equation (10) is explicitly extensive inside the square root. Relative
entropy per volume alone does not prove uniform local trace closeness for
arbitrary simultaneous volume and coupling limits, or phase uniqueness.
An infinite-volume Kubo-Martin-Schwinger identification is not inferred
from that bound. Source Gibbs stationarity does, through (4)--(5) applied
to a single commutator, make any weak spin limit stationary under Hspin.

## Faithfulness, clocks and admitted state selection

The exact source Gibbs state has the explicit strict lower bound

```
lambda_min(rho_beta)>=4^(-V) exp[-beta V(nu+6c)]>0.       (11)
```

Indeed `-3cV<=G<=nu V+3cV` as quadratic forms. Its modular automorphism
on the complete finite source algebra is exactly

```
sigma_s(A)=rho_beta^(is) A rho_beta^(-is)
          =exp(-i beta s G) A exp(i beta s G)=tau_G^(-beta s)(A). (12)
```

It is nontrivial because G is not scalar; the same source time appears
in its Hamiltonian and Gibbs reference. Finite spin gamma_beta has the
corresponding identity with Hspin, is strictly positive, and has
noncommuting neighboring bond interactions. These identities do not
identify a proper gravitational time or an empirical temperature unit.

Finite G and its Gibbs matrix are computable from the admitted simulator
parameters. This result uses Gibbs and ground selection as specified input
rules. It does not derive those selection rules from the unread unital
repair process or assert that a finite-c Gibbs state is stationary under
continuing repair. The continuing-repair theorem can consume its bounded
energy and symmetry premises, but its finite-c semigroup has its own
reference law. No exponentially large preparation circuit is claimed
executed by the finite matrix controls in this package.
