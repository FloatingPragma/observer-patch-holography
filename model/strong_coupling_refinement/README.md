# Controlled spin exchange and the growing-contact free boundary

The centered source has a controlled interacting **spin** limit on every
fixed finite graph. The result includes an isometry into its actual Fock
space, a uniform-in-input real-time error, and a charged-sector gap. The
same source's continuum charged fields remain free for the larger coupling
class `U(c)=o((log c)^(3/2))`, `c=v/a -> infinity`. These are different
limits. Neither result constructs a massive interacting relativistic
continuum sector or identifies a measured particle.

## Source and exact finite carrier

In the [source's onsite spin frame](../field_symmetry/README.md), the
centered Hamiltonian on a finite simple graph is

```
K = T + U D - (U/2) N,
T = sum_{edges xy,spin} t_xy (c_xspin^* c_yspin+c_yspin^* c_xspin),
D = sum_x n_xup n_xdown,                  t_xy = +/- c/2.
```

For the open cubic source the sign is
`(-1)^(sum_{j>axis} x_j)`, verified by the exact onsite frame. Write `V`
for the number of sites, `E` for its edge set, and `tau=c |E|`. Each
one-spin hopping term has norm `|t_xy|`; thus `||T||<=tau` on the entire
Fock space. This deliberately extensive bound requires no ground-state
hypothesis. The statements below hold for arbitrary real edge amplitudes
with `tau=2 sum |t_xy|` as well; the calculators specialize to the source.

Restrict to `N=V`, add the scalar `UV/2`, and set `H=T+UD`.
Let `P` project onto exactly one particle at every site, and `Q=I-P`.
Then `rank P=2^V`, `PTP=0`, `QUDQ>=UQ`, and

```
(UD) T P = U T P.                                      (1)
```

Indeed one allowed hop out of a singly occupied configuration creates
exactly one empty site and one double occupancy. A second hop returns
to P only by repairing this same pair. Applying the two CAR hops with
their actual fermionic signs gives the operator identity

```
Hspin = -(1/U) P T^2 P
      = sum_{xy} (4 t_xy^2/U) (S_x.S_y-1/4).             (2)
```

Here `S` is the hidden onsite spin in this frame; it is not an assertion
about physical Lorentz spin. The hopping signs cancel in (2), but need
not cancel at higher orders. `algebra.py` checks (1)--(2) with integer
matrices in the complete half-filled sectors of the dimer, three-site
chain, square, and eight-site cube (dimensions 6,20,70,12870). It does
not project the hopping before multiplying it.

## Low spectrum with a finite error

Assume `U>tau`, put `B=QTP`, `C=QTQ`, and `A=QUDQ+C`.
Then `A>=U-tau>0`. Min--max on P and the norm perturbation of UD give
exactly `m=2^V` eigenvalues of H in `(-infinity,0]`; all others are at
least `U-tau`. For `z<=0` their Feshbach operator is

```
F(z)=-B^*(A-z)^(-1)B.
```

Each low eigenvector has a nonzero P component satisfying `F(z)p=zp`.
Conversely every such solution reconstructs an eigenvector with Q
component `-(A-z)^(-1)Bp`. In particular
`-tau^2/(U-tau)<=z<=0`. Since `(QUDQ)^(-1)B=B/U`, the resolvent identity
implies, throughout this interval,

```
||F(z)-Hspin||
 <= tau^2(tau+|z|)/(U(U-tau))
 <= tau^3/(U-tau)^2.                                   (3)
```

Ordered eigenvalues of F(z) are continuous nonincreasing functions of z,
since `F'(z)=-B^*(A-z)^(-2)B<=0`. Their differences from z are strictly
decreasing, positive for sufficiently negative z, and nonpositive at
zero. Their m roots, counted with multiplicity, are precisely the
ordered low spectrum. Applying the variational eigenvalue bound to (3)
therefore pairs every low eigenvalue with the corresponding eigenvalue
of (2) within `tau^3/(U-tau)^2`. No spectral gap between the individual
spin levels is required.

## An isometry and the complete real-time error

Assume the stronger explicit condition `U>=4tau`; put `r=tau/U<=1/4`.
Let `Ptilde` be the low spectral projection, and choose any isometry
`W0:P -> ran(Ptilde)`, with low generator `L0=W0^*H W0<=0`.
The Q part obeys a Sylvester equation. Its norm-convergent solution is

```
QW0=-int_0^infinity exp(-sA) B PW0 exp(sL0) ds,
||QW0|| <= rho := r/(1-r).
```

Thus the low subspace is the graph of an operator `X:P -> Q` with
`||X||<=rho/sqrt(1-rho^2)=r/sqrt(1-2r)=:x<=sqrt(2)r`.
The canonical graph isometry and its positive P compression are

```
W=(P+X) a0,          a0=(I+X^*X)^(-1/2),
||W-P|| <= rho+x^2/2 <= 2r.                            (4)
```

Invariance of the graph gives
`QUDQ X+B=X B^*X-CX`. Using (1),
`||X+B/U||<=(tau/U)(x+x^2)`. The exact, self-adjoint low generator is
`Hlow=W^*HW=a0^(-1) B^*X a0`. Bounds
`||a0^(-1)||<=9/8`,
`||a0-I||,||a0^(-1)-I||<=x^2/2` give

```
||Hlow-Hspin||
 <= (tau^2/U)[(9/8)x+(17/8)x^2]
 <= 3 tau^3/U^2.                                      (5)
```

For the last inequality use `x<=3r/2`, `x^2<=2r^2`, and `r<=1/4`;
the bracket is at most `(11/4)r`. The intertwining
`exp(-itH)W=W exp(-itHlow)`, two applications of (4), and Duhamel yield

```
|| exp(-itH)P-P exp(-itHspin) ||
 <= min(2, 4tau/U+3|t|tau^3/U^2).                       (6)
```

This controls every input spin state and its full Fock-space output,
including leakage. There is no successful-outcome conditioning or
unmeasured-branch deletion. Restoring K supplies the known common
phase `exp(itUV/2)` within N=V.

On the entire Fock space,

```
K+UV/2 = T+(U/2) sum_x(n_x-1)^2.
```

Every sector `N!=V` has at least one occupancy defect and hence energy
at least `U/2-tau` in this normalization, whereas a P trial vector has
energy zero. Therefore `U>2tau` puts every global ground vector in
`N=V`, with charged-sector separation at least `U/2-tau`. This is a
finite-graph bound; it gives no quantitative neutral spin gap.

## A genuinely interacting limit and its volume boundary

At fixed finite graph, choose `U=c^2/J`, with fixed `J>0`, and let c
increase. Then (2) is the interacting Heisenberg generator
`J sum(S_x.S_y-1/4)`, and (6) becomes

```
4 |E| J/c+3 |t| |E|^3 J^2/c -> 0.
```

Adjacent bond terms do not commute. More directly, a dimer evolves the
product `|up,down>` into a maximally entangled spin state at
`t=pi/(2J)`: the singlet acquires phase i and the triplet phase 1.
At `c=128,U=16384,J=1`, (6) is strictly less than `5/64` at this time.
The actual unconditioned Fock output consequently has fidelity greater
than `1/2` with this Bell state, which certifies site entanglement even
allowing empty/doubly occupied outputs. The exact integer controls and
independent full-Fock exponential tests retain these branches.

For the eight-site cube, `c=65536,U=c^2,J=1` gives the unit-time bound
`327/4096<0.08`. These are controlled model examples, not measured mass
fits. The original cube parameters `c=1/8,U=1` have `tau=3/2` and do
not satisfy the new expansion premise. Their earlier exact results are
neither replaced nor used to certify (6) outside its domain.

Refining a fixed physical open cube is different. For L sites per side,
`|E|=3L^2(L-1)`, `c=vL`, and `tau=O(L^4)`. If `U` has power L^p, the
two errors in (6) scale as `L^(4-p)` and `L^(12-2p)`. These global
bounds prove vanishing error for `p>6`; the exchange then scales as
`J=c^2/U=O(L^(2-p))` and tends to zero. The fixed-J choice `p=2`
fails the hypothesis. Failure of this extensive estimate is not a
no-go theorem for a local expansion or another continuum trajectory.
It does mean the fixed-graph interacting limit cannot be presented as
a controlled interacting 3+1-dimensional field theory.

## A quantitative growing-U extension of the free continuum theorem

The [centered-ground construction](../centered_source_ground/README.md)
exists at every finite positive c,U. Its pull-through estimate is
`Delta<=11(U/c)^2/8`. The [free-limit proof](../continuum_free_limit/README.md)
uses the canonical composite `d=(2n_opposite-1)c`, not the uncentered
contact field. Its microtime, phase-space, and Gaussian estimates hold
for each finite U; their derivation does not require a common U bound.
For either charge branch, uniformly in every unit lattice smearing,

```
||1_[0,R](L_a) d(f)Omega_a||^2 <= M_a(R),
M_a(R) <= min(1, exp(R^2/(2c^2 delta^2))
 [1536 delta^6+3(U/c)/delta
  +624 sqrt(11/8)(U/c) exp(81/(2delta^2))]).             (7)
```

L_a is the positive generator of the **centered** ground. Negative
response-frequency labels on the removal branch are not negative
physical excitation energies. For kinetic windows bounded by P, the
same source equation and output-energy split give, for both adjoints,

```
||(C_a(f,t)-C_a(exp(it h_a)f,0))Omega_a||
 <= (U/2)[|t|sqrt(M_a(R))+(2+|t|P)/R] ||f||.           (8)
```

On a positive kinetic band bounded below by b>0, the corresponding
initial annihilator bound is
`(U/2)[sqrt(M_a(R))/b+1/(R+b)]`; the negative-band creator has the same
bound. These follow by integrating the exact Sylvester equation on
the conjugate one-particle input space. Their dependence on U is
explicit; mere convergence of M to zero would not suffice.

Let `ell=log c`, choose `delta=9/sqrt(ell)` and `R=sqrt(c)` for c>1.
The prefactor in (7) is `exp(ell/(162c))`. Its three terms are

```
1536*9^6/ell^3,
(U/c)*sqrt(ell)/3,
624 sqrt(11/8) U/sqrt(c).
```

Consequently

```
U=o((log c)^(3/2))  ==>  U^2 M_a(sqrt(c))->0,
                         U/sqrt(c)->0.               (9)
```

This includes unbounded U. Equations (8)--(9) kill both vacuum
Duhamel differences and both initial wrong-band components. They
establish the same pure free filled-sea initial CAR state and free
vacuum field evolution as in the bounded-U theorem.

### Why unbounded U does not silently reuse a bounded-U compactness premise

The earlier field construction used a time-Lipschitz bound containing U;
it cannot be directly reused uniformly here. Instead, for all real
times and smooth compactly supported physical smearings in all eight
valleys, consider the formal bounded generators whose evaluations are
the corresponding cell-averaged lattice fields. Their norms are bounded
by the smearing norms, independent of a and U. The universal unital
C*-algebra generated by these bounded symbols exists (equivalently use
the free product of their universal norm-bounded single-generator
algebras). Each lattice ground defines a state on it. Weak-star
compactness supplies a convergent subnet of any refinement subnet;
no temporal continuity or countability of this generator set is assumed.
All finite-word correlators thus have cluster states simultaneously.

Equal-time CAR pass to each cluster state by the source's cell-average
and valley-orthogonality estimates. Equations (8)--(9), first on bounded
kinetic windows away from zero and then by norm density, give both
initial free-band annihilation identities and both free vacuum time
identities there. In comparing input groups, the band-isometry and
central-difference estimates involve c=v/a and smooth fixed physical
packets, not U. Compactly supported packets follow by L2 approximation.

The source spatial locality estimate is also uniform in onsite U: the
onsite evolution can be removed exactly in the interaction picture;
it changes no support or bond norm. The full-CAR bound recorded in the
[thermodynamic-ground proof](../thermodynamic_ground/README.md) and its
[cell-averaged physical cone](../continuum_field_net/README.md) therefore
pass to these cluster states. In particular a field at (t,f), with f
compactly supported, graded-commutes with every initial field outside a
bounded region containing the `96v|t|` enlargement of supp f.

For completeness, the final promotion from vacuum vectors to the whole
hierarchy does not assume a Reeh--Schlieder axiom. In each free Weyl band,
projections of smooth functions supported outside any bounded region
are dense: a vector orthogonal to them would have compact spatial
support and entire Fourier transform F, with
`M_b(p)F(p)=+/-|p|F(p)`. Fix real nonzero transverse momentum and continue
in the remaining complex variable around one simple branch point of
`sqrt(z^2+r^2)`. F and M_b are single-valued whereas the root changes
sign, forcing F=0. This applies to all eight independent valleys.
Wick recursion then makes the outside initial CAR algebra cyclic on
the free filled-sea subspace H0, by density of its particle/hole wedge
vectors.

The free Weyl equation has finite propagation speed v, so the free
evolution of a compactly supported smearing also stays within a bounded
spatial enlargement. The difference between a time field and this free
initial expression is therefore bounded and odd, annihilates the vacuum
with its adjoint by (8), and graded-commutes with the same outside algebra. It therefore vanishes with
its adjoint on the dense outside-polynomial set in H0. Thus H0 reduces
every time field. The full cluster GNS representation is generated by
those time fields from the vacuum, so H0 is the whole representation.
Every cluster hierarchy is the identical free hierarchy. Compactness
then implies convergence of each original finite-word correlator,
without extracting a privileged subnet. Strong continuity, positive
energy, and free Weyl propagation follow from this identified limit.
This argument supplies the missing compactness step for unbounded U;
it does not assume the old uniform time-Lipschitz estimate survives.

The result concerns the canonical charged CAR fields: all eight Weyl
valleys have energies `v|p|` and zero rest mass. It neither sets the
composite equal-time CAR field to zero nor retains a nonzero continuum
contact interaction. No physical species assignment or gravity
conclusion follows from the free charged hierarchy alone.

### An explicit growing coupling and resolution schedule

The asymptotic assertion has a finite rational certificate. For integer
`k>=1`, set

```
c_k=2^(256k^2),  U_k=k^2,  delta_k=1/k,  R_k=k^3.
```

Use `exp(x)<=2^ceil(3x/2)` for x>=0 and
`sqrt(11/8)<=3/2`. The outside exponential in (7) is at most 2, and the
inside one is at most `2^(61k^2)`. Therefore

```
U_k^2 M_k(R_k)
 <= 3072/k^2+6k^7 2^(-256k^2)+1872k^6 2^(-195k^2)
 <= 3073/k^2,
U_k sqrt(M_k(R_k)) < 56/k,        U_k/R_k=1/k.          (10)
```

To check the last inequality exactly, use `k<=2^k`: the last two terms,
after multiplying by k^2, are bounded by
`2^(3+9k-256k^2)` and `2^(11+8k-195k^2)`, each less than 1/2.
For `||f||<=1`, `||h_a f||<=1`, and `|t|<=1`, (8) is at most
`min(2,30/k)`. These extreme scale examples demonstrate a constructive
bound, not a practical claim about accessible lattice sizes. The
calculator keeps log2(c) as an integer and proves the exponent
inequalities; it does not create a false zero by floating underflow.

## Evidence and remaining obligations

Run `python -W error -m pytest -q model/strong_coupling_refinement`, then
`python -W error -m model.strong_coupling_refinement --write` and
`--check`. Integer identities and rational constants are the proof
controls. Full small-carrier diagonalization and time evolution are
independent diagnostic tests of the bounds, not substitutes for the
proof for arbitrary finite graphs.

The receipt binds current proof/code sources and records whether
inherited historical receipts match those sources. A dependency whose
receipt is stale is reported as such; the new replay does not claim
that prior expensive experiments were rerun. The growing-U theorem uses
the stated analytic premises and their source proofs. Final integrated
custody requires the owning project's dependency refresh.

An interacting relativistic charged continuum requires a trajectory
outside the proved free class, or a changed source/field construction,
and a new volume-uniform theorem. The finite spin limit supplies actual
interaction and charged separation but cannot supply that theorem by
changing its name. It also supplies no Einstein dynamics or measured
particle dictionary.
