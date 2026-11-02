# Geometric modular flow and stress of the proved Weyl vacuum

The centered bounded-contact continuum established in
[the free-limit theorem](../continuum_free_limit/README.md) has a symmetric,
conserved Lorentz stress tensor and a geometrically normalized modular group
for every round double cone. They belong to the **same eight-valley field,
vacuum and clock**. This is a matter-sector attachment to the existing Einstein
theorem, not a derivation of an interacting theory or a dynamical metric.

The new mathematical content includes the explicit source-frame map, the
identification of the fermionic modular generator with the appropriate Wick
stress charge, its required spinor correction, and a controlled small-ball
stress expansion. The modular-group existence and its Cauchy-data generator
are imported from the applicable free-fermion theorem, not inferred from a
finite numerical fit.

## Source, state and units

The input theorem is the complete hierarchy result at scientific checkpoint
`6a66874bb`, recorded as FA82. At every spacing the reference is the
particle-hole-invariant ground of `K_a=H_a-U(a)N/2`, with `c(a)=v/a`, fixed
`v>0`, and `0<=U(a)<=Umax<infinity`. The field limit keeps all eight
`b in {0,1}^3` and identifies their joint vacuum hierarchy with the free
Weyl hierarchy. In particular its one-particle symbol, in `x0=v*t`, is

```
h_b(p)=-sum_j (-1)^b_j sigma_j p_j.
```

Write `chi_b=-product_j(-1)^b_j` and
`R_b=diag(-chi_b*(-1)^b_j)`. This is a proper rotation. The constant lift
`U_b` is the identity when `R_b=I`, and `i sigma_j` when its positive axis
is `j`. Exactly,

```
U_b^* h_b(p) U_b=chi_b sigma.p,     U_b^*U_b=I, det U_b=1.
```

These transformations act on spin components at the same point: they
preserve the CAR, localization, reference covariance and time evolution.
Pair the four negative with the four positive chiralities to obtain four
independent massless Dirac fields. This is a unitary reorganization of all
sixteen complex Cauchy components. It neither deletes valleys nor identifies
four measured species. The vacuum is the tensor/graded-product vacuum fixed
by the FA82 hierarchy; it is not a newly chosen reference.

Use `alpha_j=diag(-sigma_j,sigma_j)`, the off-diagonal identity `beta`, and
`gamma^0=beta, gamma^j=beta alpha_j`. Then
`{gamma^mu,gamma^nu}=2 eta^{mu nu}` for `eta=diag(1,-1,-1,-1)`.
The geometry stack uses `g=-eta`. These conventions have the same Lorentz
group, causal cone and covariant two-index stress: lowering both indices
changes no sign, and the homogeneous trace and Ward equations remain zero.
We use dimensionless source action units, `hbar=1`. The model-time energy is
`H_model=v P0`. No conversion to seconds, joules, or a measured Newton
constant is introduced.

## Symmetric stress from the same propagation law

For each of the four Dirac fields let `bar(Psi)=Psi^* gamma^0`. On the
smooth finite-particle domain define the vacuum-normal-ordered distribution

```
T^{mu nu}=(i/4) sum_A :bar(Psi_A)
  [gamma^mu <->partial^nu+gamma^nu <->partial^mu] Psi_A: .       (1)
```

The vacuum subtraction is taken in the FA82 reference. Products in (1) are
Wick distributions, not unrenormalized products at a point. Smearing them
and taking matrix elements between finite-particle Schwartz wave packets
gives well-defined quadratic forms. The field equations give

```
T^{mu nu}=T^{nu mu},       partial_mu T^{mu nu}=0,
eta_{mu nu} T^{mu nu}=0.                                    (2)
```

Here is an algebraic proof which also fixes the normalization. For cross
waves of null momenta `p,q`, the coefficient between `bar(u(p))` and `u(q)`
is

```
t^{mu nu}(p,q)=[gamma^mu(p+q)^nu+gamma^nu(p+q)^mu]/4.
```

Contraction with `(p-q)_mu` gives
`[(slash(p)-slash(q))(p+q)^nu+gamma^nu(p^2-q^2)]/4`.
The left and right Dirac equations kill this expression. Its trace is
`(slash(p)+slash(q))/2` and vanishes for the same reason. Fourier expansion,
normal ordering, and distributional smearing extend these identities to the
operator distributions. The exact calculator proves the corresponding
matrix polynomial identities before imposing `p^2=q^2=0`. It therefore
checks more than a list of numerical on-shell samples.

At time zero,
`T00=(i/2):Psi^*<->partial_0 Psi:`. Integrating gives the normal-ordered
second quantization of `h=alpha.p`, with the filled-sea hole convention.
Thus `integral T00=P0`, and the spatial charges generate the same
translations as the source limit. This fixes the coefficient `1/4` in (1);
an arbitrary multiple would not generate those translations. The Ward
identity uses the field equation alone, independently of any Einstein
equation. Flat-space tracelessness is asserted here; this does not discard
the curved-background conformal anomaly.

## The applicable imported modular theorem

Theorem 3.4 and Theorem 4.1 of F. La Piana and G. Morsella,
*The Fermionic Massless Modular Hamiltonian*, Commun. Math. Phys. 406,
81 (2025), [arXiv:2405.07888v2](https://arxiv.org/html/2405.07888v2),
[DOI](https://doi.org/10.1007/s00220-025-05253-1), compute the vacuum
double-cone modular group for the Weyl and massless Dirac field algebras.
The domain is the ordinary free-field vacuum representation, not a finite
CAR truncation. The theorem uses the free-field conformal representation
and its vacuum standardness. Both are applicable after the complete
unitary identification of the preceding section. The same theorem applies
to each chirality and a finite direct sum, hence to all eight valleys.

With `D_R={|x0|+|x|<R}`, set

```
w_R(x)=(R^2-|x|^2)/(2R),
q_R=pi {w_R,h}=-i*pi/R [(R^2-|x|^2) alpha.grad-alpha.x].     (3)
```

Here `q_R` is the signed-CAR Cauchy-data expression for the complete modular charge,
with sign convention `Q_R=-log Delta_R`; the modular automorphism is
`sigma_s(A)=exp(-is Q_R) A exp(is Q_R)`. Formula (3) on all of space has
weight negative outside the ball. It is not a positive Hamiltonian and is
not the logarithm of a trace-class reduced density matrix.

The cited paper's wave Hilbert space uses the complex structure
`iota=i*h/|h|`, not ordinary multiplication by `i` on the full signed CAR
Cauchy space. With its `log Delta=-iota*K_D`, equation (3) becomes `q_R`
on the particle block and `-conjugate(q_R)` on the hole block. This is
exactly the filled-sea normal-ordering convention below; it is not an
identification of the signed generator with a positive-energy one-particle
generator using the wrong complex structure.

The term `-i alpha.grad(w_R)/2` in `q_R/(2*pi)` is necessary. Integration
by parts shows that `w_R h` alone is not symmetric. The exact control
deletes this term and detects its nonzero adjoint defect. The formula is
the radius-rescaled version of the cited unit-double-cone theorem.

## Identifying the modular charge with Wick stress

An identification of a conformal generator with an unrestricted stress
integral must be proved; the cited fermion paper itself distinguishes this
question from computing the modular group. The following argument supplies
the identification on its natural common form domain.

In momentum space write `r=|p|`, `P_+=(I+h/r)/2`. Fourier transformation of
(3) gives

```
q_R=pi [R h+(h Delta_p+sum_j alpha_j partial_{p_j})/R].       (4)
```

On smooth momentum functions supported away from zero, direct Clifford
algebra gives

```
2h partial_j P_+ + [alpha_j,P_+]=0,
h Delta_p P_+ +sum_j alpha_j partial_j P_+=0.                (5)
```

Indeed `partial_j P_+=(alpha_j/r-h p_j/r^3)/2` and
`Delta_p P_+=-h/r^3`. Equations (5) prove `[q_R,P_+]=0`, and also
`[q_R,P_-]=0`. No filled-sea particle/hole pair-creation term occurs in this
complete conformal charge. All eight signed Clifford blocks satisfy the
identities separately. The formulas extend distributionally through
`p=0`: the displayed singular derivatives are locally integrable and
their boundary terms vanish. Removing a ball of radius epsilon from a
Schwartz test function gives an error of order `sqrt(epsilon)` in the
graph norm of the differential expression. Thus the origin supplies no
extra point mode or boundary term.

Define the normal-ordered charge by second quantizing these particle and
hole blocks with the usual opposite hole sign. Its vacuum value is zero.
On finite-particle smooth vectors, integration by parts in (1) yields

```
Q_R = 2*pi integral_{R^3} w_R(x) T00(0,x) d^3x              (6)
```

as a quadratic-form identity on the common domain. One can obtain this
domain first with momentum test functions away from zero and then by the
graph approximation just described. The imported modular group supplies
the self-adjoint generator and its extension; no claim of essential
self-adjointness of a sharply cut off ball-only operator is needed.
The charge (6) and that group act by the same transformation on every CAR
field and both fix the vacuum. Vacuum cyclicity then makes their unitary
implementers identical. This fixes the scalar ambiguity as well as the
noncentral action.

For fields supported strictly inside `B_R`, the contribution outside the
ball commutes with them. Consequently the local modular derivation and its
finite vacuum-subtracted energy pairings are described by

```
K_B = 2*pi integral_{B_R} w_R T00
    = (2*pi/v) integral_{B_R} w_R energy_density_model.      (7)
```

`K_B` in (7) denotes this local derivation/form and pairing. We do not
assert that either separated term `K_B` or `K_outside` is a globally
defined self-adjoint operator, although their full difference (6) is the
modular generator. In particular we never use a fictitious finite vacuum
entropy for the continuum type-III algebra. The coefficient `2*pi` is
common to all eight valleys; summing their stress does not multiply it
by eight a second time. This agrees with the geometric ball normalization
in [Casini--Huerta--Myers, section 2.1](https://arxiv.org/html/1102.0440v2).

## Geometry of the same modular action

The future-oriented inverse modular flow has vector field `2*pi*xi`, with

```
xi^0=(R^2-(x0)^2-|x|^2)/(2R),      xi^j=-x0*x^j/R.
```

It satisfies the conformal Killing equation
`partial_mu xi_nu+partial_nu xi_mu=(partial.xi) eta_mu nu/2`.
In the radial null coordinates `u=x0+|x|` and `v=x0-|x|`, each coordinate
evolves by

```
z(s)=R[(R+z)-exp(-2*pi*s)(R-z)]/
       [(R+z)+exp(-2*pi*s)(R-z)].                           (8)
```

Equation (8) preserves the double cone. Equations (1)--(2) make the
conformal current `T^{mu nu} xi_nu` conserved. Translations and Lorentz
boosts carry this result to every round double cone with any unit timelike
normal. Thus the directional charges use one tensor, not independently
fitted densities. All future timelike and null directions exist in this
sector's mathematical Minkowski domain. Identification with the native
event/count geometry and availability of physical observer readouts are
separate questions.

## A finite vacuum-relative first law

For a normal state rho with finite Araki relative entropy and finite
pairing (7), use the vacuum-subtracted entropy difference

```
Delta S_B(rho)=Delta<K_B>_rho-S(rho_B || omega_B).            (9)
```

This is a relative-entropy definition, not an assertion that two divergent
von Neumann entropies may be subtracted without a prescription. No
particular lattice entropy subtraction is proved to converge to (9).

There is an explicit nonempty class of two-sided perturbations for which
the first law follows with a quadratic remainder. Choose a self-adjoint
involution `B'` in the commutant of the ball algebra with `omega(B')=0`.
Its projections `P_+,P_-` have vacuum weights `1/2`. The normal states
`rho_±(A)=2<Omega,P_± A Omega>` satisfy
`omega=(rho_++rho_-)/2`. Let
`rho_epsilon=((1+epsilon)rho_++(1-epsilon)rho_-)/2`.
Monotonicity of Araki relative entropy under the normal classical-to-quantum
preparation channel gives

```
0 <= S(rho_epsilon||omega)
  <= [(1+epsilon)log(1+epsilon)+(1-epsilon)log(1-epsilon)]/2
  <= epsilon^2/[2(1-r0^2)],       |epsilon|<=r0<1.          (10)
```

The last inequality follows from the exact second derivative
`1/(1-epsilon^2)` and its zero value and slope at zero.
The data-processing theorem used here is the von Neumann algebra
relative-entropy monotonicity theorem; see D. Petz, *Quasi-entropies for
states of a von Neumann algebra*, Publ. RIMS 21 (1985), 787--800,
[primary text](https://ems.press/content/serial-article-files/42275).
It is not a statement about a guessed density matrix.
For example, take `B'=2 c^*(g)c(g)-I`, with normalized real smooth
scalar smearing times one fixed unit spinor, supported outside the ball.
The even scalar Fourier modulus and odd kinetic symbol give
`omega(c^*(g)c(g))=1/2`. This even operator commutes with the entire ball
CAR algebra. The conditioned vectors have smooth stress expectations on
the ball; compact smooth smearing and positive separation make their Wick
contractions and spatial derivatives finite. Thus the pairings in (9)
exist for this class.

Writing the linear energy variation as `epsilon*k`, (9)--(10) imply
`|Delta S_B-epsilon*k|<=epsilon^2/[2(1-r0^2)]`, so
`delta S_B=delta<K_B>`. This supplies an actual matter tangent first law
on stated state directions. It does not assert generalized-entropy
stationarity, a MaxEnt geometry equation, or the area term of G4.

## Controlled small-ball matter coefficient

Let `e(x)` be the varied model-time energy density at the central time,
with spatial Hessian operator norm at most M on the ball. Taylor's theorem,
radial symmetry and (7) give

```
delta<K_B> = [8*pi^2 R^4/(15v)] e(0) + remainder,
|remainder| <= 4*pi^2 M R^6/(35v).                          (11)
```

The exact integrals are `integral_B w_R=4*pi R^4/15` and
`integral_B w_R |x|^2=4*pi R^6/35`. The linear Taylor term integrates
to zero. For a fixed smooth state direction, M is bounded on a shrinking
neighborhood and (11) is `o(R^4)`. In Lorentz-time units `T00=e/v`, this
is precisely the matter coefficient of
`Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean`.

This closes the *continuum matter kernel* coefficient and its remainder.
It is not a joint estimate in lattice spacing, collar width, and ball
radius. FA82's convergence of bounded field words does not by itself imply
convergence of unbounded lattice stress or lattice modular generators.
Those uniform refinement estimates and the geometric/edge remainders of
`ScalingTailReceipt` are not imported from (11).

## Attachment to the existing Einstein stack

The target is the G1--G6 list in the spacetime paper and
`EinsteinRegisterData` in
`Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean`.

| Requirement | Supplied on this sector | Not supplied |
| --- | --- | --- |
| G1 | Minkowski domain, Lorentz action, smooth free-field equations, causal localization | Identification with the native event/count and varying geometric readout |
| G2 | The same-vacuum round-cap modular action, exact `2*pi` normalization and the same propagation stress | Physical cap access and the full interacting/massive matter sector |
| G3 | Symmetric flat-space operator Ward identity independent of Einstein dynamics | A curved, dynamical, universally coupled matter Ward identity |
| G4 | Vacuum-relative matter first law (9)--(10) on its explicit domain | Edge/center area term, renormalized area density, generalized-entropy stationarity |
| G5 | All mathematical Lorentz directions and the controlled matter remainder (11) | Common lattice/collar/chart error limit and source-selected observer coverage |
| G6 | One common centered free vacuum and a translation-fixed matter normalization | Newton calibration, dynamical geometry and universal coupling to missing sectors |

There is no dynamical gauge constraint or boundary flux algebra in this
free Dirac CAR sector. An additive scalar in a modular generator has zero
variation at fixed normalization; it cannot stand in for the stack's
`Z` with `delta<Z>=delta A/(4G)`. The finite central-sector entropy split
in `Entropy.lean` is not a proof of such an identification here. An
area-divergent regulator entropy would likewise require a renormalization
prescription and an independently specified gravitational coupling.

The missing geometric response is substantive. In the common flat metric,
the vacuum has zero normal-ordered stress and fixes a putative cosmological
term to zero. A nonzero positive-energy one-particle wave packet in the
same field has `integral T00>0`. The same flat Einstein tensor cannot
equal `kappa*T` for both states for any `kappa>0`. The exact control records
this integrated nonzero defect. This is a counterexample to obtaining
Einstein dynamics from these matter identities alone, not an exclusion of
a responsive OPH geometry or of other source-selected branches.

## Executable scope and custody

Exact symbolic checks cover all source frames, Clifford relations, polynomial
stress Ward/trace identities, filled-sea modular identities, the spinor
correction, diamond flow, and small-ball integrals. Mutation tests reject
deleted valleys, the wrong spinor coefficient, wrong units, and promotions
of the remaining gravity premises. They certify these finite identities;
the analytic proof and named imported theorem carry the continuum claims.

The new receipt binds its own proof/code, the precise current FA82 proof
inputs, and the existing geometry theorem sources. It also reports the
staleness of the older broad FA82 custody map after maintenance changes.
`prior_full_source_replay=false` is deliberate. Neither a stale old map nor
an `all_pass` flag is substituted for the historical mathematical theorem,
and no prior receipt or global ledger is edited by this package.

`new_checks_pass` records this package's exact replay. `dependency_custody_pass`
records whether the current broad FA82 source map matches; `all_pass` is their
conjunction. A false broad-custody flag remains visible even when the new
identities pass. The historical checkpoint's 512 source bindings were separately
checked against their Git objects, with zero mismatches.

Run `python -m pytest -q -W error model/weyl_geometric_modular`, then
`python -m model.weyl_geometric_modular --write` and
`python -m model.weyl_geometric_modular --check`. A successful new replay
does not turn a false dependency-custody flag into a true one.
