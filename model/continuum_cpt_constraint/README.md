# A CPT constraint on the unchanged continuum fields and clock

The bounded-contact continuum construction has an exact obstruction to
one natural physical identification. If `U(a)->U0>0`, its canonical charged
field multiplet, original ground and original time evolution cannot jointly
have a vacuum CPT antiunitary that maps these fields to their adjoints.
The obstruction is the nonzero **signed response first moment** `U0/2`,
not a mismatch obtained by fitting individual particle masses.

This is a conditional theorem about a stated field/ground/clock dictionary.
It does not prove that Lorentz covariance alone is impossible, nor exclude
other field dictionaries, unbounded coupling trajectories, or the full
OPH architecture. No CPT theorem is silently applied to the existing net.
The usual relativistic CPT results require additional field hypotheses;
see, for historical context, [Lüders, Proof of the TCP theorem (1957)](https://doi.org/10.1016/0003-4916(57)90032-5).
The elementary argument below
assumes and tests the precise antiunitary action it needs.

## Source inputs and the averaged spectral measure

Use the full eight-valley continuum net of
[continuum_field_net](../continuum_field_net/README.md), its positive GNS
generator `L`, ground vector `Omega`, and the unshifted first-moment theorem
of [continuum_ground_constraints](../continuum_ground_constraints/README.md).
Write `A_s=C_b(phi e_sigma)` for all sixteen pairs `s=(b,sigma)`, with
`phi` a real, even, unit-normalized scalar Schwartz function on `R^3`.
The sixteen components have equal-time CAR; no valley is discarded.

Let `P_L` be the spectral measure of `L>=0`. Define positive measures on
nonnegative **physical excitation energy** by

```
alpha(B) = (1/16) sum_s ||P_L(B) A_s^* Omega||^2,
beta(B)  = (1/16) sum_s ||P_L(B) A_s Omega||^2.
```

The probability measure `nu=alpha+reflection(beta)` is on signed response
frequency. Its total mass is one by CAR. Zero-energy weight, if present,
is retained in both terms. A negative frequency labels a removal event;
it is not a negative excitation energy of `L`.

The exact continuum sum rule gives

```
integral E dnu(E) = mu = U0/2,
integral E^2 dnu(E) <= C_phi
                    = v^2 ||gradient phi||^2 + U0^2/2.             (1)
```

Here is the second-moment justification, including its limit boundary.
At finite spacing divide the corresponding averaged measure by
`m_a=||J_a phi||^2`, which is positive for all sufficiently small `a`
and tends to one. The exact identities from
[spectral_readout](../spectral_readout/README.md) give
`M1=h_a+U(a)r_a I` and
`M2=h_a^2+2U(a)r_a h_a+U(a)^2 r_a I` as quadratic forms.
Tracing both spin components annihilates every Pauli matrix. Thus the
kinetic first moment and the cross term vanish exactly. Anticommutation
of the Pauli matrices gives the averaged kinetic second moment

```
(v^2/m_a) sum_j ||J_a D_(a,j) phi||^2 -> v^2 ||gradient phi||^2.
```

The original ground theorem gives `r_a->1/2`. Uniform second moments
allow first moments to pass, with first-moment tail bounded by `C/R`.
For second moments apply weak convergence to `min(E^2,R)` and then
monotone convergence: only the upper bound in (1) passes. No equality
of limiting second moments or survival of ultraviolet memory is assumed.

## The antiunitary Ward identity and its obstruction

Suppose a proposed physical dictionary supplies an antiunitary `Theta`
and a unitary sixteen-by-sixteen matrix `Q` satisfying

```
Theta Omega = Omega,
Theta L Theta^(-1) = L,
Theta A_s Theta^(-1) = sum_t Q_st A_t^*.                         (2)
```

A harmless phase in the first equality can be removed. For the standard
spacetime-inverting action on charged fields, real even `phi` makes the
spatial inversion and conjugation preserve the chosen packet. Arbitrary
unitary mixing of spins and valleys is allowed in `Q`. No separate charge
conjugation, parity or time-reversal symmetry is assumed, and no specific
value of `Theta^2` is required.

The second equality in (2) means equality as self-adjoint operators and
therefore gives `Theta P_L(B) Theta^(-1)=P_L(B)`. Applying `Theta` to the
vectors `P_L(B)A_s^*Omega`, preserving norms, and summing over `s` gives

```
sum_s ||P_L(B) A_s^* Omega||^2
 = sum_s ||sum_t conjugate(Q_st) P_L(B) A_t Omega||^2
 = sum_t ||P_L(B) A_t Omega||^2.
```

The last step is finite-dimensional unitarity, so it does not require
orthogonal excitation vectors. Consequently `alpha=beta` as measures,
`nu(B)=nu(-B)`, and its finite first moment is zero. Equation (1) rules
out (2) for every `U0>0` on the unchanged ground and clock. This is a
complete spectral statement, including continua and degeneracies; it
does not assume stable poles or equal particle masses in advance.

In a proposed relativistic realization one must first derive the relevant
CPT theorem's hypotheses and its action on the chosen fields. If that
derivation supplies (2), the obstruction applies. The existing finite
propagation cone, positive energy, and canonical CAR by themselves are
not asserted to supply that derivation. Nonlinear/dressed fields or a
different ground need a new calculation of their spectral moments.

## A bounded observable detects the asymmetry

Use `TV(p,q)=sup_B |p(B)-q(B)|`, half the total variation norm, and let
`nu^-` be the reflection of `nu`. Cauchy--Schwarz with
`|nu-nu^-|<=nu+nu^-` proves

```
(2 mu)^2 <= (2 C_phi) (2 TV(nu,nu^-)),
TV(nu,nu^-) >= mu^2/C_phi,
W1(nu,nu^-) >= 2 mu = U0.                                      (3)
```

The last inequality follows by integrating `E`, a one-Lipschitz function;
finite first moments make the Wasserstein distance well defined.
These are bounds, not numerically computed infinite-ground distributions.

There is also a concrete bounded odd witness. For `U0>0` put
`R=C_phi/(2mu)` and `f_R(E)=max(-R,min(R,E))`. The scalar inequality
`|E-f_R(E)|<=(E^2)/(4R)` follows from `(abs(E)-2R)^2>=0` outside
`[-R,R]`. Hence

```
integral f_R dnu >= mu/2,
integral (f_R/R) dnu >= mu^2/C_phi.                             (4)
```

Every reflection-symmetric law gives zero for this odd observable.
The same witness reproduces (3), with the stated TV convention.
The bounded Lipschitz spectral readout estimates from `spectral_readout`
apply to `f_R`; actual native realization and ground-state error bounds
are still needed before claiming measured infinite-volume values.

For `phi_l=(pi l^2)^(-3/4) exp(-|x|^2/(2l^2))`,
`||gradient phi_l||^2=3/(2l^2)`. At `v=1/8,U0=1` the exact bounds are:

| Packet width `l` | Signed mean | Second moment upper bound | `TV(nu,nu^-)` lower bound | Clipping radius `R` |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1/2 | 67/128 | 32/67 | 67/128 |
| 2 | 1/2 | 259/512 | 128/259 | 259/512 |
| 4 | 1/2 | 1027/2048 | 512/1027 | 1027/2048 |
| 8 | 1/2 | 4099/8192 | 2048/4099 | 4099/8192 |

For any fixed positive `U0`, broad packets make the lower bound approach
`1/2`. This cannot be removed by relabeling a finite set of spectral lines.

## The necessary charge-clock change is not a rest mass

For a gauge-invariant ground the global charge rotation has an implementer
with generator `Qcharge` fixing `Omega` and commuting strongly with `L`.
It is normalized by `[Qcharge,A_s]=-A_s`. Rephasing the annihilators as

```
tau'_t(A_s)=exp(i kappa t) tau_t(A_s)
```

changes the generator to `L'=L-kappa Qcharge` and translates signed
response frequencies by `E->E-kappa`. The first moment is `mu-kappa`.
Therefore the unique scalar rephasing compatible with a zero CPT first
moment is

```
kappa=U0/2.                                                    (5)
```

This is only a necessary condition. A constant added to the Hamiltonian
would cancel from excitation energies and cannot perform (5). A charge
subtraction changes the charged dynamics. Moreover, positivity of the
new generator requires the joint energy/charge support to satisfy
`E-kappa q>=0` for **every** excitation; original positivity `E>=0`
does not imply this. The old ground may cease to be a ground, and the
shifted law need not become symmetric or Lorentz covariant.

An exact full-Fock counterexample uses two empty fermion modes with
`H=n0+3 n1`, its empty ground, and `kappa=2`. Its averaged signed response
is `(delta_1+delta_3)/2`, with first moment two. After centering it is
`(delta_-1+delta_1)/2`, but the original empty vector is no longer a
ground: `H-2N` has eigenvalues `0,-1,1,0`. This excludes inferring ground
positivity even when centering happens to make the response symmetric.

Thus (5) identifies an exact required change *if* one insists on this
canonical CPT dictionary. Proving that a changed dynamics and reference
give the physical theory is separate work. The scalar term is not a
derived Lorentz-invariant mass and is not removed in any source receipt.

## Free and massive controls

For `U(a)->0`, the already proved full-CAR limit is free, with
`h_b(p)=-v sum_j (-1)^b_j sigma_j p_j`. Its signed spectral matrix is

```
P_(b,+)(p) delta_(v|p|) + P_(b,-)(p) delta_(-v|p|).
```

Addition and removal have physical energy `v|p|`; all eight free sectors
have rest mass zero. The two-spin trace is reflection symmetric and
satisfies (1) with mean zero. Factoring the orientation matrix into its
determinant and a proper rotation shows four Weyl chiralities of each
sign, `chi_b=-(-1)^(b1+b2+b3)`. Pairing opposite chiralities yields four
free massless Dirac multiplets as a mathematical field reorganization;
it neither removes valleys nor identifies measured species.

The exact Clifford control uses
`H_D=Z tensor (sigma dot p)+m X tensor I`. Its square is
`(|p|^2+m^2)I`, with zero trace even when `m` is nonzero. A genuine Dirac
mass couples the opposite-chirality blocks and is compatible with zero
averaged signed first moment. A scalar energy shift instead changes that
moment. This distinguishes the obstruction from a claim that CPT
prohibits massive particles.

Finally, `algebra.py` constructs every matrix of a two-mode full-CAR
example with one occupied negative band. For energies `p+mu,-p+mu`
and `p>mu>=0`, the positive ground generator and all four Fock states
are retained. A norm-preserving particle-hole antiunitary has exactly the
action (2) on its two fields and preserves the ground. It commutes with
the generator at `mu=0` and fails at `mu>0`; the physical addition and
removal costs are respectively `p+mu` and `p-mu`. This is an independent
finite check of the Ward identity's generator requirement, not a claim
that the interacting continuum has this two-line spectrum.

## Replay and scope

```
.venv/bin/python -W error -m pytest -q model/continuum_cpt_constraint
.venv/bin/python -W error -m model.continuum_cpt_constraint --write
.venv/bin/python -W error -m model.continuum_cpt_constraint --check
```

The cold receipt binds the imported continuum, ground and spectral proofs,
the present proof and exact controls. It computes rational asymmetry bounds
and complete finite counterexamples. It does not compute the nonzero-
coupling continuum distribution, execute its native detector, establish
Lorentz covariance, identify physical species, or discharge gravity premises.
