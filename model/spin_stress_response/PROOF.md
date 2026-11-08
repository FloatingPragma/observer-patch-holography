# Common-state energy and modular response

## Source, clock and physical operations

Use the protected spin carrier and the declared finite action of
`native_spin_clock_join/energy.py`, with its scalar energy origin retained:

```
h_xy = J(S_x.S_y-I/4),
e_x = (1/2) sum_(y adjacent x) h_xy - B epsilon_x S_x^1,
H = sum_x e_x,          j_xy = i[e_x,e_y] = -j_yx.
```

The nonnegative selection field is used only on a declared bipartite graph.
Setting it to zero gives the exchange clock. The finite ground preparation
does not change when a subsequent pulse changes the evolution action.
All time parameters below use `hbar=1` and the same declared clock scale.
The native observer implements each bond exponential with 62 source pulses,
including the scalar phase, and each local X probe with 25. Repair operations
may be scheduled between complete encoded gate words under their owning
complete-channel error bounds; no immunity during a gate word is assumed.

## Driven continuity and energy response

For real differentiable source weights `a_x(t)`, define the actual action
and densities

```
H_a(t) = sum_x a_x(t)e_x,          E_x(t)=a_x(t)e_x,
J_xy(t)=a_x(t)a_y(t)j_xy,          W_x(t)=dot(a_x(t))e_x.
```

Direct commutation gives operator identities on the full finite spin carrier:

```
partial_t E_x + i[H_a,E_x] + sum_y J_xy = W_x,                 (1)
d <sum_x phi_x E_x>/dt
 = sum_(x<y) (phi_y-phi_x)<J_xy> + sum_x phi_x<W_x>.           (2)
```

Thus `d<H_a>/dt=sum_x<W_x>`. For a source perturbation `a+lambda f`,

```
V_f=sum_x f_x e_x,
partial_lambda J_xy |_(lambda=0)=(f_x a_y+a_x f_y)j_xy.        (3)
```

Equations (1)-(3) hold on every finite graph. Disjoint local densities
commute, so the current has graph range at most two. The source weights
are controllable coefficients in the bond/onsite action. Calling them a
geometric lapse requires an independently supplied metric and clock map.
A weighted action is an explicitly programmed source intervention; its
coefficients are not asserted to arise autonomously as a gravitational field.
A recorded intervention contributes its complete instrument energy jump,
including failure outcomes. Equations (1)-(2) describe the unitary intervals;
they do not assert conservation of an unmodelled apparatus energy.

For time-independent `H` and a fixed Hermitian pulse `lambda V`, keep the
initial density `rho_0` and observation time `T` fixed. Differentiating the
finite propagator, with `A_H(t)=exp(iHt)A exp(-iHt)`, gives

```
d <A>_lambda/dlambda |0
 = i integral_0^T Tr rho_0 [V_H(s),A_H(T)] ds.                 (4)
```

This retarded response uses the actual input, without a stationary-state
assumption. Duhamel gives `||U'_lambda||<=T||V||` and
`||U''_lambda||<=T^2||V||^2`, uniformly for real lambda. Consequently

```
|<A>_lambda-<A>_0-lambda d<A>_lambda/dlambda|0|
 <= 2 lambda^2 ||A|| T^2 ||V||^2,                             (5)
D(rho_lambda(T),rho_0(T)) <= min(1,|lambda|T||V||).             (6)
```

Here `D` is half trace norm and `rho_0(T)` denotes the unperturbed evolved
state. The derivatives are taken before evaluating at zero in (5). Local
`V` gives a volume-independent norm bound; approximation of an infinite
ground marginal and a continuum stress limit are additional conditions.
The inequalities hold with arbitrary spectator references. Complete native
preparation, pulse and detector errors add to (6), and their expectation
error is at most twice the summed half-trace bound times `||A||`.

## Regional first law for the actual evolved state

Fix a region and observation time, and take the unperturbed source marginal
`sigma=Tr_outside rho_0(T)`. Assume its least eigenvalue is `mu>0` on the
declared common support. Write `delta=rho-sigma`, `Tr delta=0`, and
`epsilon=||delta||_1/2<mu`. Its modular readout is `K=-log sigma`. Then

```
S(rho)-S(sigma)=Tr(delta K)-D_rel(rho||sigma),
0 <= D_rel(rho||sigma) <= epsilon^2/(mu-epsilon).              (7)
```

For completeness, trace zero implies `||delta||<=epsilon` and
`Tr delta^2<=2epsilon^2`. Along `sigma+s delta`, every eigenvalue is at
least `mu-epsilon`. The logarithm derivative has quadratic form

```
integral_0^infinity
  Tr[delta (sigma+s delta+uI)^-1 delta (sigma+s delta+uI)^-1] du
 <= Tr(delta^2)/(mu-epsilon).
```

Integrating this Hessian with weight `1-s` proves the upper bound in (7).
Its nonnegative integrand proves the lower bound. The derivative at zero
therefore equals the modular pairing for every differentiable physical
source-pulse family through this marginal. Neither a thermal replacement
nor a fitted temperature is used. The estimate applies to noncommuting
regional perturbations. A growing-region application needs its own
quantitative spectral floor.

## Constructed preparation and complete modular readout

Take two sites, `J=1`, `B=3/8`, and staggered signs `(1,-1)`. Exact source
diagonalization gives the unique selected ground

```
Omega_B=(-|00>+3|01>-3|10>+|11>)/sqrt(20),
H_B Omega_B=-(9/8)Omega_B.
```

The complete selection spectrum is `{-9/8,0,0,1/8}`, so the displayed
ground is unique.

With `s=(|01>-|10>)/sqrt(2)` and
`q=(|11>-|00>)/sqrt(2)`, this is `(3s+q)/sqrt(10)`.
The exchange action is `H_J=-|s><s|`. Evolve this same prepared vector
under `H_J`; its first-site marginal is exactly

```
sigma_A(t) = [I+(3/5)cos(t) X]/2,
lambda_min(sigma_A(t)) >= 1/5.                               (8)
```

This is the two-site instance of the declared source family. It is not a
marginal asserted to equal the eight-site selected state or an infinite
ordered vacuum. Its uniform floor applies to the displayed one-site
marginal, without establishing a floor for those other states.

The complete X probe has probabilities `(1+/-r)/2`, where
`r=(3/5)cos(t)`. Assigning its two records the values
`-log((1+/-r)/2)` implements the spectral readout of the fixed reference
`K_A(t)`. Both records are used. The four reserved outcomes of the six-level
pointer have zero ideal probability; native complement records remain
explicit outcomes with their own error budget.

Pulse the bond coefficient by `1+lambda` for duration `t`, keeping
`Omega_B`, the energy origin and the observation clock fixed. The resulting
marginal is (8) with cosine argument `(1+lambda)t`. At `t=pi/4`,

```
X_A = d sigma_A,lambda/dlambda |0 = -3 sqrt(2) pi X/80,
dS_A/dlambda |0 = Tr X_A K_A
 = [3 sqrt(2) pi/80] log[(10+3 sqrt(2))/(10-3 sqrt(2))] > 0.   (9)
```

Thus (7) has a nonzero realized source tangent. At the same time,
`<H_J>=-9/10` and `<e_0>=-9/20` for every lambda. This is an explicit
countercontrol to identifying this regional modular response with the
physical local energy differential. The selected state is not stationary
under `H_J`; geometric stationarity cannot be inferred from its preparation.
The assigned energy density `e_0=H_J/2` has support on both sites; it is
not an operator on the single-site regional algebra. The countercontrol
compares the two specified expectation differentials without identifying
their operator domains.

The native experiment uses rational duration `3/4` and
`lambda in {-1/16,0,1/16}`. A signed amplitude tree prepares the encoded
integer vector from blank registers in 57 pulses. Each experiment has
144 pulses including the exchange and X-pointer word. The code executes
their six-level source matrices and retains both pointer probabilities,
four reserved zero outcomes and every first-failure weight, including
padding positions. There are 510 physical first-failure slots. No success
weight is divided by the success probability.

The preparation-angle error is bounded from the actual emitted tree
coefficients, with rational enclosures of its signed two-vector rotations.
All source-generator errors use the actual binary64 coefficients and the
exact readable source matrices. Their sum bounds the complete isometry
error. Floating diagnostics compare the literal execution with the
independently exponentiated spin action; they are not substitutes for the
rational error bound. A refinable implementation uses the owning native
observer's dyadic integrated-port interface and computable preparation
angles. Fixed binary64 execution is not a zero-error limit. Entropy of an
ideal regional state is inferred through certified intervals on its two
probabilities, never through a renormalized success state.

The implementation errors also control the state premise itself. Apply a
fixed trace-preserving decoder which sends off-code and failure output to
a declared regional state, retaining the failure flag in the archive.
Contractivity bounds its regional distance from (8) by the complete source
error `nu`. For `nu<1/5`, its actual marginal has eigenvalue floor
`1/5-nu`; the resolvent identity gives
`||log(sigma_actual)-log(sigma_A)||<=nu/(1/5-nu)`.
Equation (7) therefore applies to that actual normalized decoded reference
with its smaller floor. Comparing two implemented pulses adds both
complete source errors to (6). This completion is a fixed channel on all
outcomes, without conditional renormalization. The displayed logarithmic
readout uses the ideal reference weights, with the stated modular error
when interpreted on the decoded implementation.

## Attachment to the finite gravity theorems

Equation (7), on the same source-pulse tangent as (4), supplies the first
conjunct `entropyDifferential = modularPairing` of
`OPH.EinsteinBranch.FiniteFirstLawPremises` in `Entropy.lean`.
The physical split `K=2 pi B+Z` and an edge differential equal to the
central pairing require independent geometric and algebraic identifications.
Defining a symbol `B=K/(2 pi)` would satisfy a scalar rewriting without
identifying a boost or stress.

Equations (1)-(6) supply finite energy continuity, its source derivative and
its controlled retarded response. `ContinuumEinsteinPremises.ward` requires
all components of a common symmetric Lorentz stress tensor, including
spatial momentum transport. Its other inputs include the geometric Bianchi
identity, null balance, universal coupling and a calibrated Newton scale.
The finite spin result does not supply these identifications.

`SmallBallPremises.stationarity`, `continuumKernel` and `fixedVolumeArea`
require a responsive Lorentzian geometry and admissible fixed-proper-volume
variations. Equation (9) supplies a nonzero entropy differential on an
actual source experiment; it supplies no area variation cancelling it.
The fair-source regional modular-clock result belongs to its own reference
and carrier. Its exchange-clock collar limit is not a substitution for
the selected spin state used here.
