# Common native clock, arrival and bending records

## One complete source

Use exactly the native matrix, ground-relative energy, two-state geometry,
accepted three-dimensional profile, source strength and inertia of
`openwave_universal_energy_geometry`. Add the field vacuum and one scalar
wave on every site of its seven-by-seven-by-seven cube. In geometry
coordinate `Q=+q` or `-q`, define

```
D(Q)=diag exp[(1+gamma)Q F_x],
h_field(Q)=positive sqrt[D(Q) L D(Q)],    gamma=1.
```

The field vacuum has energy zero. The wave block contains this positive
operator plus the same `Hg`, source `exp(Q)K_source` and clock
`exp(F_clock Q)K_clock`. The native clock sits at `(3,4,3)` and the source
at `(3,3,3)`. The full action retains eight source levels, eight native
clock levels, two geometry states and `1+343` field states, giving dimension
44,032 before idle references. Every block is positive; the global vacuum
is unchanged. Reference flags and every unprepared native sector tensor
with the complete action.

The exact Jordan--Wigner identity of `openwave_joint_native_clocks` factors
the native eight-level matrix into three normal fermion occupations. With
the middle occupation unexcited, the other two independent clocks have
native energy labels `(0,3,1,5)`. The complementary spectator labels
`(2,6,4,7)` belong to the same carrier. The producer rechecks that exact
factorization. It does not infer two clocks from approximately equal gaps.

The field's spatial metric response is declared through `gamma`. A uniform
lapse changes native frequencies by `exp(s)`; the coordinate field frequency
changes by `exp((1+gamma)s)` because the spatial inverse-length response is
also included. This distinction is explicit. It is not a derivation of a
universal local stress tensor or of the metric response from the rest-energy
coupling alone. The field has no physical photon identification here.

## Joint preparation and instruments

An orthogonal source/reference flag retains equal source-off and source-on
interventions. Conditional on that flag, prepare the geometry ground of the
source with the probe in its ground, the two native clocks independently
in `|+>`, and an equal field-vacuum/wave preparation with an orthogonal
reference flag. The packet is defined by all its accepted dyadic complex
coordinates. The native clock energy is added to the source Hamiltonian
before propagation; no freely evolving clock phase is substituted afterward.

There are four independent X/Y clock-setting pairs. Their coin weight is
`1/4`; the complete clock effects, including both spectator values, sum to
the native identity. Terminal geometry-coordinate and field-position
measurements commute with these clock effects. Their recorded tables are
the actual complete probabilities. The field Fourier measurement is a
separate instrument, with every three-dimensional Fourier outcome retained;
it does not measure sharp position and momentum together.

For first arrival, evolve the whole source by half a parameter unit, then
project onto receiver sites with `x>=5` or their complement. Retain all
twelve possible arrival histories and the final no-arrival history,
including the field vacuum. Terminal clocks, geometry and field positions
are read on every history. Their original probabilities sum to one for
each source flag; multiplying by `1/2` gives the original source-flag law.
The source energy account includes the coarse receiver measurement and
every terminal readout. The total postmeasurement energy equals initial
energy plus the summed measurement work, with no lost failed branch.

The finite parameter labels specify ordered source updates. Operational
clock records are the two native Ramsey phases and their visibilities in
these same histories. A phase is circle-valued, with its period retained.
No external clock is used to unwrap it. Position versus native phase is
therefore a joint record on that circle, rather than a fundamental global
time coordinate.

## Finite error control

The profile consists of exact dyadics and the native energy intervals are
retained. Exact rational exponential bounds enclose every field coefficient.
For a numerical positive root `R` and exact positive stiffness `K`, the
Sylvester identity gives

```
||R-sqrt(K)|| <= ||R^2-K||/(lambda_min(R)+sqrt(lambda_min(K))).
```

The Dirichlet inequality supplies
`lambda_min(K)>=12 exp[-2|1+gamma|q]/(n+1)^2`.
The complete eigensystem Gram residual and multiplication-roundoff bound
supply the independent lower bound on `R`. The receipt verifies the
conservative denominator floor `2/5`. Thus these are comparisons to the
declared positive-root field, not a redefinition using the numerical root.

For each complete 686-dimensional wave block, let `V,D` be the numerical
eigensystem, `eta>=||V*V-I||` and `r>=||HV-VD||`. Polar correction gives a
unitary eigensystem with generator error at most `2r+4||D||eta`, when
`eta<1/2`. Matrix-product residual bounds use the binary64
round-to-nearest `gamma_n=n u/(1-n u)` rule and absolute-product norm
bounds. Scalar coefficient uncertainty, positive-root error and preparation
rounding are included separately. The bounds require this stated floating
arithmetic model, rather than an exact-real eigensolver assumption.
Each live residual must pass the fixed outward budgets `eta<=4*10^-11`,
generator error `<=10^-9`, and source-operator error `<=5*10^-11`.
The receipt uses those verified budgets so an arbitrary change of basis in
a degenerate eigenspace cannot change its accepted uncertainty statement.
Every actual scalar propagation phase is also compared with a 45-digit
interval exponential of its exact binary energy and parameter inputs.
The verified component error is below `10^-13`; its accumulated contribution
is included in the preparation/propagation allowance.

Duhamel bounds a unitary error by time times the generator error; replacing
the polar basis adds at most `3eta`. Projective history maps are contractions,
and their direct sum is an isometry. Telescoping therefore controls the
complete history state without multiplying by the number of outcomes.
For state error `epsilon`, the displayed conservative total-variation bound
is `2epsilon+epsilon^2`. This includes the field vacuum and final failure.

A bounded arrival coordinate in `[0,6]`, conditioned on success `p`, has
error at most `12 TV/(p-TV)`. The conditional transverse position error is
at most `12 TV` because the wave probability is exactly `1/2`. Both
source interventions are included in the difference bounds. The executed
arrival delay and inward transverse shift exceed these errors. Fourier
momentum bins and the angle of their mean provide a separate bending
diagnostic, with the corresponding bounded-observable error retained.

## Native clock response and its limit

The exact dressed clocks use the same geometry/source family with field
vacuum. For a probe strength `s`, write

```
f_E(s)=h+E cosh(q)+s cosh(qF)
       -sqrt[h^2+(E sinh(q)+s sinh(qF))^2],  h=omega/2.
```

Its first derivative at zero is the source-ground expectation of
`exp(FQ)`. Its second derivative is negative and has magnitude at most
`2 sinh(qF)^2/omega`. Hence

```
|(f_E(s)-f_E(0))/s - f_E'(0)| <= s sinh(qF)^2/omega.
```

Every denominator slope is at least `exp(-qF)`. The frequency ratio between
source-on and source-off clocks therefore differs from the common
zero-strength ratio by at most
`2s sinh(qF)^2 exp(qF)/omega`. This proves the common weak-probe response
for every native energy used as a probe. The seven full-strength ratios are
strictly redshifted and unequal, with about 288 parts per million between
their largest and smallest values. They are not claimed to establish exact
finite-strength clock universality. Their dressed vacuum ratios also are
not substituted for the actual joint wave/clock Ramsey phases.

## Controls and physical boundary

The `gamma=-1,0,1` controls execute the complete three-dimensional field
records with the native clock in its ground. At `gamma=-1`, `D=I` exactly;
the field decouples from geometry, so its source-on/off spatial and arrival
records agree. The other controls change the declared spatial response
without refitting the geometry strength. They are separate interventions,
not additional independently fitted physical predictions.

The construction supplies a common finite mechanism for attractive rest
energy, clock redshift, delayed arrival and bending. It gives no physical
Newton strength, seconds, continuum geodesic, Einstein equation or quantum
Yang--Mills gap. The pure-gauge finite-volume results have a separate carrier
and vacuum and do not imply a gap in this matter/wave theory.
