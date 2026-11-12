# Metric response, regional force errors and capacity histories

## Metric readouts of the reciprocal source

The source, preparation and update law are those of
[`finite_area_response`](../finite_area_response/README.md). There are two
geometry values `r0=9/10`, `r1=11/10`, eight matter levels in each, and
supplied inverse temperature `beta=1/sqrt(14)`. The scalar action density
and its metric dictionary give

```
C(r)=diag(-1/r,1/r,r²,r²),
g=sqrt(-det C) C^-1=diag(-r²,r²,1/r,1/r).
```

Consequently the proper volume is one, the two-face area is `A=2/r`,
the lapse is `L=r`, and the radial and transverse squared coordinate
speeds of null vectors are `1` and `r³`. These are readouts of the
specified action family, without a physical spacetime identification.
For a fixed branch, `d tau=r dt` and the proper generator is `H_r/r`,
so `exp(-i H_r t)=exp(-i(H_r/r) tau)` on the complete matter carrier.
Equal coordinate durations have proper-duration ratio `9/11`. This is
not a calculation of frequency exchange between separated clocks.

Apply the complete matter randomization mixture of strength epsilon,
then one original thermal update. The owning exact calculation gives
`Delta p1=epsilon d`, with `d` strictly positive and approximately
`0.00155196518`. The intervention itself leaves geometry unchanged.
For every two-valued metric observable f,

```
Delta <f>=epsilon [f(r1)-f(r0)] d.                 (1)
```

This proves, without additional simulation or conditioning,

```
Delta <L> = -99/200 Delta <A>,
Delta <L>/epsilon approximately 0.0003103930363,
Delta <v_transverse²>/epsilon = (301/500)d,
Delta <V>=Delta <v_radial²>=0.                    (2)
```

The source spectra are enclosed by exact rational inertia counts; bounded
exponential series enclose the complete thermal rates. The receipt reports
rational intervals rather than treating the displayed decimal values as
proofs. Independent tests multiply the original sixteen-state transition
matrix at four mixture strengths. Equal conditional Hamiltonians erase
the metric response, and conditioning on accepted cross-geometry moves
changes the answer. These controls retain the source's bath, readout,
preparation and controller assumptions.

The two geometry labels do not supply a spatial potential, curvature,
light-bending path, or universal Newton strength. Their distinct lapse
readout is a finite source effect. Promoting it to observed gravitational
redshift requires a spatial source and a clock-comparison experiment.

## Regional entropy errors to a Newton-force error

The existing
[`matched_einstein_limit`](../matched_einstein_limit/PROOF.md) proves that
ten matched timelike regional measurements bound every component of

```
E=G+Lambda g-(2pi/eta)T
```

by

```
epsilon_E = 255/(32pi eta ell^4)
            * (e_bulk+eta e_area+e_stationarity).       (3)
```

Its exact decoder constant is `17/8`. This implication requires the same
source, state, proper clock, regional volume, stress and curvature
identifications. The reciprocal two-geometry process does not provide
those premises. Equation (3) is reused directly; null-tensor balance is
not imposed as an additional independent premise.

To state an observable error consequence, also supply a static spherical
weak-field domain with metric
`ds²=-(1+2Phi)dt²+(1-2Psi) dx²`, pressureless matter, and negligible
cosmological term or a consistently subtracted cosmological background.
At first order, direct differentiation of the metric gives

```
G00=2 Laplacian(Psi),
sum_i Gii=2 Laplacian(Phi-Psi),
Laplacian(Phi)=4pi G_N rho+(E00+E11+E22+E33)/2,
Laplacian(Psi)=4pi G_N rho+E00/2.                       (4)
```

The Python executor constructs the linearized Ricci tensor independently
of these target equalities and verifies all displayed identities. The
source and residual must be spherical; anisotropic residuals require a
different Green-function estimate. Each residual component bounded by
epsilon_E gives Poisson-source errors at most `2 epsilon_E` and
`epsilon_E/2`. Their factor-four difference is necessary: equal-sign
diagonal residuals attain the first bound. Energy density alone cannot
control the lapse potential.

If these residual sources vanish outside radius R, integration of the
radial Poisson equation with a regular origin and vanishing potential at
infinity gives, for `r>=R`,

```
|delta Phi'(r)| <= 2 epsilon_E R³/(3r²),
|delta Phi(r)|  <= 2 epsilon_E R³/(3r),
|delta Psi'(r)| <= epsilon_E R³/(6r²),
|delta Psi(r)|  <= epsilon_E R³/(6r).                    (5)
```

Uniform positive error sources attain the scalar radial bounds. Thus a
matched regional error certificate controls the inverse-square force
without a second tensor reconstruction. A finite box needs its own
boundary-error term; nonlinear metric remainders and pressure also need
separate bounds. Neither (4) nor (5) asserts that a native source meets
(3), fixes G_N, or supplies an observed gravitational force.

## Capacity law without a fitted amplitude

Reuse the theorem in
[`FixedCapacityWLaw.lean`](../../Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean).
On its declared branch `rho_DE=kappa/N(a)`, positive constant kappa and
the separately conserved Friedmann continuity relation imply

```
w(a)=-1+(1/3)d ln N/d ln a.                            (6)
```

Constant capacity gives `w=-1` and exact CPL coordinates `(-1,0)`;
nondecreasing differentiable capacity gives `w>=-1`. The power family
`N=N0 a^epsilon` gives `w=-1+epsilon/3`. Choosing epsilon=-1 violates
the monotonicity premise and gives `w=-4/3`; the no-phantom conclusion
cannot be transferred to that changed branch.

Supply, separately, a flat Friedmann comparison history with matter and
radiation density fractions Omega_m and Omega_r, and define
`E(a)=H(a)/H(1)`. For positive dark density the amplitude-free identities
are

```
D(a)=E(a)²-Omega_m a^-3-Omega_r a^-4,
N(a)/N(1)=Omega_DE/D(a),
N(a)D(a)=Omega_DE N(1),
rho_DE(a1)/rho_DE(a2)=N(a2)/N(a1).                     (7)
```

These are useful tests of the capacity history even when its absolute
amplitude is unselected. The input fractions and expansion law are
comparison assumptions, not outputs of the finite metric process.
The twelve rational control histories in the receipt verify the forward
and inverse maps, including growing, fixed and decreasing capacities.
They are algebraic controls, not observational data or fits.

The transformation `(kappa,N)->(c kappa,c N)` leaves every density
unchanged. Normalized histories also erase N(1). Consequently agreement
with (6) or (7) cannot determine the cosmic capacity or the numerical
cosmological constant. The conditional horizon relation
`Lambda ell_P²=3pi/N` additionally requires a source-selected global
record capacity, a same-carrier horizon identification and an independently
attached physical scale. The finite local record counts and raw carrier
dimensions are different quantities. The source counterfamilies in
[`CapacityNonidentifiability.lean`](../../Lean/ObserverPatchHolography/CapacityNonidentifiability.lean)
and [`GlobalCapacityAttachment.lean`](../../Lean/ObserverPatchHolography/GlobalCapacityAttachment.lean)
prove why local saturation or cut composition alone cannot supply
that selection; no weaker numerical fixed-point search is substituted.

## Executable and formal scope

`response.py` owns (1)-(5), `capacity.py` owns (6)-(7), and the complete
source bindings are in `results.json`. `Transfer.lean` proves the exact
lapse/area ratio, the four-component Poisson bound, and capacity scaling
identities while importing the existing regional and capacity theorem
modules. None of these equalities closes a physical premise by renaming
a free tensor, capacity or scale.
