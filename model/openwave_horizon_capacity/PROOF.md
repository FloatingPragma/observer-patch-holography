# Horizon capacity and expansion histories

The finite capacity source has mechanical pressure minus its density at fixed
capacity. The conservation-based effective equation of state is supplied by
`EinsteinBranch/FixedCapacityWLaw.lean`: w=-1+(1/3)d log N/d log a. Combining
this proved algebra with a declared area law N proportional to R^2 gives
w=-1+(2/3)d log R/d log a. Horizon identification and the Einstein/Friedmann
branch are explicit additional premises.

A fixed radius gives w=-1, without selecting its density amplitude. A radius
comoving with the expansion gives w=-1/3. With dust, the latter has positive
deceleration parameter at every finite scale factor. Neither statement needs
a fitted equation-of-state parameter.

For a Hubble radius R=1/H, density is proportional to H^2. If matter and the
effective dark fluid are separately conserved, the flat Friedmann equation
has a constant dark fraction. It reduces exactly to H^2/H0^2=a^-3, hence
w_dark=0 and q=1/2. This is the familiar Hubble-cutoff obstruction described by
[Hsu](https://arxiv.org/abs/hep-th/0403052), applied to the declared capacity
law. Maintaining mechanical pressure -rho while exchanging energy with
matter is a different model: constant dark fraction Omega gives
rho_m proportional to a^[-3(1-Omega)] and q=(1-3 Omega)/2. It can accelerate,
but has no matter-dominance-to-acceleration transition and does not preserve
separately conserved dust. No energy exchange is silently omitted.

For the future event radius, dR/dt=HR-1. Declaring
rho_dark=3 gamma^2/(8 pi G R^2) gives

```
w=-1/3-2 sqrt(Omega)/(3 gamma),
dOmega/d log a=Omega(1-Omega)(1+2 sqrt(Omega)/gamma).
```

This is the future-horizon branch of
[Li](https://arxiv.org/abs/hep-th/0403127). Its horizon choice and
initial dark fraction are additional data. The control gamma=1 is a declared
architecture choice. It does not derive the amplitude. For x=sqrt(Omega),
the exact primitive

```
F(x)=2 log x+log(1+x)-(1/3)log(1-x)-(8/3)log(1+2x)
```

satisfies F'(x)=2/[x(1-x^2)(1+2x)]>0. Cubing its exponential eliminates all
logarithms from the actual solver: for
R(x)=x^6(1+x)^3/[(1-x)(1+2x)^8], the equation is
R(x)=a^3 R(x0). This rational function is strictly increasing on (0,1).
Rational interval bisection and integer-square-root bounds enclose every
dark fraction and H prediction. An independent logarithmic implementation
and an ODE comparison test the same solution.

The future solution increases to x=1. The same rational equation gives
H_infinity^2/H0^2=3^8 x0^6(1+x0)^4/[16(1+2x0)^8], a strictly positive finite
limit. Hence R_event=1/(H x) tends to a finite radius and R_event/a tends to
zero. Integrating d(R_event/a)/dt=-1/a proves the future-event integral
condition for this gamma=1 branch; it is not independently imposed on a
numerical local solution. The amplitude retains its H0 and Omega_m inputs.

All four branch comparisons use the fifteen BC03 chronometer rows and their
complete correlated covariance. The M11 stellar-model alternative is scored
separately on the same observations. H0 and Omega_m are inherited supplied
Planck inputs; neither chronometer data nor branch comparison selects them.
The exact rank-two inverse and quadratic residual enclosure are reused from
`openwave_capacity_stress`. Every row is retained. The event-horizon control
at redshift 100 omits radiation and is not an early-Universe data prediction.

The comparisons discriminate capacity histories. They supply neither a
source-selected cosmic capacity nor a prediction of Lambda's amplitude.
