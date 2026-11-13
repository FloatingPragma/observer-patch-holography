# Native energy, autonomous geometry and attractive response

## Complete declared source

Rebuild the exact eight-level native matrix
`M0=D/2-(XXI+IXX)/7` and retain `K=M0-E0 I>=0`, with its actual ground
energy and all eight spectral intervals. Let the geometry be a complete
two-state register with

```
Q=q sigma_z,    Hg=(omega/2)(I-sigma_x),    q=1/8, omega=1.
```

For a source at the center and a probe at a site with profile value `F`, set

```
H=Hg+exp(Q)K_source+exp(F Q)K_probe.
```

Both native factors retain all eight levels. The complete carrier has dimension
128; an arbitrary idle reference tensors with it. Each term is positive.
The global vacuum has both native factors in their ground and the geometry
in the zero eigenline of `Hg`. Energy eigenlines reduce the action, but native
coherences and their complete instruments belong to the same carrier.

An external uniform lapse shift satisfies the operator identity
`d[exp(s+FQ)K]/ds=exp(s+FQ)K>=0`. Thus the native energy relative to its
actual ground is exactly its lapse charge, including reference-entangled
states. This is a declared energy coupling. The reciprocal scalar source has
positive-energy contrasts with opposite response signs, as its retained
counterexample proves. There is no equivalence map to that reciprocal action.
Positivity of all possible local kinetic stress densities is not asserted.

## Spatial profile and independent source strength

On a seven-site-per-axis three-dimensional cube use the Dirichlet Laplacian
with diagonal six and nearest-neighbour entries minus one. The mode profile
is the normalized central Green column, rounded to exact dyadics with
denominator `2^40`. Every accepted entry is positive and the central value
is exactly one. The code evaluates its integer Laplacian and exact stiffness
`S=F^T L F`; no fitted `1/r` array is inserted.

Let `e` be the noncentral Laplacian residual. The discrete maximum principle
and the one-axis quadratic barrier bound the inverse Dirichlet norm by
`(n+1)^2/8`. Applying this to the cube with its central value pinned bounds
the difference from the exact normalized harmonic profile by the reported
residual times that constant. The central-hole boundary can only decrease
the positive Green inverse bound. The accepted dyadic profile itself defines
the source, so this approximation is not silently replaced during execution.

The geometry is one declared spatial mode. Its inertia is chosen equal to
`S`, giving oscillator frequency one. The two-state harmonic compression
has coordinate size `q^2=1/(2 kappa S)`, so the displayed source determines
`kappa=1/(2q^2 S)`. Exponentiation is of this finite coordinate operator;
it is not claimed to equal compression of an infinite exponential.
No full local geometry field or representation-cutoff limit is inferred.

## Exact attraction

For source energies `Es,Et>=0` the smaller eigenvalue is

```
f(Es,Et;F)=h+Es cosh(q)+Et cosh(qF)
                 -sqrt[h^2+(Es sinh(q)+Et sinh(qF))^2],   h=omega/2.
```

Its mixed derivative is

```
partial_Es partial_Et f
 =-h^2 sinh(q)sinh(qF)/[h^2+(Es sinh(q)+Et sinh(qF))^2]^(3/2)<0.
```

Integrating over the positive source rectangle proves that the mutual energy
`f(Es,Et)-f(Es,0)-f(0,Et)+f(0,0)` is negative. The derivative of that mutual
energy with respect to `F` is strictly negative, because
`z/sqrt(h^2+z^2)` is strictly increasing. Consequently decreasing positive
profiles give an attractive finite-difference force after subtracting both
self energies. Exact rational exponential and square-root enclosures verify
all 49 excited native pairs at each of three separations, and certify a
strict inward force for the first pair. Geometry and matter complements
are retained. Position interventions are part of this static force
instrument; a freely moving body's geodesic is not identified by it.

At weak source strengths, the cross term is
`-2 Es Et q^2 F/omega=-Es Et F/(kappa S)`. For the exact harmonic profile,
`F/S` is the central Green column. Its large-distance continuum comparison
has coefficient `G_model=1/(4 pi kappa)`. This coefficient follows from the
declared source parameters, rather than a fit to a measured force. The
finite cube and finite-strength errors are not removed from the receipt.
If one subsequently supplies energy unit `u` and length unit `ell`, the
natural-unit Newton coefficient would be `ell/(4 pi kappa u)`. Choosing
stiffness to reproduce measured Newton strength is an empirical calibration,
and does not establish the strength criterion as a derived prediction.

## Autonomous observations and energy

The source executes its exact unitary from the uncoupled geometry ground,
with both a source-off and a source-on intervention and four native probe
levels. Complete geometry measurements include both coordinate outcomes.
For `Vm=H-Hg`, the identity
`i[H,Hg]+i[H,Vm]=0` proves combined energy conservation and reciprocal
energy exchange. The geometry oscillates and entangles; no unitary
relaxation to its ground is asserted. Preparing the dressed ground for
spectroscopy is a separate declared service.

The full measurement work is the postmeasurement energy minus the
premeasurement energy. An idle record may copy an orthogonal source flag,
and no unknown quantum state is copied. All remaining native sectors have
explicit zero preparation probability, rather than being deleted from the
action. The source contains no Einstein equation, physical clock conversion,
or universal continuum-gravity theorem.
