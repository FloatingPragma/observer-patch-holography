# Native lapse response and finite linearized curvature

The carrier, native energy operator, geometry coordinate and source preparation
are those of `openwave_universal_energy_geometry`. No new physical mass or
Newton-strength input is introduced. This result derives a clock response and
finite difference identities. The assignment of that response to a physical
metric remains a separate hypothesis.

## The actual clock response

For native source energy E, the geometry block, after its scalar part is
removed, is `-h sigma_x+E sinh(q) sigma_z`, with h=omega/2. Its ground has
`<sigma_z>=-E sinh(q)/sqrt(h²+E² sinh²(q))`. The weak native clock at profile
value F reads

```
g_E(F)=cosh(qF)-E sinh(q)sinh(qF)/sqrt(h²+E² sinh²(q)).
```

The reference is the actual source-off response g_0=cosh(qF), not one.
Therefore the source-energy derivative of the log lapse is exactly

```
chi(F)=d_E log(g_E/g_0)|0=-2 sinh(q)tanh(qF)/omega.
```

This is a consequence of the complete two-geometry-state source. All eight
native energies remain; a derivative in a scaled native energy direction
multiplies this coefficient by that energy. No finite-strength clock
universality or physical metric interpretation is inferred. The parent retains
both geometry outcomes, reference flags and preparation/readout work.

For 0<=F<=1, the inequalities x-x³/3<=tanh(x)<=x follow by differentiating
x-tanh(x) and tanh(x)-x+x³/3. Hence

```
0 <= chi(F)+2q sinh(q)F/omega <= 2sinh(q)(qF)³/(3omega).
```

At q=1/8, the relative profile error is at most1/192. The far-field slope,
relative to the small-q coefficient -2q²F/omega, is sinh(q)/q. The package
encloses these expressions with the parent's rational exponential bounds.
It retains every one of the343 accepted dyadic profile coordinates.

## A common finite difference geometry

Extend the accepted profile by zero to the whole integer lattice. Set
`D_i=(E_i-E_i^-1)/(2a)`. These operators commute, including at the boundary
because no difference halo is discarded. The declared first metric variation
is h00=-2phi, hij=-2gamma phi delta_ij, h0i=0. The finite linearized Ricci
formula is obtained by replacing partial derivatives by D in its ordinary
first-variation expression. For any symmetric h, direct polynomial algebra
proves `sum_a eta_aa D_a G_ab=0` and zero curvature of
`h_ab=D_a xi_b+D_b xi_a`. These identities require neither a field equation
nor a continuum limit. The source is static, so its time derivative is zero.

On the stated static metric family the calculation gives

```
G00=2gamma Delta_c phi,
Gij=(gamma-1)(D_i D_j phi-delta_ij Delta_c phi), G0i=0.
```

The source profile has a certified nonzero mixed Hessian. Thus gamma=1 is
necessary and sufficient for vanishing spatial linearized stress in this
family. Pressureless static stress is the condition here; it is not proved
for the complete native matter/wave theory. No unobserved stress tensor is
defined to be the Einstein tensor to manufacture a physical field equation.

The parent's variational profile uses the nearest-neighbour Laplacian. It
cannot be substituted for the centered-derivative curvature Laplacian:

```
Delta_c phi-Delta_nn phi=(a²/4) sum_i Delta_nn,i² phi.
```

This is an exact shift-polynomial identity, including all halos. The receipt
retains its actual native-profile correction and the dyadic Poisson residual.
It does not set either to zero. The weak clock profile chi(F) also has the
explicit nonlinear-in-F error above. These budgets are necessary before the
two source descriptions can share a null-balance claim. No proper-area law,
Lorentz stress, nonlinear Einstein equation, physical Newton constant or
uniform physical refinement is proved by these finite identities.
