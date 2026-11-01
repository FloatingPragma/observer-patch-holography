# A controlled continuum test of the repulsive contact candidate

The existing Weyl/contact action has nontrivial lattice scattering and an
isolated repulsive pair. Those facts do not yet supply an interacting
continuum theory. This package closes a precise test of that inference:
at zero total momentum, with the **same model energy zero**, `c=v/a` and
fixed `v>0`, its positive-energy two-body scattering tends to the identity
as `a->0`, **uniformly over every nonnegative contact scaling U(a)**.

This is an analytic exclusion on the existing contact family. It does not
exclude a different physical vacuum, finite density, energy counterterms,
additional interactions or the whole OPH class. It does not assert a
many-body continuum limit. The finite-spacing interaction and bound pair
remain valid results. It advances CP-06/#1142 and informs CP-07/#1143.

## Exact principal-value sign, with every node retained

Use the already established [on-shell contact measure](../contact_scattering/ONSHELL.md).
At `c=1` it is an even probability measure `mu` on `[-sqrt(12),sqrt(12)]`.
The singlet has no dark zero-energy point mass. For `0<|t|<2` its density is

```
rho(t) = t^2 B(t^2),
B(t^2) = A(t^2/4)/(4 pi^2),
A(z) = sum_(n>=0) A_n z^n,   0<A_n<=1.
```

This is the complete eight-node density, not one selected low-momentum
pocket. The inverse second moment `M=integral t^-2 dmu(t)` is finite:
the displayed density controls zero and the remaining measure is bounded
away from zero. Since the support obeys `t^2<=12` and mu has mass one,
`M>=1/12`.

Let `a_1(x)=Re g_1(x+i0)`. For `0<x<=1/2`, the even-measure identity gives

```
a_1(x)/x = -M + x^2 PV integral dmu(t)/(t^2 (x^2-t^2)).
```

The part with `|t|>=1` is nonpositive. On `(0,1)`, expand B in its uniformly
convergent differentiable series. For each monomial,

```
I_n(x) = PV integral_0^1 t^(2n)/(x^2-t^2) dt
       = x^(2n) atanh(x)/x - sum_(k=0)^(n-1) x^(2k)/(2(n-1-k)+1)
       <= x^(2n) atanh(x)/x.
```

Termwise principal-value integration is justified by convergence in C1
on this compact interval; subtracting the numerator's value at t=x
removes its singularity. Also `atanh(x)/x<=1/(1-x^2)` and
`B(x^2)<=1/[4 pi^2(1-x^2/4)]`. Therefore, using `pi>3`,

```
a_1(x)/x <= -1/12 + 2x^2/[4 pi^2(1-x^2/4)(1-x^2)]
          < -1/12 + 8/405 = -103/1620 < -1/16.
```

No extrapolation of a broadened finite grid supplies this sign. Scaling
the actual free Hamiltonian gives `g_c(E)=g_1(E/c)/c`, and hence

```
Re g_c(E+i0) <= -E/(16 c^2),
rho_c(E) <= 4 E^2/(15 pi^2 c^3),   0<E<=c/2.
```

## A bound uniform over all repulsive coupling scalings

The rank-one bright-channel eigenvalue is

```
S(E) = (1-U a_c(E)-i pi U rho_c(E))/(1-U a_c(E)+i pi U rho_c(E)).
```

For `U>=0`, the denominator's real part is at least both one and
`U |a_c(E)|`. Thus the full scattering operator on its absolutely
continuous subspace obeys

```
||S(E)-I|| <= min(2, 128 E/(45 c), 8 U E^2/(45 c^3)),  0<E<=c/2.
```

The last term requires a specified finite U; the middle term is independent
of U. The other reducing scattering channels are identities already. Point
spectrum is outside the scattering operator's domain and has not been
discarded from the model Hilbert space.

For every `Emax>0`, all `0<E<=Emax` obey the bound whenever `a Emax/v<=1/2`.
With `c=v/a`, it follows uniformly on that entire positive-energy window
that `||S-I||<=128 a Emax/(45 v)->0`, even when U(a) diverges. Fixing U
gives the stronger cubic-in-spacing estimate. A repulsive resonance cannot
cancel this denominator in the stated positive-energy domain.

The same energy convention also sends every isolated **above-band** pair
at K=0 beyond `2 sqrt(3) v/a`. Any such pair leaves every bounded energy
window as a shrinks. This statement concerns a pair that exists; it does
not assert existence for every U, exclude embedded states, or infer a
physical mass from a band-edge separation.

The free kinetic symbol still has all eight nodes. At each `pi*b`, for
`b in {0,1}^3`, its expansion in the physical momentum label p obeys

```
||h(pi*b+a*p)+v sum_j (-1)^b_j sigma_j p_j||
  <= v a^2 sum_j |p_j|^3 / 6.
```

Thus removing a while merely retuning the repulsive onsite U cannot turn
this zero-reference two-body family into the desired nontrivial continuum
scattering sector with a finite-energy version of its above-band pair.
This identifies which additional vacuum/dynamics construction is needed;
it does not identify any Weyl node as a physical massless species.

## Retained boundaries and reproduction

Negative energies can have a different principal-value sign and are not
covered. The filled interacting vacuum, alternative energy normalization,
nonzero total-momentum scaling, many-body asymptotic completeness, gauge
fields, chiral/taste removal and gravitational response remain separate.
The nonnegative U class is fixed before this calculation; choosing an
attractive term or a counterterm changes its stated scope.

The exact certificate evaluates the rational inequalities above. Tests
independently integrate principal-value monomials, compare all eight
kinetic expansions with the actual field symbol, and verify the density
bound against the earlier exact angular certificate. Separate positive-
broadening Brillouin quadrature checks the Green function normalization;
it is explicitly not an on-shell proof or a numerical continuum limit.

```sh
python -m model.contact_continuum --check model/contact_continuum/results.json
python -m pytest -q -W error model/contact_continuum
```
