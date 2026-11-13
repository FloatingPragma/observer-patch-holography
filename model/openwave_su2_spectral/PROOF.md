# Certified full character spectrum and curvature response

## Exact physical restriction

Consider a single four-edge cycle with the same coefficients and full Haar
link space as `openwave_pure_su2_source`. Successive gauge averaging along
a spanning tree leaves the conjugacy class of the loop holonomy. The full
physical Hilbert space is therefore the SU(2) class-function space, with
orthonormal basis

```
chi_n(theta)=sin((n+1)theta)/sin(theta), n=0,1,2,...,
dmu(theta)=(2/pi)sin(theta)^2 dtheta.
```

Every link in a loop character has spin `n/2`; its four electric Casimirs
sum to `n(n+2)`. The character multiplication identity gives
`W chi_n=(chi_(n+1)+chi_(n-1))/2`, with `chi_(-1)=0`. Thus the complete
Hamiltonian is the positive infinite Jacobi operator

```
H_nn=alpha n(n+2)+beta,    H_n,n+1=-beta/2.
```

This is an exact Gauss reduction, including all representations. It is a
separate one-plaquette graph; it is not a reduction of all interacting
plaquettes in a three-dimensional cube.

## Exact infinite-tail inertia

Let `P_N` retain `0<=n<=N`, let `A=P_N H P_N`, and let `D` be the
orthogonal tail. Positivity of the magnetic term gives
`D>=d_N=alpha(N+1)(N+3)`. The only off-block entry is `b=-beta/2`
between `N` and `N+1`. For a rational test energy `x<d_N`, the Schur
complement obeys

```
A-x-b^2/(d_N-x)|N><N| <= A-x-B(D-x)^-1 B* <= A-x.
```

The tail has no negative direction, so the negative index of the full
operator is the index of that Schur complement. Exact rational Sturm
recurrences count the negative eigenvalues of both bounding tridiagonal
matrices. Whenever the two counts agree, that is the exact count for the
infinite source. Bisection accepts endpoints only with these agreeing
counts; an unresolved count increases the character cutoff. The recorded
intervals for the first four levels have width below `2^-40`.
Finite floating eigensolves propose brackets only. Neither the endpoint
acceptance nor the infinite-tail certificate uses their eigenvalue accuracy.

## Ground vector and complete curvature measure

The accepted numerical ground vector on characters `0,...,64` is stored as
exact dyadics with denominator `2^48` and normalized. Its Rayleigh value and full
residual are computed as rationals. The residual includes the omitted
`65` component `b u_64`, so it is a residual of the infinite Hamiltonian.
If `mu` is its Rayleigh value and `r` the residual norm, the certified
first excited energy gives

```
distance(ground projectors) <= r/(E1_lower-mu).
```

The producer verifies a positive denominator. The curvature `W` has norm
one. A projector distance `delta` therefore bounds its mean error by
`2delta` and its variance error by `6delta`. The explicit trial variance
includes the `65` component of `W u`, and has a strict positive lower
bound after this error. The exact double-commutator identity gives the
centered energy upper bound `alpha` on the true ground, as on the full cube.

All 65 finite spectral weights are retained in the Galerkin calculation,
including any residual vacuum weight. These weights are not claimed to be
the entire infinite spectral measure. The first four full energies have
the separate exact tail certificate.

## Euclidean correlation enclosure

The ground comparison bounds the centered-vector correlation error by
`12delta`, using contraction of `exp[-t(H-E0)]` and `||W||<=1`.
For the finite evolution, round its complete eigenvector matrix to
denominator `2^52` and eigenvalues to denominator `2^44`. Exact integer
matrix products give the Frobenius Gram error `eta` and eigenvector residual
`r_basis`. For `eta<1/2`, polar correction to an orthogonal matrix bounds
the self-adjoint generator discrepancy by

```
2 r_basis + 4 max|E_j| eta.
```

The dyadic ground, spectral basis and energy witnesses are pinned scientific
inputs. Every replay recomputes their exact Gram and full-source residual
bounds; changing solver rounding does not define a different source or evade
the certificate. Independent large-cutoff eigensolves check the accepted
enclosures in the tests.

Both corrected spectral and original Hamiltonian generators are positive.
Duhamel bounds their contraction-semigroup difference by time times this
quantity. Replacing the polar basis by its rounded basis costs at most
`3eta`. Every scalar exponential in the correlation is evaluated with
exact rational Taylor bounds. For a negative exponent below `-64`, the
bound `0<=exp(-x)<=2^-92` follows from the positive Taylor lower bound
`exp(64)>2^92`. Floating numbers are outward-rounded displays.

An independent bound includes the infinite character tail in that
semigroup. Split the centered trial vector above character `m=16`.
Its tail norm is computed exactly. A diagonal-contraction Dyson series for
the tridiagonal hopping must make at least `N-m` steps to reach character
`N=64`, so its remaining boundary amplitude is at most

```
2 exp(beta t) (beta t)^(N-m)/(N-m)!.
```

Duhamel across the single boundary edge then bounds leakage by the time
times `|b|` times this amplitude plus the initial tail norm. The receipt
includes a factor two for the quadratic correlation, the initial outside
component of `W u`, the exact ground-energy uncertainty and the explicit
`exp(t E0)` normalization. No imaginary-time filter is relabeled as an
admitted real-time detector.

## Continuum boundary

All assertions concern the stated compact-group plaquette and source clock.
At fixed coupling `H(a)=H(1)/a`; a positive constant dimensionless gap
therefore diverges in that clock as `a` decreases. A selected nonzero
curvature correlator does not prove dense local decay in an infinite-volume
theory. The source supplies no joint spatial limit, four-dimensional
Yang--Mills reconstruction, or physical GeV glueball identification.
