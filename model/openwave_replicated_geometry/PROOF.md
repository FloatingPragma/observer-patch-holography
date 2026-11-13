# Replicated quantum geometry and its controlled static limit

The complete carrier is N two-state geometry registers times both full
eight-level native matter factors. Define

```
Q_N=(q/N) sum_i Z_i,
H_g=omega_N/(2N) sum_i(I-X_i),
H=H_g+exp(Q_N) K_source+exp(F Q_N) K_probe,
omega_N=omega*N*cosh(q/N)^N*tanh(q/N)/sinh(q).
```

Here q=1/8 and omega=1 are inherited. The frequency normalization preserves
the original source's finite-q far-field susceptibility exactly. It uses no
measured Newton constant and does not repair the failed physical-unit
attachment. At N=1 this is exactly the parent geometry. Every term is
positive, all native sectors remain, and the vacuum is the matter ground
times the all-X-positive geometry product. An untouched reference tensors
with the full construction.

## Complete weak response

The free geometry energies are r omega_N/N with binomial multiplicities
C(N,r), for every r=0,...,N. Write a=q/N, b=qF/N and x=tanh(a)tanh(b).
The exact matrix elements of exp(Q_N) and exp(F Q_N) between the vacuum
and a subset of r X-flipped registers are products of cosh and sinh.
Second-order spectral perturbation theory therefore gives the zero-source,
zero-probe log-lapse derivative

```
chi_N(F)=-(2N/omega_N) cosh(a)^N
          sum_(r=1)^N C(N,r) x^r/r.
```

The denominator vacuum lapse is cosh(b)^N. Every excitation multiplicity
appears; the sum is not a projection onto the symmetric carrier. The
frequency normalization gives d chi_N/dF at F=0 exactly
-2q sinh(q)/omega for every N. For F>0, divide chi_N by this slope times F.
The r=1 term, tanh(b)>=b-b³/3, and
sum C(N,r)x^r/r <= (1+x)^N-1 <= Nx exp(Nx) imply

```
1-q²F²/(3N²) <= ratio <= exp(q²F/N).
```

Thus the response approaches a common linear spatial profile uniformly on
0<=F<=1, with a retained nonzero source normalization. All coefficients are
exact rational enclosures of finite sums, using the existing exponential
bounds. Independent full-hypercube diagonalization checks small N.

## Fluctuations and complete coordinate readout

The vacuum mean lapse is cosh(qF/N)^N. Its relative variance is exactly
(1+tanh(qF/N)^2)^N-1 and is bounded by exp(q²F²/N)-1. The complete
collective-coordinate readout has N+1 outcomes q(2k-N)/N with original
probabilities C(N,k)/2^N. All multiplicities are retained. Unread readout
removes every X expectation and costs omega_N/2 in the source; it is not
free. Coherences within a coordinate outcome are retained by its Lüders
projection.

## Full finite-strength static classical limit

Fix any finite nonnegative native source/probe energies E,P and profile F.
For a geometry state write z=<sum Z_i/N> and x=<sum X_i/N>. The average of
the individual Bloch vectors has norm at most one, also for entangled
states, so x<=sqrt(1-z²). Spectral Jensen gives
<exp(q sum Z_i/N)> >= exp(qz). Therefore the full quantum ground energy is
at least the minimum on[-1,1] of

```
g_N(z)=omega_N/2*(1-sqrt(1-z²))+E exp(qz)+P exp(qFz).
```

This function is strictly convex: g_N''>=omega_N/2. For E+PF>0 its unique
minimum is inside(-1,0); exact rational interval bisection encloses it.
At a minimizer choose the complete product state with single-register
Bloch vector (sqrt(1-z²),0,z). Its geometry energy is exactly the displayed
classical value. Hoeffding's bounded-variable exponential inequality for
the N independent Z outcomes gives a trial-energy excess at most

```
Delta_N=(E exp(q)+P exp(qF))*(exp(q²/(2N))-1).
```

This proves an explicit O(1/N) ground-energy bracket without solving a
2^N matrix or discarding any sector. It is uniform on bounded E,P ranges.
The frequency numerator tends to omega*q/sinh(q); the difference between
the corresponding classical minima is at most half the frequency change.
The receipt includes this change when comparing with the limiting action.

For the exact quantum ground, strict convexity bounds its mean-coordinate
distance from the classical minimizer by sqrt(4 Delta_N/omega_N). Uniform
convexity of the exponential potential also gives

```
Var(sum Z_i/N) <= 2 Delta_N /
 [q² E exp(-q)+q² F² P exp(-qF)],
```

when the denominator is positive. The program keeps that condition. These
are static, single-mode classical-limit bounds on the complete source.
They neither prove a dynamical limit nor create a local geometry field.

The free geometry gap is omega_N/N and tends to zero in the retained clock.
With unscaled register energy the same weak susceptibility instead vanishes
as1/N. No uniform preparation-time bound or finite-strength test-probe
universality is inferred. Finite probe energy changes the classical geometry
through backreaction. The original physical Newton-strength mismatch,
continuum geometry and full Einstein requirements remain open.
