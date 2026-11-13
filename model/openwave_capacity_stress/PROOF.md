# Capacity-dependent stress and expansion records

The source is a declared finite register coupled to a volume coordinate. Its
work gives negative pressure at fixed capacity. A conditional expansion
comparison uses two supplied cosmological parameters and all fifteen
chronometer rows of each specified stellar-population table. Neither the
cosmological amplitude nor a physical horizon follows from this register.

## Finite source and complete records

Let an active register contain k binary records. Its dimension is 2^k; its
capacity is k bits, or k log(2) nats. For positive epsilon and volume V,
declare H(V,k)=epsilon V/k times the register identity. On two orthogonal
geometry states V0=1 and V1=3/2 this is a finite diagonal Hamiltonian. At fixed
k the work derivative gives

```
rho = H/V = epsilon/k,     P = -partial_V H = -epsilon/k,     w = -1.
```

The derivative is an operator identity for every register state. The coupling
is an architectural term. An equally large register with H=epsilon/k,
independent of V, has P=0 and rho=epsilon/(kV). Capacity without a stress
coupling does not imply vacuum pressure. This is consistent with the existing
[`CapacityNonidentifiability.lean`](../../Lean/ObserverPatchHolography/CapacityNonidentifiability.lean)
counterfamilies; their theorem is not extended to a physical cosmological
capacity by identifying a carrier dimension with a logarithmic record count.

The executed apparatus starts in V0 with the maximally mixed active register.
A geometry Hadamard preparation costs DeltaE/2. Evolution under the actual
diagonal source for t=theta/DeltaE supplies the relative phase theta. Complete
X-sign and register-word readout probabilities are

```
Pr(word,sign) = [1+sign*cos(theta)] / 2^(k+1).
```

The receipt contains every outcome, including zero probabilities, for k=1,...,6
and theta=0, pi/2, pi, 2pi. No word or sign is postselected. Both conditional X
states have the same mean energy, (E0+E1)/2; the mean readout work is zero.
The phase experiment measures an apparatus gap. It identifies no cosmological
clock or particle mass. Independent dense circuit tests check its probabilities
and preparation energy.

## Changing capacity is an energy exchange

Appending an active bit uses a degenerate reference bit in a maximally mixed
state. Its two outcomes are retained. Switching k to k+1 at fixed volume
changes the source energy by epsilon V[1/(k+1)-1/k]. This energy is controller
work. The reference is part of the complete apparatus; k is the active
register's capacity, not the dimension of the combined hardware.

For the four-stroke cycle (V0,k), (V1,k), (V1,k+1), (V0,k+1), (V0,k), all
four work terms telescope to zero. Omitting the two capacity strokes gives
a nonzero false work budget. The exact six cycles in the receipt retain these
terms. Disengaging the bit returns a mixed reference; no pure-memory reset or
entropy removal is asserted.

For differentiable volume and capacity, dE=-P dV-rho V dk/k. Thus the
mechanical pressure stays -rho, while capacity variation exchanges energy
with the controller at Q=-rho (dk/dt)/k. If instead one imposes a separately
conserved effective fluid, the existing
[`FixedCapacityWLaw.lean`](../../Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean)
gives w_eff=-1+(1/3)dlog(k)/dlog(a). This effective pressure absorbs the
capacity exchange. It is not automatically the source's mechanical pressure.
A metric lapse/volume dictionary and the Einstein branch are additional
premises for interpreting the register as a gravitational vacuum stress.

## Conditional expansion comparison

Supply the flat matter plus constant-density Friedmann branch and the rounded
Planck2018 inputs H0=67.4 km/s/Mpc and Omega_m=0.315. The source's fixed-k
equation of state then gives

```
H(z)^2 = H0^2 [Omega_m (1+z)^3 + 1-Omega_m].
```

This uses a calibrated dark-energy amplitude. No chronometer values enter
either input, but it is a retrospective conditional comparison, not a
source-only prediction or a prospective test. The cosmological constant's
absolute value, the physical horizon, and a common metric-source join are not
identified. OpenWave's no-fitted-Lambda requirement is therefore unsatisfied.

The pinned author tables and covariance prescription come from
[Moresco et al., ApJ898,82](https://arxiv.org/abs/2003.07362) and
[the author's repository](https://gitlab.com/mmoresco/CCcovariance), commit
881413330a7f1e1e5203607d6964db49b4c6c461. The measurements are from
[Moresco et al.2012](https://arxiv.org/abs/1201.3609),
[Moresco2015](https://arxiv.org/abs/1503.01116), and
[Moresco et al.2016](https://arxiv.org/abs/1601.01701).
The input parameters are from
[Planck2018VI](https://arxiv.org/abs/1807.06209).

Every one of the fifteen BC03 rows is used. The complete M11 table is a
stellar-model control; it is not an independent fifteen-row dataset to add
to the primary likelihood. A further BC03 control uses all stellar-population
model systematics rather than the author's omitted-outlier column. None of
these comparisons removes a measured row or fits a chronometer parameter.

The covariance is C=D+uu^T+vv^T. D contains each quoted squared uncertainty;
u and v are the measured H times the interpolated fractional initial-mass
function and stellar-population systematics. The interpolation and endpoint
clamping match the author's recipe. In particular the z=1.965 row uses the
last tabulated systematic at z=1.475; this is a declared prescription beyond
the simulated redshift range. The full covariance has nonzero off-diagonal
entries. No unsupported joint likelihood of the larger chronometer catalogue
is used.

All decimal data are exact rationals. A rank-two Woodbury inverse avoids
floating matrix inversion. Integer square roots enclose every H prediction.
For residual r=r0+e with |e_i|<=delta_i, positivity C>=D gives

```
|r^T C^-1 r - r0^T C^-1 r0|
 <= 2 sum_i |(C^-1 r0)_i| delta_i + sum_i delta_i^2/D_i.
```

The receipt reports outward decimal bounds on chi-square and all individual
predictions. A fixed-anchor chi-square tail uses fifteen degrees of freedom
as a floating diagnostic under the stated Gaussian model. It is neither a
probability that the construction is true nor an independent discovery.
Radiation in [0,0.0001] and the supplied parameter box H0 in [66.9,67.9],
Omega_m in [0.308,0.322] have explicit enclosing sensitivity intervals. These
are deterministic boxes, not a propagated Planck posterior.

[`CapacityStress.lean`](../../Lean/ObserverPatchHolography/CapacityStress.lean),
consumed by [`Transfer.lean`](Transfer.lean), proves the volume work, density, capacity
work, complete cycle and negative-pressure identities while importing the
existing capacity and conservation modules. These algebraic theorems supply
no cosmic-capacity selector.
