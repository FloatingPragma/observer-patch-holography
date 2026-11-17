# Entropy, energy and responsive finite geometry

The specified reciprocal scalar process has an exact area-response equation
with memory. Its fourteen memory modes and fifteen visible relaxation modes
are certified from the complete sixteen-state law. The same process supplies
a positive thermal entropy/area coefficient on its conditional equilibrium
family and an exact entropy-production balance with matter energy and every
bath record retained. These are finite heat-bath dynamics; deriving Einstein
dynamics requires a common spatial refinement and regional entropy/stress law.

The source is the unchanged pair of three-oscillator Hamiltonians at
`r=9/10,11/10`, with supplied `beta=1/sqrt(14)`. Its action gives the metric
`diag(-r^2,r^2,1/r,1/r)`, proper volume one and area `A=2/r`.
All sixteen energy states, all 64 cross-geometry collision edges, rejected
proposals and both native proposal laws participate. The result extends
the [one-step matter response](../finite_area_response/README.md) to every
finite update count. No Hamiltonian, reservoir gap or transition rate changes.

## A derived thermal area coefficient

Let `tau_g` be the actual conditional Gibbs state and `e_g=Tr(tau_g H_g)`.
Among all joint quantum states with geometry-1 probability p, the unique
minimum of `Tr(sigma H)-S(sigma)/beta` is
`sigma_p=(1-p)tau_0 + p tau_1`, with the two terms on their original blocks.
The source determines

```
eta=beta(e_1-e_0)/(A_0-A_1),       0.578<eta<0.579.
S(sigma_p)+eta<A>_p = constant-D_Bernoulli(p||p_equilibrium).
```

This family has a strict fixed-volume entropy/area maximum at equilibrium.
Its variation changes the probabilities of two fixed metrics.
The actual process does not preserve that family. For an arbitrary joint
state define `E_ex=Tr(sigma H)-[(1-p)e_0+p e_1]`. Its exact balance is

```
Delta S+eta Delta<A>-beta Delta E_ex = Delta S-beta Q >= 0.
```

Here Q is the total heat absorbed from the specified collision baths after
the preparation. The right side includes the relative entropy removed by
energy dephasing and subsequent population relaxation. Conditional matter
energy cannot be dropped: the full-state stationarity obstruction
applies to `S+eta<A>` without this correction. The derived eta is a thermal
energy-to-area conversion for these two geometries, not a Newton constant.

## Memory and complete records

The geometry deviation x obeys a source-derived linear equation
`x_(n+1)=a x_n+sum_(j<n) k_j x_(n-1-j)+f_n`.
The coefficient is `a=0.6391668062...`, and
`k_0=0.0217017000...>0`. The forcing is determined by the prepared matter
state. The same eliminated variables reconstruct its energy and heat at
each update. For the admitted field randomization, the first area response
is the certified `-0.0006270566...` times its mixture amplitude.
No Markov approximation of the geometry marginal is made.

Rational Gram pivots certify fourteen independent memory directions.
Reversibility and strict laziness put all memory poles on the positive
relaxation interval. A source-derived minorization gives an explicit
exponential bound for every memory tail. These poles describe equilibration,
without particle or spacetime interpretation.

All 512 proposal/bath outcomes admit an exact reversal involution for each
native law. Its likelihood ratio is `exp(-beta Q)` per collision. This gives
the full path entropy-production identity and, for faithful initial
populations, the integral fluctuation identity. The bath likelihood flow
differs from the von Neumann entropy change of the spent bath. Tapes are
retained; no memory-erasure cost is silently set to zero. Repeated use
requires the source's declared fresh thermal inputs and readout services.
Their preparation and pulse-controller work are supplied interfaces.

The [derivation](DERIVATION.md) gives the projection, positivity, variational
and full-record proofs. The update count is the source collision clock;
the Hamiltonians use their common model-coordinate clock. A physical
duration per update and a common proper temperature are not derived.
The result supplies no spatial refinement, local conserved stress tensor,
regional edge-area normalization or Einstein stationarity premise.

```sh
python -m pytest -q -W error model/reciprocal_area_dynamics
python -m model.reciprocal_area_dynamics --write
python -m model.reciprocal_area_dynamics --check
```

Replay isolates every source level by rational matrix inertia, encloses
thermal weights with a bounded Taylor series, and checks all memory pivots
with outward rational arithmetic. Independent tests use the original
quantum source channel and its complete records. The source pulse experiments
and the full reservoir-bank controller are not reexecuted by this package.
