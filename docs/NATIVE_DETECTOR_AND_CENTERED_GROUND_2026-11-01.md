# Native energy measurement and the centered continuum limit

The [continuum field construction](CONTINUUM_FIELD_AND_SPECTRAL_READOUT_2026-11-01.md)
now has a more precise reference-state choice. A signed-energy symmetry
constraint identifies a necessary charge-clock correction for the stated
canonical field dictionary. Requiring particle-hole invariance derives
that correction directly from the source and supplies a new ground.
The finite cube's corrected ground is uniquely selected, rigorously half
filled and separated from every other charge sector. Its one-site states
are faithful on the complete local algebra. A complete native finite-clock detector
and an exact regional modular construction also extend the original
ground's operational and entropy results. A further uniform spectral
estimate now determines the centered bounded-coupling continuum: its
entire canonical fermion hierarchy is free and massless.

These constructions preserve the distinction between the original and
centered dynamics. Their energies are in model units. They do not identify measured particle species or establish the common
Einstein premises. The proved free limit excludes this particular
bounded-coupling trajectory from the required interacting massive
endpoint; a different justified trajectory or action is necessary.

## Deriving the corrected clock and selecting its ground

The [CPT constraint](../model/continuum_cpt_constraint/README.md) uses all
sixteen canonical spin/valley components. In the original bounded-contact
continuum their averaged signed response has first moment `U0/2`.
Any vacuum antiunitary that preserves the ground and energy generator
and maps this same field multiplet unitarily to its adjoints forces that
measure to be reflection symmetric, with zero first moment. For `U0>0`
these requirements are incompatible. This is a theorem about that
specific field, ground and clock dictionary, not a no-go theorem for
OPH or a deduction from Lorentz covariance alone.

The asymmetry has a bounded energy witness. If
`C=v^2||gradient phi||^2+U0^2/2`, its total variation distance from its
reflection is at least `(U0/2)^2/C`. For broad normalized Gaussian
packets this lower bound approaches `1/2`. Relabeling individual
spectral lines cannot remove it.

A charge rephasing translates the signed frequency by `E -> E-kappa`.
Zero first moment therefore requires `kappa=U0/2`. This is a change of
charged dynamics, not an additive energy-origin shift. It need not
preserve the old ground's positive-energy property.

The [centered source construction](../model/centered_source_ground/README.md)
derives the corresponding exact finite law. The original bipartite
source has a canonical particle-hole transformation satisfying

```
P(H)=H+U(V-N),
P(H-kappa N)-(H-kappa N)=(U-2kappa)(V-N).
```

Thus particle-hole invariance on the full Fock algebra fixes
`K=H-(U/2)N`. Its full lowest-energy projector defines the new finite
reference. Centered finite-ground projectors, a particle-hole invariant
exterior extension, translation averaging and the established local
dynamics bounds then construct a centered thermodynamic ground. It has
positive GNS energy and exactly one fermion per site in expectation.
The same physical-smearing argument constructs its own eight-valley
continuum net. The averaged addition/removal symmetry and zero first
moment are now obtained on this new reference.

## Exact finite checks and energy readouts

On the actual cube at `c=1/8,U=1`, the old global operator bound
`H >= -sqrt(3)/4` and the new particle-hole symmetry prove that every
centered ground vector has exactly `N=8`. This is sharp charge-sector
selection, stronger than merely having mean particle number eight.
The old `N=4` ground has centered energy `-2-sqrt(3)/4 > -4`, while an
eight-particle singly occupied trial has expectation `-4`; the old
reference therefore cannot be reused.

A product of four exact dimer grounds improves the finite bounds to

```
-4-sqrt(3)/4 <= E_K <= -2-sqrt(17)/2,
charged-sector gap >= sqrt(17)/2-sqrt(3)/4-3/2 > 1/8.
```

The gap applies to every `N != 8` sector. Separately, exact checks establish
real spin-independent hopping, a connected bipartite graph with four sites
in each part, and strictly positive on-site coupling. Together with the
proved half filling, these are the hypotheses of [Lieb's ground-state
theorem](https://arxiv.org/pdf/cond-mat/9311033), Theorem 1. The centered
cube therefore has a unique global ground with hidden spin zero. This
imports a proved theorem with its model premises checked; it does not
supply a quantitative neutral excitation gap or a relativistic mass.

The actual centered cube also has a faithful one-site state, without
using a fitted eigenvector. Exact signed vertex symmetries make its local
double occupation `d` site independent. A kinetic hop takes every
zero-doublon, eight-particle configuration out of that subspace. Combining
this fact with the kinetic norm and the strict dimer-trial energy
improvement gives

```
1/15000 < d < 1/8,
rho_site = diag(d,1/2-d,1/2-d,d),
K_site = -log(1/2-d) I + log((1/2-d)/d) (N_site-1)^2.
```

All four eigenvalues are strictly positive. The entropy first law now
holds for every one-site density variation, with relative-entropy
remainder at most `7500 epsilon^2 ||delta||_HS^2/(1-epsilon)` along
`rho_site+epsilon delta`, `delta=sigma-rho_site`, `0<=epsilon<1`. This
finite result does not infer faithfulness of an infinite-volume local
state or identify a geometric modular generator.

The complete two-site control is solved exactly on all sixteen Fock
states. Its centered ground is the unique two-particle singlet with a
doublon component. At hopping `t=1/16,U=1`, each fixed-spin addition or
removal branch has the following unnormalized weights:

| Positive excitation energy | Approximate model energy | Weight in each charge branch |
| --- | ---: | --- |
| `(2sqrt(17)-1)/16` | 0.452888203202 | `(1+1/sqrt(17))/4` |
| `(2sqrt(17)+1)/16` | 0.577888203202 | `(1-1/sqrt(17))/4` |

Each branch has total weight `1/2`. The signed response has both positive
and negative versions of these frequencies, mean zero and second moment
`65/256`. The one-site density is faithful on all four local states:
`diag(d,1/2-d,1/2-d,d)`, with `d=(1-4/sqrt(17))/4 > 0`.
These are exact dimer excitation energies, not the cube's complete
charged spectrum or measured particle masses.

## Complete native readout of the original ground

The [native detector](../model/native_spectral_detector/README.md)
keeps the original Hamiltonian and its invariant five-dimensional
ground mixture. Fifteen conditional native spin-rotation words realize
an exact finite spin average: five azimuths and three polar quadrature
nodes. A native occupation-copy and Majorana operation implement both
`c^dagger` and `c` charge outcomes.

Two clock qubits then control original-action evolution for times one
and two. Their positive Fourier transform produces the complete
four-outcome Fejer energy POVM, including coherent histories through
time three. All 78 previously certified positive energy lines contribute;
four broad bins do not individually resolve those lines.

The implementation retains every field row, selector record, charge
branch, clock outcome and first export failure. Its complete analytic
output-error bound is below `3/100`, including finite action splitting
and native source errors. No sampled shot data are invented. This
detector's reference is the original `N=4` ground; it does not prepare
the distinct centered `N=8` ground.

The executed eighteen-bit experiment propagates all 262,144 quantum
rows and all thirty selector/charge columns. Its joint probabilities,
rounded here, are:

| Clock bin | Addition | Removal |
| --- | ---: | ---: |
| 0 | 0.570054120459 | 0.246358517815 |
| 1 | 0.157161295637 | 0.001616803608 |
| 2 | 0.012037531812 | 0.000723176545 |
| 3 | 0.010747052064 | 0.001301502022 |

The observed output error against the exact original-action comparator
is about `4.50e-5`, below the certified complete bound `0.02879333504`.
The numerical probability-accounting residual is `3.65e-11`; clock
export failures are retained, totaling about `1.32e-27`. These are
computed instrument probabilities, not sampled experimental records
or four inferred particle masses.

## Stronger ground bounds and the resolved force channel

The centered ground admits a nonperturbative occupation bound. Integrating
the exact ground equation against the positive-energy semigroup controls
each nonzero kinetic momentum band. Since `1/sum sin^2(k_j)` is integrable
in three dimensions,

```
D_centered <= min(2,11(U/c)^2/8),
D_original <= min(2,11(U/c)^2/4).
```

The kinetic scaling itself follows from the source hopping. In its Fourier
convention, near valley `b`,

```
h_a(b*pi+a*p) = -c(a) sum_j (-1)^b_j sigma_j sin(a*p_j).
```

Its operator-norm difference from
`-a*c(a) sum_j (-1)^b_j sigma_j p_j` is at most
`|c(a)|a^3|p|^3/6`. A finite nonzero kinetic velocity in the fixed clock
and momentum chart therefore requires `a*c(a) -> v > 0`. Choosing
`c=v/a` gives the explicit remainder `v*a^2*|p|^3/6`. This argument does
not select `U(a)`, a surviving interacting trajectory, or physical units.
Bounded coupling remains a stated hypothesis of this construction.

Both occupation bounds improve to `O(a^2)` under that scaling. Their proof
uses bounded maps on momentum windows and does not assume pointwise
plane-wave operators in the ground representation. It sharpens the
microscopic convergence estimates without asserting macroscopic
Gaussianity or an interacting particle mass shell.

The [composite force channel](../model/continuum_composite_channel/README.md)
is also resolved algebraically. On each site,
`d=(2n_opposite-1)c` is an exact canonical CAR field, and

```
[H,c(f)] = -c((h+U/2)f) - (U/2)d(f).
```

Under the original bounded-contact refinement, its limiting equal-time
CAR family anticommutes with the original continuum fields. This follows
from controlled averaged fluctuations and a uniform finite-degree Dyson
approximation, tested in all polynomial GNS contexts. No extra
microscopic degree of freedom is inserted. The new family does not inherit a continuous finite-energy time evolution
on the enlarged representation merely from its equal-time algebra.

The [centered free-limit theorem](../model/continuum_free_limit/README.md)
now resolves the physical-energy question for `0<=U(a)<=Umax` and
`c(a)=v/a`. Rescale time by `1/c` and write `u=U/c`. An exact graded
anticommutator expansion cancels its degree-six term; every remaining
term contains a connecting free propagator. Combining its absolute
summability with the quadratic ground bound gives, uniformly over
**all unit lattice wavepackets**,

```
|F_d(s/c)-F_d,free(s)| <= 3u|s|+312 sqrt(Delta) exp(9|s|).
```

This uniform estimate avoids a divergent sum of absolute wavepacket
amplitudes. In the free filled sea, `d` is a normal-ordered cubic field:
its linear term cancels, so it creates three excitations. Two independent
momenta must lie near the Weyl nodes to have small total energy. The
resulting complete signed spectral measure obeys
`mu_d,free([-epsilon,epsilon]) <= 32 epsilon^6`. Gaussian smoothing then
proves that the actual composite loses **all weight in every fixed
physical-energy window**, uniformly over unit wavepackets.

Splitting the positive GNS energy at a cutoff turns this estimate into
vanishing time-integrated contact force on both particle and hole
vacuum vectors. It also proves that the entire initial-time state is
the free filled sea. Finally, a compact-support Fourier argument proves
cyclicity of the free vacuum for fields outside any bounded region.
Together with the already established finite propagation, that upgrades
the vacuum-vector result to operator equality on the full limiting
canonical field space. No Lorentz symmetry is assumed in this step.
Thus the complete hierarchy, rather than only a two-point function,
is the free eight-valley Weyl theory.

The resulting continuum energy and rest-mass readout is therefore

| Centered continuum sector | Positive excitation energy | Relativistic rest mass |
| --- | --- | --- |
| Each of the eight Weyl valleys, particles and holes | `E(p)=v|p|` | `0` |

There are four valleys of each chirality. These are mathematical field
sectors, not assignments to measured particle species. The finite dimer
energies and cube charge-gap bound above remain valid at their stated
finite parameters; they do not become nonzero continuum rest masses on
this trajectory. A nonzero equal-time composite norm is compatible with
its energy escaping to infinity.

The required interacting massive endpoint must consequently use a
different justified coupling/refinement trajectory or microscopic action.
Its ground, clock, locality, spectral and interaction estimates must be
proved for that changed construction; the bounded-coupling theorems
cannot simply be carried over.

## Regional modular data from the particle reference

The [regional construction](../model/ground_regional_modular/README.md)
uses the original invariant `N=4` ground, rather than a separately
chosen thermal state. Its exact regional spectrum follows from the
Slater restriction and spin average. All eight sites and all 28 site
pairs are checked with signed CAR traces.

For adjacent sites the density has rank eight in a sixteen-state local
Fock space. On its positive support,

```
K_R = log(6)I + 16log(2+sqrt(3)) H_R + log(N_R+1).
```

This gives an entropy first law for every supported trace-zero
variation and an explicit relative-entropy remainder. The last term
is central on the number-preserving supported algebra. Zero eigenvalues
are retained: disturbances into the orthogonal kernel have a singular
entropy derivative and cannot use that supported first law.

The existing Einstein theorem stack remains the consumer of these
matter and entropy data. Its common geometric modular kernel, local
Lorentz stress, area stationarity, null matching, scaling tails and
independent normalization are still required. The finite central term
has not been renamed as a geometric edge, nor has the energy pairing
been assumed to be a diamond charge.
