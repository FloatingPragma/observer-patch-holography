# Continuum spin fields, finite entropy limits and joint spectroscopy

2 November 2026. This continues the [particle/modular construction](PARTICLE_SCALING_MODULAR_CONSTRUCTION.md) under the same five issue owners: [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033), [#1138](https://github.com/FloatingPragma/observer-patch-holography/issues/1138), [#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142), [#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143) and [#1145](https://github.com/FloatingPragma/observer-patch-holography/issues/1145). New derivations take priority; repository-wide CI is deferred at the user's request.

| Attachment | Established result | Scope |
| --- | --- | --- |
| [FA99](../ballistic_free_boundary/README.md) | The canonical charged continuum is free for `U=o(c^(2/7))`, improving the logarithmic sufficient bound | Same centered clock and all eight valleys; this is not a sharp interaction threshold |
| [FA100](../spin_canonical_field/README.md) | A nontrivial regular continuum net of Weyl (CCR) operators, continuous time and space translations, positive energy and a finite propagation cone | One canonical pair in the actual selected ordered spin ground; no assumed Gaussian state or ergodicity |
| [FA101](../shrinking_modular_entropy/README.md) | Actual finite-density entropy differences converge to the existing Weyl stress response with an explicit common shrinking schedule and invertible nine-frame response | `U=0`, full regional CAR algebras and observer-frame shell preparations on ordinary rest balls |
| [FA102](../native_joint_spin_spectroscopy/README.md) | Complete joint energy/symmetry instruments, all 32 ideal cube-bin enclosures, and an executed native two-site detector | 19,228 pulses with all 128 work rows and every failure retained; the cube native program is not executed |

The ledger contains **102 reviewed domains and 97 realized component attachments**. The 84 canonical premise rows, their discharge statuses, and integrated physical acceptance are unchanged. The proofs close the stated component domains; they do not close their stronger owning issues.

## What the field construction establishes

For the selected cubic spin ground, with actual order parameter `m_s>1/1000`, the scaled fields are

```
Q_p(F) = p^2 sum_x (-1)^x F(px) S_x^3,
P_p(G) = -p sum_x G(px) S_x^2.
```

Maximality of the selected order among ground states makes the macroscopic staggered observable a scalar, without an extra purity or ergodicity premise. Exact spin algebra and a uniform energy-form estimate then control every bounded multitime field history. A common subsequence gives the Weyl relation with symplectic coefficient `m_s`, strong regularity on the full history representation, continuous translations and nonnegative energy. The microscopic locality estimate yields a conservative continuum propagation cone.

This constructs a field algebra and its dynamics from the actual ordered ground. It does not yet determine the limiting dynamical equation, uniqueness, a sharp relativistic shell or surviving connected interaction. One canonical pair is proved; a second pair cannot simply be added because its cross-commutator requires another estimate. The order parameter has not been numerically calculated.

The charged-field calculation complements this construction: polynomial free propagation and a direct mixed-time CAR comparison prove that the canonical fields on the centered clock become free even for the power-growing class `U=o(c^(2/7))`, `c=v/a`. For example, `U=c^(1/4)` lies inside it. A stronger microscopic coupling outside that class is not itself a proof of continuum interaction. The protected spin limit is a different construction.

## Energies and particle identification

The charged continuum in the proved free class has energies `E=v|p|` and zero relativistic rest masses in all eight valleys. These are model predictions for that class, without a measured species assignment. For the protected spin field, no proved sharp shell yet exists from which to calculate a relativistic mass.

The [complete finite energy readout](ORDERED_SPIN_ENERGY_READOUT.md) retains all 256 certified levels of the selected eight-site Hamiltonian. The first six positive gaps, in source units with `J=1` and `B=1/4`, have enclosure midpoints

```
1.092284677687, 1.442323754769, 2.166719725510,
2.543094000283, 2.665699560844, 2.729285256694.
```

Relative to the first gap these are approximately
`1 : 1.320465062114 : 1.983658445249 : 2.328233703387 : 2.440480595672 : 2.498694078978`.
They are finite selection-Hamiltonian energies, not relativistic rest masses. The joint detector adds certified probe-dependent energy/symmetry probabilities; it does not make every spectral level visible to every probe. Open-cube graph reflections are not a low-momentum dispersion curve. All rejection and failure probabilities remain in the instrument.

## How this advances the existing Einstein stack

The new entropy proof uses natural exterior-shell vectors and ordinary finite regional density matrices. Projection-block factorization controls the endpoint singularity of the modular logarithm. A fractional Sobolev boundary estimate gives an explicit positive mesh schedule for a joint shrinking-ball limit; its extremely small scales are represented symbolically rather than rounded to zero.

The resulting finite entropy derivatives and finite differences converge to the previously derived local Weyl stress response with error `o(R^4)` at fixed shell scale. Nine observer-frame preparations have a uniformly invertible response matrix. This supplies a quantitative matter/modular input to the existing geometry theorem stack. The frames are changes of preparation in the free continuum vacuum; a finite lattice boost symmetry is not asserted.

Area and edge identification, Newton normalization, generalized-entropy stationarity and responsive geometry remain separate premises. The entropy result is for `U=0`; applying it to the ordered-spin continuum field would require another derivation. No combination of these separate sectors is silently treated as one completed physical universe.

## Verification and remaining work

Each new package includes its analytic proof, exact or controlled finite checks, mutation controls, a source-bound receipt and a fresh-process cold checker. Independent reviews check the arguments as well as the finite examples. The [campaign](continuum_field_entropy_campaign.json) and [review](continuum_field_entropy_review.json) record those scopes. These checks are not a claim of full repository or remote CI success.

The next physical construction needs the spin dynamical equation and sharp spectral component, surviving continuum interaction, common physical charges and units, and the remaining geometric area/stationarity/response join. The [mass-shell diagnostic](../spin_shell_diagnostic/README.md) identifies the stronger ordered-ground current or force cancellation sufficient for a sharp cone; the present field-algebra result does not substitute for that estimate.
