# C2 algebraic event-frame soldering contract

The C2 Lean packet closes the finite algebraic contract and nothing stronger.
It does not construct a physical event manifold.

## Modules and exact results

| Module | Result |
|---|---|
| `Geometry/LorentzOverlapCocycle.lean` | A time-oriented Lorentz linear equivalence on `Herm2`; affine chart transitions; identity, reverse-transition, and triple-overlap cocycle theorems; exact cancellation of affine translations in coordinate differences. |
| `Geometry/EventGermDisplacement.lean` | An actual coincidence setoid is required. A raw readback descends uniquely iff it is coincidence-invariant. Separately, an `EventGermAtlas` supplies a chart-coordinate family and its affine overlap law; only for that supplied family do displacements reverse, chain, transform covariantly, and have chart-independent Lorentz interval. A three-germ control proves that reflexive symmetric pairwise overlap need not be transitive. |
| `Geometry/CelestialSoldering.lean` | Oriented Lorentz maps descend to future-null rays and the C1 celestial sphere, respect identity and composition, and carry a future-null event displacement covariantly across overlaps. Nullness remains an input. |
| `Geometry/EventFrameSoldering.lean` | One base-chart frame selection and one Lorentz cocycle construct all chart-frame readbacks. The premise-reduction theorem adjoins every exact algebraic consequence while retaining eight named source/physical receipts. A finite two-event/two-chart control has a nonzero affine transition and a nonzero future-null displacement. |
| `Geometry/SpatialReadbackSoldering.lean` | Canonical rest projection, temporal-plus-rest decomposition, Lorentz-natural rest-space equivalence, exact rest-metric preservation, overlap-covariant local spatial readback, and an exact isometric linear bridge from the source `FrameQuotient` to the standard internal rest fiber. |
| `Geometry/EventPopulationChartInterface.lean` | `SourceDerivedOrderEventPopulationChartInterface` fixes the legacy finite precedence field to the reflexive authenticated informational order while leaving the event population, atlas, cone attachment, and physical interpretation supplied. |
| `Geometry/SourceOrderFrameCompatibilityPacket.lean` | Conditionally combines that order interface with supplied per-event oriented Lorentz transports and the independent rank-three quotient/rest bridge. It proves exact equivalences `SourceUnitDirection ≃ CelestialSphere` and `SourceUnitDirection ≃ FutureNullRay`, plus the transported rest-metric identity. No inhabitant is constructed. |
| `Geometry/SourceOrderEinsteinComposition.lean` | Defines the exact inverse images of the existing coordinate-fixed nine-vector null-tomography frame in `SourceUnitDirection`; the resulting distinguished family is neither source-selected nor invariant as a set. `sourceOrderEinstein_from_source_directions` derives all-null balance and an Einstein-form tensor identity on the packet's finite `Event` carrier from nine supplied balance identities plus supplied symmetric tensor fields, Ward/Bianchi, connectivity, vacuum-reference, and scale data. The generated order shares the carrier but is not used in the tensor calculation. A universal-source equality is copied into the conclusion for downstream composition but is not used in that equation. |

## What the contract discharges

- coincidence-class readback descent once a genuine setoid and invariant raw readback are supplied;
- affine translation cancellation and Lorentz covariance of displacement for a separately supplied chart-coordinate family satisfying `coordinate_overlap`;
- chart independence of the intrinsic Lorentz interval;
- identity, inverse, and triple-overlap consistency;
- induced future-null-ray and celestial transport;
- frame transport from one base selection;
- canonical three-dimensional local rest projection and metric-preserving transport;
- an internal candidate bridge from the rank-three source Gram quotient to the standard-frame rest fiber;
- exact algebraic equivalences from source-unit Gram directions to the C1
  celestial two-sphere and future-null-ray labels;
- a conditional composition theorem exposing authenticated informational order,
  supplied future-cone compatibility, supplied per-event frames, and the
  rank-three rest-fiber bridge in one interface.
- an exact same-finite-carrier reduction from nine supplied balances on the
  fixed algebraic inverse images of the pre-existing coordinate tomography
  frame to all-null balance and an Einstein-form tensor relation. The
  distinguished nine-set is neither source/poset-selected nor Lorentz- or
  `SO(3)`-invariant as a set, and the theorem works inside supplied `Mat 3`,
  `eta 3`, and `Fin 4` types rather than deriving event dimension.

## What remains external

The premise-reduction theorem keeps all of the following explicit:

1. a source-produced event-germ atlas and coincidence relation, including a
   supplied chart-coordinate/readback family satisfying the affine overlap law;
2. population/realization of completion points;
3. certified separation;
4. rank-four open, locally bi-Lipschitz charts with an interior receipt;
5. attachment of the source cone to the intrinsic `Herm2` form with a positive margin;
6. refinement naturality on one source tower;
7. semantic causal reachability;
8. an operational clock.

The first and sixth receipts pass to the common-tower construction in E2.
Population, certified separation, open-chart, cone, and causal attachment
remain with F1 and the event/Einstein continuation. The clock remains D1. The
standalone quotient-descent theorem does not identify its descended readback
with `EventGermAtlas.coordinate`; that identification and the chart overlap law
must be supplied by a source construction. The realized screen receipt is only
`1+2` dimensional, while the existing four-dimensional reconstruction script
imports synthetic ground truth. Neither is a C2 source witness.

The `FrameQuotient` equivalence is deliberately targeted at
`RestSpace standardFrame`. It does not select that frame physically and cannot
be read as absolute space, a global slice, a rod, or a position observable.
Likewise, the unit-direction equivalences identify algebraic types, not a
physical sky or signal distribution. `SourceOrderFrameCompatibilityPacket`
has no constructed inhabitant: its atlas, cone equivalence, event population,
and frame transports are data fields, and it carries no count-density law.
It therefore proves no physical causal faithfulness, dense or isotropic link
coverage, event dimension, volume calibration, faithful embedding,
manifoldlikeness, refinement convergence, or continuum limit.

The same-carrier Einstein-form theorem does not strengthen the order into a
physical causal relation. Its nine directions come from the independent
source Gram quotient, not from provenance links. No link selects a direction,
tensor field, step, or balance, and no `SourceOrderFrameCompatibilityPacket`
or `SourceIndexedEinsteinPremises` inhabitant is constructed. The theorem's
geometry and stress fields are supplied finite functions. Reading them as
smooth curvature and physical stress-energy still requires external
event-manifold, refinement, count--volume, small-ball, and continuum
geometric-identification receipts. The inherited adapter name
`toContinuumEinsteinPremises` is a type-interface name, not a continuum proof.
