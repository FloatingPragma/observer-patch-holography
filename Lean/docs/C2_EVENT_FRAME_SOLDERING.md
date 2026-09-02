# C2 algebraic event-frame soldering contract

The C2 Lean packet closes the finite algebraic contract and nothing stronger.
It does not construct a physical event manifold.

## Modules and exact results

| Module | Result |
|---|---|
| `Geometry/LorentzOverlapCocycle.lean` | A time-oriented Lorentz linear equivalence on `Herm2`; affine chart transitions; identity, reverse-transition, and triple-overlap cocycle theorems; exact cancellation of affine translations in coordinate differences. |
| `Geometry/EventGermDisplacement.lean` | An actual coincidence setoid is required. A raw readback descends uniquely iff it is coincidence-invariant. Separately, an `EventGermAtlas` supplies a chart-coordinate family and its affine overlap law; only for that supplied family do displacements reverse, chain, transform covariantly, and have chart-independent Lorentz interval. A three-germ control proves that reflexive symmetric pairwise overlap need not be transitive. |
| `Geometry/CelestialSoldering.lean` | Oriented Lorentz maps descend to future-null rays and the C1 celestial sphere, respect identity and composition, and carry a future-null event displacement covariantly across overlaps. Nullness remains an input. |
| `Geometry/EventFrameSoldering.lean` | One base-chart frame selection and one Lorentz cocycle construct all chart-frame readbacks and their exact algebraic consequences. A finite two-event/two-chart control has a nonzero affine transition and a nonzero future-null displacement. |
| `Geometry/SpatialReadbackSoldering.lean` | Canonical rest projection, temporal-plus-rest decomposition, Lorentz-natural rest-space equivalence, exact rest-metric preservation, overlap-covariant local spatial readback, and an exact isometric linear bridge from the source `FrameQuotient` to the standard internal rest fiber. |
| `Geometry/EventPopulationChartInterface.lean` | `SourceDerivedCausalChartInterface` uses semantic events as the population and generated ancestry as the order, with no separate population map or free precedence relation. The finite atlas, total visibility, and exact order/cone attachment remain explicit data; base-chart separation follows from order/cone exactness and source-order antisymmetry and then propagates through the overlap cocycle. |
| `Geometry/SourceDerivedSpacetimeCarrier.lean` | Constructs `SourceSpacetimeCarrier := ℝ × FrameQuotient` from an independent real axis and the exact rank-three source quotient; proves real dimension four, exact equivalence to `Herm2`, and a displayed `(+---)` Gram matrix; maps every source-unit direction to `(1,q)`, which is future-null. Positive `timeScale` times canonical source height gives only the temporal coordinate of a supplied event placement in that independently defined target, not a physical event time. A supplied spatial readback and authenticated-edge speed bound imply generated ancestry is future-causal; an explicit cone-support converse gives an exact two-way finite order embedding, preserves every source interval exactly in the placed cone, and constructs a one-chart source-native interface. The Boolean response diamond is a fully checked non-chain control: all parent edges are null and its two responses are spacelike. |
| `Geometry/SourceOrderFrameCompatibilityPacket.lean` | Combines the source-native chart interface with the quotient/rest bridge and optional event-frame transports. An order-faithful placement constructs the chart packet; choosing the standard-frame gauge removes the frame-transport existence field without claiming physical event-dependent frames. The explicit Boolean diamond inhabits the full finite packet and its advertised consequences. Exact equivalences `SourceUnitDirection ≃ CelestialSphere` and `SourceUnitDirection ≃ FutureNullRay` and the transported rest-metric identity remain algebraic. |
| `Geometry/SourceOrderEinsteinComposition.lean` | Defines the exact inverse images of the existing coordinate-fixed nine-vector null-tomography frame in `SourceUnitDirection`; the resulting distinguished family is neither source-selected nor invariant as a set. `SourceDirectionEinsteinShapePremises` removes six normalization or unification fields from the Einstein-shape argument. Nine supplied balances plus supplied symmetric tensors, Ward/Bianchi identities, and connectivity imply all-null balance and one constant metric ambiguity. The physically normalized equation still requires the separate vacuum and Newton-scale data; the full conclusion also carries a source-identification equality that the equation does not use. The generated order shares only the finite event type with the separately supplied tensor interface; it does not construct the tensor fields, identify coordinate differences, or enter the tensor calculation. |

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
- the exact algebraic carrier `ℝ × FrameQuotient`, of dimension `1+3`, with
  one positive and three negative displayed axes;
- a canonical future-null vector `(1,q)` for every source-unit spatial
  direction;
- a forward order-to-cone theorem for the event placement whose temporal
  coordinate is scaled canonical source height, using a supplied spatial
  readback and one authenticated-edge squared-speed bound;
- an exact one-chart finite order embedding, event separation, and source-native
  chart interface under the explicit cone-support converse;
- a fully constructed non-chain Boolean-diamond packet whose parent edges are
  null and whose independent responses are spacelike;
- a conditional general composition theorem exposing authenticated
  informational order, optional per-event frames, and the rank-three
  rest-fiber bridge in one interface;
- an exact same-finite-event-type reduction, in a separately supplied `3+1`
  tensor interface, from nine supplied balances on the fixed algebraic inverse
  images of the pre-existing coordinate tomography frame to all-null balance
  and an Einstein-form tensor relation. The
  distinguished nine-set is neither source/poset-selected nor Lorentz- or
  `SO(3)`-invariant as a set. The `1+3` carrier is constructed independently;
  identifying its differences with the supplied `Mat 3`, `eta 3`, and
  `Fin 4` tensor/step data remains open.

## What remains external

The finite algebraic interfaces do not supply the source-causal continuum
certificate required for a physical open-manifold claim. One common cofinal
family must still provide:

1. a physical interpretation of retained events and authenticated links, with
   an adequacy argument relating informational support to physical signals;
2. source-selected placements with exact order/cone agreement, together with
   dense and isotropic source-`S²` link directions rather than a fitted cone;
3. refinement maps preserving events, generated order, directions, and
   placements, plus a specified convergence notion on the full cofinal family;
4. a locally calibrated Poisson/count--volume density law, separate from the
   finite exact-order theorem;
5. mutually consistent dimension, interval-profile, and manifoldlikeness
   tests, with stable topology or homology;
6. a distinguishing Lorentzian limit with open local four-dimensional charts,
   regularity, and uniqueness control;
7. an operational clock calibration kept separate from ordinal source height.

Einstein promotion then additionally needs either same-family tensor-curvature
reconstruction converging to the Einstein tensor or an independently premised
continuum small-ball/null-balance identification, as well as physical stress,
conservation, coupling, and controlled remainders. Scalar-curvature convergence
alone is only a diagnostic. A supplied `EventGermAtlas` remains a valid
multi-chart adapter for a continuum already obtained this way; it is not the
producer or first premise of the source-derived route.

At the finite semantic-event level, the source-native route needs no separate
population map, freely declared precedence, or pre-existing rank-four carrier.
An order-faithful placement constructs separation,
one global chart, and exact order/cone agreement. Generic placements still
require a spatial readback and edge-speed certificate; exact order reflection still
requires the cone-support converse. None of these finite results supplies
open chart images, completion points, a topology, local finiteness,
count--volume calibration, manifoldlikeness, refinement convergence, a
physical signal interpretation, or an operational clock. The standalone
quotient-descent theorem also does not identify its descended readback with a
nontrivial multi-chart `EventGermAtlas.coordinate` family.

The `FrameQuotient` equivalence is deliberately targeted at
`RestSpace standardFrame`. It does not select that frame physically and cannot
be read as absolute space, a global slice, a rod, or a position observable.
Likewise, the unit-direction equivalences identify algebraic types, not a
physical sky or signal distribution. The explicit Boolean diamond now
inhabits `SourceOrderFrameCompatibilityPacket` in a standard-frame gauge, but
that gauge does not derive physical event-dependent frames and the four-event
control carries no count-density law. The generic exact-order-embedding theorem
is conditional on its spatial and converse certificates. It therefore proves
no dense or isotropic link coverage, physical causal faithfulness, volume
calibration, manifoldlikeness, refinement convergence, or continuum limit.

The same-event-type Einstein-form theorem, in a separately supplied tensor
interface, does not strengthen the order into a
physical causal relation. Its nine directions come from the independent
source Gram quotient, not from provenance links. No link selects a direction,
tensor field, step, or balance, and no `SourceDirectionEinsteinShapePremises`
or `SourceIndexedEinsteinPremises` inhabitant is constructed. The theorem's
geometry and stress fields are supplied finite functions, and their coordinate
differences are not identified with those of the constructed carrier. Reading
them as smooth curvature and physical stress-energy still requires the
source-causal continuum certificate: source-selected refinement,
manifoldlikeness, count--volume, small-ball, curvature, and stress-identification
receipts. The inherited adapter name
`toContinuumEinsteinPremises` is a type-interface name, not a continuum proof.
