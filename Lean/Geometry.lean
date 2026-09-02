import Geometry.CanonicalLorentzModule
import Geometry.CelestialNullCone
import Geometry.ObserverFrameHyperboloid
import Geometry.ObserverRestSpace
import Geometry.EinsteinTensorBridge
import Geometry.LorentzOverlapCocycle
import Geometry.EventGermDisplacement
import Geometry.CelestialSoldering
import Geometry.EventFrameSoldering
import Geometry.SpatialReadbackSoldering
import Geometry.InverseSquareShellLaw
import Geometry.MassShellKinematics
import Geometry.InternalClockRestFrequency
import Geometry.FreeEvolutionPersistence
import Geometry.CausalOrderComposition
import Geometry.EventPopulationChartInterface
import Geometry.SourceDerivedSpacetimeCarrier
import Geometry.SourceOrderFrameCompatibilityPacket
import Geometry.SourceOrderEinsteinComposition
import Geometry.CommonWorldKinematicsWitness
import Geometry.RestFiberShellTransport
import Geometry.CommonWorldIslandBridge
import Geometry.EinsteinBranchInhabitant
import Geometry.CommonWorldInstrumentJoin
import Geometry.CommonWorldMaxwellClockJoin
import Geometry.IntegerKCombInvariance
import Geometry.CommonWorldJointAction
import Geometry.ScreenCarrierMapCandidate
import Geometry.CarrierDynamicsCompatibility
import Geometry.ChargeFixedInteraction
import Geometry.InternalEnergyInertia
import Geometry.CompositeMomentumCovariance
import Geometry.PortChargeMinimalCoupling
import Geometry.ProperTimeInternalAction
import Geometry.WorldlineHopTransport
import Geometry.TransportedChargeForceLaw
import Geometry.SeamStepSpeedLimit
import Geometry.SourceClockRateAlongWorldlines
import Geometry.ProperLengthClockedChain
import Geometry.TimelikeClassForceLaw

/-!
# OPH construction geometry umbrella

The C1 modules construct the intrinsic four-dimensional Hermitian Lorentz
module, its future-null-ray celestial sphere, the future-unit-timelike frame
hyperboloid, three-dimensional positive-definite rest spaces, and the exact
sign-convention bridge into the existing Einstein tensor coordinates.  The C2
modules add the bounded algebraic soldering contract: coincidence-invariant
readbacks descend to event classes; separately, one affine Lorentz cocycle and
a supplied chart-coordinate family satisfying its overlap law induce
displacement, interval, celestial, frame, and local-rest-space covariance; and
the source Gram quotient is linearly and isometrically identified with the
standard internal rest fiber as a candidate readback. The quotient-descent
result does not construct or identify the atlas coordinate family. The
contract does not construct its source event atlas, population, open charts,
physical cone,
causal reachability, refinement tower, rods, clocks, or physical spacetime.
`Geometry.SourceDerivedSpacetimeCarrier` now constructs the finite algebraic
carrier `ℝ × FrameQuotient` from an independent real axis and the exact
rank-three source quotient, proves dimension four and one-positive/
three-negative Lorentz signature, and maps source-unit directions to
future-null vectors. On every finite authenticated log, canonical source
height is the attained maximum authenticated-parent-chain length. The
separate finite-relation compiler realizes every supplied finite decidable
strict partial order exactly in an abstract semantic log; for a transitive
input it authenticates every strict-order pair as a parent edge rather than
reducing the relation to Hasse covers, and it does not construct a threaded
execution. Positive-scale source height enters only in an event placement.
Every finite log has an
explicit enumeration-dependent injective forward-causal placement, but that
placement is not order-reflecting or source-selected. Given a general supplied
spatial readback, its explicit edge-speed bound proves every generated
ancestry relation future-causal. Equal-height spatial separation together
with strict spacelikeness of every increasing-height nonancestor derives the
converse cone-support clause; exact two-way agreement then derives event
separation and a one-chart exact finite order embedding. The physical spatial readback,
geometric conditions, and continuum family remain inputs, and the result is
not a topology, open atlas, physical-causality theorem, count--volume law, or
refinement/manifold limit.
The finite `SourceOrderFrameCompatibilityPacket` keeps those boundaries
explicit. An order-faithful placement constructs its one-chart order/cone attachment;
event-dependent Lorentz-frame transports remain optional supplied data, while
the standard-frame specialization removes that field for a purely algebraic
existence result. The Boolean-diamond control constructs a complete non-chain
packet. The independent rank-three quotient/rest-fiber metric bridge and its
exact algebraic celestial/future-null-ray equivalences are reused. None of
these results identifies a physical sky or signal relation or proves physical
causality, a volume law, manifoldlikeness, or a continuum limit.
`Geometry.SourceOrderEinsteinComposition` then takes the fixed algebraic
inverse images, in the source unit sphere, of the Einstein branch's
pre-existing nine-vector coordinate tomography frame.  The distinguished
nine-set is neither selected by the generated poset/source dynamics nor
invariant, as a set, under Lorentz or `SO(3)` transformations; only the
`(1,n)` representative of each already chosen direction is canonical.
Balance on those nine algebraic source directions, together with explicit
symmetry, Ward, Bianchi, and connectivity on the same finite event carrier,
derives the all-null balance and one constant metric ambiguity. This minimal
Einstein-shape theorem omits six normalization or source-unification fields.
The physically normalized Einstein-form equation additionally consumes the
vacuum reference and Newton-scale identification. The full packet also carries
a universal-source equality for later composition, but that equality is not
used in the tensor equation itself. The nine directions are not observed
provenance links or a physical sky, and this finite conditional theorem
supplies no count--volume law, manifoldlikeness, smooth limit, or physical
spacetime.  Its `V 3`, `Mat 3`, `eta 3`, and `Fin 4` tensor/step types remain
the pre-existing Einstein algebra: the new source-carrier theorem proves an
exact `1+3` algebraic carrier but does not yet identify its finite differences
with those tensor fields.  The generated order is not used to select the
tensor fields, balance law, or nine algebraic directions.
The inverse-square shell law carries the Newtonian-limit row: under the
declared radial readout and shell-flux normalization the strength falls as
the inverse square, with the exponent supplied by the carrier dimension.
`Geometry.EinsteinBranchInhabitant` inhabits the Einstein-branch register
surface for every declared repair law and horizon record with a synthetic
demo-tier value whose modular flow supply is the thermodynamic first-law
datum on the simplex tangent space; the Einstein clause holds at the value by
definition of the chart geometry and selects nothing.
`Geometry.CommonWorldInstrumentJoin` extends the bridged common-world record
by the certified step-scaled evolution on the shared screen carrier, the
declared Lüders instrument with its typed pinching transport to the committed
diagonal partition, a declared step-to-clock dictionary, and the first-law
join of the run-state diagonal; the missing joins are named exactly and the
common-world row stays owed.
-/
