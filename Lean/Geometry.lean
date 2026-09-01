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
The finite `SourceOrderFrameCompatibilityPacket` keeps those boundaries explicit:
it composes authenticated informational order with a supplied cone/chart
attachment and supplied per-event Lorentz-frame transports, while reusing the
independent rank-three source quotient/rest-fiber metric bridge. The unit
Gram directions of that quotient are exactly equivalent to the algebraic
celestial two-sphere and future-null rays. This is not a physical sky or
signal identification; no inhabitant of the packet is constructed, and it
proves no physical causality, volume
law, manifoldlikeness, or continuum limit.
`Geometry.SourceOrderEinsteinComposition` then takes the fixed algebraic
inverse images, in the source unit sphere, of the Einstein branch's
pre-existing nine-vector coordinate tomography frame.  The distinguished
nine-set is neither selected by the generated poset/source dynamics nor
invariant, as a set, under Lorentz or `SO(3)` transformations; only the
`(1,n)` representative of each already chosen direction is canonical.
Balance on those nine algebraic source directions, together with explicit
symmetry, Ward, Bianchi, connectivity, vacuum-reference, and scale premises
on the same finite event carrier, derives the all-null balance and the
Einstein-form tensor equation.  The universal-source equality is carried in
the premise and returned conclusion for later composition, but it is not used
to derive that displayed equation.  The nine directions are not observed
provenance links or a physical sky, and this finite conditional theorem
supplies no count--volume law, manifoldlikeness, smooth limit, or physical
spacetime.  Its `V 3`, `Mat 3`, `eta 3`, and `Fin 4` coordinate/step types are
the pre-existing `3+1` Einstein algebra, not a dimension derived from the
poset.  The generated order is not used to select the tensor fields, balance
law, or nine algebraic directions.
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
