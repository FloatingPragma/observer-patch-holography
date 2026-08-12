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
The inverse-square shell law carries the Newtonian-limit row: under the
declared radial readout and shell-flux normalization the strength falls as
the inverse square, with the exponent supplied by the carrier dimension.
-/
