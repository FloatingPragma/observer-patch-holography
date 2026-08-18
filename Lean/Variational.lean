import Variational.DiscreteEulerLagrange
import Variational.DiscreteNoether
import Variational.FiniteHistoryBridge
import Variational.FiniteRealTransfer
import Variational.LegendreBridge
import Variational.RealizedHistoryLegendreNoGo
import Variational.StationarySaddleCoverage
import Variational.SourceToHamiltonianComposed
import Variational.ModeExtremalEnrichment
import Variational.EnrichmentCharacterization
import Variational.TranslationInvariantComposedInstance
import Variational.MechanicsAdequacySurface
import Variational.PolynomialEnrichmentSelection

/-!
# Variational umbrella root

The discrete Euler-Lagrange equation for local path actions, the
discrete Noether theorem for invariant Lagrangians, and the certified
interface obstruction between finite-state Gibbs histories and universal
real single-site variations.
`Variational.LegendreBridge` supplies the finite Legendre
transform with involutivity and Fenchel-Young receipts, the equivalence
of discrete Euler-Lagrange transport with the discrete Hamilton step for
the strictly convex class, the two-faces theorem tying log-transition
minimizers to most probable paths of the realized chain, and the Noether
correspondence carrying the constant chain current to a conserved
quantity of the Hamilton flow.
`Variational.RealizedHistoryLegendreNoGo` proves the complementary boundary:
the binary source-history law fixes only the four endpoint values of a real
two-point Lagrangian.  Infinitely many strictly convex real enrichments agree
on every realized history while giving different Legendre transforms, so a
source or physical curvature receipt is required to select one.
`Variational.SourceToHamiltonianComposed` composes the two halves under
one typed antecedent bundle: from a supplied source law, the registered
embedding and enrichment, and the registered derivative data, one
extremal history threads the derived action, the exponential-tilt path
law, the interior mode clause, the discrete Hamilton equations, and,
with supplied symmetry data, a constant Noether segment momentum.
`Variational.ModeExtremalEnrichment` selects one representative inside a
declared one-parameter ansatz by the mode-extremality principle, with
curvature `2 * log (362055879 / 271780)`, and discharges the composed
hypothesis at the transition-weight global mode.  Its machine-checked
two-parameter counterfamily proves that the same stationarity and fixed-endpoint
single-site variation minimum do
not select a unique real enrichment beyond that ansatz.
`Variational.TranslationInvariantComposedInstance` discharges the
symmetry-extended bundle on a declared symmetric kernel, with the
constant Noether current exactly zero.
-/
