import Variational.DiscreteEulerLagrange
import Variational.DiscreteNoether
import Variational.FiniteHistoryBridge
import Variational.FiniteRealTransfer
import Variational.LegendreBridge

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
-/
