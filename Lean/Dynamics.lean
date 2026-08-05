import Dynamics.ProtectedCharge
import Dynamics.WardLimitManifest
import Dynamics.ConditionalExpectationGenerator
import Dynamics.PublicMarkov
import Dynamics.PrivateInner

/-!
# OPH construction dynamics umbrella

Cross-track dynamical interfaces live under `OPH.Dynamics`. Concrete finite
realizations remain in their issue-owned libraries and import these interfaces
when needed.

`ConditionalExpectationGenerator.lean` contains the finite matrix identity
for the generator of partition pinching, its projector-rate GKSL expression,
the single-collar fixed algebra, and the multi-collar intersection under an
explicit no-cancellation premise. It supplies no formal channel object or
physical rate attachment.

`PublicMarkov.lean` and `PrivateInner.lean` carry the bounded B3 dichotomy:
positive unital maps of finite record functions are row-stochastic, continuous
label-permutation flows are trivial, full private matrix-block automorphisms
are unitarily inner, and a fixed self-adjoint Hamiltonian generates a unitary
real-parameter von Neumann flow. Classification of arbitrary public star
automorphisms as label permutations, arbitrary central-block algebras, and the
converse continuous-group-to-Hamiltonian theorem remain outside this result.
-/
