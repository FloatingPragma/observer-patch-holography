# Source-law and operational selection controls

This package checks finite witnesses for
`paper/tex_fragments/PAULI_SOURCE_SELECTION_BOUNDARY.tex`. It tests whether the
supplied source hopping, positive metric, finite capacity, covariance
and exact number-phase law select exclusion. They do not suffice: a symmetric
four-mode Fock space capped at total occupation two has all these mathematical
properties and retains double occupation.

From the research repository root:

```sh
python3 code/pauli_source_selection/build.py --write
python3 code/pauli_source_selection/verify.py
python3 -m pytest -q code/pauli_source_selection/test_pauli_source_selection.py
```

The producer uses occupation polynomials and factorial Gram weights. The
independent verifier embeds the 15 states into a labeled tensor space with
21 coordinates, checking adjoints, quadratic source dynamics and its Cayley
lift. Four normalized mode controls retain square-creation norm two. The
analytic proof establishes all-mode covariance and phase identities; these
finite controls are not an exhaustive test of infinitely many modes.
The exact source edge has energy spectrum `-2,-1,0,1,2` with multiplicities
`3,2,5,2,3`. Its one-particle restriction is the supplied current action.
The global occupation cap is a new declared model and does not reproduce the
parent's fifteen-particle preparation, quantum Gauss or local relativistic fields.

A constructive conditional route avoids the scalar-deformation ansatz: on a
finite Hilbert space, a creator that is both a contraction and an exact unit
phase raiser is square-zero. Applying this to linear mode superpositions gives
creation anticommutation. Mixed CAR additionally follows if creator and adjoint
are the complete pair of Kraus outcomes. The package checks a finite CAR
realization and a contrasting three-outcome capped-boson instrument. Scaling
the latter's branches by 1/2 makes them contractions but changes their phase
coefficient to 1/4, so the same operator does not meet both premises.
Under the exact phase law, contraction is equivalent to exclusion; identifying
the source creator as that contraction is the substantial missing condition.

Basis-only checks are insufficient: commuting hard-core bit raisers and signed
CAR raisers have the same individual occupation probabilities. Their normalized
`3/5,4/5` superposition distinguishes them; the commuting model's double-created
vacuum amplitude is `24/25`. The signed model passes exact real and imaginary
pair acceptance checks. These are operator identities, not conclusions inferred
from a finite sample of detector outcomes.

The source-history control reads the eight committed repair-count actions.
Their eigenspace dimensions are `2,4,2`. Treating that particular count operator
as `a†a` with exact unit phase would require an isometry from dimension four to
dimension two, which is impossible. Enlarged carriers and different number
operators are not excluded.

The next OPH attachment must construct the same physical creation map as a
contractive operation and a unit phase raiser on bounded self-reading patches,
with local state, ports, readback, records, feedback or repair and a public
evidence bundle. An arbitrary measurement branch, an independently supplied
field or classical occupation records do not establish that identification.
This package does not derive those premises or physical spin--statistics,
and it is not a countermodel to the full observer axioms.
