# Assumptions and source-specific theorem use

The [inventory](assumption_closure.json) classifies the assumptions consumed
by a fixed, named dependency set. Its canonical premise statuses are copied
without promotion. A declared architecture is a legitimate construction
choice. An implementation claim about a narrower source requires an
operation and record comparison; its declaration alone supplies no such
comparison. Empirical inputs retain their calibration ancestry.

## Audited set

The manifest fixes revision `44d42ff4e6f754c67fd6e67d1b3156d76ff3bf34` for
all 84 canonical premise rows, the conditional inputs of the joint field
and particle dependency graph, and 16 named claims. Every named assumption, imported-result edge
and each consumed, open or boundary premise edge in those claim rows is
represented. The missing identifier between PR-82 and PR-84 is preserved.
The theorem-premise declarations are consumed by reference to
[`claims/axiom_registry.yaml`](../../claims/axiom_registry.yaml).

The native execution headline has 48 named assumptions and declares no
direct canonical-premise consumption. This inventory preserves both facts.
Sharing terminology with a canonical premise does not create a dependency.
Conversely, the gravity graph node's empty conditional-input array does not
remove the gravity claims' explicit modular, stress and geometry inputs.
Canonical rows absent from the selected claims are marked as non-consumed
within this audited set, with their standing owners and statuses retained.

Named imported results retain an explicit applicability obligation unless
a separately reviewed finite application matches their hypotheses. Empirical
comparison imports retain their data and calibration role. A theorem name
alone provides no source-specific discharge.

The selected claims' wider evidence files are unreviewed pointers in this
inventory. Eight finite applications have separately reviewed mathematical
witnesses, code and independent tests in the source manifest. Their owning
packages retain transitive source-custody and mathematical replay checks.
The inventory checker verifies its direct manifest; it does not replay
all owning packages or certify every claim in the repository. Coverage of
this declared dependency set is distinct from truth or completeness of the
full theory.

The three additional batch admissions bind the lepton comparison, the
linear electroweak covariance budget and the species-identifiability
calculation to their producer, proof and independent checker. They inherit
explicit service, calibration and physical-map obligations. Live claim
prose, claim evidence additions and completion-program pointers are excluded
from live hash dependencies, avoiding a claim/inventory citation cycle.
The baseline input graph and canonical contract are immutable within this
inventory; changing that graph requires an explicit review and revision.

## Coverage and interface proposition

Let `P` be the baseline canonical identifier list, `C` the selected claims,
and `J` the joint graph. Define the canonical edge set

```
E = {(p,c,r) : p occurs in c.premise_dependencies[r]},
r in {consumed, open, boundary}.
```

Let `A={(c,a): a occurs in c.assumptions}`,
`I={(c,i): i occurs in c.imported_results}` and
`G={(a,n): a occurs in n.conditional_assumptions}`. The validator enforces
ordered equality of the canonical row list with `P`, exact equality of all
canonical contracts/statuses, equality of the recorded edge set with `E`,
and bijective coverage of `A`, `I` and `G`. It separately checks the reviewed
service set and source-specific applications. Missing rows, duplicated rows
and unclassified inputs cannot satisfy these equalities.

Every admitted source-specific application has a tuple

```
D = (state, algebra, operation, reference, clock, domain).
```

The validator requires all six nonempty coordinates and equality to the
reviewed source interface. It also requires a mathematical proof anchor,
a distinct independent checker, explicit finite hypotheses, the reviewed
hypothesis match and hashed source/proof/checker contents. Thus every
accepted source-specific application names all six coordinates and its
actual premises. This is a structural proposition about the admission
contract. It does not prove a mathematical statement merely because a
path exists or a hash matches.

The imported conjugation theorem is a concrete example. The finite state
and involution, native branch equivariance, even preparation and odd
observable are specified independently. The native checker exhausts the
local and closed-cell branches; the generic Lean theorem propagates
finite-kernel evenness. This combination is an executable hypothesis check
and a formal finite theorem. It is not a Lean extraction of the Python
program, nor a physical baryon-number identification.

## Finite premises and physical transfers

`hypothesis_inputs` contains only the inputs needed by the finite result.
`downstream_obligations` contains additional conditions for wider uses.
For example, the closed-cell absorption theorem does not require a physical
species decoder or a federation transfer. The fixed charged cube does not
require the primary program's physical vacuum. Those transfers are separate
obligations with an owner, a smallest witness and a counterexample test.

The finite lepton theorem retains every outcome and arbitrary untouched
reference. Its weak generator has spectrum `{-1,0,1}`, diameter two, and
its half-diamond bound is `14 t^2/(81 N)` for dimensionless time `t`.
The mass-round bound is `349 t^2/(1296 N)`. Freshness means joint
independence from the continuing system, reference and archive; correct
source marginals do not suffice. A fixed-analyzer trace-distance bound
controls effect eigenvalues and conditional effect-rate ratios. It does
not establish a physical pole-mass error.

The population law concerns a closed finite cell and a nonincreasing
load-plus-register cost. The unresolved probability is at most
`25^(-floor(n/32))` for every initial distribution. The limiting occupied
share `1007/13824` instead consumes independent uniform loads `0..5`
and zero initial flux. This class occupancy is not a thermal species count.
In the declared fermion-only dictionaries the class census identifies the
quark and total lepton aggregates. With bosons admitted it does not
separate the lepton, Higgs and gauge shares. Species-specific and physical
counting maps therefore retain distinct obligations.

The electroweak variance bounds are sharp for the supplied marginal errors
and fixed linear coefficients. Independence of separate run sets is an
explicit additional assumption. The unrestricted positive-semidefinite
covariance range removes that assumption; neither range bounds nonlinear
remainder or supplies confidence coverage. The ideal-gas benchmark is also
external: `T=1 TeV`, zero chemical potentials, massless states, three active
neutrino families, no populated right-handed neutrinos, 90 fermionic and
28 bosonic helicities. Number and energy density conventions are distinct.
No native thermal history follows from that declaration.

All listed source joins have `join_established=false`. A declaration of
physical interpretation closes disclosure only; its implementation or
physical-map obligation is retained. No canonical premise status or
integrated gate acceptance changes, and no independent physical comparison
is added by this inventory.

## Reproduction

```sh
python3 -m model.completion.assumption_closure
python3 -m pytest -q model/completion/test_assumption_closure.py
python3 tools/check_claim_registry.py
```

`--write` records the reviewed classification policy and source bindings.
It is a maintainer action after review, not a proof-search operation.
The tests remove each input family, erase interface coordinates, substitute
receipts for proofs, promote calibrations, invent source joins and mix
physical transfers into finite premises. They also independently enumerate
the canonical, selected-claim and joint-graph coverage sets.
