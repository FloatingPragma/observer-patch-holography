# B5 Finite Continuity and the Ward-Premise Boundary

## Machine-checked finite layer

`Screen/RegionalContinuity.lean` fixes the sign convention of the exact
thirty-seam incidence map. The existing boundary is net inflow, with `-n` at
the smaller endpoint and `+n` at the larger endpoint. For

```text
qNext - q = source + seamBoundary current
```

the theorem `regional_continuity` gives

```text
regional load change = regional source - signed outward flux.
```

Internal seams cancel. On the complete twelve-port graph, total-load change
equals the total source, and a zero-total source preserves total stored load.

`Screen/DiscreteGauss.lean` packages the exact rational incidence results:

- a port load admits a seam-flux solution exactly when its total load is zero;
- the declared port-zero spanning-tree section constructs one explicit
  solution; no naturality or canonical-choice claim is made;
- every other solution differs by a cycle current;
- the cycle kernel has rational dimension nineteen.

Its displayed Gauss equation uses the repository boundary convention
`net inflow = charge`. A physical convention written with outward divergence
would insert the corresponding minus sign. No electromagnetic identification
is made in either convention.

`Dynamics/ProtectedCharge.lean` proves the reusable channel criterion. For a
declared Schrödinger/Heisenberg dual pair, a charge expectation is preserved
on every finite signed state exactly when the dual channel fixes the charge
observable.

The following notions are deliberately separate.

- **Antisymmetric transfer** is an orientation convention for a current:
  reversing a seam reverses its signed value. It gives internal incidence
  cancellation.
- **Action symmetry** says weights or an operator are invariant under a
  declared relabelling action.
- **Channel covariance** says a state channel commutes with that action.
- **Conservation** says the dual channel fixes a specified charge observable.

None of the first three conditions implies the fourth. The theorem
`channel_covariance_does_not_imply_charge_conservation` is an exact two-state
negative control: uniform averaging commutes with the nontrivial swap while
destroying the nonzero swap-odd charge. The identity-channel theorem supplies
a nonzero positive control. Thus B5 neither obtains conservation from a
symmetry slogan nor hides the charge in a zero-observable case.

All audited B5 declarations depend only on `propext`, `Classical.choice`, and
`Quot.sound`. They contain no `sorry`, `admit`, project-level `axiom`, or
`native_decide` proof.

## Einstein-branch boundary

B5 supplies the finite incidence precursor of a Ward identity.
`OPH.EinsteinBranch.ContinuumEinsteinPremises.ward` remains an explicit input.
The finite theorem does not construct a continuum stress tensor or identify
the seam load and current with physical stress-energy components.

Discharging the Ward field requires the following typed attachments and
limit statements:

1. a source-derived finite datum at every regulator;
2. refinement-natural transport and a positive regulator scale tending to
   zero along one injective monotone cofinal sequence;
3. a local-chart realization into the temporal and spatial components of one
   stress tensor;
4. an admissible separating test-field class and weak convergence of the
   stress pairing on that class;
5. a finite residual defined from the actual continuity witness, its exact
   vanishing, and convergence to a continuum residual operator that is not
   globally zero;
6. identification of that continuum residual with physical distributional
   divergence;
7. covariance of the realization under chart changes;
8. an inhabited common-domain tower carrying the same source through the
   finite and continuum layers.

`Dynamics/WardLimitManifest.lean` records the mathematical fields and guards
in this list. The physical distributional-divergence identification in item 6
and the concrete meanings of the abstract relations remain F1 obligations.
Every regulator carries an actual `FiniteContinuityWitness`, so the manifest cannot
replace B5 with a detached source label. It also requires a not-all-zero
finite datum, a positive scale tending to zero, admissible tests that separate
stress tensors, a nonzero finite stress response, and a residual that is not
globally zero. The declared residual equality on every admissible test at the
candidate limit is a theorem: exact
finite-residual vanishing and residual convergence force it by uniqueness of
real limits. The candidate limit may itself be zero. Source derivation,
transport compatibility, the meaning of the finite and continuum residuals,
chart realization, and common-tower carriage remain relations whose concrete
physical meanings and evidence must be supplied and audited by F1. Defining
or generically inhabiting the structure does not establish that its residual
is physical distributional divergence or that `AdmissibleTest` is the intended
compactly supported/domain test class.

`residualConverges` is the decisive analytic premise. F1 must prove it from an
independent finite-to-continuum compatibility estimate; deriving it from the
desired residual equality would be circular and does not discharge the
manifest. The physical Ward reading additionally requires the separate
distributional-divergence identification.

F1 may consume B5 only through a receipt that supplies these fields. Reusing
the finite theorem by name without the realization and limit receipts does not
discharge the Ward premise.
