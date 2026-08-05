# B10 finite quantum limitations and matter corollaries

This note indexes the exact finite package. Every construction takes its
state, observable, partition, exterior carrier, or component-to-sector map as
an explicit input. The package supplies no source selection or laboratory
attachment.

## Supplied-state Robertson inequality

`EventAlgebra/Robertson.lean` proves:

- `supplied_pairing_robertson`: Cauchy--Schwarz for a supplied complex
  seminormed pairing;
- `finite_state_pairing_robertson`: the corresponding inequality for a
  supplied finite density matrix and two supplied Hermitian observables;
- `neg_I_mul_commutator_expectation_eq_readout`: the complex identity between
  the real readout and `-i` times the ordinary commutator expectation;
- `finite_state_robertson_commutator`: the conventional bound using the
  squared complex norm of that commutator expectation;
- `stateVariance_nonnegative` and the two abstract zero-variance lemmas;
- `pauli_xy_noncommuting_control`: an exact noncommuting saturation control;
- `pauliZ_pauliX_ne_pauliX_pauliZ`: the second Pauli pair is explicitly
  noncommuting;
- `pauli_z_zero_variance_control`: its exact zero-variance and zero-readout
  values.

The theorem does not infer a state, observable, Hamiltonian, uncertainty
scale, or instrument.

## Partition-relative superselection

`EventAlgebra/Superselection.lean` proves:

- `partitionOperationallyEquivalent_iff_pinching_eq`: two matrices give the
  same trace statistic against every matrix in the partition commutant exactly
  when their partition pinchings agree;
- `trace_mul_eq_zero_of_partitionOffDiagonal`: the kernel of pinching is
  invisible to the whole sector-preserving commutant;
- `partitionPinching_partitionCorner_eq_zero`: every cross-sector corner is
  in that kernel;
- the corresponding partition-average and projector-span corollaries.

This is superselection relative to a supplied projective partition and its
declared readout. It does not construct an edge center or identify a physical
charge sector.

## Exterior basis and algebraic exclusion

`Screen/ExteriorComponentBridge.lean` uses Mathlib's
`Module.Basis.ExteriorAlgebra` on the supplied five-coordinate carrier. It
proves:

- `exterior_basis_label_count`: the exterior basis has 32 subset labels;
- `basis_bidegree_exhaustive` and `bidegree_count_table`: those labels split
  into the twelve color/weak bidegrees with exact multiplicities;
- `componentDegree_exact_nontrivial_menu`: removing the vacuum and top line
  leaves exactly the ten component bidegrees in `ExteriorSelection.lean`;
- `component_dimension_binding`, `component_charge_binding`,
  `component_parity_binding`, and `component_conjugation_binding`: every
  frozen table column is bound to the exterior bidegree construction;
- `creation_square_zero`, `creation_generators_anticommute`, and
  `creation_actions_anticommute`: exact algebraic exclusion and
  anticommutation.

These statements concern the declared exterior module. A physical
spin--statistics theorem, continuum locality, and matter attachment are not
consequences.

## Typed finite composition

`Screen/QuantumMatterIntegration.lean` defines
`FiniteQuantumMatterBridge`. Its fields explicitly identify:

1. a projective partition and a declared block readout equal to pinching;
2. each exterior component row with a nonzero partition sector;
3. each mapped sector with its exterior-derived central weight;
4. a separate nonempty chiral anomaly-free selection mask on the same ten
   rows.

The transported conclusions are:

- `declaredBlockReadout_eq_iff_operationallyEquivalent`;
- `mapped_cross_component_corner_invisible`;
- `even_component_weights_eq_matterWeights`;
- `kernel_on_exterior_component_weights` and
  `kernel_on_mapped_component_weights`;
- `bridge_selection_is_parity_sector`.

The block readout lands in the generally noncommutative partition commutant.
It is not the commutative public-record readout. Partition averaging supplies
the latter on the projector span, and `Superselection.lean` proves that both
readouts erase cross-sector corners. The weight-kernel theorem is arithmetic
on supplied labels; it does not define a group action on the projectors. The
selection mask is constrained on the component rows but is not attached to a
source-selected physical sector action.

`coordinateBridgeControl` is an exact synthetic witness that the interface is
inhabited. `coordinate_nonzero_offDiagonal_control` supplies an explicit
nonzero element of the pinching kernel, while
`coordinate_diagonal_not_partitionOffDiagonal` checks that a sector-diagonal
datum is retained rather than erased. These are not source evidence. The
weight mutation
`fractional_singlet_mutation_collapses_component_kernel` reduces the common
central kernel to the identity and checks that the kernel result is sensitive
to an inadmissible extra weight.

## Trust and physical boundary

The declarations compile with standard Lean/Mathlib axioms only:
`propext`, `Classical.choice`, and `Quot.sound`. Finite enumerations use kernel
evaluation and no `native_decide`.

The bounded package supports a conditional finite structural claim. It does
not supply a source-selected state, observable, sector partition, matter
action, parity representative, global gauge-group form, continuum Spin
attachment, laboratory charge, numerical prediction, or postdiction.
