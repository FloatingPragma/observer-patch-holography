"""Conditional theorem navigation with checked source anchors and premise IDs.

This checker validates the declared interface and its references. It neither
recompiles the Lean proofs nor supplies physical instances of their premises.
Receipt generation and source hashing belong to the owning package runner.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REGISTER = 'tracking/premise_register.json'


def _reference(path, *anchors):
    return {'path': path, 'anchors': list(anchors)}


def _row(identifier, premises, references, conclusion, boundary):
    return {'id': identifier, 'canonical_premises': premises,
            'references': references, 'conclusion': conclusion, 'boundary': boundary}


def _interface(identifier, applies_to, premises, requirements, references=()):
    return {'id': identifier, 'applies_to': applies_to, 'canonical_premises': premises,
            'requirements': requirements, 'references': list(references),
            'physical_instance_supplied_by_contract': False}


def _specification():
    matter = 'model/common_source/native_matter_attachment.py'
    structure = 'Lean/Screen/SMStructureAdequacySurface.lean'
    chain = [
        _row('coarse_gauge_and_matter_selection', ['PR-09', 'PR-10', 'PR-11', 'PR-59'], [
            _reference('Lean/Screen/SMStructureComposition.lean',
                       'theorem composed_family_unique', 'theorem smStructure_two_coarse_masks'),
            _reference(structure, 'structure SMStructurePremiseData', 'theorem forced_gauge_type')],
            'The supplied weighting, compact-classification and exterior-mask premises force '
            'centre dimension1, simple dimensions3 and8, the coarse family G and either of '
            'two conjugate15-component masks.',
            'A coarse family label does not determine a continuous bracket parameter, '
            'physical current, global group or laboratory species.'),
        _row('same_mask_anomalies_and_family_band', ['PR-11', 'PR-36', 'PR-59'], [
            _reference(structure, 'theorem masked_generation_count',
                       'theorem masked_anomaly_forms', 'theorem masked_z6_kernel',
                       'theorem band_identified', 'theorem three_family_receipt'),
            _reference('Lean/Screen/ExteriorSelection.lean', 'theorem witten_automatic'),
            _reference('model/common_source/chiral_candidate/action.py',
                       'def _identities(', '"SU3_Cartan_112"')],
            'The same mask has15 internal components and vanishing mixed/U1 anomalies; '
            'the actual color cubic traces vanish and each generation has four weak doublets. '
            'The complete faithful cost-minimal band premises select dimension3 and hence45 '
            'internal components.',
            'Band-to-family identification is a supplied attachment. Anomaly arithmetic does '
            'not construct a chiral quantum regulator or count spin, color or Nambu rows as flavors.'),
        _row('same_current_higgs_lines', ['PR-48', 'PR-50', 'PR-54', 'PR-59'], [
            _reference(matter, 'def invariant_higgs_lines(', 'def relative_spectrum_freedom(')],
            'The complete current invariance equations in the stated symmetric matter-pair '
            'and one-doublet grammar have three one-dimensional Higgs lines. With trivial '
            'gauge action on the three-family factor they allow three arbitrary complex3x3 '
            'Yukawa matrices.',
            'Exhaustion is within that field/interaction grammar. The currents do not fix '
            'the27 complex entries, a mass hierarchy or loop-corrected poles.'),
        _row('weyl_and_conjugation_conventions', ['PR-47', 'PR-48', 'PR-50'], [
            _reference('Lean/Screen/WeylYukawaConventions.lean',
                       'theorem weyl_up_neutral', 'theorem weyl_down_neutral',
                       'theorem weyl_lepton_neutral', 'theorem mixed_convention_not_neutral',
                       'theorem epsilon_SL2_invariant'),
            _reference('model/common_source/chiral_candidate/kinetic.py',
                       'def weyl_kinetic(', 'def bdg_symbol(')],
            'All-left and four-component charge conventions agree after conjugating the '
            'right-handed singlets. The Weyl epsilon contraction and complete(p,-p) Nambu '
            'symbol give the stated quadratic mass-shell identity with one supplied speed.',
            'The Weyl principal symbol and physical Spin/CAR attachment are inputs; a finite '
            'matrix identity is not a continuum chirality or interacting-pole construction.'),
        _row('canonical_electroweak_breaking', ['PR-48', 'PR-49', 'PR-67', 'PR-68', 'PR-69'], [
            _reference('Lean/Screen/ElectroweakBreakingComposition.lean',
                       'theorem minimizer_iff', 'theorem gaugeMassForm_exact',
                       'theorem stabilizer_exactly_one_line',
                       'theorem goldstone_exactly_three_parameters', 'theorem photon_massless')],
            'The positive one-doublet quartic potential has the declared vacuum shell, '
            'one electromagnetic stabilizer and three broken electroweak directions. '
            'With canonical real-vector normalization and w=vSM/sqrt2, the tree coefficients '
            'are mW^2=gW^2*vSM^2/4, mZ^2=(gW^2+gY^2)*vSM^2/4 and mH^2=2*lambda*vSM^2.',
            'The Lean gauge mass form is the coefficient of|D H|^2; comparison with '
            '(1/2)A^T M^2 A supplies the factor2. Potential minimization does not select '
            'the interacting quantum phase or identify running parameters with pole masses.'),
        _row('local_gauge_jet_covariance', ['PR-45', 'PR-48', 'PR-54'], [
            _reference('Lean/Screen/LocalGaugeJetAction.lean',
                       'theorem covariant_first_jet', 'theorem curvature_second_jet',
                       'theorem quadratic_trace_ward', 'theorem missing_parameter_derivative_control')],
            'The inhomogeneous connection variation cancels the derivative of the local '
            'gauge parameter; curvature transforms covariantly and its quadratic cyclic '
            'trace obeys the finite Ward identity.',
            'These algebraic jet identities do not select a source connection, prove a '
            'quantum Ward identity or construct a dynamical gauge-field measure.'),
        _row('native_current_kinetic_higgs_intertwiners', ['PR-47', 'PR-48', 'PR-50', 'PR-54', 'PR-59'], [
            _reference(matter, 'def representation(', 'def action_intertwiner(',
                       'def field_intertwiner(', 'def kinetic_attachment(')],
            'The explicit exterior-basis map intertwines all twelve native currents, '
            'both Higgs components and all three Weyl kinetic axes for arbitrary complex '
            'family matrices, preserving the physical CAR pairing and tensor signs.',
            'The spin lift, exterior choice and physical CAR functor are supplied. '
            'Internal representation matching does not identify source spatial propagation '
            'with the Weyl symbol.'),
        _row('complete_static_quasiparticle_attachment', ['PR-15', 'PR-47', 'PR-50', 'PR-53', 'PR-54'], [
            _reference('model/defect_field_benchmark/quasiparticles.py',
                       'def quasiparticles(', 'def validate_chart(', 'def current_probe_data('),
            _reference('model/defect_field_bridge/colored.py',
                       'def ordered_cover_embedding(', 'def complete_native_certificate('),
            _reference('model/defect_field_benchmark/README.md',
                       'All 534,600 hard-core states are present.',
                       'The 48,060 same-port collision', 'There are 84 positive charged directions'),
            _reference('model/defect_field_bridge/REPRESENTATION.md',
                       'K_h = (1/sqrt(2))', 'E_col=[I;J]/sqrt(2)')],
            'The specified calibrated static Weyl-Higgs action has a complete canonical '
            'quasiparticle chart, and its conserved two-quasiparticle sector attaches to '
            'the actual colored repair law with crossing action-SWAP and retained complements.',
            'This benchmark imports its calibration. Its fixed two-particle law does not '
            'include collision sectors, other particle numbers or broken-current/Higgs '
            'pair creation as a full interacting continuum action.'),
        _row('full_rectangular_spectral_transfer', ['PR-15', 'PR-53', 'PR-71', 'PR-72'], [
            _reference('model/conditional_standard_model/matching.py',
                       'def check_intertwiner(', 'residual = hs*j - j*he'),
            _reference('model/conditional_standard_model/PROOF.md',
                       '## Transport to a microscopic defect sector',
                       'dist(E,spectrum(Hs))<=epsilon')],
            'For self-adjoint finite Hs,Ht, isometric J and full rectangular error '
            '||Hs J-J Ht||<=epsilon, every target eigenvalue lies within epsilon of '
            'the source spectrum. The off-real resolvent error is at most '
            'epsilon/|Im z|^2 in the same clock and energy reference.',
            'Compression Jdagger Hs J=Ht alone omits leakage. Multiplicity, isolated '
            'residues, defect topology and continuum limits need additional hypotheses.')]

    interfaces = [
        _interface('finite_premise_bundle', 'conditional_representation_and_quadratic_action',
                   ['PR-09', 'PR-10', 'PR-11', 'PR-36', 'PR-59'], [
            'Supply SMStructurePremiseData and the complete faithful minimal-cost family-band '
            'data, plus its isometry to the family factor of the same-current action.',
            'Either conjugate matter mask is admissible; the second mask is not another '
            'generation. Reading the simulator laws may supply these specifications '
            'without deriving their selection.']),
        _interface('full_fock_spin_chiral_multiplicity', 'continuum_interpretation',
                   ['PR-47', 'PR-53', 'PR-54', 'PR-59'], [
            'Specify graded CAR, spacetime Spin covariance and chiral gauge representations '
            'on the complete field algebra, including pair-changing operations and collision sectors.',
            'Control every surviving light taste or mirror in the same limit; do not delete '
            'unwanted bands or identify weak SU2 with spacetime spin.',
            'Distinguish45 internal Weyl components,90 spin-CAR components and180 Nambu '
            'rows from flavor multiplicity. The minimal grammar has three massless '
            'neutrino fields; a Weinberg term or right-handed neutrinos is a separate extension.']),
        _interface('same_source_continuum_matching', 'continuum_interpretation',
                   ['PR-16', 'PR-52', 'PR-53', 'PR-54', 'PR-58', 'PR-71', 'PR-72', 'PR-74'], [
            'Use one indexed source law, state/reference, local observable algebra, clock '
            'and spatial refinement for the current, kinetic, Higgs and propagation maps.',
            'Supply controlled convergence of local correlations or resolvents, positivity, '
            'unitarity, locality, chiral Ward identities and renormalization on that family; '
            'the hypothesis is not the assertion that its limit is the Standard Model.',
            'Bound source-word, regulator and discarded-sector errors in the claimed '
            'finite-time, volume and continuum domains, retaining complete outcomes.']),
        _interface('effective_action_and_normalization', 'conditional_representation_and_quadratic_action',
                   ['PR-45', 'PR-48', 'PR-49', 'PR-50', 'PR-54', 'PR-67', 'PR-68', 'PR-69'], [
            'State the flat-spacetime local Lorentz/gauge-invariant polynomial field '
            'and operator grammar, the dimension-four truncation or controlled '
            'higher-operator remainder, and '
            'positive kinetic matrices whitened to a common canonical normalization.',
            'Supply gauge couplings, Higgs quadratic/quartic coefficients and the three '
            'complex Yukawa matrices. Allowed topological theta and vacuum-energy '
            'coefficients are not fixed or set to zero by the current/Higgs identities.',
            'The three-line exhaustion is within the stated pair-Higgs grammar; no '
            'classification of every possible effective operator is inferred.']),
        _interface('vacuum_and_phase', 'continuum_interpretation', ['PR-07', 'PR-08', 'PR-49'], [
            'Supply the reference state and an interacting quantum phase with the '
            'required positivity, spectrum and symmetry properties on the same algebra.',
            'An algebraic Higgs shell or a stationary finite repair distribution does '
            'not itself instantiate that physical vacuum. Source selection may be '
            'investigated separately from the conditional implication.']),
        _interface('defect_sector_intertwiner', 'microscopic_defect_identification',
                   ['PR-15', 'PR-47', 'PR-53', 'PR-54', 'PR-71', 'PR-72'], [
            'Provide an isometric graded sector map J, an actual localized or transported '
            'defect sector and its charge/topology dictionary, and same-reference observable maps.',
            'Intertwine currents, Higgs tensors, spatial dynamics and the generator with '
            'the full residual Hs J-J Ht, including leakage outside the image.',
            'The radial zero doublet and finite interacting charged poles have no supplied '
            'map to the45-component field target; matching degeneracies or an index is insufficient.']),
        _interface('physical_global_group', 'physical_global_group_identification',
                   ['PR-35', 'PR-46', 'PR-47', 'PR-54'], [
            'Identify the actual source loop/character and admissible line sectors with '
            'the representation kernel when selecting a physical global quotient.',
            'Local matter descends through all four quotients by subgroups of Z6; '
            'the compatible faithful representation quotient alone selects no physical '
            'global form or fermion-parity/Spin attachment.'], [
            _reference('Lean/Screen/Z6Descent.lean', 'theorem four_admissible_global_forms',
                       'theorem no_universal_fermion_minus_one'),
            _reference(matter, 'def global_cover_contract(')]),
        _interface('parameter_and_unit_provenance', 'numerical_or_physical_comparison',
                   ['PR-15', 'PR-17', 'PR-50'], [
            'Record whether every coupling, texture, phase, scale and state is supplied, '
            'source-derived or empirically calibrated; conditional symbolic formulas '
            'do not require a numerical parameter-selection theorem.',
            'Calibrated_matter supplies running MSbar common-MZ inputs and CKM data. '
            'Physical units require a same-generator clock anchor, and independent '
            'mass predictions cannot reuse fitted target masses as derived outputs.'], [
            _reference('model/common_source/calibrated_matter.py',
                       'def parameters(', 'def mass_coordinates(', 'def bosonic_tree_control(')]),
        _interface('interacting_poles_and_observables', 'physical_particle_spectrum',
                   ['PR-15', 'PR-28', 'PR-53', 'PR-54', 'PR-71', 'PR-72'], [
            'Match source and target correlation functions, states and observable '
            'interpolators with a common energy origin and clock; use gauge-invariant '
            'or explicitly dressed charged observables.',
            'Specify spectral separation and residues for stable isolated poles, or '
            'the corresponding resonance/infrared prescription where such poles are absent.',
            'Separate quadratic w*sigma(Y), running parameters, finite-volume energy '
            'gaps and physical interacting pole or resonance masses. A resolvent bound '
            'off the real axis alone proves neither pole multiplicity nor a stable particle.']),
        _interface('optional_same_source_einstein_join', 'optional_gravity_composition',
                   ['PR-16', 'PR-23', 'PR-24', 'PR-25', 'PR-26', 'PR-27', 'PR-52', 'PR-54'], [
            'Use the same limiting matter action, Hilbert stress, state, metric/tetrad '
            'and proper clock in the ten regional entropy/area/stress readouts.',
            'Require the matched bulk, fixed-volume area and generalized-entropy '
            'stationarity errors to be o(ell^4), with eta tending to a positive '
            'spacetime constant and sufficient metric/tensor convergence.',
            'Independent Ward conservation and Bianchi conservation make Lambda '
            'constant. Actual finite-ball deficits yield nonlinear Einstein; a '
            'first variation alone yields its linearized equation.',
            'Include nonminimal Higgs-curvature and curvature-counterterm contributions '
            'in the stress/entropy functional or a controlled truncation. Finite '
            'responsive Gram data do not instantiate these regional hypotheses.'], [
            _reference('model/matched_einstein_limit/PROOF.md',
                       '## The common source interface', '## Quantitative Einstein residual')])]
    premise_ids = sorted({p for row in chain + interfaces for p in row['canonical_premises']})
    return {
        'schema': 'oph.conditional_standard_model.contract.v1',
        'conclusion_kind': 'CONDITIONAL_REPRESENTATION_AND_CANONICAL_FIELD_ACTION',
        'conditional_derivation_requires_source_selection': False,
        'conditional_derivation_requires_all_physical_gates_accepted': False,
        'source_selection_derived': False,
        'physical_parameter_values_predicted': False,
        'interacting_continuum_constructed': False,
        'microscopic_defect_matching_instantiated': False,
        'einstein_regional_premises_instantiated': False,
        'canonical_premises_discharged': [],
        'physical_comparison_count_added': 0,
        'canonical_premise_register': {'path': REGISTER, 'ids': premise_ids},
        'theorem_chain': chain,
        'conditional_interfaces': interfaces,
        'verification': {
            'mode': 'STRUCTURAL_REFERENCE_AND_SCOPE_CHECK',
            'paths_anchors_and_canonical_ids_checked': True,
            'lean_modules_recompiled_by_contract': False,
            'python_scientific_proofs_reexecuted_by_contract': False,
            'source_hashes_owned_by_package_receipt': True}}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _local_file(root, name):
    _require(type(name) is str and bool(name), 'a referenced path must be a nonempty string')
    relative = Path(name)
    _require(not relative.is_absolute() and '..' not in relative.parts,
             'referenced path must stay inside the repository')
    path = (root / relative).resolve()
    _require(path.is_relative_to(root) and path.is_file(), 'missing or external referenced path: '+name)
    return path


def validate_contract(value, root=ROOT):
    """Validate references and the exact scoped interface, without proof replay."""
    root = Path(root).resolve()
    _require(type(value) is dict, 'contract must be a JSON object')
    registration = value.get('canonical_premise_register')
    _require(type(registration) is dict and registration.get('path') == REGISTER,
             'canonical premise register changed')
    register = json.loads(_local_file(root, registration['path']).read_text())
    rows = register.get('rows')
    _require(type(rows) is list and all(type(row) is dict for row in rows),
             'canonical premise register needs typed rows')
    ids = [row.get('id') for row in rows]
    _require(all(type(identifier) is str for identifier in ids) and len(ids) == len(set(ids)),
             'canonical premise IDs must be unique strings')
    cited = set()
    reference_count = 0
    for key in ('theorem_chain', 'conditional_interfaces'):
        entries = value.get(key)
        _require(type(entries) is list and entries, 'missing theorem or interface rows')
        labels = []
        for row in entries:
            _require(type(row) is dict and type(row.get('id')) is str and bool(row['id']),
                     'theorem/interface row needs an identifier')
            labels.append(row['id'])
            premises = row.get('canonical_premises')
            _require(type(premises) is list and premises and
                     all(type(p) is str and p in ids for p in premises) and
                     len(premises) == len(set(premises)), 'unknown or duplicate canonical premise ID')
            cited.update(premises)
            references = row.get('references')
            _require(type(references) is list, 'references must be a list')
            for ref in references:
                _require(type(ref) is dict and set(ref) == {'path', 'anchors'},
                         'reference needs exactly path and anchors')
                text = _local_file(root, ref['path']).read_text()
                anchors = ref['anchors']
                _require(type(anchors) is list and anchors and
                         all(type(a) is str and bool(a) and a in text for a in anchors),
                         'missing or invalid theorem anchor: '+ref['path'])
                reference_count += 1
        _require(len(labels) == len(set(labels)), 'duplicate theorem/interface identifier')
    _require(registration.get('ids') == sorted(cited), 'canonical premise reference index differs')
    try:
        actual = json.dumps(value, sort_keys=True, allow_nan=False)
    except (TypeError, ValueError) as error:
        raise ValueError('contract must contain finite JSON values') from error
    _require(actual == json.dumps(_specification(), sort_keys=True, allow_nan=False),
             'conditional theorem scope or interface changed')
    return {'theorem_rows': len(value['theorem_chain']),
            'conditional_interfaces': len(value['conditional_interfaces']),
            'canonical_premise_references': len(cited), 'anchored_references': reference_count,
            'proof_reexecution_claimed': False, 'canonical_discharges_added': 0}


def theorem_contract(root=ROOT):
    """Return a detached JSON contract after checking all cited source anchors."""
    value = _specification()
    validate_contract(value, root)
    return value
