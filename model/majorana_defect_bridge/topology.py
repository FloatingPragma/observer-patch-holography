"""Universal additive charges with Majorana masses and a declared mixing graph."""
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

from ..conditional_standard_model.defect_constraints import SPECIES as CHARGED

SPECIES = CHARGED+('nu_1', 'nu_2', 'nu_3')
FULL_PMNS = tuple((i, j) for i in range(3) for j in range(3))


def _domain(massive, pmns_edges, exclude_massless, electroweak_instanton):
    if type(massive) is not tuple or len(massive) != 3 or any(type(x) is not bool for x in massive):
        raise ValueError('three Boolean massive-neutrino flags required')
    if any(type(x) is not bool for x in (exclude_massless, electroweak_instanton)):
        raise ValueError('explicit Boolean charge hypotheses required')
    if (type(pmns_edges) is not tuple
            or any(type(edge) is not tuple or len(edge) != 2
                   or any(type(x) is not int or not 0 <= x < 3 for x in edge) for edge in pmns_edges)):
        raise ValueError('distinct exact charged-lepton/neutrino index pairs required')
    if len(set(pmns_edges)) != len(pmns_edges):
        raise ValueError('distinct exact charged-lepton/neutrino index pairs required')


def relation_matrix(massive=(False, True, True), *, pmns_edges=FULL_PMNS,
                    exclude_massless=True, electroweak_instanton=False):
    """Relations retain each field before using mixing or Majorana constraints.

    A mixing edge is an independently admitted nonzero charged-current
    entry. It is not inferred from an arbitrary symmetric mass matrix.
    """
    _domain(massive, pmns_edges, exclude_massless, electroweak_instanton)
    rows, names = [], []

    def add(name, coefficients):
        rows.append([coefficients.get(field, 0) for field in SPECIES])
        names.append(name)

    for up in SPECIES[:3]:
        for down in SPECIES[3:6]:
            add('quark charged current: '+up+','+down, {up: 1, down: -1, 'W-': 1})
    for lepton, neutrino in pmns_edges:
        add('lepton charged current: '+SPECIES[6+lepton]+','+SPECIES[12+neutrino],
            {'W-': 1, SPECIES[6+lepton]: -1, SPECIES[12+neutrino]: 1})
    for i, flag in enumerate(massive):
        if flag:
            add('Majorana self-conjugation: nu_'+str(i+1), {'nu_'+str(i+1): 2})
        elif exclude_massless:
            add('massless field excluded from defect assignment: nu_'+str(i+1), {'nu_'+str(i+1): 1})
    add('neutral conjugate-pair Z vertex', {'Z': 1})
    add('neutral conjugate-pair Higgs vertex', {'Higgs': 1})
    if electroweak_instanton:
        add('three-generation electroweak (udd nu)^3 process',
            {'up': 1, 'charm': 1, 'top': 1, 'down': 2, 'strange': 2, 'bottom': 2,
             'nu_1': 1, 'nu_2': 1, 'nu_3': 1})
    return s.ImmutableMatrix(rows), names


def universal_group(massive=(False, True, True), **hypotheses):
    matrix, names = relation_matrix(massive, **hypotheses)
    diagonal = smith_normal_form(matrix, domain=s.ZZ)
    values = [abs(int(diagonal[i, i])) for i in range(min(diagonal.shape)) if diagonal[i, i] != 0]
    return {'species_order': list(SPECIES), 'massive_neutrino_mask': list(massive),
            'relation_names': names, 'relation_matrix': [list(map(int, matrix.row(i))) for i in range(matrix.rows)],
            'Smith_nonzero_diagonal': values, 'free_rank': len(SPECIES)-len(values),
            'torsion_invariant_factors': [x for x in values if x > 1]}


def allowed_assignment(values, massive=(False, True, True), *, modulus=None, **hypotheses):
    if (type(values) is not dict or set(values) != set(SPECIES)
            or any(type(x) is not int for x in values.values())):
        raise ValueError('complete exact integer charge dictionary required')
    if modulus is not None and (type(modulus) is not int or modulus < 2):
        raise ValueError('integer modulus at least two required')
    matrix, _ = relation_matrix(massive, **hypotheses)
    residual = matrix*s.Matrix([values[name] for name in SPECIES])
    return all(x == 0 if modulus is None else x % modulus == 0 for x in residual)


def certificate():
    d, e, n = s.symbols('D E N', integer=True)
    general = [d-e+n]*3+[d]*3+[e]*3+[e-n, 0, 0]+[n]*3
    matrix, names = relation_matrix((True, True, True), electroweak_instanton=True)
    residual = matrix*s.Matrix(general)
    for label, value in zip(names, residual, strict=True):
        expected = 2*n if label.startswith('Majorana') else (9*d-3*e+6*n if label.startswith('three-generation') else 0)
        if s.expand(value-expected) != 0:
            raise ArithmeticError('the full field-charge reduction changed')
    massive = universal_group((True, True, True))
    massive_instanton = universal_group((True, True, True), electroweak_instanton=True)
    massless = universal_group()
    massless_instanton = universal_group(electroweak_instanton=True)
    if ((massive['free_rank'], massive['torsion_invariant_factors']) != (2, [2])
            or (massive_instanton['free_rank'], massive_instanton['torsion_invariant_factors']) != (1, [6])
            or (massless['free_rank'], massless['torsion_invariant_factors']) != (2, [])
            or (massless_instanton['free_rank'], massless_instanton['torsion_invariant_factors']) != (1, [3])):
        raise ArithmeticError('the Majorana additive charge presentation changed')
    return {'hypotheses': {'abelian_additive_charge': True, 'antiparticle_charge_is_inverse': True,
                           'each_listed_mass_eigenfield_has_a_definite_additive_charge': True,
                           'all_nine_quark_charged_current_entries_nonzero': True,
                           'full_PMNS_support_or_sufficient_connected_bipartite_graph': True,
                           'massive_Majorana_field_is_self_conjugate': True,
                           'massless_fields_excluded_from_defect_assignment': True,
                           'these_hypotheses_are_proved_for_OPH_defect_dynamics': False},
            'general_connected_graph_assignment': dict(zip(SPECIES, map(str, general), strict=True)),
            'Majorana_relation': '2*N=0', 'instanton_relation_before_reduction': '9*D-3*E+6*N=0',
            'all_three_positive_Majorana_masses': massive,
            'all_three_positive_with_optional_instanton': massive_instanton,
            'calibrated_one_massless_excluded': massless,
            'calibrated_one_massless_with_optional_instanton': massless_instanton,
            'integer_valued_neutrino_winding_for_any_positive_Majorana_mass_connected_graph': 0,
            'one_excluded_massless_neutrino_kills_all_neutrino_charges_in_connected_graph': True,
            'distinct_generation_additive_charges_under_connected_mixing': False,
            'disconnected_graphs_classified_by_same_presented_matrix_not_same_group_claim': True,
            'torsion_discarded': False, 'native_instanton_or_defect_correspondence_constructed': False,
            'neutral_composite_nonadditive_or_approximate_mechanisms_excluded': False}
