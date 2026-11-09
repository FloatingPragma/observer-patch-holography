"""Independent representation, Hessian, neutrino and actual-texture controls."""
from collections import Counter
from itertools import combinations
import json

import numpy as np
import pytest
import sympy as s

from . import conditional_standard_model as m


def test_charges_from_independent_exterior_representation():
    # wedge^2(5) + wedge^4(5) is 10 + conjugate5, without a spin doubling.
    fundamental_y = [s.Rational(-1, 3)]*3+[s.Rational(1, 2)]*2
    fundamental_t = [0, 0, 0, s.Rational(1, 2), s.Rational(-1, 2)]
    actual_y, actual_q = [], []
    for degree in (2, 4):
        for word in combinations(range(5), degree):
            actual_y.append(sum(fundamental_y[j] for j in word))
            actual_q.append(sum(fundamental_y[j]+fundamental_t[j] for j in word))
    declared_y, declared_q = [], []
    for r in m.representations():
        declared_y += [s.Rational(r['hypercharge'])]*r['color_dimension']*r['weak_dimension']
        declared_q += [s.Rational(q) for q in r['electric_weights']]*r['color_dimension']
    assert Counter(actual_y) == Counter(declared_y)
    assert Counter(actual_q) == Counter(declared_q)
    assert len(actual_q) == 15


def test_anomalies_are_signed_and_conjugate_singlet_matters():
    a = m.anomaly_coefficients()
    assert all(v == 0 for k, v in a.items() if 'doublets' not in k)
    assert a['weak_doublets_three_families'] % 2 == 0
    # Treating e^c as e without conjugation fails both tests by two units.
    rows = m.representations()
    weights = [s.Rational(r['hypercharge']) for r in rows
               for _ in range(r['color_dimension']*r['weak_dimension'])]
    assert weights[-1] == 1
    weights[-1] = -1
    assert sum(weights) == -2 and sum(y**3 for y in weights) == -2


def test_all_left_yukawa_contractions_and_weak_epsilon():
    h_y = s.Rational(1, 2)
    y = {r['name']: s.Rational(r['hypercharge']) for r in m.representations()}
    assert y['Q']+h_y+y['u_c'] == 0
    assert y['Q']-h_y+y['d_c'] == 0
    assert y['L']-h_y+y['e_c'] == 0
    eps = s.Matrix([[0, 1], [-1, 0]])
    for generator in (s.Matrix([[0, 1], [1, 0]])/2,
                      s.Matrix([[0, -s.I], [s.I, 0]])/2, s.diag(1, -1)/2):
        assert generator.T*eps+eps*generator == s.zeros(2)


def test_vector_mass_characteristic_polynomial_and_neutral_rotation():
    q = m.quadratic_action()
    g, gp, v, _ = q['symbols']
    matrix = q['gauge_mass_squared']
    x = s.Symbol('x')
    expected = x*(x-g*g*v*v/4)**2*(x-(g*g+gp*gp)*v*v/4)
    assert s.factor(matrix.charpoly(x).as_expr()-expected) == 0
    # Test a nonspecial weak angle with normalized, orthogonal A and Z.
    numeric = matrix.subs({g: 3, gp: 4, v: 2})
    photon = s.Matrix([0, 0, s.Rational(4, 5), s.Rational(3, 5)])
    massive = s.Matrix([0, 0, s.Rational(3, 5), -s.Rational(4, 5)])
    assert photon.dot(massive) == 0
    assert numeric*photon == s.zeros(4, 1)
    assert numeric*massive == 25*massive


def test_scalar_hessian_and_higgs_normalization():
    q = m.quadratic_action()
    _, _, v, lam = q['symbols']
    assert q['scalar_mass_squared'] == s.diag(0, 0, 0, 2*lam*v*v)
    assert q['Higgs_electric_charge'] == s.zeros(2, 1)
    # The physical radial coordinate has kinetic coefficient 1/2.
    h, dh = s.symbols('h dh', real=True)
    doublet = s.Matrix([0, (v+h)/s.sqrt(2)])
    derivative = doublet.diff(h)*dh
    assert (derivative.H*derivative)[0] == dh*dh/2


def test_complex_yukawa_masses_use_singular_values_not_signed_eigenvalues():
    a = s.Matrix([[0, s.I, 0], [2, 0, 0], [0, 0, 3]])
    assert a != a.H and a*a.H != a.H*a
    assert sorted((a.H*a).eigenvals()) == [1, 4, 9]
    left = s.Matrix([[s.Rational(3, 5), s.Rational(4, 5), 0],
                     [-s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, 1]])
    right = s.diag(1, s.I, -1)
    transformed = left*a*right.H
    assert sorted((transformed.H*transformed).eigenvals()) == [1, 4, 9]


@pytest.mark.parametrize('option,weyl,nu', [('minimal_massless',45,6),
                                         ('majorana_eft',45,6), ('dirac',48,12)])
def test_neutrino_alternatives_and_no_double_count(option, weyl, nu):
    c = m.field_census(option)
    assert c['left_Weyl_components'] == weyl
    assert c['two_component_CAR_coordinates'] == 2*weyl
    assert c['neutrino_on_shell_degrees'] == nu
    assert c['charged_fermion_on_shell_degrees'] == 84
    assert c['boson_on_shell_degrees_unbroken'] == c['boson_on_shell_degrees_broken'] == 28
    assert c['neutrino_mass_values'] is None


def test_majorana_takagi_requires_transpose_not_adjoint():
    # A supplied symmetric complex coefficient has nonnegative masses.
    u = s.diag(s.I, 1, s.exp(s.I*s.pi/4))
    diagonal = s.diag(1, 2, 3)
    mass = u.conjugate()*diagonal*u.H
    assert mass == mass.T
    assert s.simplify(u.T*mass*u) == diagonal
    assert u.H*mass*u != diagonal


def test_no_undeclared_neutrino_completion():
    with pytest.raises(ValueError):
        m.field_census('three_measured_flavor_masses')


def test_eighteen_target_rows_do_not_invent_neutrino_masses():
    rows = m.benchmark_rows()
    assert len(rows) == len({r['target'] for r in rows}) == 18
    assert {r['target'] for r in rows} == {
        'electron','muon','tau','up','down','strange','charm','bottom','top',
        'electron_neutrino','muon_neutrino','tau_neutrino',
        'photon','gluon','W_plus','W_minus','Z','Higgs'}
    neutrinos = [r for r in rows if 'neutrino' in r['target']]
    assert len(neutrinos) == 3 and all(r['tree_mass'] is None for r in neutrinos)
    assert all(r['source_predicted_physical_mass'] is None for r in rows)


def test_actual_686_row_texture_has_nonzero_interior_translation_commutator():
    from model.compact_defect_pair import spectral
    source = spectral.source(radius=3, mass=64, half_separation=1)
    h, points = source['h'], source['points']
    assert h.shape == (686, 686)
    index = {p:i for i,p in enumerate(points)}
    zero, ez = index[(0,0,0)], index[(0,0,1)]
    sigma_z = np.diag([1, -1])
    for spin in range(2):
        image = h[:,2*ez+spin].toarray().ravel()
        shifted = np.zeros(686, dtype=complex)
        col = h[:,2*zero+spin].toarray().ravel()
        for point,i in index.items():
            q = (point[0],point[1],point[2]+1)
            if q in index:
                shifted[2*index[q]:2*index[q]+2] = col[2*i:2*i+2]
        expected = np.zeros(686, dtype=complex)
        expected[2*ez:2*ez+2] = 64*sigma_z[:,spin]
        assert np.array_equal(image-shifted, expected)


def test_same_pauli_kinetics_admit_no_scalar_mass_matrix():
    b0, b1, b2, b3 = variables = s.symbols('b0 b1 b2 b3', real=True)
    sigma = (s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    beta = b0*s.eye(2)+sum((b*t for b,t in zip(variables[1:],sigma)),s.zeros(2))
    equations = [entry for t in sigma for entry in beta*t+t*beta]
    matrix,_ = s.linear_eq_to_matrix(equations,variables)
    assert matrix.rank() == 4
    assert matrix.nullspace() == []


def test_certificate_json_and_immutable_replay_scope():
    report = m.certificate()
    before = json.dumps(report, sort_keys=True)
    assert report['physical_OPH_assignments'] == 0
    for key in ('source_selected_renormalized_parameters','source_selected_vacuum_and_clock',
                'naive_eight_valleys_removed','interacting_continuum_existence_proved',
                'finite_scalar_identified_with_continuum_Higgs',
                'all_reference_masses_inverted_as_tree_parameters'):
        assert report[key] is False
    report['representations_one_family'][0]['electric_weights'][0] = '0'
    report['primary_reviews'].clear()
    assert json.dumps(m.certificate(),sort_keys=True) == before
