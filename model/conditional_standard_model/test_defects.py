"""Independent finite controls of the scoped defect-family interface."""
from dataclasses import replace
from fractions import Fraction as F
import json

import pytest
import sympy as s

from ..native_mass_channels import algebra as mass
from ..topological_defect_spectrum import algebra as chiral
from ..bounded_defect_index import finite
from . import defects, spectrum


def test_full_edge_nullspace_and_component_endpoints():
    value = defects.component_certificate()
    assert value['edge_nullspace']['constraint_rank'] == 28
    assert value['edge_nullspace']['solution_dimension'] == 4
    assert [row['homotopy_type'] for row in value['components']] == ['point', 'point', 'S2']
    for initial in ((10, -2, 1, 3), (-8, -1, 2, -2), (1, 1, 1, 1), (0, -7, 0, 0)):
        for numerator in range(21):
            point = mass.component_retraction(initial, F(numerator, 20))
            assert mass.phase(point) == mass.phase(initial)
            assert mass.bulk_gap_squared(point).is_positive
    with pytest.raises(ValueError):
        mass.component_retraction((1, 1, 0, 0), F(1, 2))


def test_octahedral_degree_orientation_and_nonvanishing_margin():
    value = defects.degree_certificate()
    assert value['hedgehog']['degree'] == 1
    assert value['antihedgehog']['degree'] == -1
    assert value['stability']['face_separation_lower'] == '7/10'
    # A proper cyclic rotation keeps degree; a single coordinate reflection reverses it.
    rotated = tuple((z, x, y) for x, y, z in mass.VERTICES)
    reflected = tuple((-x, y, z) for x, y, z in mass.VERTICES)
    assert mass.texture_degree(rotated, (2, 3, 5))['degree'] == 1
    assert mass.texture_degree(reflected, (2, 3, 5))['degree'] == -1
    with pytest.raises(ValueError):
        mass.hedgehog_stability(tuple((F(x)+F(1, 3), y, z) for x, y, z in mass.VERTICES))


def test_internal_chirality_and_spatial_inversion_have_different_orientation_rules():
    import numpy as np
    gamma = chiral.generators()[-1]
    for sign in (-1, 1):
        columns = chiral.zero_basis(sign)
        assert s.trace(columns.H*gamma*columns) == 2*sign
        # The complete atomic finite cube has two even origin zeros for BOTH signs.
        h = finite.source_matrix(1, mass=4, sign=sign, hopping_fraction=0)
        even, odd = finite.inversion_bases(1)
        assert np.linalg.matrix_rank(h.toarray(), tol=1e-12) == h.shape[0]-2
        assert even.shape[1]-odd.shape[1] == 2
    value = defects.operator_certificate()
    assert value['strong_mass_exact_kernel_example']['one_body_gap_lower'] == '62'
    assert value['lattice_spatial_inversion']['atomic']['reflection_index_antihedgehog'] == 2
    assert value['continuum_internal_chiral']['algebra']['negative_texture_chiral_index'] == -2


def test_hopf_integral_and_the_static_energy_domain():
    value = defects.ordered_texture_certificate()
    f = s.symbols('f', real=True)
    assert s.integrate(-2*s.sin(f)**2/s.pi, (f, s.pi, 0)) == 1
    assert value['compact_Hopf_charge_absolute'] == 1
    assert value['macroscopic_finite_mass_particle_in_these_families'] is False
    assert value['microscopic_conserved_topological_operator_proved'] is False


def test_canonical_doublet_fixed_norm_is_a_transitive_three_sphere():
    value = defects.electroweak_orbit_certificate()
    assert value['pi1'] == value['pi2'] == '0'
    z1, z2 = (1+s.I)/2, (1-s.I)/2
    u = s.Matrix([[s.conjugate(z2), z1], [-s.conjugate(z1), z2]])
    assert s.simplify(u.det()) == 1
    assert (u.H*u).applyfunc(s.simplify) == s.eye(2)
    assert u*s.Matrix([0, 1]) == s.Matrix([z1, z2])
    assert value['other_global_groups_or_stabilizer_embeddings_classified'] is False


def test_zero_yukawas_cannot_enter_massive_defect_assignment_candidates():
    inputs = replace(spectrum.calibrated_inputs(), Yu=s.diag(0, 1, 2),
                     Yd=s.zeros(3), Ye=s.diag(0, 0, 3))
    result = defects.family_mapping_contract(inputs)
    by_name = {row['target']: row for row in result['candidates']}
    assert by_name['up']['massive_family_mass_orders'] == [2, 3]
    assert by_name['down']['eligible_massive_target'] is False
    assert by_name['charged_lepton']['massive_family_mass_orders'] == [3]
    assert len(result['excluded_massless_targets']) == 11
    assert by_name['neutrino_mass_extension']['eligible_massive_target'] is False
    assert all(row['assigned_defect_type'] is None for row in result['candidates'])


def test_catalogue_is_json_with_separate_invariants_and_null_assignments():
    result = defects.defect_catalogue()
    assert json.loads(json.dumps(result)) == result
    assert result['onsite_mass_components']['point_defect_type_count_in_spin_component'] == 1
    assert result['onsite_mass_components']['three_components_identified_as_three_generations'] is False
    assert result['sphere_degree']['zero_mode_index_for_arbitrary_integer_degree_proved'] is False
    assert all(value is False for value in result['scope'].values())
    rows = result['family_mapping']['candidates']
    assert all(row['assigned_topological_charge'] is None for row in rows)
    assert all(row['SM_representation_intertwiner'] is None for row in rows)
    assert {'photon', 'eight_gluons', 'nu_1', 'nu_2', 'nu_3'} == set(
        result['family_mapping']['excluded_massless_targets'])
