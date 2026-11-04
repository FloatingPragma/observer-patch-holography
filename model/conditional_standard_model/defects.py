"""Conditional defect classes and unassigned massive-field interfaces.

Homotopy groups classify the stated coefficient or ordered-state spaces.
Operator indices have their own symmetry domains. Neither supplies a
Standard Model particle assignment without a field/state intertwiner.
"""
from fractions import Fraction as F

import sympy as s

from ..native_mass_channels import algebra as mass
from ..bounded_defect_index import algebra as reflection
from ..topological_defect_spectrum import algebra as chiral
from ..radial_defect_gap import bounds as radial
from ..spin_texture_energy import energy as spin
from . import spectrum


SOURCE_PATHS = (
    'model/native_mass_channels/algebra.py',
    'model/native_mass_channels/README.md',
    'model/topological_defect_spectrum/algebra.py',
    'model/topological_defect_spectrum/PROOF.md',
    'model/bounded_defect_index/algebra.py',
    'model/bounded_defect_index/PROOF.md',
    'model/radial_defect_gap/bounds.py',
    'model/radial_defect_gap/PROOF.md',
    'model/spin_texture_energy/energy.py',
    'model/spin_texture_energy/PROOF.md',
)


def component_certificate():
    """Replay finite identities supporting the proved three-component retract."""
    classification = mass.exhaustion_certificate()
    rows = []
    examples = ((5, 1, 2, 1), (-5, 1, 2, 1), (1, 2, -3, 1))
    types = {'positive_scalar': 'point', 'negative_scalar': 'point', 'spin_texture': 'S2'}
    for initial in examples:
        component = mass.phase(initial)
        path = [mass.component_retraction(initial, t) for t in (F(0), F(1, 3), F(1))]
        if any(mass.phase(point) != component for point in path):
            raise ArithmeticError('the source component retraction crossed the gap boundary')
        endpoint = path[-1]
        if component == 'spin_texture':
            valid = endpoint[0] == 0 and sum(abs(x) for x in endpoint[1:]) == 1
        else:
            valid = endpoint == (F(1 if component == 'positive_scalar' else -1), F(0), F(0), F(0))
        if not valid or mass.component_retraction(endpoint, F(2, 5)) != endpoint:
            raise ArithmeticError('the source deformation retraction endpoint changed')
        rows.append({'component': component, 'homotopy_type': types[component],
                     'pi1': '0', 'pi2': 'Z' if component == 'spin_texture' else '0',
                     'pi3': 'Z' if component == 'spin_texture' else '0',
                     'initial': list(map(str, initial)),
                     'retraction_witness': [list(map(str, point)) for point in path]})
    walls = [mass.wall_rotation(t) for t in (F(0), F(1, 2), F(1), F(3, 2), F(2))]
    if any(mass.bulk_gap_squared(point) != 1 for point in walls):
        raise ArithmeticError('the uniform spin sign path lost its unit gap')
    return {'domain': 'Hermitian onsite masses anticommuting with the fixed two-spin Weyl hopping',
            'edge_nullspace': classification, 'gapped_component_count': 3, 'components': rows,
            'scalar_to_spin_or_opposite_scalar_wall_requires_leaving_gapped_component': True,
            'three_components_identified_as_three_generations': False,
            'spin_sign_wall_constant_gap_path': [list(map(str, point)) for point in walls],
            'uniform_coefficient_path_proves_interface_operator_gap': False,
            'point_defect_type_count_in_spin_component': 1,
            'point_defect_charge_set': 'Z', 'nontrivial_line_defect_pi1': False,
            'classification_exhausts_all_OPH_defects': False}


def degree_certificate():
    values = mass.VERTICES
    reverse = tuple(tuple(-x for x in row) for row in values)
    plus, minus = mass.texture_degree(values), mass.texture_degree(reverse)
    perturbed = tuple(tuple(F(x)+F((j+axis) % 3-1, 10) for axis, x in enumerate(row))
                      for j, row in enumerate(values))
    stability = mass.hedgehog_stability(perturbed)
    if (plus['degree'], minus['degree'], mass.texture_degree(perturbed)['degree']) != (1, -1, 1):
        raise ArithmeticError('the continuous octahedral degree witness changed')
    return {'hedgehog': plus, 'antihedgehog': minus,
            'perturbed_hedgehog': mass.texture_degree(perturbed), 'stability': stability,
            'domain': 'nonvanishing continuous radial piecewise-linear interpolation on the octahedral sphere',
            'bare_lattice_vertex_values_define_degree_without_interpolation': False,
            'zero_mode_index_for_arbitrary_integer_degree_proved': False}


def operator_certificate():
    internal = chiral.clifford_certificate()
    lattice = reflection.atomic_certificate()
    strong = radial.encode(radial.spectral_window(mass=64))
    if (internal['positive_texture_chiral_index'], internal['negative_texture_chiral_index'],
            lattice['reflection_index_plus_hedgehog'], lattice['reflection_index_antihedgehog']) != (2, -2, 2, 2):
        raise ArithmeticError('the internal and spatial-inversion indices were conflated')
    return {
        'continuum_internal_chiral': {'algebra': internal,
            'grading': 'Gamma=i*gamma1*gamma2*gamma3*beta1*beta2*beta3',
            'profile_domain': 'the explicit bounded radial hedgehog/antihedgehog or strictly convex smooth radial core',
            'protection': 'continuous bounded chiral perturbation path with zero outside essential spectrum',
            'index_is_mass_or_electric_charge': False},
        'lattice_spatial_inversion': {'atomic': lattice,
            'grading': '(R psi)(n)=psi(-n)',
            'protection': 'bounded self-adjoint reflection-odd compact perturbations, or reflection-odd perturbations of norm less than m',
            'kernel_dimension_lower': 2,
            'reflection_index_is_hedgehog_degree': False},
        'strong_mass_exact_kernel_example': strong,
        'two_zero_modes_identified_with_three_generations': False,
        'zero_operator_energy_identified_with_massless_propagating_particle': False,
        'native_repair_preserves_index_symmetry_proved': False,
        'autonomous_texture_formation_energy_or_stable_particle_proved': False,
    }


def ordered_texture_certificate():
    identities = spin.exact_identities()
    if identities['compact_Hopf_lift_absolute_degree'] != 1:
        raise ArithmeticError('the explicit compact Hopf lift changed')
    return {'ordered_orbit': 'SO(3)/SO(2)=S2 in the selected ordered-spin ground orbit',
            'point_boundary_homotopy': 'pi2(S2)=Z', 'compact_texture_homotopy': 'pi3(S2)=Z',
            'exact_identities': identities,
            'compact_Hopf_charge_absolute': 1,
            'compact_Hopf_energy': 'DeltaE(R(x/L))=L*E_cont+O(1), E_cont>0, fixed J clock',
            'isolated_radial_hedgehog_energy': 'liminf DeltaE/L>=4*pi*J*kappa>0',
            'macroscopic_finite_mass_particle_in_these_families': False,
            'finite_lattice_unconstrained_rotation_assignment_is_contractible': True,
            'microscopic_conserved_topological_operator_proved': False,
            'all_quantum_defects_or_finite_energy_spectral_components_excluded': False}


def electroweak_orbit_certificate():
    """Exact transitive SU(2) action on a canonical doublet of fixed norm."""
    x1, x2, x3, x4 = s.symbols('x1 x2 x3 x4', real=True)
    w = s.symbols('w', positive=True)
    z1, z2 = x1+s.I*x2, x3+s.I*x4
    norm2 = x1*x1+x2*x2+x3*x3+x4*x4
    unitary = s.Matrix([[s.conjugate(z2), z1], [-s.conjugate(z1), z2]])/w
    if ((unitary.H*unitary-norm2/w**2*s.eye(2)).applyfunc(s.expand) != s.zeros(2)
            or s.expand(unitary.det()-norm2/w**2) != 0
            or unitary*s.Matrix([0, w]) != s.Matrix([z1, z2])):
        raise ArithmeticError('the fixed-norm canonical Higgs orbit changed')
    return {'Higgs_representation': 'one canonical complex SU(2) doublet',
            'minimum_constraint': 'H^dagger H=vEW^2/2, vEW>0',
            'vacuum_orbit_before_quotient_by_all_local_gauge_transformations': 'S3',
            'exact_SU2_transitive_orbit_identity': True,
            'pi1': '0', 'pi2': '0', 'pi3': 'Z',
            'topologically_required_string_or_point_monopole_in_this_orbit': False,
            'pi3_proves_stable_particle': False,
            'global_group_condition': 'the canonical doublet action; a central quotient acting trivially on it leaves this orbit unchanged',
            'other_global_groups_or_stabilizer_embeddings_classified': False,
            'embedded_strings_singular_monopoles_or_gauge_saddles_excluded': False,
            'primary_reference': 'https://arxiv.org/abs/1610.05623'}


def family_mapping_contract(inputs=None):
    """Candidate massive targets under a supplied family-to-defect premise.

    No topological assignment is invented. Zero Yukawa singular values are
    excluded individually, including within a charged family channel.
    """
    inputs = spectrum.calibrated_inputs() if inputs is None else inputs
    if not isinstance(inputs, spectrum.ActionInputs):
        raise ValueError('ActionInputs required')
    candidates, excluded = [], ['photon', 'eight_gluons', 'nu_1', 'nu_2', 'nu_3']
    for name, matrix in inputs.yukawas().items():
        _, roots = spectrum.singular_spectrum(matrix)
        massive = []
        for j, root in enumerate(roots):
            if root == 0:
                excluded.append(name+':family_mass_order_'+str(j+1))
            else:
                massive.append(j+1)
        candidates.append({'target': name, 'kind': 'charged fermion family channel',
                           'massive_family_mass_orders': massive,
                           'eligible_massive_target': bool(massive),
                           'assigned_defect_type': None, 'assigned_topological_charge': None,
                           'SM_representation_intertwiner': None})
    for target in ('W_plus', 'W_minus', 'Z', 'Higgs'):
        candidates.append({'target': target, 'kind': 'massive bosonic field',
                           'eligible_massive_target': True,
                           'assigned_defect_type': None, 'assigned_topological_charge': None,
                           'SM_representation_intertwiner': None})
    candidates.append({'target': 'neutrino_mass_extension', 'kind': 'separate action extension',
                       'eligible_massive_target': False,
                       'condition': 'a supplied nonzero neutrino mass action and its field intertwiner',
                       'assigned_defect_type': None, 'assigned_topological_charge': None,
                       'SM_representation_intertwiner': None})
    return {'family_correspondence_is_an_extra_conditional_premise': True,
            'candidates': candidates, 'excluded_massless_targets': excluded,
            'required_witnesses': [
                'same source, vacuum sector and admissible texture configuration space',
                'continuous invariant and full dynamical/repair preservation or quantified lifetime',
                'gauge-current representation and electric/color charges on defect states',
                'Lorentz spin and statistics of quantized collective modes',
                'finite formation energy, localization and translational dynamics',
                'Hamiltonian, field-observable and state intertwiners with controlled errors',
                'mass and decay amplitudes in one calibrated clock',
                'compatibility with the interaction and decay charge constraints'],
            'massless_fields_assigned_defects': False,
            'assignment_derived_from_family_dimension_or_integer_charge_count': False,
            'all_candidate_targets_admit_nonzero_additive_topology': False,
            'additive_charge_constraints_owner': 'model/conditional_standard_model/defect_constraints.py',
            'exact_charge_decay_rule': 'a conserved additive charge of an unstable parent equals the total charge of its decay products',
            'exact_topology_implies_every_massive_species_stable': False}


def defect_catalogue(inputs=None):
    return {'schema': 'oph.conditional_standard_model.defects.v1',
            'onsite_mass_components': component_certificate(),
            'sphere_degree': degree_certificate(), 'operator_indices': operator_certificate(),
            'ordered_spin_textures': ordered_texture_certificate(),
            'canonical_electroweak_orbit': electroweak_orbit_certificate(),
            'family_mapping': family_mapping_contract(inputs),
            'scope': {'all_possible_OPH_defects_classified': False,
                      'topological_families_identified_with_SM_families': False,
                      'protected_operator_index_proves_stable_massive_particle': False,
                      'native_repair_preserves_topological_charge_proved': False,
                      'massless_species_assigned_defects': False}}
