"""Conditional tree spectrum from the actual retained Weyl/Higgs tensors.

The default parameter packet is an empirical common-scale calibration.
General exact finite Yukawa matrices are accepted, with their full family
Gram spectra and the actual45-row all-left mass kernel retained.
"""
from dataclasses import dataclass
import sympy as s

from model.common_source import calibrated_matter as calibrated
from model.common_source.chiral_candidate import action
CHANNELS = ('up', 'down', 'charged_lepton')
SPECIES = {'up': ('up', 'charm', 'top'), 'down': ('down', 'strange', 'bottom'),
           'charged_lepton': ('electron', 'muon', 'tau')}


def _exact(value, *, nonnegative=False):
    if isinstance(value, (bool, float, complex)):
        raise ValueError('exact finite numeric action coefficients required')
    value = s.sympify(value)
    if (value.free_symbols or value.atoms(s.Float) or value.is_finite is not True
            or value.is_number is not True):
        raise ValueError('exact finite numeric action coefficients required')
    if nonnegative and (value.is_real is not True or value.is_nonnegative is not True):
        raise ValueError('real nonnegative bosonic action coefficients required')
    return value


def _matrix(value):
    try:
        matrix = s.ImmutableMatrix(value)
    except (TypeError, ValueError) as error:
        raise ValueError('an exact three-by-three complex Yukawa matrix required') from error
    if matrix.shape != (3, 3):
        raise ValueError('an exact three-by-three complex Yukawa matrix required')
    return s.ImmutableMatrix(3, 3, [_exact(x) for x in matrix])


def _zero(matrix):
    return all(s.expand(x) == 0 or s.simplify(x) == 0 for x in matrix)


def _clean(matrix):
    return s.ImmutableSparseMatrix(matrix.applyfunc(lambda x: s.simplify(s.expand(x))))


def _serialize(matrix):
    value = s.ImmutableSparseMatrix(matrix)
    return {'shape': list(value.shape),
            'entries': [[i, j, str(x)] for (i, j), x in sorted(value.todok().items())]}


def _display(value):
    """A decimal display of an exact expression, not a measurement interval."""
    return str(s.N(value, 15))


@dataclass(frozen=True)
class ActionInputs:
    """Canonical kinetic metric, H=(0,vEW/sqrt2), exact supplied couplings.

    Bosonic coefficients are strictly positive, fixing the regular broken
    electroweak domain. Every finite3x3 complex Yukawa is allowed, including repeated and zero
    singular values. Flavor labels denote sorted mass eigenvalues; a
    degenerate eigenspace does not select its own flavor basis.
    """
    vEW: object
    gY: object
    gW: object
    g3: object
    lambda_SM: object
    Yu: object
    Yd: object
    Ye: object
    provenance: tuple
    energy_unit: str = 'supplied model energy unit'

    def __post_init__(self):
        for key in ('vEW', 'gY', 'gW', 'g3', 'lambda_SM'):
            value = _exact(getattr(self, key), nonnegative=True)
            if value.is_positive is not True:
                raise ValueError('strictly positive broken-electroweak bosonic inputs required')
            object.__setattr__(self, key, value)
        for key in ('Yu', 'Yd', 'Ye'):
            object.__setattr__(self, key, _matrix(getattr(self, key)))
        if (type(self.provenance) is not tuple or not self.provenance
                or any(type(x) is not str or not x.strip() for x in self.provenance)
                or type(self.energy_unit) is not str or not self.energy_unit.strip()):
            raise ValueError('explicit action provenance and energy-unit convention required')

    def yukawas(self):
        return dict(zip(CHANNELS, (self.Yu, self.Yd, self.Ye), strict=True))


def calibrated_inputs():
    """Retain the existing empirical M_Z action and its conventions."""
    bosons, yukawas = calibrated.bosonic_parameters(), calibrated.family_yukawas()
    return ActionInputs(
        bosons['v_SM_GeV'], bosons['gY'], bosons['gW'], bosons['g3'], bosons['lambda_SM'],
        yukawas['up'], yukawas['down'], yukawas['charged_lepton'],
        ('model/common_source/calibrated_matter_parameters.json',
         'empirical MSbar M_Z parameters and supplied CKM boundary',
         'tree comparison m2=-lambda_SM*vEW^2; distinct from loop-matched m2'),
        energy_unit='GeV in the supplied empirical MSbar M_Z calibration')


def _inputs(value):
    if value is None:
        return calibrated_inputs()
    if type(value) is not ActionInputs:
        raise ValueError('ActionInputs required')
    return value


def mass_kernel(inputs=None):
    """Actual K45: directed Y* and its transpose, in the retained source order."""
    inputs = _inputs(inputs)
    result = s.zeros(45)
    for name, yukawa in inputs.yukawas().items():
        oriented = calibrated._oriented_components(name)[1]
        directed = s.kronecker_product(oriented, yukawa.conjugate())
        result += directed+directed.T
    return s.ImmutableSparseMatrix(inputs.vEW*result/s.sqrt(2))


def singular_spectrum(yukawa):
    """Exact squared singular eigenvalues, including algebraic RootOf output."""
    y = _matrix(yukawa)
    gram = _clean(y.H*y)
    variable = s.Symbol('z')
    polynomial = s.Poly(gram.charpoly(variable).as_expr(), variable)
    if gram.is_diagonal():
        eigenvalues = list(gram.diagonal())
    elif all(x.is_Rational for x in polynomial.all_coeffs()):
        eigenvalues = polynomial.all_roots(radicals=False)
    else:
        found = gram.eigenvals(error_when_incomplete=True)
        eigenvalues = [value for value, multiplicity in found.items() for _ in range(multiplicity)]
    if len(eigenvalues) != 3:
        raise ArithmeticError('complete family Gram spectrum required')
    # Hermitian Gram positivity is exact; numeric evaluation orders the real
    # roots for conventional family labels without changing their expressions.
    eigenvalues.sort(key=lambda x: float(s.re(s.N(x, 50))))
    roots = tuple(s.sqrt(s.simplify(x)) for x in eigenvalues)
    return {'matrix': _serialize(y), 'Gram': _serialize(gram),
            'Gram_characteristic_polynomial': str(polynomial.as_expr()),
            'squared_singular_values_exact': list(map(str, eigenvalues)),
            'singular_values_exact': list(map(str, roots)),
            'singular_values_decimal_display': list(map(_display, roots)),
            'rank': int(y.rank()),
            'degenerate_flavor_basis_selected': False}, roots


def bosonic_spectrum(inputs=None):
    inputs = _inputs(inputs)
    v, gy, gw, lam = inputs.vEW, inputs.gY, inputs.gW, inputs.lambda_SM
    # |D H|^2 is one half of the real-vector quadratic mass term.
    mass2 = s.Matrix([[gw*gw, -gw*gy], [-gw*gy, gy*gy]])*v*v/4
    photon, z = s.Matrix([gy, gw]), s.Matrix([gw, -gy])
    if mass2*photon != s.zeros(2, 1) or not _zero(mass2*z-(gw*gw+gy*gy)*v*v*z/4):
        raise ArithmeticError('neutral gauge spectrum normalization failed')
    h = s.symbols('h0:4', real=True)
    radius = sum(x*x for x in h)
    potential = -lam*v*v*radius/2+lam*radius*radius/4
    hessian = s.hessian(potential, h).subs(dict(zip(h, (0, 0, v, 0), strict=True)))
    if hessian != s.diag(0, 0, 2*lam*v*v, 0):
        raise ArithmeticError('declared Higgs tree-background Hessian failed')
    rows = []
    for name, square, charge, multiplicity, spin in (
        ('photon', s.S.Zero, '0', 1, 1), ('gluon', s.S.Zero, '0', 8, 1),
        ('W+', gw*gw*v*v/4, '1', 1, 1), ('W-', gw*gw*v*v/4, '-1', 1, 1),
        ('Z', (gw*gw+gy*gy)*v*v/4, '0', 1, 1),
        ('Higgs_radial', 2*lam*v*v, '0', 1, 0)):
        square = s.simplify(square)
        rows.append({'name': name, 'spin_under_conditional_Lorentz_dictionary': spin,
                     'electric_charge': charge, 'multiplicity': multiplicity,
                     'tree_mass_squared_exact': str(square), 'tree_mass_exact': str(s.sqrt(square)),
                     'tree_mass_decimal_display': _display(s.sqrt(square)),
                     'physical_pole_mass': None,
                     'free_polarizations_per_species': 1 if spin == 0 else 2 if square == 0 else 3})
    return {'Higgs_background': 'H=(0,vEW/sqrt(2))',
            'potential': 'V=-lambda_SM*vEW^2*(Hdagger H)+lambda_SM*(Hdagger H)^2',
            'tree_m2_exact': str(-lam*v*v),
            'stationary_tree_potential_is_supplied_comparison': True,
            'tree_m2_is_inferred_from_SMDR_loop_minimization': False,
            'neutral_real_vector_mass_squared_matrix': _serialize(mass2),
            'Higgs_real_coordinate_Hessian': _serialize(hessian),
            'photon_Z_basis_selected_when_gY_and_gW_both_zero': False,
            'rows': rows,
            'Goldstone_directions_in_broken_regime': 3 if v > 0 and gw > 0 else None,
            'gauge_invariant_nonperturbative_particle_poles_constructed': False}


def kernel_certificate(inputs=None):
    inputs = _inputs(inputs)
    kernel = mass_kernel(inputs)
    charge = action.electric_charge()
    if kernel.T != kernel or not _zero(charge.T*kernel+kernel*charge):
        raise ArithmeticError('actual K45 violates all-left symmetry or electric neutrality')
    pairing = s.kronecker_product(action.SPIN_EPSILON, kernel.conjugate())
    if pairing.T != -pairing:
        raise ArithmeticError('actual90-CAR creation pairing is not antisymmetric')
    pairs, rows = set(), []
    for name, y in inputs.yukawas().items():
        oriented = calibrated._oriented_components(name)[1]
        entries = oriented.todok()
        if len(entries) != (1 if name == 'charged_lepton' else 3):
            raise ArithmeticError('source color multiplicity changed')
        for (left, right), sign in entries.items():
            if sign not in (-1, 1) or left in pairs or right in pairs:
                raise ArithmeticError('neutral Higgs source pairs are not disjoint unit couplings')
            pairs.update((left, right))
            block = kernel.extract([3*left+i for i in range(3)], [3*right+i for i in range(3)])
            if block != sign*inputs.vEW*y.conjugate()/s.sqrt(2):
                raise ArithmeticError('actual family kernel is not the declared Yukawa action')
        rows.append({'channel': name, 'neutral_source_pairs': len(entries),
                     'family_rank': int(y.rank()), 'K45_rank_contribution': 0 if inputs.vEW == 0 else 2*len(entries)*int(y.rank())})
    neutral = set(range(15))-pairs
    if len(neutral) != 1:
        raise ArithmeticError('minimal source must retain one unpaired left neutrino per family')
    rank = sum(row['K45_rank_contribution'] for row in rows)
    return {'K45': _serialize(kernel), 'K45_symmetric_exact': True,
            'electric_charge_neutrality_exact': True,
            'creation_Delta90_antisymmetric_exact': True,
            'source_pair_blocks': rows, 'K45_rank_exact': rank,
            'K45_nullity_exact': 45-rank, 'unpaired_source_indices': sorted(neutral),
            'canonical_family_kinetic_metric': 'I3',
            'native_current_intertwiner': 'model/common_source/native_matter_attachment.py:kinetic_attachment',
            'native_control_execution_repeated': False}


def representation_certificate():
    """Read the actual exact charge matrices and neutral Higgs source pairs."""
    labels = action.ordering()['all_left_components']
    groups = {}
    for row in labels:
        group = groups.setdefault(row['internal_type'], {'internal_type': row['internal_type'],
                                  'hypercharge': str(s.Rational(row['q6Y'])/6),
                                  'source_indices': [], 'electric_charges': []})
        group['source_indices'].append(row['index'])
        group['electric_charges'].append(row['electric_Q'])
    for group in groups.values():
        group['dimension_per_family'] = len(group['source_indices'])
    if sorted((name, len(row['source_indices'])) for name, row in groups.items()) != [
            ('L', 2), ('Q', 6), ('d^c', 3), ('e^c', 1), ('u^c', 3)]:
        raise ArithmeticError('retained finite source does not supply the five M15 multiplets')
    return {'internal_dimension_per_family': 15,
            'multiplets_in_all_left_convention': [groups[name] for name in ('Q', 'u^c', 'd^c', 'L', 'e^c')],
            'gauge_current_generators': len(action.generators(families=False)),
            'charged_invariant_Higgs_channels': list(CHANNELS),
            'family_dimension_is_conditional_selected_band': 3,
            'current_and_Higgs_native_intertwiner': 'model/common_source/native_matter_attachment.py:action_intertwiner',
            'physical_family_count_selected_by_current_algebra_alone': False}


def build_catalogue(inputs=None):
    default = inputs is None
    inputs = _inputs(inputs)
    families, rows = {}, []
    for name, y in inputs.yukawas().items():
        record, roots = singular_spectrum(y)
        record['ratios_to_lightest_in_channel_exact'] = (
            [str(s.simplify(value/roots[0])) for value in roots] if roots[0] != 0 else None)
        record['ratios_are_supplied_action_mass_ratios_not_pole_predictions'] = True
        families[name] = record
        colors, charge = (1, '-1') if name == 'charged_lepton' else (3, '2/3' if name == 'up' else '-1/3')
        for family, root in enumerate(roots):
            mass = s.simplify(inputs.vEW*root/s.sqrt(2))
            rows.append({'name': SPECIES[name][family], 'family_mass_order': family+1,
                         'channel': name, 'color_multiplicity': colors,
                         'electric_charge': charge, 'antiparticle_electric_charge': str(-s.Rational(charge)),
                         'conditional_Lorentz_type': 'Dirac charged pair',
                         'tree_running_mass_exact': str(mass), 'tree_running_mass_decimal_display': _display(mass),
                         'physical_pole_mass': None, 'internal_Weyl_components': 2*colors,
                         'continuum_on_shell_states_including_antiparticles': 4*colors,
                         'confined_quark_parameter': name != 'charged_lepton'})
    for family in range(3):
        rows.append({'name': 'nu_'+str(family+1), 'family_basis': 'unselected degenerate massless neutral basis',
                     'electric_charge': '0', 'color_multiplicity': 1,
                     'conditional_Lorentz_type': 'massless left Weyl field and right-helicity antiparticle',
                     'tree_running_mass_exact': '0', 'tree_running_mass_decimal_display': '0',
                     'physical_pole_mass': None, 'internal_Weyl_components': 1,
                     'continuum_on_shell_states_including_antiparticles': 2})
    internal = sum(row['internal_Weyl_components'] for row in rows)
    states = sum(row['continuum_on_shell_states_including_antiparticles'] for row in rows)
    if internal != 45 or states != 90:
        raise ArithmeticError('the conditional three-family continuum count changed')
    result = {
        'schema': 'oph.conditional_standard_model.catalogue.v1',
        'inputs': {'vEW': str(inputs.vEW), 'gY': str(inputs.gY), 'gW': str(inputs.gW),
                   'g3': str(inputs.g3), 'lambda_SM': str(inputs.lambda_SM),
                   'provenance': list(inputs.provenance), 'energy_unit': inputs.energy_unit},
        'families': families, 'fermions': rows, 'bosons': bosonic_spectrum(inputs),
        'finite_representation': representation_certificate(),
        'kernel_checks': kernel_certificate(inputs),
        'counts': {'families': 3, 'internal_Weyl_components': internal,
                   'continuum_fermion_on_shell_states_including_antiparticles': states,
                   'static_spin_times_internal_CAR_annihilation_modes': 90,
                   'static_CAR_mode_count_is_identified_with_on_shell_degrees': False},
        'dispersion': 'E_f(p)^2=|p|^2+m_f^2 in the supplied canonical continuum kinetic action, c=hbar=1',
        'scope': {'conditional_on_given_action_and_electroweak_background': True,
                  'couplings_and_Yukawas_are_inputs': True,
                  'default_uses_existing_empirical_calibration': default,
                  'physical_parameters_selected_by_native_repair': False,
                  'interacting_physical_pole_masses_or_widths': False,
                  'three_massless_neutrinos_fit_observed_oscillations': False,
                  'optional_Weinberg_neutrino_extension_included': False,
                  'these_fields_are_identified_as_radial_topological_defects': False,
                  'vacuum_or_defect_sector_selected': False,
                  'constructive_interacting_relativistic_continuum_proved': False},
    }
    if default:
        result['calibration'] = calibrated.parameters()
        result['calibrated_CKM_exact'] = _serialize(calibrated.ckm_matrix())
        if not _zero(mass_kernel(inputs)-calibrated.mass_kernel()):
            raise ArithmeticError('conditional catalogue differs from the existing calibrated action')
        result['exact_same_calibrated_K45'] = True
    return result
