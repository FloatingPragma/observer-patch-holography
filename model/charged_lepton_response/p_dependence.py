"""Bounded evaluation of retained capacity-dependent source laws, without fits."""
from fractions import Fraction as F
import hashlib
from pathlib import Path

from tools import strict_json
from model.native_flavor_orbit.intervals import I, pi_interval, sqrt_rational, rounded_out
from model.mcpr_native_memory_bridge import construction, certificate as shape
from model.mcpr_native_memory_bridge import packet as bridge_owner
from model.pixel_mass_closure import packet as closure_owner

ROOT = construction.ROOT
REFERENCE = 'model/charged_lepton_response/calibration_reference.json'
SOURCE_INPUTS = ('model/mcpr_native_memory_bridge/results.json',
                 'model/pixel_mass_closure/results.json', construction.FORWARD,
                 construction.MEASURED, construction.MCPR, REFERENCE)


def _interval(value):
    if isinstance(value, I):
        return I(F(value.lo.numerator, value.lo.denominator),
                 F(value.hi.numerator, value.hi.denominator))
    if type(value) not in (str, int, F):
        raise ValueError('exact rational source coordinates required')
    return I(value)


def evaluate(pixel, unified_coupling):
    """Enclose the source shape throughout a rectangular coordinate box."""
    p, u = _interval(pixel), _interval(unified_coupling)
    if not 1 < p.lo <= p.hi < 2 or not 0 < u.lo <= u.hi < 1:
        raise ValueError('the stated source domain is 1<P<2, 0<u<1')
    rounded = lambda v: rounded_out(v, 28)
    phi = (1+sqrt_rational(5))/2
    chi = rounded((u*u/512)/(1-u/77))
    zeta = rounded(((p-phi)**2/(21*pi_interval()))/(1-u/27-u*u/135))
    delta = rounded(F(2, 9)+zeta)
    rho = rounded(sqrt_rational(2)*shape.exponential(-chi))
    roots = [rounded(1+rho*shape.cosine(rounded(delta+2*pi_interval()*j/3)))
             for j in range(3)]
    roots.sort(key=lambda x: x.lo)
    if roots[0].lo <= 0 or any(a.hi >= b.lo for a, b in zip(roots, roots[1:])):
        raise ValueError('three ordered positive roots are not resolved')
    squares = [rounded(x*x) for x in roots]
    ratios = [I(1)]+[rounded(v/squares[0]) for v in squares[1:]]
    q = rounded((1+shape.exponential(-2*chi))/3)
    if q.hi >= F(2, 3):
        raise ValueError('strict deformation from balanced Koide is unresolved')
    return {'P': p, 'u': u, 'chi': chi, 'zeta': zeta, 'delta': delta,
            'rho': rho, 'positive_roots': roots, 'shape_squares': squares,
            'charged_action_ratios': ratios, 'positive_root_Koide_ratio': q}


def comparison(evaluation, reference):
    """Read observations only after the source calculation; no parameter update."""
    result = {}
    for index, name in ((1, 'muon_electron_ratio'), (2, 'tau_electron_ratio')):
        row = reference[name]
        center, error = F(row['value']), F(row['standard_uncertainty'])
        if center <= 0 or error <= 0:
            raise ValueError('positive reference and standard uncertainty required')
        value = evaluation['charged_action_ratios'][index]
        result[name] = {
            'source_minus_reference': rounded_out(value-center, 14),
            'relative_residual_ppm': rounded_out((value/center-1)*10**6, 10),
            'residual_over_quoted_standard_uncertainty': rounded_out((value-center)/error, 10),
            'source_interval_inside_reference_one_sigma':
                center-error <= value.lo <= value.hi <= center+error,
            'source_interval_disjoint_from_reference_one_sigma':
                value.hi < center-error or value.lo > center+error}
    return result


def authenticate(path, expected_bindings):
    """Authenticate a retained receipt and every complete official source pin."""
    value = strict_json.load(path)
    if value.get('all_pass') is not True or value.get('source_bindings') != expected_bindings:
        raise ValueError('owning receipt is incomplete or not admitted')
    for name, digest in expected_bindings.items():
        source = (ROOT/name).resolve()
        if not source.is_relative_to(ROOT) or hashlib.sha256(source.read_bytes()).hexdigest() != digest:
            raise ValueError('stale owning source: '+name)
    return value


def _box(row):
    return I(row['lo'], row['hi'])


def _inverse_count(memory):
    count = len(memory['inverse_branches'])
    if (count != 2 or memory['quadratic_root_count'] != 2
            or memory['cubic_root_count'] != 3
            or memory['all_complete_record_probabilities_retained'] is not True
            or memory['root_domain'] != 'entire certified implicit-coordinate enclosure'):
        raise ValueError('complete root-box memory inverse is not certified')
    return count


def _json(value):
    if isinstance(value, I):
        return value.strings()
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    return value


def certificate():
    bridge_pins, closure_pins = bridge_owner.source_bindings(), closure_owner.source_bindings()
    bridge = authenticate(bridge_owner.OUTPUT, bridge_pins)
    closure = authenticate(closure_owner.OUTPUT, closure_pins)
    reference_path = ROOT/REFERENCE
    reference = strict_json.load(reference_path)
    rows = []
    for branch in construction.BRANCHES:
        p, u = construction.source(branch)
        retained = [x for x in bridge['designs'] if x['source_branch'] == branch]
        if (len(retained) != 2 or any(type(x['inverse_branch']) is not int for x in retained)
                or {x['inverse_branch'] for x in retained} != {0, 1}
                or any((x['P'], x['alpha_U']) != (p, u) for x in retained)):
            raise ValueError('declared endpoint differs from the retained complete inverse')
        value = evaluate(p, u)
        # The owning rational shape calculation is independent of this readout.
        old = shape.source_shape_boxes(branch)
        if any(max(a.lo, b.lo) > min(a.hi, b.hi)
               for a, b in zip(value['shape_squares'], old)):
            raise ValueError('source law disagrees with the owning rational calculation')
        rows.append({'branch': branch,
                     'ancestry': ('supplied forward decimal coordinates; no complete endpoint selection theorem'
                                  if branch == 'forward' else
                                  'supplied endpoint inferred from measured alpha; no charged masses fitted here'),
                     'source': value, 'comparison': comparison(value, reference),
                     'optional_electron_unit_mass_intervals_MeV':
                         [rounded_out(v*F(reference['electron_MeV']['value']), 16)
                          for v in value['charged_action_ratios']],
                     'admitted_complete_memory_inverse_branches': len(retained)})
    for root in closure['closures']:
        if not root['strict_self_map'] or not root['conditional_unique_root_in_domain']:
            raise ValueError('conditional closure root is not certified')
        memory = root['memory_admissibility']
        inverse_count = _inverse_count(memory)
        p = _box(root['fixed_point_P'])
        u = _box(root['implicit_and_mass_enclosures']['alpha_u'])
        value = evaluate(p, u)
        masses = [_box(x) for x in root['implicit_and_mass_enclosures']['masses_over_v']]
        for actual, old in zip(value['charged_action_ratios'], [I(1)]+[m/masses[0] for m in masses[1:]]):
            if max(actual.lo, old.lo) > min(actual.hi, old.hi):
                raise ValueError('capacity readout disagrees with the retained mass boxes')
        rows.append({'branch': root['design'],
                     'ancestry': 'unique root within the admitted interval of the declared census, scale and response law',
                     'source': value, 'comparison': comparison(value, reference),
                     'optional_electron_unit_mass_intervals_MeV':
                         [rounded_out(v*F(reference['electron_MeV']['value']), 16)
                          for v in value['charged_action_ratios']],
                     'admitted_complete_memory_inverse_branches': inverse_count,
                     'closure_alpha_inverse': root['fixed_point_alpha_inv'],
                     'v_over_declared_E_star': root['implicit_and_mass_enclosures']['v']})
    if bridge_owner.source_bindings() != bridge_pins or closure_owner.source_bindings() != closure_pins:
        raise RuntimeError('owning sources changed during bounded evaluation')
    receipt_paths = (bridge_owner.OUTPUT, closure_owner.OUTPUT, reference_path)
    return _json({'schema': 'oph.charged_lepton_response.capacity_dependence.v1',
                  'authenticated_receipts_and_reference': {
                      str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in receipt_paths},
                  'authenticated_source_pin_counts': [len(bridge_pins), len(closure_pins)],
                  'branches': rows,
                  'source_law': 'MCPR declared response; delta=2/9+zeta, rho=sqrt(2) exp(-chi)',
                  'source_parameters_fitted_by_this_certificate': False,
                  'optional_MeV_values_use_the_supplied_electron_unit': True,
                  'optional_MeV_unit_calibration_is_source_derived': False,
                  'experimental_values_consumed_by_source_evaluation': False,
                  'cold_root_or_native_memory_producers_replayed_here': False,
                  'existing_receipts_and_source_bindings_authenticated': True,
                  'primitive_finite_assumptions_force_this_capacity_law': False,
                  'inverse_memory_existence_implies_source_selection': False,
                  'physical_interacting_mass_predictions': [],
                  'joint_statistical_or_theory_error_claim': False})
