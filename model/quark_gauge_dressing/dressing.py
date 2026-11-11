"""Positive electric energy from exact Gauss constraints, with native probes.

The complete joined source occupies one vertex of a separate committed
finite cone. This is a charge-island extension, not local gauging of the
343 internal sites. The field is its longitudinal electrostatic sector.
"""
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from model.affine_gauss_source import model as gauge
from model.defect_mass_join.source import car, one_cell_polynomial, representation
from model.particle_common_contract.provenance import import_closure
from . import symmetry

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = 'model/native_quark_response/results.json'


def serialize(value):
    if isinstance(value, s.Basic):
        return str(value)
    if isinstance(value, dict):
        return {k: serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(v) for v in value]
    return value


@lru_cache(maxsize=1)
def geometry():
    _, _, d, c, _ = gauge.cone.load_geometry()
    maps = gauge.cone.cone_maps(d, c)
    gradient, curl = maps['D'], maps['C']
    laplacian = gradient.T*gradient
    green = (laplacian+s.ones(13)).inv()-s.ones(13)/169
    assert gradient.shape == (42, 13) and curl.shape == (50, 42)
    assert gradient.rank() == 12 and curl*gradient == s.zeros(50, 13)
    assert laplacian*green == s.eye(13)-s.ones(13)/13
    assert green*s.ones(13, 1) == s.zeros(13, 1)
    rows = []
    for site in range(13):
        charge = s.eye(13)[:, site]-s.eye(13)[:, gauge.REFERENCE_VERTEX]
        field = gradient*green*charge
        coefficient = (field.T*field)[0]/2
        assert gradient.T*field == charge and curl*field == s.zeros(50, 1)
        assert coefficient == (charge.T*green*charge)[0]/2
        assert coefficient >= 0 and (coefficient > 0) == (site != gauge.REFERENCE_VERTEX)
        rows.append({'source_vertex': site, 'reference_vertex': gauge.REFERENCE_VERTEX,
                     'unit_charge_vector': list(charge), 'unit_longitudinal_field': list(field),
                     'positive_electric_coefficient': coefficient})
    return rows


def validate_parent():
    value = json.loads((ROOT/PARENT).read_text())
    for name, expected in value['input_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('native quark source ancestry differs: '+name)
    assert value['same_source']['local'] == 'full native thirty-CAR plus five-state scalar, all Yukawa coefficients one'
    return value


def shifted_records(site=0, coupling=s.Integer(1)):
    if type(site) is not int or site not in range(13) or site == gauge.REFERENCE_VERTEX:
        raise ValueError('one declared cone vertex separated from the reference required')
    coupling = s.sympify(coupling)
    if coupling.has(s.Float) or coupling.is_positive is not True:
        raise ValueError('one exact positive common electric coefficient required')
    parent = validate_parent()
    coefficient = geometry()[site]['positive_electric_coefficient']*coupling
    charges = {'uc': s.Rational(2, 3), 'dc': s.Rational(-1, 3), 'Q': s.Rational(-1, 6)}
    rows = {}
    for family, charge in charges.items():
        inherited = parent['exact_local_response'][family]
        shift = coefficient*charge**2
        records = []
        for record in inherited['complete_positive_records']:
            records.append({'target_energy_interval': [s.Rational(x)+shift for x in record['source_energy_interval']],
                'energy_transfer_interval': [s.Rational(x)+shift for x in record['source_excitation_interval']],
                'original_probability_interval': record['probability_interval'],
                'exactly_zero_residue': record['exactly_zero_residue']})
        assert len(records) == s.Poly(inherited['complete_target_annihilator'], s.Symbol('x')).degree()
        for original, shifted in zip(inherited['complete_positive_records'], records, strict=True):
            assert shifted['original_probability_interval'] == original['probability_interval']
            assert shifted['exactly_zero_residue'] == original['exactly_zero_residue']
        rows[family] = {'source_hypercharge_after_removal': charge, 'reference_hypercharge': -charge,
            'positive_electric_energy': shift, 'complete_energy_records': records,
            'pump_success': inherited['pump_success'], 'pump_failure': inherited['pump_failure'],
            'outside_support_probability': inherited['outside_support'],
            'reference_intrinsic_energy': '0',
            'parent_reference_cost_four_imported': False,
            'complete_target_annihilator': inherited['complete_target_annihilator'],
            'dressed_resolvent': 'R_dressed(z)=R_native(z-electric_energy)',
            'source_reference_prep': 'same source physical ground doublet, retained singlet purification; additional compensator blank and flux zero'}
    assert parent['computed_equal_coupling_uc_dc_resolvent_equality'] is True
    assert rows['uc']['positive_electric_energy']-rows['dc']['positive_electric_energy'] == coefficient/3
    return {'source_vertex': site, 'common_electric_coupling': coupling,
            'geometry_coefficient': geometry()[site]['positive_electric_coefficient'],
            'uc_minus_dc_corresponding_line_shift': coefficient/3, 'channels': rows}


def column_charge_controls():
    h = one_cell_polynomial()
    y = car.gauge_current(11)
    assert not car.commutator(y, h)
    weights = representation()[2][11].diagonal()
    for mode in range(30):
        c = car.monomial(((mode, False),))
        assert car.commutator(y, c) == car.scale(c, -weights[mode % 15])
    return {'all_thirty_elementary_charge_identities_exact': True,
        'full_native_quantum_H_commutes_with_Y': True,
        'full_spatial_extension': 'h is identity on native labels, Omega is identity on doublet indices, local Ward identity holds at every site; hence [H,Y_total]=0',
        'spectral_shift_proof': '[H,Y]=0 and Y a=-y a on a Y=0 reference imply (H+g C Y²) a P0=(H+g C y²) a P0; every spectral projector and zero residue translates with unchanged weight',
        'Gauss_isometry': 'J|psi_q>_source|ref_-q> = |psi_q>|ref_-q>|ell=q>; G_v J=0, J† H_electric J=g C Y²',
        'finite_flux_register': 'ell=k/6, |k|<=36V+3; every charge eigenvalue of CAR(30V) and the global scalar cap is retained',
        'electric_operator': 'E_e=(D G (delta_source-delta_reference))_e ell; H_electric=g sum_e E_e²/2 >=0',
        'reference_carrier': 'blank plus full direct/dual source representation sectors sufficient for every charge; intrinsic Hamiltonian zero; all charge and nonabelian probe coordinates retained',
        'physical_constraints': 'G_v=(D^T E)_v-delta_v,source Y_source-delta_v,reference Y_reference=0 for all13vertices',
        'Gauss_preserving_probe': 'native charged probe with opposite quantum reference and flux-register shift by its exact hypercharge',
        'clock': 'exp(-it(H_native tensor I+H_electric)), one hbar=1 clock',
        'complete_observer_records': 'every shifted energy with original residue, exact zero residues, pump failure and outside-support outcome retained',
        'positive_electric_energy_implies_full_new_physical_ground_proved': False,
        'selected_neutral_source_reference_stationary_and_energy_unchanged': True}


def bindings():
    seeds = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py', '.md')}
    seeds.update({PARENT, 'model/defect_mass_join/PROOF.md', 'model/affine_gauss_source/PROOF.md',
        'code/electromagnetism/cone_whitney_bridge.py', 'code/electromagnetism/discrete_coulomb_green.py'})
    seeds.update(gauge.pins())
    seeds.update(gauge.cone.PIN_PATHS)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def certificate():
    pins = bindings()
    result = serialize({'schema': 'oph.quark_gauge_dressing.v1', 'source_bindings': pins,
        'full_joined_quark_symmetry': symmetry.certificate(),
        'native_charge_and_Gauss_transfer': column_charge_controls(),
        'complete_cone_electric_geometry': geometry(),
        'literal_one_cell_complete_energy_example': shifted_records(),
        'scope': {'positive_Gauss_electric_energy_derived': True,
            'full_native_source_as_one_charge_island': True,
            'compact_343_internal_sites_locally_gauged': False,
            'transverse_photon_oscillator_added': False,
            'scalar_or_Yukawa_coefficients_fitted': False,
            'one_common_primitive_electric_coefficient': '1',
            'empirical_inputs': [], 'source_has_dynamic_fermion_propagation': 'the joined native law retains its own h; this electrostatic extension supplies no spatial gauge hopping',
            'complete_343_interacting_ground_computed': False,
            'physical_rest_mass_or_GeV_output': False,
            'observed_generation_count_derived': False}})
    if bindings() != pins:
        raise ValueError('gauge dressing source changed during calculation')
    return result
