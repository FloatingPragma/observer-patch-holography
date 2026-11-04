"""Receipt and independent full-record checks for the chosen memory design."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from . import construction as c
from .certificate import admissibility

ROOT = c.ROOT
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENT = 'model/native_flavor_orbit/results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    inherited = strict_json.load(ROOT/PARENT)['source_bindings']
    for path,digest in inherited.items():
        if not (ROOT/path).resolve().is_relative_to(ROOT.resolve()):
            raise ValueError('source escapes repository')
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError('stale inherited source: '+path)
    paths = set(inherited)|{PARENT,c.FORWARD,c.MEASURED,c.MCPR}
    paths.update(str(p.relative_to(ROOT)) for p in PACKAGE.iterdir()
                 if p.is_file() and p.suffix in ('.py','.md'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def _serialize(row):
    return {key:(value.tolist() if isinstance(value,np.ndarray) else
                 list(value) if isinstance(value,tuple) else value)
            for key,value in row.items()}


def independent_check(row):
    """Check a supplied design through full-record dynamics and the old law.

    This function does not call the inverse-design routine. The twelve-input
    stationary solve is independent of the compact inverse coefficient.
    """
    expected = c.target(row['source_branch'])
    if (row['P'],row['alpha_U']) != c.source(row['source_branch']):
        raise ValueError('source coordinates changed')
    p,q = row['p'],row['q']
    if not 0 < p < 1/12 or not 0 < q < 1/6:
        raise ValueError('native gate domain changed')
    matrix = np.asarray(row['mass_matrix_over_v'],float)
    coefficient = np.asarray(row['coefficient_matrix'],float)
    normalized = np.asarray(row['normalized_matrix'],float)
    if any(x.shape != (3,3) or not np.isfinite(x).all() for x in (matrix,coefficient,normalized)):
        raise ValueError('complete three-family matrices required')
    direct = c.direct_memory(row)
    old = np.asarray([float(x) for x in expected['masses_over_v']])
    errors = {
        'full_record_coefficient':float(np.linalg.norm(direct['coefficient_matrix']-coefficient,2)),
        'unconditional_duty_trace':abs(float(np.trace(coefficient))-p),
        'normalized_matrix':float(np.linalg.norm(coefficient-p*normalized,2)),
        'feedback_amplitude_identity':float(np.linalg.norm(matrix-row['gain_over_v']*coefficient,2)),
        'source_spectrum':float(np.max(np.abs(np.linalg.eigvalsh(matrix)-old))),
        'source_gain':abs(row['gain_over_v']-sum(old)/p),
        'left_gate':abs(np.sin(row['theta_left']*np.sqrt(9/8))**2/12-p),
        'right_gate':abs(np.sin(row['theta_right']*np.sqrt(3))**2/6-q),
        'full_stationarity':float(np.max(np.abs(direct['transition']@direct['stationary']-direct['stationary']))),
        'all_input_normalization':float(np.max(np.abs(direct['transition'].sum(axis=0)-1))),
        'complete_current_instruments':max(c.instrument_errors(row)),
    }
    if not all(np.isfinite(x) for x in errors.values()) or max(errors.values()) > 2e-12:
        raise ValueError('independent native-memory identity failed: '+str(errors))
    if np.min(direct['stationary']) < -2e-13:
        raise ValueError('negative stationary probability')
    # Display numbers are checked against source values too; no hidden target.
    if not np.allclose(row['masses_over_v'],old,rtol=2e-14,atol=0):
        raise ValueError('source mass display changed')
    return {'all_checks_pass':True,'numerical_tolerance':2e-12,
            'max_residual_below_tolerance':True,
            'retained_binary_native_records_per_input':512,
            'retained_old_memory_inputs':12,'current_outcomes_per_leg':13,
            'fine_current_native_records_per_input':21632}


def report():
    pins = source_bindings()
    rows = []
    for branch in c.BRANCHES:
        for inverse in (0,1):
            row = _serialize(c.design(branch,inverse))
            row['independent_check'] = independent_check(row)
            rows.append(row)
    result = {'schema':'oph.mcpr_native_memory_bridge.v1','all_pass':True,
              'issues':[1141,1142,1143],
              'generic_algebra':c.exact_identities(),
              'weak_pulse_limit':{
                  'generator':'X tensor (gain_over_v * Z_C)',
                  'arbitrary_CQ_reference_trace_error':'4*g*T/(N*(1-p))+(8/(1-p)+4)*g^2*T^2/N',
                  'stationary_product_preparation':'the first term is zero',
                  'domain':'finite untouched reference, N>=1, T>=0, 0<=p<1, g>=0',
                  'finite_step_channel_equals_limit_unitary':False},
              'admissibility':{branch:admissibility(branch) for branch in c.BRANCHES},
              'designs':rows,
              'source_design':{
                  'repair_law':'unchanged native C/E repair',
                  'memory_law':'same twelve-input reset2 half-turn controller',
                  'chosen_parameters':'independent current durations and C-duty feedback gain from MCPR source law',
                  'gain':'sum(m_i/v)/p, applied to every enabled pulse without conditioning',
                  'family_frame':'canonical kinetic I3; unitary diagonalization commutes with all native gauge currents',
                  'full_operator_identity':'gain*Z_C is unitarily equivalent to diag(m_i/v)',
                  'source_endpoints':'supplied decimal coordinates; forward branch has no full physical endpoint proof',
                  'measured_endpoint_control':'uses the endpoint inferred from measured alpha, not a forward prediction'},
              'scope':{
                  'old_P_mass_law_realized_in_declared_native_mean_feedback':True,
                  'every_native_record_and_memory_input_retained':True,
                  'algebraic_inverse_and_admissibility_certified':True,
                  'correlated_C_duty_weak_pulse_generator_derived':True,
                  'runtime_measured_charged_masses_consumed':False,
                  'MCPR_design_historically_target_informed':True,
                  'unique_selection_by_bare_OPH_proved':False,
                  'static_mean_equals_finite_correlated_channel':False,
                  'compiled_feedback_service_in_native_pulse_words':False,
                  'physical_interacting_rest_mass_or_clock_derived':False,
                  'all_quark_and_neutrino_parameters_derived':False,
                  'particle_equals_created_topological_defect_proved':False,
                  'continuum_or_Einstein_limit_proved':False,
                  'integrated_gate_accepted':False},
              'canonical_physical_premises_discharged':[],
              'source_bindings':pins}
    if pins != source_bindings():
        raise ValueError('source changed during calculation')
    return result


def validate(value):
    if value.get('source_bindings') != source_bindings():
        raise ValueError('stale or incomplete source manifest')
    for row in value.get('designs',[]):
        independent_check(row)
    if canonical(value) != canonical(report()):
        raise ValueError('receipt differs from full replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',nargs='?',const=str(OUTPUT))
    group.add_argument('--check',nargs='?',const=str(OUTPUT))
    args = parser.parse_args()
    if args.write:
        Path(args.write).write_text(canonical(report()),encoding='utf-8')
        print('wrote '+args.write)
    else:
        validate(strict_json.load(Path(args.check)))
        print('MCPR native-memory source and full-record replay passed')


if __name__ == '__main__':
    main()
