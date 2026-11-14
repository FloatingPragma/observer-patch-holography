"""Bind the reviewed force result and typed consumer contract."""
import json
from pathlib import Path

from force_measurement import sha

HERE = Path(__file__).resolve().parent
result = json.loads((HERE/'RESULT.json').read_text())
review = json.loads((HERE/'REVIEW.json').read_text())
assert review['result_sha256'] == sha(HERE/'RESULT.json')
assert all(result['checks'].values()) and all(review['mutation_rejections'].values())
primary = result['measurements'][3]
contract = {
    'schema': 'oph.classical_native_force.contract.v1',
    'source_identity': 'frozen-reciprocal-two-mobile-literal-C-side2-96ports',
    'source_manifest_sha256': '4597dc9c226d203d890bc3d3a929bbe1be80e2dc57ec03898c77da4fe5b954e1',
    'field_manifest_sha256': '612aedd40f13faea2b80b5c31d4afbb19cb2332ff22db5b459ac2ae936879008',
    'source_domain': {'copies': 2, 'ports_per_copy': 96, 'seams_per_copy': 288,
                      'joint_position_states': 9216, 'quantum': False, 'C10_source_match': False},
    'preparation': {'positive_copy': {'0': '1/2', '1': '1/2'},
                    'negative_copy': 9, 'preparation_label_retained': True,
                    'reference': 'known initial pair and current-record reconstruction'},
    'supplied_inputs': {'beta': 4, 'kappa_primary': 1, 'nu_primary': 1,
                        'field_rounds': 512, 'energy': 'kappa*rho^T G_K rho/2',
                        'scheduled_rate': 'nu*exp(-beta*Delta_U/2)',
                        'native_coin_probability': '1/2',
                        'signs_energy_controller_are_architecture': True},
    'clock': {'parameter': 'one common continuous Poisson source time t',
              'nu_units': 'inverse source time per baseline seam',
              'uncoupled_total_attempt_rate': '576*nu',
              'uncoupled_total_transfer_rate': '6*nu',
              'controller_arithmetic_rounds': 'separate work count, no supplied latency',
              'physical_seconds': None, 'C10_clock_transfer': False},
    'instrument': {'horizon': '1/1000', 'attempt_cap': 1,
                   'outcomes_per_preparation': 1153,
                   'outcome': 'preparation, no_event OR copy/seam/coin/current/terminal_positions/cap_exhausted',
                   'cap_probability_retained': primary['cap_probabilities'],
                   'waiting_law': 'no_event exp(-Lambda*T); mark exp(-Lambda*t)*r_e/2 dt for 0<t<T',
                   'success_conditioning': False},
    'observable': {'seam': [0, 1], 'score': 'j_e/g_s, g_s=(1-exp(-Lambda_s*T))/Lambda_s; zero on no event',
                   'score_absolute_bound': primary['score_absolute_bound'],
                   'force': 'F=U_K(0,y)-U_K(1,y), declared unit seam length',
                   'drift': 'J=E(score)=(nu/2)*sinh(beta*F/2)',
                   'linear_response': 'J=beta*D0*F+R, D0=nu/4',
                   'remainder_bound': 'nu*abs(beta*F/2)^3*cosh(abs(beta*F/2))/12',
                   'measured_force': primary['force_minus_energy_difference'],
                   'measured_drift': primary['measured_drift_per_source_time']},
    'controlled_field_to_rate_transfer': result['fixedpoint_variant_transfer'],
    'rational_force_score_transfer': review['rational_force_transfer'],
    'rejecting_controls': review['mutation_rejections'],
    'boundaries': {'native_integer_field_energy_equals_real_energy': False,
                   'uniform_mean_current_Green_identity_used_under_feedback': False,
                   'static_field_refinement_implies_dynamic_force_refinement': False,
                   'mechanical_acceleration_or_inertia': False, 'Maxwell_or_Lorentz_limit': False,
                   'alpha_derived': False},
    'artifact_sha256': {name: sha(HERE/name) for name in ('RESULT.json', 'FIRST_EVENT_OUTCOMES.npz', 'REVIEW.json', 'PROOF.md')},
    'official_rows_changed': []}
(HERE/'FORCE_CONTRACT.json').write_text(json.dumps(contract, indent=2)+'\n')
names = ['README.md', 'PROOF.md', 'force_measurement.py', 'verify.py', 'finalize.py',
         'RESULT.json', 'FIRST_EVENT_OUTCOMES.npz', 'REVIEW.json', 'FORCE_CONTRACT.json']
manifest = {'schema': 'oph.classical_native_force.manifest.v1',
            'files': {name: sha(HERE/name) for name in names},
            'scope': 'Finite classical native seam force, complete cap-one law, bounded field-controller variant; no C10 transfer.',
            'official_rows_changed': []}
(HERE/'MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
print(json.dumps({'manifest_sha256': sha(HERE/'MANIFEST.json'),
                  'contract_sha256': sha(HERE/'FORCE_CONTRACT.json')}))
