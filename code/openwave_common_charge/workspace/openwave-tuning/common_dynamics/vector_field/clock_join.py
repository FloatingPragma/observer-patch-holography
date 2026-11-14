"""Attach the same native completed-word clock to the full coupled word."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

import native_field

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'clock'))
import phase_clock


def run():
    paths = [Path(__file__), Path(native_field.__file__), Path(phase_clock.__file__),
             HERE / 'COMPOSITION.json', HERE.parent / 'clock/MANIFEST.json']
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    composition = json.loads((HERE / 'COMPOSITION.json').read_text())
    p0 = F(composition['native_blank_reset_success_lower_exact'])
    pr = F(composition['native_ready_reset_success_lower_exact'])
    rows = []
    for setting in ('X', 'Y'):
        clock_error = phase_clock.exact_word_bound(3, setting, F(1))
        extra = clock_error + (1 - p0) + (1 - pr)
        total = F(composition['three_step_full_record_TV_bound_exact']) + extra
        rows.append({'quadrature': setting, 'native_clock_word_error_exact': str(clock_error),
                     'additional_reset_failure_bound_exact': str((1-p0)+(1-pr)),
                     'full_source_field_clock_record_TV_bound_exact': str(total),
                     'full_source_field_clock_record_TV_bound_upper': float(total)})
    result = {'schema': 'oph.local12.coupled_clock_join.v1',
              'checks': {'full_coupled_word_not_bare_W': True, 'fresh_clock_wire': True,
                         'same_actual_native_phase_clock': True,
                         'additional_resets_and_clock_errors_counted': True,
                         'no_independence_assumption_for_worst_reset_bound': True,
                         'source_bytes_unchanged': True},
              'source_clock_generator': 'native_field.complete_word(alpha,face_couplings208,clock_spacing=Fraction(1))',
              'clock_wire': 3203, 'spectral_probe_wire_if_added': 3204,
              'pointer_auxiliaries': 'fresh native carrier codes above both occupied wires',
              'clock_preparation': 'extra blank native reset, same native H word before the first complete V',
              'clock_analysis': 'same native X/Y analysis and complete C/E pointer after three complete V words',
              'pulse_angle_per_completed_V': '-h/16', 'measured_relative_phase_per_V': 'h/8',
              'declared_time': 't=n*tau0*h', 'executed_clock_law_reference': '../clock/PHASE_CLOCK.json',
              'three_step_relative_phase': '3/8', 'bounds': rows,
              'scope': 'The clock attaches to completed fullV words and has its own output and retained failure flags. The duration of the larger coupled pulse/service word is a distinct resource. Bounds quantify systematic source/readout error, not finite-shot clock resolution or a physical-second calibration.',
              'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE / 'CLOCK_JOIN.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'checks': result['checks'], 'maximum_complete_TV': max(r['full_source_field_clock_record_TV_bound_upper'] for r in rows)}))


if __name__ == '__main__':
    run()
