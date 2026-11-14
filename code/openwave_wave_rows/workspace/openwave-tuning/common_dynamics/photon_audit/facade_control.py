"""Adversarial fixed-data replay without any eigensolver availability."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys

sys.dont_write_bytecode = True
import numpy as np
import frozen

HERE = Path(__file__).resolve().parent


def forbidden(*args, **kwargs):
    raise AssertionError('frozen replay called a forbidden eigensolver')


def run():
    before = frozen.verify_inputs()
    original = {name:getattr(np.linalg, name) for name in ('eigh', 'eigvalsh', 'eig')}
    try:
        for name in original: setattr(np.linalg, name, forbidden)
        first = frozen.spectral_certificate()
        initial = first['frame_integers'][0][0]
        first['frame_integers'][0][0] += 1
        assert frozen.spectral_certificate()['frame_integers'][0][0] == initial
        program, bound, phases, vacuum = frozen.word('X')
        target = program[0]['targets'][0]
        program[0]['targets'][0] = 999
        assert frozen.word('X')[0][0]['targets'][0] == target
        diagonal, phase_bound, phase_values, phase_vacuum = frozen.photon_phase()
        assert np.array_equal(phases, phase_values) and vacuum.json() == phase_vacuum.json()
        assert len(phases) == 256 and len(diagonal) > 0 and phase_bound > 0
        assert len(frozen.phase_targets(F(1, 2))[0]) == 256
        rejected = 0
        for value in (1., F(1, 2), F(2)):
            try: frozen.photon_phase(value)
            except ValueError: rejected += 1
        assert rejected == 3
    finally:
        for name, value in original.items(): setattr(np.linalg, name, value)
    assert before == frozen.verify_inputs()
    result = {'schema':'oph.quadratic_photon.frozen_facade_control.v1',
              'checks':{'all_eigensolvers_disabled_during_replay':True,
                        'caller_frame_mutation_isolated':True, 'caller_word_mutation_isolated':True,
                        'frozen_phase_word_and_targets_agree':True,
                        'analytic_time_does_not_authorize_native_time':True,
                        'all_input_bytes_unchanged':True},
              'provider_manifest':frozen.MANIFEST_SHA256, 'native_time_rejections':rejected,
              'diagonal_pulses':len(diagonal), 'official_rows_changed':[]}
    (HERE/'FACADE_CONTROL.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__': run()
