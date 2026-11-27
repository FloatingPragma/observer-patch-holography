"""Calibrated native nondemolition signed sum on the two far receiver rails."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from model.executable_model.carrier import federation
from model.experiments.material_inertia.pair_spectrum.native_source_probe import physical
from model.experiments.material_inertia.pair_spectrum.native_source_probe.measure import routes
from model.experiments.interval_echo_consensus.signed_effects import fixed_pointer as pointer
from ...row_20_electric import experiment as prior


def calibration():
    ids = [[int(np.array((a,b,1,2,0)) @ [81,27,9,3,1]), int(np.array((a,b,1,0,2)) @ [81,27,9,3,1])]
           for a in range(3) for b in range(3)]
    initial = np.eye(243)[:, [row[0] for row in ids]]
    result = {}
    basis = pointer.neutral_basis(3, 4)
    for setting in ('parity', 'sign'):
        ratio = Fraction(2 if setting == 'parity' else 1, 65)
        word = basis + pointer.diagonal_word(0, 4, ratio) + pointer.diagonal_word(1, 4, ratio)
        if setting == 'sign':
            word += [pointer.pulse('reference', 4, 4, Fraction(1,64))]
        full = physical.primitive(word + pointer.inverse(basis), initial)
        allowed = np.zeros_like(full)
        blocks = []
        for col, targets in enumerate(ids):
            vec = full[targets, col]
            allowed[targets, col] = vec
            blocks.append(vec)
        leakage = float(np.linalg.norm(full-allowed, 2))
        equality = float(np.linalg.norm(blocks[2]-blocks[6]))
        if max(leakage, equality) > 3e-12:
            raise ArithmeticError('native read disturbed an orientation or its phase')
        coherent = []
        for z in (1, -1, 1j, -1j):
            v = (full[:, 2] + z*full[:, 6])/np.sqrt(2)
            expected = np.zeros(243, complex)
            expected[ids[2]] = blocks[2]/np.sqrt(2)
            expected[ids[6]] = z*blocks[2]/np.sqrt(2)
            coherent.append(float(np.linalg.norm(v-expected)))
        result[setting] = dict(head_20_probabilities=[float(abs(z[0])**2) for z in blocks],
                               maximum_column_leakage=leakage, guarded_orientation_amplitude_difference=equality,
                               coherent_test_errors=coherent, native_pulses=len(word)+len(basis))
    board = federation(2)
    _, path, _ = prior.guide(board)
    layouts = []
    for n in (23, 45):
        ends = path[2*n-1:2*n+1]
        geo = routes(board, {p//12 for p in ends})
        paths = [geo['endpoint_transposition_paths'][geo['target_ports'].index(port)] for port in ends]
        # Both target coordinates act at the same physical interaction port on successive restored visits.
        lens = [physical.validate_route(board, route, geo['source_park'], geo['source_cross']) for route in paths]
        layouts.append(dict(n=n, target=ends, head=[geo['source_park'],geo['source_cross']],
                            station=geo['interaction_port'], paths=paths, roundtrip_swap_counts=[2*v for v in lens]))
    return dict(source_order=['00','01','02','10','11','12','20','21','22'],
                settings=result, receiver_routing=layouts,
                scope='Exact columns of a nondemolition local sum; charged replacements and coherent neutral orientations.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.parse_args().output.write_text(json.dumps(calibration(), indent=1, allow_nan=False)+'\n')
