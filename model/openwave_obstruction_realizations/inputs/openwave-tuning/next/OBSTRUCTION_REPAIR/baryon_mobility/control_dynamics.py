"""Literal three-quark CAR control and all four placement/flux outcomes.

This g=0,m=0 control does not replace the full shared-auxiliary vacuum.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import cmath
import json
import math
import numpy as np
from scipy.linalg import expm
from run import annihilate, create

HERE = Path(__file__).resolve().parent


def hopping(vector, forward_only=False):
    output = defaultdict(int)
    directions = [(i, i+12) for i in range(12)]
    if not forward_only:
        directions += [(j, i) for i, j in directions]
    for word, coefficient in vector.items():
        for target, origin in directions:
            middle, a = annihilate(origin, word)
            if middle is None:
                continue
            end, b = create(target, middle)
            if end is not None:
                output[end] += coefficient*a*b
    return {word: c for word, c in output.items() if c}


def main():
    initial = sum(1 << i for i in (12, 13, 14))
    final = sum(1 << i for i in (0, 1, 2))
    frame = [{initial: 1}]
    for _ in range(3):
        frame.append(hopping(frame[-1], True))
    assert frame[-1] == {final: 6}
    norms = [sum(c*c for c in v.values()) for v in frame]
    assert norms == [1, 3, 12, 36]
    for n, v in enumerate(frame):
        expected = defaultdict(int)
        if n < 3:
            for w, c in frame[n+1].items(): expected[w] += c
        if n:
            for w, c in frame[n-1].items(): expected[w] += n*(4-n)*c
        assert hopping(v) == dict(expected)
    # A two-baryon fixed-frame control detects the inter-cell fermion sign
    # that the empty-reference one-baryon test alone cannot detect.
    pair_start = 7 | (56 << 12)
    pair = {pair_start: 1}
    for _ in range(3):
        pair = hopping(pair, True)
    actual_pair_element = pair.get(63, 0)
    local_state, local_sign = 7, 1
    for mode in (5, 4, 3):
        local_state, sign = create(mode, local_state)
        local_sign *= sign
    assert local_state == 63 and actual_pair_element == 6
    commuting_tensor_element = 6*local_sign
    assert commuting_tensor_element == -6 != actual_pair_element
    result = []
    t = 1/8
    max_error = 0.0
    for kap in (1/2, 1, 2):
        Delta = 4*kap/3
        H = np.diag([0., Delta, Delta, 0.])
        for n in range(3):
            H[n, n+1] = H[n+1, n] = -t*math.sqrt((n+1)*(3-n))
        third = 6*t**3/Delta**2
        times = [0., 1., 8., 32., 128., math.pi/(2*third)]
        for tau in times:
            vector = expm(-1j*tau*H)[:, 0]
            probabilities = list(map(float, abs(vector)**2))
            assert abs(sum(probabilities)-1) < 1e-12
            amplitudes = []
            for middle in (Delta-2*t, Delta+2*t):
                frequency = math.sqrt(middle*middle+12*t*t)
                amplitudes.append(cmath.exp(-.5j*middle*tau)*
                    (math.cos(.5*frequency*tau)+1j*middle/frequency*math.sin(.5*frequency*tau)))
            transferred = abs((amplitudes[0]-amplitudes[1])/2)**2
            held = abs((amplitudes[0]+amplitudes[1])/2)**2
            error = max(abs(transferred-probabilities[3]), abs(held-probabilities[0]))
            assert error < 1e-12
            max_error = max(max_error, error)
            result.append({'kappa': kap, 'source_time': tau,
                'left_quark_count_outcomes': [0, 1, 2, 3], 'probabilities': probabilities,
                'complete_colour_flux_outcomes_retained': True,
                'source_derived_third_order_half_transfer_time': math.pi/(2*third)})
    out = {'schema': 'oph.baryon_mobility.full_CAR_control.v1',
           'raw24_CAR_triple_transfer_coefficient': 6,
           'fixed_frame_two_baryon_matrix_element': actual_pair_element,
           'dropped_parity_control_matrix_element': commuting_tensor_element,
           'dropping_parity_rejected_on_fixed_preparation_and_probe': True,
           'invariant_four_state_frame_norms': norms,
           'full_Hamiltonian_control_mass': 0, 'full_Hamiltonian_control_auxiliary_coupling': 0,
           'actual_half_filled_auxiliary_vacuum': False,
           'native_program_execution': False, 'all_placement_outcomes': result,
           'independent_even_odd_amplitude_max_error': max_error,
           'physical_halo_lifetime': False}
    (HERE/'CONTROL_DYNAMICS.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'raw_triple_coefficient': 6, 'complete_probability_laws': len(result),
                      'max_amplitude_comparison_error': max_error}), flush=True)


if __name__ == '__main__':
    main()
