"""An exact unitary history representation of the declared native instrument.

The completion and its Hermitian logarithms are representations at native
event boundaries. Their intermediate parameter is not an observer clock.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
from itertools import product
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer import native
from ..common_source import slow_geometry
from ..observer.scalar_response.experiment import canonical, compare_receipt
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    'model/native_action/history_action.py',
    'model/observer/native.py',
    'model/observer/scalar_response/experiment.py',
    'tools/strict_json.py',
    'model/common_source/slow_geometry.py',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
)


def carrier_embedding():
    """Embed local tensor ports 0,1,2 into an actual native carrier triangle."""
    ports = (0, 1, 5)
    seams = tuple(slow_geometry.source().seams())
    local_seams = ((0, 1), (0, 2), (1, 2))
    if any((ports[a], ports[b]) not in seams for a, b in local_seams):
        raise ArithmeticError('the history support is not a native carrier triangle')
    branches = 0
    for mask in range(8):
        local = tuple((mask >> p) & 1 for p in range(3))
        loads = tuple(local[ports.index(p)] if p in ports else 0 for p in range(12))
        for seam, coin in product(local_seams, (0, 1)):
            expected, record = native.native(local, seam, coin)
            output, source_record = native.native(loads, tuple(ports[p] for p in seam), coin)
            if (tuple(output[p] for p in ports) != expected or source_record != record
                    or any(output[p] for p in range(12) if p not in ports)):
                raise ArithmeticError('carrier restriction changed an endpoint or record')
            branches += 1
    return {'local_tensor_ports': [0, 1, 2], 'native_carrier_ports': list(ports),
            'native_seams': [[ports[a], ports[b]] for a, b in local_seams],
            'endpoint_and_record_restrictions_checked': branches,
            'other_ports_are_untouched': True}


def _domain(law, seam, ports):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('native C or E instrument required')
    if type(ports) is not int or not 2 <= ports <= 4:
        raise ValueError('two through four binary ports required')
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(v) is not int for v in seam)
            or not 0 <= seam[0] < seam[1] < ports):
        raise ValueError('two ordered distinct port addresses required')


def kraus_numerators(law, seam=(0, 1), ports=2):
    """Integer W_alpha with K_alpha=W_alpha/sqrt(2), from native endpoints.

    The label is (coin, visible current, candidate record). Only the current
    is accessible to the control below; the other label entries are traced.
    """
    _domain(law, seam, ports)
    dim = 1 << ports
    groups = {}
    for x in range(dim):
        loads = tuple((x >> p) & 1 for p in range(ports))
        for coin in (0, 1):
            output, (_, current) = native.native(loads, seam, coin)
            y = sum(value << p for p, value in enumerate(output))
            record = current if law == 'C' else loads[seam[0]]-loads[seam[1]]
            label = (coin, current, record)
            groups.setdefault(label, np.zeros((dim, dim), dtype=object))[y, x] = 1
    return tuple(sorted(groups.items()))


def completion_blocks(law):
    """Exact blocks of U=[[D,W/sqrt(2)],[W.T/sqrt(2),0]].

    W.T W=2I, D=I-W W.T/2, D W=0 and D^2=D prove U*=U=U^-1.
    This avoids approximate square-root arithmetic in the certificate.
    """
    groups = kraus_numerators(law)
    w = np.vstack([matrix for _, matrix in groups])
    identity = np.eye(w.shape[0], dtype=object)
    d = identity-F(1, 2)*(w@w.T)
    return groups, w, d


def controlled_history_numerators(law):
    """Two overlapping repairs and local current-controlled feedback.

    The first current chooses the second seam and a phase on port zero.
    Fresh hidden labels are retained for both repairs. No federation is used.
    The returned branch maps are divided by 2, not independently normalized.
    """
    first = kraus_numerators(law, (0, 1), 3)
    branches = []
    identity = np.eye(8, dtype=object)
    phase = np.diag([(-1)**(x & 1) for x in range(8)]).astype(object)
    for label, left in first:
        current = label[1]
        seam = (0, 2) if current > 0 else (1, 2)
        control = phase if current == 0 else identity
        for following, right in kraus_numerators(law, seam, 3):
            branches.append(((label, following), right@control@left))
    return tuple(branches)


def history_unit(law, left, right):
    """Complete terminal operator for every two-current word on |left><right|."""
    if (type(left) is not int or type(right) is not int
            or not 0 <= left < 8 or not 0 <= right < 8):
        raise ValueError('three-port matrix unit required')
    output = defaultdict(F)
    for labels, matrix in controlled_history_numerators(law):
        word = tuple(label[1] for label in labels)
        for i, j in product(range(8), repeat=2):
            value = F(matrix[i, left]*matrix[j, right], 4)
            if value:
                output[word, i, j] += value
    return {key: value for key, value in output.items() if value}


def direct_history_unit(law, left, right):
    """Independent floor/ceiling instrument recursion; no dilation matrices.

    It keeps every current branch and applies the same feedback policy.
    Equality on all matrix units proves arbitrary-reference equality.
    """
    _domain(law, (0, 1), 3)
    if (type(left) is not int or type(right) is not int
            or not 0 <= left < 8 or not 0 <= right < 8):
        raise ValueError('three-port matrix unit required')

    def image(mask, seam, coin):
        a, b = seam
        x, y = (mask >> a) & 1, (mask >> b) & 1
        lo, hi = (x+y)//2, (x+y+1)//2
        out_a, out_b = (lo, hi) if coin == 0 else (hi, lo)
        out = mask & ~((1 << a) | (1 << b))
        return out | (out_a << a) | (out_b << b), x-out_a, x-y

    states = {((), left, right): F(1)}
    for step in range(2):
        following = defaultdict(F)
        for (word, a, b), weight in states.items():
            seam = (0, 1) if step == 0 else ((0, 2) if word[0] > 0 else (1, 2))
            for coin in (0, 1):
                first, second = image(a, seam, coin), image(b, seam, coin)
                label_index = 1 if law == 'C' else 2
                if first[label_index] != second[label_index]:
                    continue
                assert first[1] == second[1]
                i, j, current = first[0], second[0], first[1]
                phase = (-1)**((i & 1)+(j & 1)) if step == 0 and current == 0 else 1
                following[word+(current,), i, j] += weight*F(phase, 2)
        states = {key: value for key, value in following.items() if value}
    return states


def _digest(rows):
    return hashlib.sha256(canonical(rows).encode()).hexdigest()


def run():
    rows = []
    for law in ('C', 'E'):
        groups, w, d = completion_blocks(law)
        assert np.array_equal(w.T@w, 2*np.eye(4, dtype=object))
        assert np.array_equal(d.T, d) and np.array_equal(d@d, d)
        assert not np.any(d@w)
        histories = controlled_history_numerators(law)
        gram = sum((matrix.T@matrix for _, matrix in histories), np.zeros((8, 8), dtype=object))
        assert np.array_equal(gram, 4*np.eye(8, dtype=object))
        exact_history = []
        for left, right in product(range(8), repeat=2):
            actual = history_unit(law, left, right)
            assert actual == direct_history_unit(law, left, right)
            exact_history.append([left, right, [list(word)+[i, j, str(value)]
                for (word, i, j), value in sorted(actual.items())]])
        rows.append({
            'candidate': law,
            'local_source_dimension': 4,
            'hidden_kraus_labels': [list(label) for label, _ in groups],
            'completion_dimension': w.shape[0]+4,
            'exact_isometry_and_involutive_completion': True,
            'single_attempt_action': 'A_0=pi P_minus, exp(-i A_0)=U',
            'same_endpoint_actions': 'A_m=pi P_minus+2 pi m P_plus, m in Z',
            'two_repair_history_hidden_branches': len(histories),
            'three_port_input_matrix_units_checked': 64,
            'all_history_maps_equal': True,
            'arbitrary_reference_extension': True,
            'exact_history_operator_sha256': _digest(exact_history),
            'completion_integer_isometry_sha256': _digest(w.tolist()),
        })
    return {
        'schema': 'oph.native_action.history_action.v1',
        'native_carrier_embedding': carrier_embedding(),
        'native_instrument_history_actions': rows,
        'interpolation_nonuniqueness': {
            'comparison': 'A_0 and A_2 are positive and have identical native endpoints',
            'fractional_completion_parameter': '1/3',
            'input_sector_return_probability_A0': '3/4',
            'input_sector_return_probability_A2': '0',
            'readout_requires_access_to_auxiliary_completion_sector': True,
            'native_event_observations_distinguish_completions': False,
        },
        'scope': {
            'single_carrier_local_derivation': True,
            'federation_required': False,
            'additional_native_physics_law_asserted': False,
            'finite_history_environment_kept': True,
            'arbitrary_finite_adaptive_composition_follows_by_instrument_equality': True,
            'fresh_record_resources_are_not_free_physical_services': True,
            'unitary_completion_and_logarithm_selected_by_repair': False,
            'autonomous_source_hamiltonian_derived': False,
            'Lorentzian_continuum_action_or_lab_clock_established': False,
        },
        'source_bindings': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                            for path in SOURCE_PATHS}),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native unitary history action and endpoint ambiguity replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
