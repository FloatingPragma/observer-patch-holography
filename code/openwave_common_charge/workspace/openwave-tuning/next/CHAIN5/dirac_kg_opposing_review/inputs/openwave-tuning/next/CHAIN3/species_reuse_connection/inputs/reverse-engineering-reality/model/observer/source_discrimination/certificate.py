"""Exact encoded-channel checks and a conditional apparatus error contract.

The channel comparison independently derives floor/ceiling repair. It does
not call the prediction oracle or the native update. Exact permutation and
projector identities concern the mathematical source construction; they do
not certify floating compilation, numerical execution or physical services.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement

from . import circuit


APPARATUS_ALLOWANCE = F(1, 64)
ERROR_TERMS = (
    "pulse_operator", "load_channel", "export_channel", "preparation",
    "reset", "coin", "diagonal_randomizer", "record", "storage", "readout",
    "compilation", "numerical_execution",
)
RESOURCE_TERMS = ("native_pulses", "proper_code_load_events", "proper_code_export_events")


def _states():
    """Independent enumeration; no imported carrier or endpoint oracle."""
    states = []
    for first, second in combinations_with_replacement(range(12), 2):
        state = [0] * 12
        state[first] += 1
        state[second] += 1
        states.append(tuple(state))
    return tuple(sorted(states))


def _endpoint_image(state, candidate, coin):
    total = state[0] + state[1]
    low, high = total // 2, (total + 1) // 2
    endpoints = (low, high) if coin == 0 else (high, low)
    output = endpoints + state[2:]
    value = state[0] - endpoints[0] if candidate == "C" else state[0] - state[1]
    return (coin, value), output


def _permutation_image(candidate, coin, index):
    # E uses the same complete five-bit C permutation, tensored with
    # identity on its extra record bit, followed by one parity CNOT.
    local = index & 31
    for left, right in circuit.REPAIR_SWAPS["C", coin]:
        if local == left:
            local = right
        elif local == right:
            local = left
    index = (index & ~31) | local
    if candidate == "E" and (index >> coin) & 1:
        index ^= 32
    return index


def _matrix_unit_image(left, right):
    result = {}
    for coin in (0, 1):
        record_left, target_left = left[coin]
        record_right, target_right = right[coin]
        if record_left == record_right:
            result[record_left, target_left, target_right] = F(1, 2)
    return result


def channel_certificate(candidate):
    """Intertwining on every matrix unit, with all committed records intact.

    For the nine-qubit embedding J, the encoded recorded channel satisfies
    E_code(J rho J*) = (id_record tensor J) E_native(rho) (...)*.
    Matrix-unit equality implies the same identity with an arbitrary
    reference system. It is not just a comparison of basis populations.
    """
    if candidate not in ("C", "E"):
        raise ValueError("candidate C or E required")
    states = _states()
    codes = tuple(circuit.encode(state) for state in states)
    if len(states) != 78 or len(set(codes)) != 78:
        raise AssertionError("the full 78-state carrier must embed injectively")
    if any(circuit.decode(code) != state for code, state in zip(codes, states, strict=True)):
        raise AssertionError("encoding must preserve every native configuration")
    # Independently spelled out labels: E's raw code is j_code+4*parity,
    # and its contrast meaning depends on the already committed fair coin.
    decoders = {
        ("C", 0): {0: 0, 1: 1, 2: -1}, ("C", 1): {0: 0, 1: 1, 2: -1},
        ("E", 0): {0: 0, 1: 2, 2: -2, 4: -1, 5: 1},
        ("E", 1): {0: 0, 1: 2, 2: -2, 4: 1, 6: -1},
    }
    expected, encoded = [], []
    for state, code in zip(states, codes, strict=True):
        native_column, encoded_column = [], []
        for coin in (0, 1):
            record, output = _endpoint_image(state, candidate, coin)
            native_column.append((record, circuit.encode(output)))
            local = _permutation_image(candidate, coin, code % 8)
            value = decoders[candidate, coin].get(local // 8)
            if value is None:
                raise AssertionError("an admitted source column reached a reserved record")
            encoded_column.append(((coin, value), (code // 8) * 8 + local % 8))
        expected.append(tuple(native_column))
        encoded.append(tuple(encoded_column))
    if encoded != expected:
        raise AssertionError("the encoded dilation does not implement floor/ceiling repair")

    # Each joint (coin, record, output) column contains two amplitudes
    # 1/sqrt(2); distinct columns have disjoint support. This proves W*W=I.
    occupied = set()
    for column in encoded:
        support = set(column)
        if len(support) != 2 or occupied.intersection(support):
            raise AssertionError("joint record/source columns are not isometric")
        occupied.update(support)
    units = 0
    for left in range(len(states)):
        for right in range(len(states)):
            actual = _matrix_unit_image(encoded[left], encoded[right])
            target = _matrix_unit_image(expected[left], expected[right])
            if actual != target:
                raise AssertionError("recorded channels disagree on a matrix unit")
            units += 1
    return {
        "candidate": candidate, "native_dimension": len(states),
        "encoded_qubits": circuit.SOURCE_QUBITS, "matrix_units_checked": units,
        "kraus_entries": len(occupied),
        "record_count": len({record for column in encoded for record, _ in column}),
        "column_amplitude_squared": "1/2", "gram_identity_exact": True,
        "floor_ceiling_intertwining_exact": True,
        "arbitrary_reference_extension": True,
        "scope": "exact encoded target channel on the full native carrier, retaining every record; floating compilation and physical identification are separate",
    }


def parity_refinement_certificate():
    """E is exactly C followed by one retained seam-parity measurement.

    The required extra classical refinement has alphabet size two. This
    is a statement about the fixed recorded instruments, not minimal
    physical memory, native gate count or a rule selecting E in nature.
    """
    states = _states()
    current_columns, contrast_columns, refined_columns = [], [], []
    labels = {}
    for state in states:
        current_column, contrast_column, refined_column = [], [], []
        parity = (state[0] + state[1]) % 2
        for coin in (0, 1):
            current_record, output = _endpoint_image(state, "C", coin)
            contrast_record, direct_output = _endpoint_image(state, "E", coin)
            j = current_record[1]
            d = 2 * j + (2 * coin - 1) * parity
            if contrast_record != (coin, d) or output != direct_output:
                raise AssertionError("current and parity do not reconstruct the contrast")
            output_parity = (output[0] + output[1]) % 2
            c_output = _permutation_image("C", coin, circuit.encode(state) % 8)
            encoded_parity = (c_output >> coin) & 1
            if output_parity != parity or encoded_parity != parity:
                raise AssertionError("the one-bit output probe does not read seam parity")
            key = (coin, j, parity)
            labels[key] = (coin, d)
            inverse_parity = abs(d) % 2
            inverse_current = (d - (2 * coin - 1) * inverse_parity) // 2
            if (coin, inverse_current, inverse_parity) != key:
                raise AssertionError("contrast relabeling is not invertible")
            current_column.append((current_record, output))
            contrast_column.append((contrast_record, output))
            refined_column.append(((coin, d), output))
        current_columns.append(current_column)
        contrast_columns.append(contrast_column)
        refined_columns.append(refined_column)
    if len(labels) != 10 or len(set(labels.values())) != 10:
        raise AssertionError("all ten admitted contrast records must be retained")

    units = cross_shell_units = 0
    for left, x in enumerate(states):
        for right, z in enumerate(states):
            refined = _matrix_unit_image(refined_columns[left], refined_columns[right])
            contrast = _matrix_unit_image(contrast_columns[left], contrast_columns[right])
            if refined != contrast:
                raise AssertionError("parity refinement differs from E on a matrix unit")
            current = _matrix_unit_image(current_columns[left], current_columns[right])
            dephased_current = {
                key: value for key, value in current.items()
                if (sum(key[1][:2]) - sum(key[2][:2])) % 2 == 0
            }
            forgotten = {}
            for ((coin, d), a, b), value in contrast.items():
                parity = abs(d) % 2
                j = (d - (2 * coin - 1) * parity) // 2
                forgotten[(coin, j), a, b] = value
            if forgotten != dephased_current:
                raise AssertionError("forgetting parity must leave its source dephasing")
            units += 1
            cross_shell_units += sum(x[:2]) != sum(z[:2])

    fibers = Counter((coin, j) for coin, j, _ in labels)
    largest_fiber = max(fibers.values())
    # X and Z have the same C record (0,0), but E records contrast 0 and -1.
    x = tuple(int(port in (2, 3)) for port in range(12))
    z = tuple(int(port in (1, 2)) for port in range(12))
    if (_endpoint_image(x, "C", 0)[0] != _endpoint_image(z, "C", 0)[0]
            or _endpoint_image(x, "E", 0)[0] == _endpoint_image(z, "E", 0)[0]
            or largest_fiber != 2):
        raise AssertionError("the binary refinement must be necessary as well as sufficient")
    return {
        "native_dimension": 78, "matrix_units_checked": units,
        "cross_spectator_shell_matrix_units_checked": cross_shell_units,
        "record_relabeling_bijective": True, "admitted_refined_records": len(labels),
        "contrast_from_current_and_parity": "d=2*j+(2*c-1)*p",
        "parity_from_contrast": "p=abs(d)%2",
        "minimum_extra_classical_refinement_alphabet": largest_fiber,
        "extra_binary_record_sufficient": True, "nontrivial_refinement_necessary": True,
        "minimum_scope": "classical refinement of fixed C records; no claim of minimal hardware, environment, or gate count",
        "extra_native_cnot_words": 1,
        "E_is_recorded_parity_after_C": True,
        "forgetting_parity_leaves_dephased_C": True,
        "arbitrary_reference_extension": True,
        "scope": "identity of the declared finite C/E instruments, not an extra dynamical hypothesis or physical selection of parity recording",
    }


def _encoded_readout_column(address, outcome):
    """Exact U |pointer-sector><pointer-sector| U* column.

    U is the bit-zero Y preparation; two factors 1/sqrt(2) combine to 1/2.
    The pointer flips only at the ALL-zero nine-bit address, not at every
    even address. The inverse basis change is followed by its restoration.
    """
    base = address & ~1
    inverse_column = ((base, 1), (base + 1, 1 if address & 1 else -1))
    result = {}
    for middle, before in inverse_column:
        if (middle == 0) != (outcome == 1):
            continue
        forward_column = ((base, -1 if middle & 1 else 1), (base + 1, 1))
        for target, after in forward_column:
            result[target] = result.get(target, F()) + F(before * after, 2)
    return {target: value for target, value in result.items() if value}


def readout_certificate():
    """Exact full complement, including coherences among all minus states."""
    states = _states()
    codes = tuple(circuit.encode(state) for state in states)
    x = tuple(int(port in (2, 3)) for port in range(12))
    z = tuple(int(port in (1, 2)) for port in range(12))
    if (circuit.encode(x), circuit.encode(z)) != (0, 1):
        raise AssertionError("the fixed witness states must occupy addresses zero and one")
    columns = {}
    for address in range(1 << circuit.SOURCE_QUBITS):
        plus = {0: F(1, 2), 1: F(1, 2)} if address in (0, 1) else {}
        minus = {target: -value for target, value in plus.items()}
        minus[address] = minus.get(address, F()) + 1
        minus = {target: value for target, value in minus.items() if value}
        expected = {1: plus, -1: minus}
        for outcome in (1, -1):
            actual = _encoded_readout_column(address, outcome)
            if actual != expected[outcome]:
                raise AssertionError("the pointer does not implement the full Lueders projector")
            columns[address, outcome] = actual
    # Since each full-register Kraus column equals P+ or I-P+, the matrix
    # unit images are their exact outer products. Check their trace sum too.
    units = 0
    for left in codes:
        for right in codes:
            trace = sum((a * columns[right, outcome].get(target, F())
                         for outcome in (1, -1)
                         for target, a in columns[left, outcome].items()), F())
            if trace != int(left == right):
                raise AssertionError("the complete readout channel is not trace preserving")
            units += 1
    return {
        "encoded_basis_columns_checked": 1 << circuit.SOURCE_QUBITS,
        "native_matrix_units_checked": units, "plus_rank": 1,
        "native_minus_rank": 77, "encoded_minus_rank": (1 << circuit.SOURCE_QUBITS) - 1,
        "both_unnormalized_branches_retained": True,
        "minus_subspace_coherence_preserved": True,
        "old_record_blocks_unchanged": True, "arbitrary_reference_extension": True,
    }


def apparatus_error_contract(resources, error_bounds):
    """Sufficient conditional bound for the unchanged eta=1/64 protocol.

    Resources are maxima over COMPLETE attempted-trial branches, including
    the recorded fair coin and diagonal-control randomizer. Pulse errors
    are operator-norm bounds. Load/export errors are half-diamond bounds
    per invocation, on complete channels including failures. All other
    terms are separately supplied aggregate whole-trial trace/channel
    bounds, covering additional errors not already charged to these words.
    Every term must be explicitly supplied as an exact nonnegative rational.
    This arithmetic never establishes any of those physical premises.
    """
    if not isinstance(resources, dict) or any(name not in resources for name in RESOURCE_TERMS):
        raise ValueError("complete pulse, load and export resource counts required")
    if any(type(resources[name]) is not int or resources[name] < 0 for name in RESOURCE_TERMS):
        raise ValueError("nonnegative integer resource counts required")
    if not isinstance(error_bounds, dict) or set(error_bounds) != set(ERROR_TERMS):
        raise ValueError("every declared error term must be explicitly supplied")
    if any(not isinstance(value, F) or value < 0 for value in error_bounds.values()):
        raise ValueError("nonnegative exact rational error bounds required")
    contributions = {
        "pulse_operator": resources["native_pulses"] * error_bounds["pulse_operator"],
        "load_channel": resources["proper_code_load_events"] * error_bounds["load_channel"],
        "export_channel": resources["proper_code_export_events"] * error_bounds["export_channel"],
    }
    contributions.update({name: error_bounds[name] for name in ERROR_TERMS[3:]})
    additive = sum(contributions.values(), F())
    bound = min(F(1), additive)
    return {
        "contributions": {name: str(value) for name, value in contributions.items()},
        "additive_error_upper": str(additive),
        "complete_instrument_trace_distance_upper_if_premises_hold": str(bound),
        "every_joint_plus_and_abort_probability_error_upper_if_premises_hold": str(bound),
        "frozen_protocol_apparatus_allowance": str(APPARATUS_ALLOWANCE),
        "sufficient_for_frozen_protocol_if_premises_hold": additive <= APPARATUS_ALLOWANCE,
        "physical_error_certificate": False, "certified_numerical_execution": False,
        "scope": "conditional telescoping implication only; all supplied upper bounds require independent justification; all failed attempts remain abort outcomes",
    }
