"""Independent endpoint, tensor-trace and hostile tests for localized charge."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import core


def independent_local(mask, law, coin):
    """Enumerate the occupied/vacant endpoint rule without calling its source."""
    left, right = mask & 3, mask >> 2
    if left == 3 or right == 3:
        return 8, mask
    na, nb = int(bool(left)), int(bool(right))
    if na != nb and na != coin:
        left, right = right, left
    current = na - int(bool(left))
    raw = {0: 0, 1: 1, -1: 2}[current]
    if law == "E":
        raw += 4 * (na != nb)
    return raw, left + 4 * right


def independent_configuration(index):
    orientation, pair = divmod(index, 66)
    a, b = core.PAIRS[pair]
    result = [0] * 12
    result[a], result[b] = ((1, 2) if not orientation else (2, 1))
    return tuple(result)


def independent_branch(law, seam, coin, raw):
    configurations = tuple(independent_configuration(i) for i in range(132))
    indices = {state: i for i, state in enumerate(configurations)}
    result = []
    a, b = seam
    for i, source in enumerate(configurations):
        flag, output = independent_local(source[a] + 4 * source[b], law, coin)
        if flag == raw:
            target = list(source)
            target[a], target[b] = output & 3, output >> 2
            result.append((i, indices[tuple(target)]))
    return result


def matrix_branch(matrix, law, seam, coin, raw):
    pairs = independent_branch(law, seam, coin, raw)
    result = np.zeros((132, 132), complex)
    if pairs:
        inputs, outputs = zip(*pairs)
        result[np.ix_(outputs, outputs)] = matrix[np.ix_(inputs, inputs)] / 2
    return result


def trace_distance(matrix):
    return np.linalg.norm(matrix, ord="nuc") / 2


def test_whole_local_completion_is_the_unchanged_colorblind_endpoint_law():
    for law in ("C", "E"):
        rows = core.local_completion(law)
        assert set(rows) == set(product(range(2), range(16)))
        assert sum(len(row) for row in rows.values()) == 32
        assert sum(len(rows[coin, 8]) for coin in (0, 1)) == 14
        for coin in (0, 1):
            seen = []
            for raw in range(16):
                for origin, output, sign in rows[coin, raw]:
                    assert (raw, output) == independent_local(origin, law, coin)
                    assert sign == 1
                    seen.append(origin)
            assert sorted(seen) == list(range(16))


def test_complete_local_tensor_intertwiner_and_chart_coordinates():
    result = core.intertwiner_certificate()
    assert result["native_raw_slots"] == 720
    assert result["supported_cover_columns"] == 15840
    assert result["input_columns_checked_per_slot"] == 132
    assert len(result["local_maps_sha256"]) == 64
    for index in range(132):
        configuration = core.cover_configuration(index)
        assert configuration == independent_configuration(index)
        assert core.cover_index(configuration) == index
    for law in ("C", "E"):
        for (seam, coin, raw), row in core.ordered_cover_branches(law).items():
            assert row == tuple(independent_branch(law, seam, coin, raw))


def test_pair_kernel_has_independent_exclusion_walk_entries():
    kernel = core.pair_kernel()
    for i, pair in enumerate(core.PAIRS):
        expected = [F(0)] * 66
        for seam in core.SEAMS:
            occupied = set(pair).intersection(seam)
            if len(occupied) != 1:
                expected[i] += F(1, 30)
            else:
                old = next(iter(occupied))
                new = next(port for port in seam if port != old)
                moved = tuple(sorted((set(pair) - {old}) | {new}))
                expected[i] += F(1, 60)
                expected[core.PAIRS.index(moved)] += F(1, 60)
        assert tuple(expected) == kernel[i]
    mixing = core.mixing_certificate()
    assert mixing["states"] == 66
    assert len(mixing["spanning_tree"]) == 65
    assert mixing["minimum_holding_probability"] == "5/6"
    assert mixing["maximum_holding_probability"] == "13/15"
    constant = (F(1, 66),) * 66
    assert core.advance(constant, 7) == constant


def test_exact_localization_transport_and_asymptotic_charge_accounting():
    initial = core.initial_coefficients()
    expected = ((F(1), F(0), F(0)),
                (F(13, 15), F(0), F(2, 15)),
                (F(113, 150), F(7, 900), F(43, 180)),
                (F(5911, 9000), F(8, 375), F(2897, 9000)))
    for step, weights in enumerate(expected):
        coefficients = core.advance(initial, step)
        report = core.region_report(coefficients, core.CORE)
        assert tuple(report[key] for key in ("inside", "outside", "split")) == weights
        assert report["charge_expectation"] == F(1, 66)
        assert report["global_trace_distance"] == F(1, 132)
        assert sum(coefficients) == 1
    for size in range(13):
        values = core.initial_coefficients(tuple(range(size)))
        assert sum(values) == size * (size - 1) // 2
        for epsilon in (F(-1), F(-1, 3), F(0), F(1, 2), F(1)):
            rho = core.density_matrix(values, epsilon)
            assert abs(np.trace(rho) - 1) < 1e-14
            # Each orientation block has exact eigenvalues (1 +/- eps f)/132.
            assert min(F(1) - abs(epsilon) * value for value in values) >= 0
            assert np.linalg.eigvalsh(rho).min() >= -1e-14
    limit = core.region_report((F(1, 66),) * 66, core.CORE)
    assert limit["inside"] == F(1, 66)
    assert limit["outside"] == F(45, 66)
    assert limit["split"] == F(20, 66)


def test_exact_pair_current_divergence_and_regional_return_flux():
    initial = core.initial_coefficients()
    first = core.current_certificate(initial, core.CORE)
    assert first["configuration_vertices"] == 66
    assert first["movement_edges"] == 300
    assert len(first["nonzero_currents"]) == 8
    assert first["regional_balance"] == {
        "inside_before": "1", "inside_after": "13/15", "gross_outflow": "2/15",
        "gross_inflow": "0", "net_outflow": "2/15", "inside_change": "-2/15"}
    second = core.current_certificate(core.advance(initial, 1), core.CORE)
    assert second["regional_balance"]["gross_inflow"] == "1/450"
    assert second["regional_balance"]["gross_outflow"] == "26/225"
    assert second["regional_balance"]["inside_change"] == "-17/150"
    # Independent nonuniform coefficients exercise currents in both directions.
    coefficients = tuple(F((7 * i + 5) % 17, 17) for i in range(66))
    changes = [F(0)] * 66
    kernel = core.pair_kernel()
    for i in range(66):
        for j in range(66):
            changes[i] -= coefficients[i] * kernel[i][j]
            changes[j] += coefficients[i] * kernel[i][j]
    for region in ((), (0,), (0, 1), (0, 5, 11), tuple(range(12))):
        report = core.current_certificate(coefficients, region)
        assert tuple(map(F, report["one_attempt_change"])) == tuple(changes)
        assert tuple(map(F, report["divergence"])) == tuple(-x for x in changes)
        selected = [i for i, pair in enumerate(core.PAIRS) if set(pair) <= set(region)]
        assert F(report["regional_balance"]["inside_change"]) == sum(changes[i] for i in selected)
    assert all(current == 0 for _, _, current in core.pair_currents((F(1, 66),) * 66))


def test_partial_trace_against_independent_ambient_qutrit_tensor():
    # Embed an arbitrary entangled vector into the literal 3^12 tensor space.
    rng = np.random.default_rng(497)
    vector = rng.normal(size=132) + 1j * rng.normal(size=132)
    vector /= np.linalg.norm(vector)
    ambient = np.zeros((3,) * 12, complex)
    for index, amplitude in enumerate(vector):
        ambient[independent_configuration(index)] = amplitude
    density = np.outer(vector, vector.conj())
    for region in ((), (1,), (0, 1), (0, 5, 11)):
        outside = tuple(port for port in range(12) if port not in region)
        rows = ambient.transpose(region + outside).reshape(3 ** len(region), -1)
        expected = rows @ rows.conj().T
        basis, actual = core.region_density(density, region)
        padded = np.zeros_like(expected)
        addresses = [sum(value * 3 ** (len(region) - position - 1)
                         for position, value in enumerate(state)) for state in basis]
        padded[np.ix_(addresses, addresses)] = actual
        np.testing.assert_allclose(padded, expected, atol=2e-15, rtol=0)


def test_regional_distance_is_not_a_claim_that_split_charge_disappeared():
    zero = core.density_matrix((F(0),) * 66)
    regions = (core.CORE, tuple(range(2, 12)), (0,), (), tuple(range(12)))
    for step in (0, 1, 2):
        coefficients = core.advance(core.initial_coefficients(), step)
        delta = core.density_matrix(coefficients) - zero
        assert abs(trace_distance(delta) - 1 / 132) < 1e-14
        for region in regions:
            _, actual = core.region_density(delta, region)
            expected = core.region_report(coefficients, region)["regional_trace_distance"]
            assert abs(trace_distance(actual) - float(expected)) < 1e-14
    # An individual split pair has no contrast on either separate side.
    split = tuple(F(pair == (0, 2)) for pair in core.PAIRS)
    delta = core.density_matrix(split) - zero
    for region in (core.CORE, tuple(range(2, 12))):
        assert np.count_nonzero(core.region_density(delta, region)[1]) == 0
    assert trace_distance(delta) > 0


def test_all_complete_two_attempt_histories_preserve_tape_probability_and_contrast():
    initial = core.initial_coefficients()
    reference = core.density_matrix((F(0),) * 66)
    defect = core.density_matrix(initial)
    summed = [F(0)] * 66
    total_mass, total_distance, exterior_distance = 0.0, 0.0, 0.0
    raw_slots = 0
    for coin1, raw1, coin2, raw2 in product(range(2), range(4), range(2), range(8)):
        history = (("C", (0, 5), coin1, raw1), ("E", (1, 7), coin2, raw2))
        sigma, rho = reference, defect
        for event in history:
            sigma, rho = matrix_branch(sigma, *event), matrix_branch(rho, *event)
        assert abs(np.trace(rho) - np.trace(sigma)) < 1e-15
        coefficients = core.history_transport(initial, history)
        expected_delta = core.density_matrix(coefficients) - reference
        np.testing.assert_allclose(rho - sigma, expected_delta, atol=2e-15, rtol=0)
        total_mass += np.trace(rho).real
        total_distance += trace_distance(rho - sigma)
        exterior_distance += trace_distance(core.region_density(rho - sigma, tuple(range(2, 12)))[1])
        summed = [a + b for a, b in zip(summed, coefficients, strict=True)]
        raw_slots += 1
    assert raw_slots == 128
    assert abs(total_mass - 1) < 2e-14
    assert abs(total_distance - 1 / 132) < 2e-14
    report = core.region_report(tuple(summed), core.CORE)
    assert abs(exterior_distance - float(report["exterior_trace_distance"])) < 2e-14
    assert core.history_transport(initial, (("C", (0, 1), 0, 3),)) == (F(0),) * 66
    assert core.history_transport(initial, ()) == initial


def test_actual_odd_exchange_acts_on_two_labels_and_not_an_added_chosen_phase():
    for law in ("C", "E"):
        loop = core.protection.exchange_loop(law)
        history = tuple((event[0], tuple(event[1]), event[2], event[3]) for event in loop["history"])
        pair = core.PAIRS.index(tuple(loop["initial_pair"]))
        for initial in (pair, pair + 66):
            current = initial
            for event in history:
                flag, current, sign, _ = core.cover_action(current, *event[:3])
                assert flag == event[3] and sign == 1
            assert current == (initial + 66) % 132
        coefficients = tuple(F(i == pair) for i in range(66))
        final = core.history_transport(coefficients, history)
        assert final[pair] == F(1, 8) and sum(final) == F(1, 8)


def test_local_label_phase_and_dephasing_break_the_charge_while_tape_stays_blind():
    rho = core.density_matrix(core.initial_coefficients())
    q = np.zeros((132, 132), complex)
    for i in range(66):
        q[i, i + 66] = q[i + 66, i] = 1
    assert abs(np.trace(rho @ q) - 1 / 66) < 1e-15
    for port in core.CORE:
        phased = core.label_phase(rho, port)
        dephased = core.label_phase(rho, port, dephase=True)
        assert abs(np.trace(phased @ q) + 1 / 66) < 1e-15
        assert np.trace(dephased @ q) == 0
        np.testing.assert_array_equal(np.diag(phased), np.diag(rho))
        np.testing.assert_array_equal(np.diag(dephased), np.diag(rho))
    np.testing.assert_array_equal(core.label_phase(rho, 4), rho)


def test_sharp_sector_exterior_obstruction_uses_the_actual_local_charge():
    rho = core.density_matrix((F(1),) * 66)
    reference = core.density_matrix((F(0),) * 66)
    for size in (0, 1, 2, 5, 10, 11, 12):
        region = tuple(range(size, 12))
        _, delta = core.region_density(rho - reference, region)
        bound = F((12 - size) * (11 - size), 2 * 132)
        assert abs(trace_distance(delta) - float(bound)) < 2e-14


def test_missing_predecessor_raw_branch_is_not_silently_accepted(monkeypatch):
    original = core.protection.branches
    monkeypatch.setattr(core.protection, "branches", lambda law: original(law)[:-1])
    with pytest.raises(ArithmeticError, match="missing or duplicated"):
        core.intertwiner_certificate()


def test_changed_local_collision_completion_is_rejected(monkeypatch):
    original = core.full_domain.branch_action
    def wrong(mask, law, coin, **kwargs):
        raw, output, sign, current = original(mask, law, coin, **kwargs)
        return (0 if raw == 8 else raw), output, sign, current
    monkeypatch.setattr(core.full_domain, "branch_action", wrong)
    with pytest.raises(ArithmeticError, match="collision complement"):
        core.intertwiner_certificate()


def test_wrong_global_car_extension_is_detected(monkeypatch):
    original = core.cover_action
    def wrong(index, law, seam, coin):
        raw, target, _, current = original(index, law, seam, coin)
        mask = sum(value << (2 * port) for port, value in enumerate(core.cover_configuration(index)))
        _, _, sign, _ = core.full_domain.branch_action(mask, law, coin, sites=12, flavors=2, seam=seam)
        return raw, target, sign, current
    monkeypatch.setattr(core, "cover_action", wrong)
    with pytest.raises(ArithmeticError, match="ordered cover"):
        core.intertwiner_certificate()


@pytest.mark.parametrize("bad", [True, 0.0, float("nan"), -1, 132, "0", None])
def test_bad_cover_index_rejected(bad):
    with pytest.raises(ValueError):
        core.cover_configuration(bad)


@pytest.mark.parametrize("bad", [(0, True), (0, 1.0), (1, 0), (0, 0), (0, 12), [0, 1], None])
def test_bad_region_rejected(bad):
    with pytest.raises(ValueError):
        core.initial_coefficients(bad)


@pytest.mark.parametrize("bad", [True, 1.0, float("nan"), F(-1), F(3, 2), "1", None])
def test_bad_coefficient_rejected(bad):
    values = (bad,) + (F(0),) * 65
    for operation in (core.density_matrix, lambda x: core.advance(x, 0),
                      lambda x: core.history_transport(x, ()),
                      lambda x: core.region_report(x, core.CORE)):
        with pytest.raises(ValueError):
            operation(values)


@pytest.mark.parametrize("history", [[], (("C", (0, 1), False, 0),),
                                    (("C", (0, 1), 0, 0.0),),
                                    (("C", (0, 1), 0, 4),),
                                    (("E", (0, 2), 0, 0),),
                                    (("X", (0, 1), 0, 0),), (("C",),)])
def test_bad_history_rejected(history):
    with pytest.raises(ValueError):
        core.history_transport(core.initial_coefficients(), history)


def test_bad_matrix_epsilon_and_auxiliary_controls_rejected():
    coefficients = core.initial_coefficients()
    for epsilon in (True, 1.0, float("nan"), F(2), "1"):
        with pytest.raises(ValueError):
            core.density_matrix(coefficients, epsilon)
    for steps in (True, 1.0, -1, float("nan")):
        with pytest.raises(ValueError):
            core.advance(coefficients, steps)
    for matrix in (np.eye(131), np.ones((132, 132), bool),
                   np.full((132, 132), np.nan), np.full((132, 132), np.inf)):
        with pytest.raises(ValueError):
            core.region_density(matrix, core.CORE)
    with pytest.raises(ValueError):
        core.label_phase(core.density_matrix(coefficients), True)
    with pytest.raises(ValueError):
        core.label_phase(core.density_matrix(coefficients), 0, dephase=1)
    with pytest.raises(ValueError):
        core.cover_index((0,) * 12)
    with pytest.raises(ValueError):
        core.local_branch_action(True, "C", 0)
    with pytest.raises(ValueError):
        core.local_branch_action(0, "C", 0.0)


@pytest.mark.parametrize("numerator,denominator", [(1, 0), (-1, -2), (2, 4),
                                                   (0, 2), (True, 2), (1, True),
                                                   (1.0, 2), (1, 2.0)])
def test_malformed_fraction_storage_is_rejected_before_arithmetic(numerator, denominator):
    bad = F(0)
    object.__setattr__(bad, "_numerator", numerator)
    object.__setattr__(bad, "_denominator", denominator)
    coefficients = (bad,) + (F(0),) * 65
    for operation in (core.density_matrix, core.pair_currents,
                      lambda x: core.advance(x, 0),
                      lambda x: core.history_transport(x, ()),
                      lambda x: core.region_flux(x, core.CORE)):
        with pytest.raises(ValueError, match="Fraction storage"):
            operation(coefficients)
    with pytest.raises(ValueError, match="Fraction storage"):
        core.density_matrix(core.initial_coefficients(), bad)


def test_valid_integer_and_canonical_fraction_inputs_agree():
    integers = tuple(int(pair == core.CORE) for pair in core.PAIRS)
    assert core.advance(integers, 1) == core.advance(core.initial_coefficients(), 1)
    np.testing.assert_array_equal(core.density_matrix(integers, 1),
                                  core.density_matrix(core.initial_coefficients(), F(1)))


def test_matrix_helper_refuses_masking_subclasses_loss_and_extended_precision():
    class ArraySubclass(np.ndarray):
        pass
    candidates = [np.ma.array(np.eye(132), mask=False),
                  np.eye(132).view(ArraySubclass),
                  np.full((132, 132), "0"), np.full((132, 132), object()),
                  np.full((132, 132), 2**53 + 1, dtype=np.int64),
                  np.full((132, 132), 2**64 - 1, dtype=np.uint64)]
    if np.dtype(np.longdouble).itemsize > 8:
        candidates.extend([np.eye(132, dtype=np.longdouble),
                           np.full((132, 132), np.finfo(np.longdouble).max, dtype=np.longdouble)])
    if np.dtype(np.clongdouble).itemsize > 16:
        candidates.append(np.eye(132, dtype=np.clongdouble))
    for value in candidates:
        with pytest.raises(ValueError):
            core.region_density(value, core.CORE)
        with pytest.raises(ValueError):
            core.label_phase(value, 0)
    for dtype in (np.int32, np.uint64, np.float32, np.float64, np.complex64, np.complex128):
        basis, result = core.region_density(np.eye(132, dtype=dtype), ())
        assert basis == ((),) and result[0, 0] == 132
