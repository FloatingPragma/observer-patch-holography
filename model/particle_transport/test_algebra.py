"""Independent CAR, exterior-algebra and complete-spectrum controls."""
import itertools
import math
from dataclasses import replace

import numpy as np
import pytest

from . import algebra as a


def jw(modes):
    """Tensor Jordan--Wigner oracle, independent of FullFock address loops."""
    one, z = np.eye(2), np.diag([1., -1.])
    lower = np.array([[0., 1.], [0., 0.]])
    operators = []
    for mode in range(modes):
        result = np.ones((1, 1))
        for bit in reversed(range(modes)):
            result = np.kron(result, lower if bit == mode else z if bit < mode else one)
        operators.append(result)
    return operators


def independent_dimer(axis):
    c = jw(4)
    sigma = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.diag([1, -1]))[axis]
    h = np.zeros((16, 16), complex)
    for left, right in itertools.product(range(2), repeat=2):
        term = 1j*sigma[left, right]/16*(c[left].conj().T@c[2+right])
        h += term+term.conj().T
    numbers = [x.conj().T@x for x in c]
    return h+numbers[0]@numbers[1]+numbers[2]@numbers[3]-sum(numbers)/2


def wedge_address(address, permutation, width):
    """Permute each occupied orbital and sort its wedge product explicitly."""
    images = [width*permutation[j//width]+j % width
              for j in range(width*len(permutation)) if address >> j & 1]
    inversions = sum(x > y for j, x in enumerate(images) for y in images[j+1:])
    return sum(1 << j for j in images), (-1)**inversions


@pytest.mark.parametrize("axis", range(3))
def test_analytic_columns_are_the_actual_full_Fock_candidates(axis):
    k = independent_dimer(axis)
    np.testing.assert_allclose(k, a.setup(axis).generator, atol=0, rtol=0)
    g = a.ground_column(axis)
    np.testing.assert_allclose(k@g, a.GROUND_ENERGY*g, atol=2e-15)
    assert np.vdot(g, g).real == pytest.approx(1, abs=1e-15)
    numbers = np.array([n.bit_count() for n in range(16)])
    for kind, number in zip(a.KINDS, (1, 3, 0, 4), strict=True):
        for band in (("lower", "upper") if number in (1, 3) else ("lower",)):
            v = a.candidate_columns(kind, axis, band)
            np.testing.assert_allclose(v.conj().T@v, np.eye(v.shape[1]), atol=5e-16)
            np.testing.assert_array_equal(v[numbers != number], 0)
            np.testing.assert_allclose(k@v, (a.GROUND_ENERGY+a.candidate_gap(kind, band))*v, atol=2e-15)
    # Complete odd-number sector, including the upper pole, remains present.
    odd = sum(a.candidate_projector(kind, axis, band)
              for kind in ("hole", "particle") for band in ("lower", "upper"))
    np.testing.assert_allclose(odd, np.diag(numbers % 2), atol=5e-16)


def test_full_graded_block_maps_and_nontrivial_three_block_composition():
    permutations = tuple(itertools.permutations(range(3)))
    # Every basis state of three complete dimers; no parity-sector sampling.
    for address in range(4096):
        for permutation in permutations:
            assert a.block_permutation_address(address, permutation) == wedge_address(address, permutation, 4)
        for first, second in (((1, 0, 2), (0, 2, 1)), ((2, 0, 1), (1, 2, 0)),
                              ((2, 1, 0), (2, 1, 0))):
            middle, phase1 = a.block_permutation_address(address, first)
            target, phase2 = a.block_permutation_address(middle, second)
            composed = tuple(second[first[j]] for j in range(3))
            assert (target, phase1*phase2) == a.block_permutation_address(address, composed)
    # The two braid words agree on the entire three-block space.
    for address in range(4096):
        outputs = []
        for word in (((1, 0, 2), (0, 2, 1), (1, 0, 2)),
                     ((0, 2, 1), (1, 0, 2), (0, 2, 1))):
            target, phase = address, 1
            for permutation in word:
                target, factor = a.block_permutation_address(target, permutation)
                phase *= factor
            outputs.append((target, phase))
        assert outputs[0] == outputs[1]


def test_neighbor_FSWAP_route_retains_double_occupancy_sign_on_every_column():
    route = a.adjacent_swap_route((1, 0))
    assert len(route) == 16
    fswap = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, -1]])
    result = np.eye(256, dtype=int)
    for left in route:
        result = np.kron(np.eye(1 << (6-left), dtype=int),
                         np.kron(fswap, np.eye(1 << left, dtype=int)))@result
    np.testing.assert_array_equal(result, a.block_swap())
    assert result[17, 17] == -1  # One occupied mode in each exchanged block.


@pytest.mark.parametrize("axis", range(3))
@pytest.mark.parametrize("kind", a.KINDS)
def test_localized_candidate_moves_with_unknown_spin_and_entangled_reference(kind, axis):
    source = a.localized_columns(kind, 0, axis)
    destination = a.localized_columns(kind, 1, axis)
    swap = a.block_swap()
    np.testing.assert_allclose(swap@source, destination, atol=1e-15)
    # Deliberately complex columns, with a nontrivial untouched reference.
    rng = np.random.default_rng(1907)
    state = rng.normal(size=(source.shape[1], 3))+1j*rng.normal(size=(source.shape[1], 3))
    state /= np.linalg.norm(state)
    actual, target = swap@(source@state), destination@state
    np.testing.assert_allclose(actual, target, atol=1e-15)
    # The whole exterior ground density is unchanged, not only its number.
    before = (source@state).reshape(16, 16, 3)
    exterior = np.einsum("air,bir->ab", before, before.conj())
    ground = a.ground_column(axis)
    np.testing.assert_allclose(exterior, np.outer(ground, ground.conj()), atol=1e-15)


def test_two_odd_candidates_exchange_internal_states_with_minus_sign():
    swap = a.block_swap()
    for left_kind, right_kind in itertools.product(a.KINDS, repeat=2):
        left, right = a.candidate_columns(left_kind), a.candidate_columns(right_kind)
        sign = -1 if left.shape[1] == right.shape[1] == 2 else 1
        for x, y in itertools.product(left.T, right.T):
            np.testing.assert_allclose(swap@np.kron(y, x), sign*np.kron(x, y), atol=1e-15)
    # Mixed spins do not have a universal -1 scalar Ramsey expectation.
    hole = a.candidate_columns("hole")
    opposite = np.kron(hole[:, 1], hole[:, 0])
    assert abs(np.vdot(opposite, swap@opposite)) < 1e-15
    mixed = a.candidate_projector("hole")/2
    assert np.trace(swap@np.kron(mixed, mixed)).real == pytest.approx(-.5, abs=1e-15)


@pytest.mark.parametrize("axis", range(3))
def test_neutral_current_complete_measure_matches_independent_CAR_and_all_three_bands(axis):
    c = jw(8)
    generator = np.kron(np.eye(16), independent_dimer(axis))+np.kron(independent_dimer(axis), np.eye(16))
    energy, vectors = np.linalg.eigh(generator)
    g = a.ground_column(axis)
    reference = np.kron(g, g)
    number = np.diag([n.bit_count() for n in range(256)])
    for mode in range(4):
        transfer = c[mode].conj().T@c[4+mode]
        for name, operator in (("directed", transfer), ("real", transfer+transfer.conj().T),
                               ("current", -1j*(transfer-transfer.conj().T))):
            np.testing.assert_allclose(a.transfer_operator(axis, mode, name), operator, atol=0, rtol=0)
            np.testing.assert_array_equal(number@operator-operator@number, 0)
            weights = abs(vectors.conj().T@operator@reference)**2
            frequencies = energy-2*a.GROUND_ENERGY
            data = a.neutral_pair_measure(axis, mode, name)
            assert len(data["weights"]) == 256
            # Moments compare complete arrays; no visibility threshold removes dark states.
            for order in range(5):
                direct = np.sum(weights*frequencies**order)
                actual = np.sum(data["weights"]*data["frequencies"]**order)
                expected = sum(w*f**order for f, w in a.pair_line_data(name))
                assert actual == pytest.approx(direct, abs=2e-13)
                assert actual == pytest.approx(expected, abs=2e-13)
            assert np.sum(data["weights"][data["numbers"] != 4]) < 1e-25


@pytest.mark.parametrize("quadrature,parity", (("current", -1), ("real", 1)))
def test_pair_exchange_sector_and_same_frequency_are_distinct_observables(quadrature, parity):
    swap = a.block_swap()
    generator = a.pair_generator()
    for band, (frequency, weight) in zip(("lower", "mixed", "upper"), a.pair_line_data(quadrature), strict=True):
        phi = a.exchange_pair_column(band=band, quadrature=quadrature)
        assert np.vdot(phi, phi).real == pytest.approx(1, abs=3e-15)
        np.testing.assert_allclose(swap@phi, parity*phi, atol=2e-15)
        np.testing.assert_allclose(generator@phi, (2*a.GROUND_ENERGY+frequency)*phi, atol=3e-15)
        assert weight > 0


def test_complete_number_read_disrupts_lower_candidate_despite_charge_conservation():
    n = np.diag([int(address & 1 != 0) for address in range(16)])
    for kind in ("hole", "particle"):
        p = a.candidate_projector(kind)
        rho = p/2
        changed = n@rho@n+(np.eye(16)-n)@rho@(np.eye(16)-n)
        assert np.trace(changed).real == pytest.approx(1, abs=1e-15)
        assert np.trace(p@changed).real == pytest.approx(.75, abs=1e-15)


@pytest.mark.parametrize("bad", (True, np.int64(0), -1, 3, 0., "0"))
def test_axis_guards_apply_to_all_candidate_constructors(bad):
    for function in (a.setup, a.ground_column, lambda x: a.candidate_columns(axis=x),
                     a.pair_generator, lambda x: a.transfer_operator(axis=x)):
        with pytest.raises(ValueError):
            function(bad)


@pytest.mark.parametrize("arguments", ((True, (1, 0), 4), (0, (True, 0), 4),
    (256, (1, 0), 4), (-1, (1, 0), 4), (0, (0, 0), 4), (0, [1, 0], 4), (0, (1, 0), True)))
def test_integer_permutation_contract_rejects_ambiguous_inputs(arguments):
    with pytest.raises(ValueError):
        a.block_permutation_address(*arguments)


def test_certificate_rejects_an_ungraded_exchange(monkeypatch):
    ordinary = np.zeros((256, 256), int)
    for address in range(256):
        ordinary[((address & 15) << 4)+(address >> 4), address] = 1
    monkeypatch.setattr(a, "block_swap", lambda: ordinary)
    with pytest.raises(ArithmeticError, match="routing"):
        a.certificate()


def test_full_native_negative_lift_dictionary_binds_every_raw_flag_and_rotation():
    result = a.native_intertwiner_certificate()
    assert result["exact_directed_movement_columns"] == 600
    assert result["complete_C_and_E_raw_slots"] == 720
    assert result["complete_C_and_E_supported_columns"] == 7920
    assert result["exact_rotation_columns"] == 3960
    assert result["relative_CAR_number_per_occupied_block"] == {"hole": -1, "particle": 1}
    assert result["48_mode_source_execution_claimed"] is False


def test_actual_three_dimer_isometry_and_local_charge_intertwine_all_permutations():
    # Materialize complete12-mode vectors of the actual interacting ground
    # and a complex internal state. No block-parity shortcut builds W or its
    # transported vectors; the independent occupied-orbital oracle acts on
    # each of4096 ambient Fock addresses.
    pairs = ((0, 1), (0, 2), (1, 2))
    permutations = tuple(itertools.permutations(range(3)))
    actions = {}
    for permutation in permutations:
        images = [wedge_address(address, permutation, 4) for address in range(4096)]
        actions[permutation] = (np.array([row[0] for row in images]),
                                np.array([row[1] for row in images]))
    local_numbers = [np.array([((address >> (4*cell)) & 15).bit_count()
                              for address in range(4096)]) for cell in range(3)]
    for axis in range(3):
        ground = a.ground_column(axis)
        for kind, charge in (("hole", -1), ("particle", 1)):
            candidate = a.candidate_columns(kind, axis)@np.array([1, 2j])/math.sqrt(5)
            columns = []
            for pair in pairs:
                column = np.ones(1, complex)
                for cell in reversed(range(3)):
                    column = np.kron(column, candidate if cell in pair else ground)
                columns.append(column)
            w = np.column_stack(columns)
            np.testing.assert_allclose(w.conj().T@w, np.eye(3), atol=2e-15)
            for cell, number in enumerate(local_numbers):
                actual = w.conj().T@((number-2)[:, None]*w)
                expected = np.diag([charge*int(cell in pair) for pair in pairs])
                np.testing.assert_allclose(actual, expected, atol=2e-15)
            for permutation, (targets, signs) in actions.items():
                actual = np.empty_like(w)
                actual[targets] = signs[:, None]*w
                wanted = []
                for pair in pairs:
                    moved = tuple(permutation[cell] for cell in pair)
                    index = pairs.index(tuple(sorted(moved)))
                    wanted.append((-1 if moved[0] > moved[1] else 1)*w[:, index])
                np.testing.assert_allclose(actual, np.column_stack(wanted), atol=1e-15)


def test_dictionary_rejects_positive_lift_in_place_of_fermionic_branches(monkeypatch):
    from ..exchange_sector import statistics
    original = statistics.branch

    def wrong(*args):
        branch = original(*args)
        return replace(branch, signs=tuple(1 for _ in branch.signs))

    monkeypatch.setattr(statistics, "branch", wrong)
    with pytest.raises(ArithmeticError, match="signed branch"):
        a.native_intertwiner_certificate()


def test_dictionary_rejects_dropped_supported_column_and_corrupted_empty_flag(monkeypatch):
    from ..exchange_sector import statistics
    original = statistics.branch

    def missing(*args):
        branch = original(*args)
        if branch.inputs:
            return replace(branch, inputs=branch.inputs[1:], outputs=branch.outputs[1:], signs=branch.signs[1:])
        return branch

    monkeypatch.setattr(statistics, "branch", missing)
    with pytest.raises(ArithmeticError, match="signed branch"):
        a.native_intertwiner_certificate()

    def false_empty(*args):
        branch = original(*args)
        return replace(branch, abort=not branch.abort) if not branch.inputs else branch

    monkeypatch.setattr(statistics, "branch", false_empty)
    with pytest.raises(ArithmeticError, match="signed branch"):
        a.native_intertwiner_certificate()


def test_certificate_rejects_a_nonfinite_later_projector_residual(monkeypatch):
    original = a.candidate_projector

    def nonfinite(kind="hole", axis=0, band="lower"):
        projector = original(kind, axis, band)
        if kind == "particle" and axis == 2 and band == "upper":
            projector[0, 0] = np.nan
        return projector

    monkeypatch.setattr(a, "candidate_projector", nonfinite)
    with pytest.raises(ArithmeticError, match="numerical identity"):
        a.certificate()


def test_same_candidate_certificate_and_no_passive_protection_claim():
    result = a.certificate()
    assert result["numerical_identity_admission"] == "Frobenius residual <= 2e-11"
    assert result["fresh_numerical_identity_checks_passed"] is True
    assert result["complete_exact_permutation_columns"] == 256
    assert result["neutral_current_pair_exchange_parity"] == -1
    assert result["number_read_lower_doublet_survival"] == .75
    assert result["passive_topological_protection"] is False


@pytest.mark.parametrize("corruption", ["duplicate", "improper", "nongeometry", "bool", "float"])
def test_native_dictionary_rejects_corrupted_warm_rotation_catalogue(monkeypatch, corruption):
    from ..common_source import slow_geometry
    carrier = slow_geometry.source()
    original = tuple(carrier.rotations())
    rows = list(original)
    if corruption == "duplicate":
        rows = [tuple(range(12))]*60
    elif corruption == "improper":
        # Antipodal inversion preserves all edges but reverses face orientation.
        replacement = tuple(carrier.antipode())
        assert replacement not in set(original)
        rows[0] = replacement
    elif corruption == "nongeometry":
        replacement = (1, 0, *range(2, 12))
        assert replacement not in set(original)
        rows[0] = replacement
    else:
        row = list(rows[0])
        index = row.index(0)
        row[index] = False if corruption == "bool" else 0.0
        rows[0] = tuple(row)
    monkeypatch.setattr(carrier, "rotations", lambda: tuple(rows))
    with pytest.raises(ArithmeticError, match="proper native rotation catalogue"):
        a.native_intertwiner_certificate()


def test_native_dictionary_accepts_all_proper_rotations_in_any_iteration_order(monkeypatch):
    from ..common_source import slow_geometry
    carrier = slow_geometry.source()
    expected = a.native_intertwiner_certificate()
    rotations = tuple(reversed(carrier.rotations()))
    monkeypatch.setattr(carrier, "rotations", lambda: rotations)
    assert a.native_intertwiner_certificate() == expected


@pytest.mark.parametrize("field", ["faces", "seams"])
def test_native_dictionary_binds_warm_face_and_seam_chart(monkeypatch, field):
    from ..common_source import slow_geometry
    carrier = slow_geometry.source()
    if field == "faces":
        changed = tuple((a, c, b) for a, b, c in carrier.oriented_faces())
        monkeypatch.setattr(carrier, "oriented_faces", lambda: changed)
    else:
        changed = tuple(reversed(carrier.seams()))
        monkeypatch.setattr(carrier, "seams", lambda: changed)
    with pytest.raises(ArithmeticError, match="face/seam chart"):
        a.native_intertwiner_certificate()


@pytest.mark.parametrize("field", ["native_states", "states", "pairs", "index", "coherent_reorder"])
def test_native_dictionary_binds_literal_integer_state_chart(monkeypatch, field):
    from ..exchange_sector import statistics
    from ..observer import native
    if field == "coherent_reorder":
        states = tuple(reversed(statistics.STATES))
        pairs = tuple(reversed(statistics.PAIRS))
        monkeypatch.setattr(statistics, "STATES", states)
        monkeypatch.setattr(statistics, "PAIRS", pairs)
        monkeypatch.setattr(statistics, "INDEX", {pair: i for i, pair in enumerate(pairs)})
    elif field == "index":
        index = dict(statistics.INDEX)
        index[statistics.PAIRS[0]] = False
        monkeypatch.setattr(statistics, "INDEX", index)
    elif field == "pairs":
        pairs = list(statistics.PAIRS)
        pairs[0] = tuple(float(port) for port in pairs[0])
        monkeypatch.setattr(statistics, "PAIRS", tuple(pairs))
    else:
        states = list(statistics.STATES)
        states[0] = tuple(bool(value) for value in states[0])
        if field == "native_states":
            monkeypatch.setattr(native, "states", lambda binary=False: tuple(states))
        else:
            monkeypatch.setattr(statistics, "STATES", tuple(states))
    with pytest.raises(ArithmeticError, match="binary state/pair/index chart"):
        a.native_intertwiner_certificate()


@pytest.mark.parametrize("field", ["coin", "raw", "abort", "law", "seam", "inputs", "outputs", "signs"])
def test_native_dictionary_rejects_equal_valued_wrong_branch_types(monkeypatch, field):
    from ..exchange_sector import statistics
    original = statistics.branch

    class StringAlias(str):
        pass

    def changed(*args):
        branch = original(*args)
        value = getattr(branch, field)
        if field == "coin":
            wrong = bool(value)
        elif field == "raw":
            wrong = float(value)
        elif field == "abort":
            wrong = int(value)
        elif field == "law":
            wrong = StringAlias(value)
        else:
            wrong = tuple(float(item) for item in value)
        return replace(branch, **{field: wrong})

    monkeypatch.setattr(statistics, "branch", changed)
    with pytest.raises(ArithmeticError, match="strict integer native signed branch"):
        a.native_intertwiner_certificate()


@pytest.mark.parametrize("corruption", ["admit_empty", "wrong_contrast", "bool_key"])
def test_native_dictionary_independently_binds_shared_record_decoder(monkeypatch, corruption):
    from ..observer.source_discrimination import circuit
    original = circuit.record_values

    def changed(law, coin):
        decoder = dict(original(law, coin))
        if corruption == "admit_empty":
            decoder[3] = 0
        elif corruption == "wrong_contrast":
            decoder[0] = 2
        else:
            value = decoder.pop(0)
            decoder[False] = value
        return decoder

    monkeypatch.setattr(circuit, "record_values", changed)
    with pytest.raises(ArithmeticError, match="integer native record decoder"):
        a.native_intertwiner_certificate()
