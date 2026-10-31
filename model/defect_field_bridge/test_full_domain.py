"""Independent CAR, complete-outcome and native-restriction controls."""
import numpy as np
import pytest

from ..defect_protection import protection
from ..exchange_sector import statistics as native
from . import full_domain as full


def _car_hop(mask, origin, target):
    if not mask >> origin & 1 or mask >> target & 1:
        return None
    sign = (-1)**(mask & ((1 << origin)-1)).bit_count()
    middle = mask ^ (1 << origin)
    sign *= (-1)**(middle & ((1 << target)-1)).bit_count()
    return middle | (1 << target), sign


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_fock_sectors_and_car_signs(law):
    sites, flavors, seam = 4, 2, (0, 3)
    rows = full.finite_branches(sites=sites, flavors=flavors, law=law, seam=seam)
    assert len(rows) == 32
    counts = np.zeros(256, dtype=int)
    output_counts = np.zeros(256, dtype=int)
    for (coin, raw), row in rows.items():
        assert len({output for _, output, _ in row}) == len(row)
        for mask, output, sign in row:
            counts[mask] += 1
            output_counts[output] += 1
            assert output.bit_count() == mask.bit_count()
            if mask != output:
                removed = (mask & ~output).bit_length()-1
                added = (output & ~mask).bit_length()-1
                assert _car_hop(mask, removed, added) == (output, sign)
                assert removed % flavors == added % flavors
            else:
                assert sign == 1
            if raw == full.COLLISION_RAW:
                assert (mask & 3).bit_count() > 1 or ((mask >> 6) & 3).bit_count() > 1
                assert output == mask
    assert set(counts) == {2}
    assert set(output_counts) == {2}
    # A retained spectator introduces a negative hopping coefficient.
    assert any(sign == -1 for row in rows.values() for _, _, sign in row)


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('flavors', (2, 3))
def test_exact_existing_colored_restriction_on_every_native_seam(law, flavors):
    checked = 0
    for row in protection.branches(law):
        for origin, target, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            x, y = native.PAIRS[origin]
            tx, ty = native.PAIRS[target]
            for a in range(flavors):
                for b in range(flavors):
                    mask = (1 << (x*flavors+a)) | (1 << (y*flavors+b))
                    out_a, out_b = (b, a) if sign == -1 else (a, b)
                    expected = (1 << (tx*flavors+out_a)) | (1 << (ty*flavors+out_b))
                    raw, output, actual_sign, _ = full.branch_action(
                        mask, law, row.coin, sites=12, flavors=flavors, seam=row.seam)
                    assert (raw, output, actual_sign) == (row.raw, expected, sign)
                    checked += 1
    assert checked == 30*2*66*flavors**2


@pytest.mark.parametrize('law', ('C', 'E'))
def test_complete_map_with_entangled_reference(law):
    # Independent Stinespring assembly, complete field-reference input.
    rows = full.finite_branches(sites=2, flavors=2, law=law, seam=(0, 1))
    operators = []
    for row in rows.values():
        k = np.zeros((16, 16))
        for source, target, sign in row:
            k[target, source] = sign/np.sqrt(2)
        operators.append(k)
    assert np.allclose(sum(k.T@k for k in operators), np.eye(16), atol=1e-15)
    psi = np.zeros((16, 3), complex)
    psi[0, 0], psi[3, 1], psi[6, 2] = 1, 1j, -2
    psi /= np.linalg.norm(psi)
    result = sum(np.outer((k@psi).ravel(), (k@psi).ravel().conj()) for k in operators)
    assert abs(np.trace(result)-1) < 1e-15
    assert np.linalg.eigvalsh(result).min() > -1e-15
    before_reference = psi.T@psi.conj()
    after_reference = np.trace(result.reshape(16, 3, 16, 3), axis1=0, axis2=2)
    assert np.allclose(after_reference, before_reference, atol=1e-15)


def test_collision_hold_reports_and_dephases_the_added_domain():
    psi = np.zeros(16)
    psi[0] = psi[3] = 1/np.sqrt(2)
    result = full.channel(np.outer(psi, psi), sites=2, flavors=2, law='C', seam=(0, 1))
    assert result[0, 3] == 0
    assert np.isclose(result[0, 0], .5) and np.isclose(result[3, 3], .5)
    assert np.isclose(np.trace(result@result), .5)
    assert full.branch_action(3, 'C', 0, sites=2, flavors=2, seam=(0, 1))[0] == 8


@pytest.mark.parametrize('law', ('C', 'E'))
@pytest.mark.parametrize('sites,flavors', ((2, 2), (3, 1)))
def test_literal_source_complete_input_isometry(sites, flavors, law):
    actual, target, resources = full.source_isometry(sites, flavors, law)
    assert np.linalg.norm(actual-target) < full.SOURCE_BOUND
    assert np.linalg.norm(actual.conj().T@actual-np.eye(target.shape[1])) < 2*full.SOURCE_BOUND
    assert resources['proper_code_service_events'] == 2*resources['native_pulses']
    assert resources['proper_code_complement_outcome_slots'] > 0


@pytest.mark.parametrize('mask,coin', ((True, 0), (-1, 0), (0, True), (0, 2)))
def test_invalid_masks_and_coins_rejected(mask, coin):
    with pytest.raises(ValueError):
        full.branch_action(mask, 'C', coin, sites=2, flavors=2, seam=(0, 1))


def test_general_large_domain_has_no_integer_overflow_or_number_cutoff():
    mask = (1 << 1079) | (1 << 10) | (1 << 11)
    raw, output, sign, current = full.branch_action(mask, 'E', 1, seam=(0, 11))
    assert (raw, output, sign, current) == (8, mask, 1, 0)
    with pytest.raises(ValueError):
        full.finite_branches(sites=12, flavors=90, law='C', seam=(0, 5))
