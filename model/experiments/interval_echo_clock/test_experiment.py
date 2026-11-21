import hashlib
import json

import numpy as np
from scipy.linalg import expm

from . import experiment as ex


def test_full_native_pointer_matches_reduced_dimer_and_its_zero_control():
    law = ex.clock.law
    parts = law.local_term_parts()
    full = parts[0] + ex.clock.KAPPA * parts[1]
    birth = np.eye(law.DIM)[:, law.SRC.index((2, 0))]
    p = ex.native_probe()
    for setting, analyzer in enumerate((0, np.pi / 2)):
        state = expm(-1j * (3 * ex.PULSE + analyzer) * full) @ birth
        assert abs(abs(state[law.SRC.index((2, 0))])**2 - p[setting, 3]) < 1e-12
    assert np.allclose(p[:, 0], (1, .5), rtol=0, atol=1e-12)


def test_uncoupled_primary_clock_and_one_opportunity_closing():
    words, weights = ex.joint_marks((False, True, False, True))
    assert np.isclose(weights @ words[:, 0], 2.)
    home_phase = 4 * ex.PULSE / 2
    assert abs(np.angle(weights @ np.exp(1j * ex.PULSE * words[:, 0])) - home_phase) < 1e-12
    assert weights @ words[:, 1] < 2
    words, weights = ex.joint_marks((True, True))
    assert weights @ words[:, 1] == 0


def test_receipt_can_be_read_from_local_words_and_complementary_pointer_counts():
    archived = json.loads((ex.HERE / 'results.json').read_text())
    assert archived['protocol_sha256'] == hashlib.sha256((ex.HERE / 'PROTOCOL.md').read_bytes()).hexdigest()
    for study in archived['studies']:
        for arm in study['readouts'].values():
            sums = np.array(arm['block_sums'])
            assert sums.shape == (2, 32, 4)
            counts = sums[:, :, 2:].sum(axis=1)
            total = archived['preparations_per_setting']
            z = 2 * counts[0] / total - 1
            y = 1 - 2 * counts[1] / total
            angles = np.arctan2(y, z)
            assert abs(angles[1] / angles[0] - arm['readings']['dimer_ratio']['value']) < 1e-9
            assert abs(sums[:, :, 1].sum() / sums[:, :, 0].sum()
                       - arm['readings']['registered_ratio']['value']) < 1e-9



def test_independent_replication_uses_frozen_settings_and_preserves_both_verdicts():
    receipt = json.loads((ex.HERE / 'replication.json').read_text())
    protocol = (ex.HERE / 'REPLICATION.md').read_bytes()
    assert receipt['protocol_sha256'] == hashlib.sha256(protocol).hexdigest()
    assert receipt['seed'] != ex.SEED
    assert receipt['preparations_per_analyzer_and_rule'] == 8 * ex.BLOCKS * ex.BLOCK_SIZE
    for name, row in receipt['readouts'].items():
        sums = np.asarray(row['block_sums'])
        block_size = receipt['preparations_per_analyzer_and_rule'] // ex.BLOCKS
        reconstructed = ex.read_ratios(sums.mean(axis=1) / block_size)
        endpoint = row['readings']['dimer_minus_registered']
        assert abs(reconstructed[2] - endpoint['value']) < 1e-9
        if name == 'accepted_tick':
            assert abs(endpoint['value']) + 4 * endpoint['se'] < ex.TOLERANCE
        else:
            assert name == 'positive_primary'
            assert endpoint['value'] - 4 * endpoint['se'] > ex.TOLERANCE
