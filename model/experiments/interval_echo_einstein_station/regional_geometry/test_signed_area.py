import json

import numpy as np

from . import measure, signed_counts
from ..collar_scaling import signed_area as protected


def test_same_marginals_both_volumes_and_increasing_areas():
    rows = json.loads((measure.HERE/'signed_area.json').read_text())['boards']
    original = json.loads((measure.HERE/'results.json').read_text())
    ref = np.array(original['reference'])
    for row,board in zip(rows, original['screens'],strict=True):
        laws = np.asarray(board['word_moments'])
        assert row['fixed_volume_area_cosine'] < -.99
        for name, case in row['witnesses'].items():
            finite = case['finite']
            assert finite['max_fractional_volume_error'] < 3e-50
            assert max(abs(v) for v in case['transferred_volume_slopes']) < 6e-12
            assert case['area_slope_enclosures'][0][0] > 0
            if name == 'both_increase':
                assert case['area_slope_enclosures'][1][0] > 0
                assert case['charge_slope_enclosure'][0] > .16
            else:
                assert case['area_slope_enclosures'][1][1] < 0
            pair = [np.array(p).reshape(3,3) for p in finite['probabilities']]
            assert min(p.min() for p in pair) > 0
            for p in pair:
                np.testing.assert_allclose(p.sum(axis=0),ref.sum(axis=0),atol=2e-16)
                np.testing.assert_allclose(p.sum(axis=1),ref.sum(axis=1),atol=2e-16)
            for i in range(2):
                reads = [measure.geometry(laws[:,:,:,i,:],p,np.zeros_like(p)) for p in pair]
                np.testing.assert_allclose(reads[0]['volume'],reads[1]['volume'],rtol=1e-12)
                secant = 16*(reads[1]['surface_area']-reads[0]['surface_area'])
                np.testing.assert_allclose(secant,finite['area_secants'][i],atol=7e-12,rtol=0)
            if 'sharp_measurement' in case:
                sharp = case['sharp_measurement']
                assert abs(sharp['charge_secant']-finite['fixed_reference_charge_secant']) < 3*sharp['counting_standard_error']


def test_rare_clock_bins_and_counting_resolution():
    rng = np.random.default_rng(7124)
    shots = 2**60
    p = np.array([4e-17, 1e-15, .125, 1-.125-1.04e-15])
    counts = signed_counts.clock_counts(rng,shots,p)
    assert sum(counts) == shots
    for n,q in zip(counts[:2], p[:2],strict=True):
        assert n > 0 and abs(n-shots*q) < 6*np.sqrt(shots*q)
    archived = json.loads((measure.HERE/'signed_counts.json').read_text())
    for name, row in archived['menus'].items():
        case = row['both_increase']
        assert case['projected_shots_per_edge_arm_for_five_se'] > 1e25
        low, high = case['acquisitions']
        if name=='first_return':
            assert low['observed_area_change'] is None
        assert abs(high['observed_area_change']-case['population_area_change']) < 3*high['asymptotic_counting_se']
        for arm in high['counts']:
            assert all(sum(edge)==high['shots_per_edge_arm'] for edge in arm)


def test_protected_annulus_charge_at_identical_routed_ranges():
    result = protected.run()
    assert abs(result['masked_ordinary_contact']['contact_gap_mean']-8.5) < 2e-14
    for row in result['boards']:
        assert row['global_tag']['per_coupling_lower'] > 0
        assert row['ancestor_tag']['per_coupling_lower'] > 0
    last = result['boards'][-1]
    assert last['global_tag']['per_coupling_lower'] > .00129
    assert last['ancestor_tag']['per_coupling_lower'] > .000245
