import json
from collections import Counter
from itertools import combinations
from math import comb, log

import numpy as np

from model.executable_model import federation
from . import coincidence, measure


def test_exact_sharp_word_reduction_and_isolated_tag():
    v, m, n = 12, 5, 6
    sharp = Counter(sum(i < m for i in word) for word in combinations(range(v), n))
    k, prob, logs = measure.class_law(v, m, n)
    np.testing.assert_allclose(prob, [sharp[j]/comb(v,n) for j in k])
    np.testing.assert_allclose(logs, [log(comb(m,j)) for j in k], atol=1e-14)
    p = measure.tag_probability()
    for j, flags in zip(k, coincidence.pair_laws(k, m, p), strict=True):
        read = np.zeros(4)
        for a,b in combinations(range(m), 2):
            for i,l in ((a,b),(b,a)):
                x,y = p*(i >= j), p*(l >= j)
                read += [(1-x)*(1-y), (1-x)*y, x*(1-y), x*y]
        np.testing.assert_allclose(read/(m*(m-1)), flags, atol=2e-14)


def test_native_shape_count_and_recorded_hypotheses():
    archived = json.loads((measure.HERE/'results.json').read_text())
    for panel in archived['panels']:
        parent, depth = panel['parent_level'], panel['microscopic_subdivisions']
        for g in panel['geometries']:
            sizes = {'cycle5':5, 'path5':7, 'cycle5_spur':6, 'path6':8}
            m, v = sizes[g['name']]*2**depth, len(g['parent_cells'])*12*4**depth
            if parent+depth <= 2:
                board = federation(parent+depth)
                edges = np.array(board.external)
                members = np.isin(edges//(12*4**depth),g['parent_cells'])
                assert int((members[:,0] ^ members[:,1]).sum()) == m
            for r in g['sources']:
                assert (r['m'],r['v']) == (m,v)
                assert sum(r['class_counts']) == measure.INFO_SHOTS
            matched = g['sources'][2]['matched_variation']
            assert matched['exact_tag_slope'] == 0
            assert matched['exact_bulk_slope'] > 0
        if (parent,depth)==(3,4):
            for r in panel['shape_transfers']:
                if r['filling']=='half':
                    assert r['frozen_half_fill_relative_error'] < .05
                elif r['filling']=='three_quarters':
                    assert r['frozen_half_fill_relative_error'] > .5


def test_two_tag_complete_flags_and_independent_held_counts():
    observed = coincidence.run()
    assert observed == json.loads((measure.HERE/'coincidence.json').read_text())
    for row in observed['panels']:
        h = row['variance_transfer']
        assert h['exact_variance_slope'] > 0
        assert abs(h['flag_slopes'][1] + h['flag_slopes'][3]) < 1e-14
        blind = row['fixed_two_ray_law']
        np.testing.assert_allclose(blind['flag_slopes'], 0, atol=2e-15)
        assert blind['exact_bulk_slope'] > 0
        assert blind['exact_mean_energy_slope'] == '0'
