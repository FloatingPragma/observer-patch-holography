import hashlib
import json

import numpy as np

from . import distance, free_selection, measure


def receipt(name):
    return json.loads((measure.HERE / name).read_text())


def test_sector_counts_and_regional_entropy():
    record = receipt('results.json')
    patch = measure.search.Patch()
    words = np.array(patch.field_component())
    for outline, data in zip(measure.shapes(patch),record['regions'],strict=True):
        region = measure.Region(patch,words,outline)
        assert region.complete and region.exterior_reads_center
        np.testing.assert_allclose(np.sort(region.schmidt),np.sort(region.center_counts/len(words)),
                                   atol=1e-14,rtol=0)
        assert data['contact'] == measure.contact_survey(patch,words,outline)
        np.testing.assert_allclose(data['derivative_by_head']['center'],data['derivative_by_head']['modular'],
                                   atol=2e-13,rtol=0)
        np.testing.assert_allclose(np.array(data['derivative_by_head']['interior_word'])
                                   + data['derivative_by_head']['exterior_completion'],0,atol=2e-13)
        directions=np.array(data['free_crossing_vacuum_and_derivatives'])[:,1:]
        assert np.max(abs(directions@data['fitted_phase_direction']))<2e-13
        internal=np.array([set(edge).issubset(outline) for edge in record['head_order']])
        assert np.linalg.matrix_rank(directions[:,internal],tol=1e-12) == internal.sum()
        assert data['blind_entropy_slope']>.03
    assert record['observable_phase_rank'] == 19


def test_held_flight_and_independently_scored_endpoints():
    experiment = receipt('results.json')
    for key,filename in [('protocol_sha256','PROTOCOL.md'), ('endpoint_protocol_sha256','ENDPOINT_PROTOCOL.md'),
                         ('distance_protocol_sha256','DISTANCE_PROTOCOL.md')]:
        assert experiment[key]==hashlib.sha256((measure.HERE/filename).read_bytes()).hexdigest()
    reader=experiment['endpoint_calibration']
    sign_index=measure.SIGNS.index(reader['sign'])
    time_index=measure.TIMES.tolist().index(reader['exposure'])
    score=np.isin(np.arange(12),reader['decoder_ports']).astype(float)
    for h in experiment['held']:
        counts=np.array(h['endpoints']['full_endpoint_counts'])
        assert (counts.sum(axis=3)==experiment['shots']).all()
        table=counts[sign_index,reader['origin'],time_index].T
        result=measure.signed_read(table,score,h['angle'])
        assert abs(result['measured_secant']/result['standard_error'])>200
        assert result['measured_secant']==h['endpoints']['measured_secant']
        for reader_name in ('modular','center'):
            result=h['entropy'][reader_name]
            assert abs(result['measured_secant']-.159127)<4*result['standard_error']
    for chart in experiment['distance_charts']['readings']:
        assert min(w['counts']['largest_held_relative_error'] for w in chart['perturbations'])>.34
    coordinates=np.random.default_rng(10).normal(size=(12,3))
    assert distance.chart(distance.square_distances(coordinates))['largest_held_relative_error']<1e-11


def test_matched_conserved_moments_do_not_fix_local_long_time_read():
    record=free_selection.run()
    archive=receipt('free_selection.json')
    assert record['ground_dimension']==1
    assert record['ground_ring_probability']<1e-24
    assert record['mixed_stationary_ring_probability']>.16
    local=[]
    for row,old in zip(record['records'],archive['records'],strict=True):
        assert row['word']==old['word'] and row['mean_energy']==4 and row['second_energy_moment']==18
        assert row['ground_probability']==1/458
        assert abs(row['limit_probability']-row['finite_probability'])<2e-5
        np.testing.assert_allclose(row['finite_probability'],old['finite_probability'],rtol=0,atol=1e-11)
        local.append(row['limit_probability'])
    assert local[1]-local[0]>.068
